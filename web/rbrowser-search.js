/**
 * RBrowser Index (.rbi) static search — frontend reference implementation.
 *
 * Mirrors refbox.sqlite_index.search() (refbox v0.5.5) exactly, against a
 * read-only *.rbi (SQLite+FTS5) served as a static file and queried with sql.js-httpvfs
 * (SQLite WASM + HTTP Range VFS). No backend.
 *
 * Schema it expects (v0.5.5):
 *   feature(id, feature_type, gene_id, gene_name, transcript_id, transcript_name,
 *           chrom, start, end, chrom_start0, chrom_end0, strand, biotype, source,
 *           exon_count, exon_starts, exon_ends, cds_start, cds_end,
 *           utr5_start, utr5_end, utr3_start, utr3_end, payload_json)
 *      start/end = 1-based inclusive; chrom_start0/end0 = 0-based half-open.
 *   alias(feature_id, alias_norm, alias_type)        -- only alias_norm is matched
 *      idx_alias_norm(alias_norm, alias_type, feature_id)  -- covering
 *   feature_fts      -- FTS5 prefix/autocomplete (names + IDs + synonyms)
 *   feature_trigram  -- FTS5 trigram substring (NAMES + synonyms only, NO IDs)
 *
 * ── SAFETY RULES (violating any one turns a query into a full scan, i.e. a
 *    full-DB download over HTTP Range) ─────────────────────────────────────
 *   1. Exact column compares MUST use `= ? COLLATE NOCASE` (matches the index).
 *      Never wrap a column in lower()/upper(); never compare without COLLATE.
 *   2. Substring/prefix only via FTS `MATCH`. Never `LIKE '%x%'`.
 *   3. Always parameterized — never string-concat user input into SQL.
 *   4. A byte budget on the VFS (maxBytesToRead) aborts a runaway query instead
 *      of downloading the whole file.
 */

import { createDbWorker } from "sql.js-httpvfs";

const DB_URL =
  "https://data.rbrowser.org/datahub/reference/hg38.gencode.v45.transcript.rbi";

// ── one-time worker setup ──────────────────────────────────────────────────
let _worker = null;
export async function initSearch(url = DB_URL) {
  if (_worker) return _worker;
  _worker = await createDbWorker(
    [{
      from: "inline",
      config: {
        serverMode: "full",
        url,
        requestChunkSize: 4096,        // = the DB page_size
      },
    }],
    // adjust these paths to where you host the sql.js-httpvfs assets:
    new URL("sqlite.worker.js", import.meta.url).toString(),
    new URL("sql-wasm.wasm", import.meta.url).toString(),
    { maxBytesToRead: 10 * 1024 * 1024 } // safety budget per query (~10 MB)
  );
  return _worker;
}

// ── normalization (must match Python normalize() byte-for-byte) ────────────
//   lowercase → strip a trailing Ensembl ".<version>" → drop separators _ - . space
export function normalize(v) {
  if (!v) return "";
  return v.trim().toLowerCase()
    .replace(/\.\d+$/, "")
    .replace(/[\s_.\-]+/g, "");
}

const COLS = `id, feature_type, gene_name, gene_id, transcript_name,
  transcript_id, chrom, start, end, strand, biotype, payload_json`;

// alias_type → [rank score, matched_field].  Lower score = better.
// Anything not listed (gene_synonym, hgnc, havana, ccds, protein_id,
// rnacentral_id, rnacentral_db, dbxref, refseq, name, alias) → generic alias_exact.
const ALIAS_FIELD = {
  transcript_id:   [1, "transcript_id_exact"],
  transcript_name: [2, "transcript_name_exact"],
  gene_name:       [3, "gene_name_exact"],
  gene_id:         [4, "gene_id_exact"],
};

// sanitize a user string into ONE safe FTS token (strip quotes, collapse spaces)
function ftsToken(q) {
  return q.replace(/["']/g, " ").trim().replace(/\s+/g, " ");
}

// rerank fuzzy (prefix/trigram) candidates: startswith → contains → shortest name
function rerank(rows, nq) {
  const key = (r) => {
    const name = r.transcript_name || r.gene_name || r.transcript_id || r.gene_id || "";
    const nn = normalize(name);
    return [nn.startsWith(nq) ? 0 : 1, nq && nn.includes(nq) ? 0 : 1, name.length, name];
  };
  return rows.slice().sort((a, b) => {
    const ka = key(a), kb = key(b);
    for (let i = 0; i < ka.length; i++) if (ka[i] !== kb[i]) return ka[i] < kb[i] ? -1 : 1;
    return 0;
  });
}

/**
 * Ranked search. Returns up to `limit` result objects, each with rank_score and
 * matched_field. Priority:
 *   exact transcript_id → transcript_name → gene_name → gene_id → alias →
 *   prefix (FTS) → trigram substring (names/synonyms only).
 */
export async function search(query, limit = 10) {
  const w = await initSearch();
  const db = w.db;
  query = (query || "").trim();
  if (!query) return [];
  const nq = normalize(query);

  // ── Tiers 1–5: ONE index-only seek into idx_alias_norm ───────────────────
  // Resolves any exact match (id / name / synonym / RNAcentral URS / …) across
  // all species in a single lookup; alias_type gives rank + matched_field.
  const cand = await db.query(
    "SELECT alias_type, feature_id FROM alias WHERE alias_norm = ? LIMIT ?",
    [nq, Math.max(limit * 8, 50)]
  );
  if (cand.length) {
    const best = new Map(); // feature_id -> [score, field]
    for (const { alias_type, feature_id } of cand) {
      const [score, field] = ALIAS_FIELD[alias_type] || [5, "alias_exact"];
      const cur = best.get(feature_id);
      if (!cur || score < cur[0]) best.set(feature_id, [score, field]);
    }
    const ordered = [...best.entries()]
      .sort((a, b) => a[1][0] - b[1][0] || a[0] - b[0])
      .slice(0, limit);
    const ids = ordered.map(([fid]) => fid);
    const ph = ids.map(() => "?").join(",");
    const rows = await db.query(`SELECT ${COLS} FROM feature WHERE id IN (${ph})`, ids);
    const byId = new Map(rows.map((r) => [r.id, r]));
    const out = [];
    for (const [fid, [score, field]] of ordered) {
      const r = byId.get(fid);
      if (r) out.push({ ...r, rank_score: score, matched_field: field });
    }
    return out; // exact hit: do NOT dilute/slow with fuzzy tiers
  }

  const seen = new Set();
  const out = [];
  const take = (rows, field, score) => {
    for (const r of rows) {
      if (seen.has(r.id)) continue;
      seen.add(r.id);
      out.push({ ...r, rank_score: score, matched_field: field });
      if (out.length >= limit) return true;
    }
    return false;
  };

  // ── Tier 6: prefix / autocomplete via feature_fts (names + IDs + synonyms) ─
  const token = ftsToken(query);
  if (token) {
    try {
      const rows = await db.query(
        `SELECT ${COLS} FROM feature WHERE id IN
           (SELECT rowid FROM feature_fts WHERE feature_fts MATCH ?) LIMIT ?`,
        [`"${token}"*`, limit * 4]
      );
      if (take(rerank(rows, nq), "prefix", 6)) return out;
    } catch (_) { /* FTS phrase error on odd input → ignore, fall through */ }
  }

  // ── Tier 7: trigram substring (NAMES + synonyms only) ────────────────────
  // Skip for letter-less queries: the trigram has no IDs, so a pure digit/symbol
  // fragment like "000003351" can never match a name — running it would just
  // scan postings for nothing (the old pathological multi-second case).
  if (query.length >= 3 && /[a-z]/i.test(query)) {
    try {
      const rows = await db.query(
        `SELECT ${COLS} FROM feature WHERE id IN
           (SELECT rowid FROM feature_trigram WHERE feature_trigram MATCH ?) LIMIT ?`,
        [`"${query.replace(/"/g, "")}"`, limit * 4]
      );
      if (take(rerank(rows, nq), "trigram", 7)) return out;
    } catch (_) { /* ignore */ }
  }

  return out;
}

/**
 * List a gene's transcripts. A gene-name / gene-id query returns the GENE record
 * (v0.5.5 behavior); call this to expand it. Index-backed (idx_feature_gene_id),
 * COLLATE NOCASE so the lookup is a seek, not a scan.
 */
export async function transcriptsOfGene(geneId, limit = 200) {
  const w = await initSearch();
  return w.db.query(
    `SELECT ${COLS} FROM feature
       WHERE gene_id = ? COLLATE NOCASE AND feature_type <> 'gene' LIMIT ?`,
    [geneId, limit]
  );
}

// ── UI helper: debounced search box ────────────────────────────────────────
// Debounce avoids firing a query on every keystroke; require ≥2 chars.
export function attachSearchBox(inputEl, renderResults, { delay = 180, limit = 10 } = {}) {
  let timer = null, seq = 0;
  inputEl.addEventListener("input", () => {
    clearTimeout(timer);
    const q = inputEl.value.trim();
    if (q.length < 2) { renderResults([], q); return; }
    timer = setTimeout(async () => {
      const mine = ++seq;
      const res = await search(q, limit);
      if (mine === seq) renderResults(res, q); // ignore out-of-order responses
    }, delay);
  });
}

/* ── usage ──────────────────────────────────────────────────────────────────
import { initSearch, search, transcriptsOfGene, attachSearchBox } from "./rbrowser-search.js";

await initSearch();
console.log(await search("TP53"));        // → gene TP53            [gene_name_exact]
console.log(await search("oct4"));        // → POU5F1               [alias_exact]
console.log(await search("ENST00000269305")); // versionless OK     [transcript_id_exact]
console.log(await search("URS000035F234"));   // ncRNA              [alias_exact]
console.log(await transcriptsOfGene("ENSG00000141510.18")); // TP53 transcripts

attachSearchBox(document.querySelector("#q"), (rows, q) => { /* render rows */ });
──────────────────────────────────────────────────────────────────────────── */
