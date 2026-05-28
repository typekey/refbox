#!/usr/bin/env bash
# release.sh — build sdist + wheel and (optionally) upload to PyPI
#
# Usage:
#   ./release.sh build           # only build dist/
#   ./release.sh test-upload     # upload to TestPyPI
#   ./release.sh upload          # upload to PyPI (production)
#
# Auth: uses ~/.pypirc, or set TWINE_USERNAME=__token__ and TWINE_PASSWORD=pypi-AgEI...

set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

cmd="${1:-build}"

ensure_tools() {
    pip install --upgrade build twine >/dev/null
}

case "$cmd" in
    build)
        ensure_tools
        rm -rf dist/ build/ src/*.egg-info
        python -m build
        ls -lh dist/
        twine check dist/*
        ;;
    test-upload)
        ensure_tools
        twine upload --repository testpypi dist/*
        echo "Try: pip install --index-url https://test.pypi.org/simple/ refbox"
        ;;
    upload)
        ensure_tools
        twine upload dist/*
        echo "Try: pip install refbox"
        ;;
    *)
        echo "usage: $0 {build|test-upload|upload}" >&2
        exit 2
        ;;
esac
