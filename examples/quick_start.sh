#!/usr/bin/env bash
set -euo pipefail

python -m pip install -e '.[test]'
moon-packing-figures --figures 2 3 4 8 9 10
python -m pytest
