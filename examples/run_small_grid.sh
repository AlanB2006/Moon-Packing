#!/usr/bin/env bash
set -euo pipefail

moon-packing-run-grid \
  --moon-type Luna --n-moons 2 --tau 100 \
  --beta-start 6.0 --beta-stop 6.2 --beta-step 0.1 \
  --max-orbits 1000 --workers 2 \
  --output data/generated/luna_demo.csv
