#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/generated/timeseries

for beta in 4.4 8.2 8.55; do
  moon-packing-run-timeseries --moon-type Luna --n-moons 2 --tau 698 \
    --beta "$beta" --years 25 \
    --output "data/generated/timeseries/luna_beta${beta}_tau698.npz"
done

for beta in 4.3 9 9.7; do
  moon-packing-run-timeseries --moon-type Pluto --n-moons 3 --tau 698 \
    --beta "$beta" --years 25 \
    --output "data/generated/timeseries/pluto_beta${beta}_tau698.npz"
done

for beta in 6 13.2 13.55; do
  moon-packing-run-timeseries --moon-type Ceres --n-moons 5 --tau 698 \
    --beta "$beta" --years 25 \
    --output "data/generated/timeseries/ceres_beta${beta}_tau698.npz"
done

moon-packing-figures --figures 5 6 7
