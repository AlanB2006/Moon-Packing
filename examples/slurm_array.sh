#!/usr/bin/env bash
#SBATCH --job-name=moon-pack
#SBATCH --array=0-30
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --output=logs/%x-%A_%a.out

set -euo pipefail
mkdir -p logs data/generated/tiles
source ~/.bashrc
conda activate moon-packing

# One eccentricity row per array task. Each task spans the full beta range.
E0=$(python -c 'import os; print(f"{int(os.environ[\"SLURM_ARRAY_TASK_ID\"])*0.01:.2f}")')

moon-packing-run-grid \
  --moon-type Luna --n-moons 2 --tau 698 \
  --beta-start 3.5 --beta-stop 9.0 --beta-step 0.01 \
  --e-start "$E0" --max-orbits 1e5 --workers "$SLURM_CPUS_PER_TASK" \
  --output "data/generated/tiles/luna_tau698_e${E0}.csv"

# After all array tasks finish, merge the rows (keeping one header) and sort
# using a small Python script or the repository I/O helpers.
