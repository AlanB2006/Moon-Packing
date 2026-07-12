# Manuscript figure map

| Figure | Description | Inputs | Reproduction command |
|---:|---|---|---|
| 1 | REBOUNDx vs secular tidal evolution | Generated validation CSVs | `moon-packing-validate-tides && moon-packing-figures --figures 1` |
| 2 | 2 Luna-mass moons: lifetime and maximum eccentricity | `data/sweeps/Luna/` | `moon-packing-figures --figures 2` |
| 3 | 3 Pluto-mass moons: lifetime and maximum eccentricity | `data/sweeps/Pluto/` | `moon-packing-figures --figures 3` |
| 4 | 5 Ceres-mass moons: lifetime and maximum eccentricity | `data/sweeps/Ceres/` | `moon-packing-figures --figures 4` |
| 5 | Luna time series at \(\beta=4.4,8.2,8.55\) | Generated NPZ files | Run the three `run-timeseries` commands, then `moon-packing-figures --figures 5` |
| 6 | Pluto time series at \(\beta=4.3,9.0,9.7\) | Generated NPZ files | Same workflow with Pluto, then `--figures 6` |
| 7 | Ceres time series at \(\beta=6.0,13.2,13.55\) | Generated NPZ files | Same workflow with Ceres, then `--figures 7` |
| 8 | Luna eccentricity-variation maps | `data/eccentricity_maps/Luna/` | `moon-packing-figures --figures 8` |
| 9 | Pluto eccentricity-variation maps | `data/eccentricity_maps/Pluto/` | `moon-packing-figures --figures 9` |
| 10 | Ceres eccentricity-variation maps | `data/eccentricity_maps/Ceres/` | `moon-packing-figures --figures 10` |

The `figures/reference/` directory contains the four PDFs supplied with the original bundle. Their descriptive filenames state the manuscript figure to which each corresponds.

## Legacy script dependency

The original `Fig02.py` workflow additionally uses REBOUND SimulationArchive `*.bin` files from `Data/Simulation_Archives/` (normalized here as `data/Simulation_Archives/`). This is a dependency of the legacy script, not of the table-based figure commands in the map above. Preserve the archive filenames when adding them to the GitHub and Zenodo release.

