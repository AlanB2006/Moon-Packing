# Reproducibility notes

## Numerical model

The repository encodes the revised manuscript setup:

- Sun mass: \(1\,M_\odot\).
- Planet mass: \(3.003\times10^{-6}\,M_\odot\), semimajor axis 1 au, eccentricity 0.01671022.
- Planet radius 6371 km, density 5.515 g cm\(^{-3}\), normalized moment of inertia 0.3308, Love number 0.298, obliquity 23.44°, initial spin period 6 hr.
- Inner moon begins at twice the fluid Roche radius.
- Initial phases follow the golden-ratio prescription.
- REBOUND IAS15, with the initial timestep based on 5% of the period at 1.8 Roche radii.
- REBOUNDx `tides_spin` constant-time-lag model with \(\tau=0,100,698\) s.

## REBOUND SimulationArchive inputs for legacy `Fig02.py`

The legacy `Fig02.py` workflow requires binary REBOUND SimulationArchive files (`*.bin`) in `Data/Simulation_Archives/`. The normalized path used by this repository is `data/Simulation_Archives/`. Preserve the original filenames because the plotting workflow may refer to them directly.  

These archives contain saved simulation states and are a separate data product from:

- the lifetime and maximum-eccentricity tables in `data/sweeps/`; and
- the eccentricity-grid tables in `data/eccentricity_maps/`.

The SimulationArchive files are too big to host on GitHub and are therefore hosted on [Zenodo](https://doi.org/10.5281/zenodo.20408204) instead.

## Production durations

- Lifetime sweeps: \(10^7 P_1\).
- Eccentricity maps: \(10^5 P_1\).
- Time-series examples: 25 yr with output every 5 \(P_1\).

The complete production grids are expensive. Run the laptop examples first and estimate resource use from `wall_time` before submitting a cluster campaign.


## Randomness and determinism

There is no random seed in the initial conditions. Moon phases are deterministic. Parallel jobs may finish in a different order, but `write_results` sorts rows by initial eccentricity and beta before writing.

## Exact environment

The manuscript reports REBOUND 4.4.6 and REBOUNDx 4.4.1. `environment.yml` and the simulation extra pin these versions. The plotting and unit tests do not require these packages.
