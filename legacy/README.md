# Legacy scripts

This directory contains the original Python scripts from the supplied archive,
unchanged. They are retained for provenance, not as the recommended interface.
Several assume Google Colab paths or interactive input and some contain missing
imports or typographical errors. Use the package in `src/moon_packing/` and the
commands in `scripts/` for reproducible runs.

## Binary archive dependency

The original `Fig02.py` workflow uses REBOUND SimulationArchive `*.bin` files stored in `Data/Simulation_Archives/`. In this repository, place them in `data/Simulation_Archives/` and preserve their original filenames. These binary archives are not interchangeable with the text tables under `data/sweeps/` or `data/eccentricity_maps/`.
