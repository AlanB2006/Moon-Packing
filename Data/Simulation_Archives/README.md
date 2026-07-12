# REBOUND SimulationArchive files

This directory is reserved for the binary REBOUND SimulationArchive files (`*.bin`) used by the original `Fig02.py` plotting workflow.

The original author data bundle refers to this location as `Data/Simulation_Archives/`. In this normalized repository layout, the corresponding location is `data/Simulation_Archives/`.

## Important

- Keep the archive filenames exactly as generated because the legacy plotting script expects specific names.
- Treat these files as binary data. The repository `.gitattributes` already marks `*.bin` as binary.
- These archives are separate from the text tables in `data/sweeps/` and `data/eccentricity_maps/`.
- They are required to run the original `Fig02.py` workflow directly.
- They are not required for the table-based reproduction command:

  ```bash
  moon-packing-figures --figures 2 3 4 8 9 10
  ```

For a public release, include the archive files in the Zenodo deposition even if they are distributed separately from the Git history. Record their filenames, sizes, and SHA-256 checksums in the release manifest.

The `.bin` files themselves were not present in the uploaded bundle used to assemble this repository, so this directory currently contains documentation only.
