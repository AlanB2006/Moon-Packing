# Data dictionary

All published tables are comma-separated text files with a commented header.

| Column | Meaning |
|---|---|
| `dH` | Initial mutual-Hill spacing parameter, \(\beta\) |
| `e0` | Initial eccentricity of every moon |
| `eMi_max` | Maximum eccentricity recorded for moon `i` |
| `eMi_min` | Minimum eccentricity recorded for moon `i` |
| `dAi_over_A0` | Fractional semimajor-axis change for moon `i` |
| `t_end(PM)` | End time in periods of the innermost moon |
| `wall_time` | Wall-clock runtime in seconds |
| `fate` | Termination condition |

Observed fate labels include:

- `tmax`: reached the requested integration duration.
- `collision`: REBOUND direct collision.
- `d<R_roche`: moon entered the fluid Roche radius.
- `a>0.4R_H`: semimajor axis exceeded 0.4 planetary Hill radii in the archived implementation.
- `e>1`: osculating eccentricity exceeded unity.
- `a<0`: unbound/hyperbolic osculating semimajor axis.

## Directory split

- `data/sweeps/`: \(e_0\approx0\) integrations to \(10^7 P_1\), used for Figures 2–4.
- `data/eccentricity_maps/`: \((\beta,e_0)\) grids to \(10^5 P_1\), used for Figures 8–10.
- `data/generated/`: outputs created by new runs. These are ignored by Git unless explicitly added.
