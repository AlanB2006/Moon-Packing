from __future__ import annotations

from dataclasses import dataclass

AU_METERS = 149_597_870_700.0
YEAR_SECONDS = 365.25 * 24.0 * 3600.0
M_SUN = 1.0
M_EARTH = 3.003e-6
A_EARTH = 1.0
E_EARTH = 0.01671022
R_EARTH = 6_371_000.0 / AU_METERS
RHO_EARTH = 5.515
K2_EARTH = 0.298
C_EARTH = 0.3308
EARTH_SPIN_PERIOD_HOURS = 6.0
EARTH_OBLIQUITY_DEG = 23.44


@dataclass(frozen=True)
class MoonModel:
    name: str
    code: str
    mass_solar: float
    mass_earth: float
    radius_km: float
    density_g_cm3: float
    manuscript_count: int
    beta_stop: float

    @property
    def radius_au(self) -> float:
        return self.radius_km * 1_000.0 / AU_METERS


MOON_MODELS = {
    "luna": MoonModel("Luna", "L", M_EARTH / 81.0, 0.0123, 1737.0, 3.3, 2, 9.0),
    "pluto": MoonModel("Pluto", "P", 0.002183 * M_EARTH, 0.0022, 1190.0, 1.86, 3, 11.0),
    "ceres": MoonModel("Ceres", "C", 4.72e-10, 0.00015, 473.0, 2.2, 5, 14.0),
}


def get_moon_model(value: str) -> MoonModel:
    key = value.strip().lower()
    aliases = {"l": "luna", "moon": "luna", "p": "pluto", "c": "ceres"}
    key = aliases.get(key, key)
    if key not in MOON_MODELS:
        choices = ", ".join(m.name for m in MOON_MODELS.values())
        raise ValueError(f"Unknown moon type {value!r}; choose one of {choices}.")
    return MOON_MODELS[key]
