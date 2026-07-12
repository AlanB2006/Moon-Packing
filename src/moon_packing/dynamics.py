from __future__ import annotations

import math
import numpy as np

from .constants import A_EARTH, M_EARTH, R_EARTH, RHO_EARTH, MoonModel


def hill_radius(planet_semimajor_axis: float = A_EARTH, planet_mass: float = M_EARTH,
                star_mass: float = 1.0) -> float:
    """Planetary Hill radius in the same length units as ``planet_semimajor_axis``."""
    return planet_semimajor_axis * (planet_mass / (3.0 * star_mass)) ** (1.0 / 3.0)


def roche_radius(model: MoonModel, planet_radius: float = R_EARTH,
                 planet_density: float = RHO_EARTH) -> float:
    """Fluid Roche radius in the same units as ``planet_radius``."""
    return 2.44 * planet_radius * (planet_density / model.density_g_cm3) ** (1.0 / 3.0)


def mutual_hill_x(model: MoonModel, index: int, planet_mass: float = M_EARTH) -> float:
    """Mass factor X used in the mutual-Hill spacing relation."""
    return 0.5 * ((2.0 * model.mass_solar) /
                  (3.0 * (planet_mass + index * model.mass_solar))) ** (1.0 / 3.0)


def satellite_semimajor_axes(beta: float, n_moons: int, model: MoonModel,
                             first_semimajor_axis: float | None = None,
                             mode: str = "archived") -> np.ndarray:
    """Return initial moon semimajor axes.

    ``mode='archived'`` exactly reproduces the expression in the supplied
    simulation scripts. ``mode='recursive'`` applies the recursion described
    in the manuscript equation one adjacent pair at a time.
    """
    if n_moons < 1:
        raise ValueError("n_moons must be at least 1")
    a0 = first_semimajor_axis or 2.0 * roche_radius(model)
    axes = [a0]
    if mode == "archived":
        for i in range(1, n_moons):
            x = mutual_hill_x(model, i)
            ratio = (1.0 + beta * x) / (1.0 - beta * x)
            if ratio <= 0:
                raise ValueError("beta produces a non-physical orbital spacing")
            axes.append(a0 * ratio**i)
    elif mode == "recursive":
        for i in range(1, n_moons):
            x = mutual_hill_x(model, i)
            ratio = (1.0 + beta * x) / (1.0 - beta * x)
            if ratio <= 0:
                raise ValueError("beta produces a non-physical orbital spacing")
            axes.append(axes[-1] * ratio)
    else:
        raise ValueError("mode must be 'archived' or 'recursive'")
    return np.asarray(axes, dtype=float)


def golden_phase(index: int) -> float:
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    return (index * phi * 2.0 * math.pi) % (2.0 * math.pi)


def beta_max(n_moons: int, model: MoonModel, a_outer: float | None = None,
             a_inner: float | None = None) -> float:
    """Idealized maximum beta from equation 5 in the manuscript."""
    if n_moons < 2:
        raise ValueError("n_moons must be at least 2")
    a1 = a_inner or 2.0 * roche_radius(model)
    a_n = a_outer or 0.4 * hill_radius()
    r = (a_n / a1) ** (1.0 / (n_moons - 1))
    return ((r - 1.0) / (r + 1.0)) * (12.0 * M_EARTH / model.mass_solar) ** (1.0 / 3.0)


def resonance_beta(period_ratio: float, model: MoonModel) -> float:
    """Spacing beta corresponding to an adjacent period ratio (equation 15)."""
    if period_ratio <= 1.0:
        raise ValueError("period_ratio must exceed 1")
    x = mutual_hill_x(model, 1)
    y = period_ratio ** (2.0 / 3.0)
    return (y - 1.0) / (x * (y + 1.0))
