from __future__ import annotations

import time
from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

from .constants import (
    A_EARTH, C_EARTH, EARTH_OBLIQUITY_DEG, EARTH_SPIN_PERIOD_HOURS,
    E_EARTH, K2_EARTH, M_EARTH, M_SUN, R_EARTH, YEAR_SECONDS,
    get_moon_model,
)
from .dynamics import golden_phase, hill_radius, roche_radius, satellite_semimajor_axes


@dataclass(frozen=True)
class SimulationConfig:
    moon_type: str
    n_moons: int
    tau_seconds: float = 0.0
    max_orbits: float = 1e7
    spacing_mode: str = "archived"
    outer_boundary: str = "semimajor-axis"
    output_interval_orbits: float | None = None


def _require_rebound() -> tuple[Any, Any]:
    try:
        import rebound
        import reboundx
    except ImportError as exc:
        raise RuntimeError(
            "REBOUND and REBOUNDx are required for integrations. "
            "Install with: python -m pip install -e '.[simulation]'"
        ) from exc
    return rebound, reboundx


def _relative_elements(rebound: Any, primary: Any, satellite: Any, central_mass: float):
    diff = satellite - primary
    temp = rebound.Simulation()
    temp.units = ("AU", "yr", "Msun")
    temp.add(m=central_mass)
    temp.add(m=0.0, x=diff.x, y=diff.y, z=diff.z,
             vx=diff.vx, vy=diff.vy, vz=diff.vz)
    temp.move_to_com()
    p = temp.particles[1]
    return p.a, p.e, p.d, p.P


def _setup_tides(sim: Any, tau_seconds: float, rebound: Any, reboundx: Any) -> Any:
    particles = sim.particles
    extras = reboundx.Extras(sim)
    force = extras.load_force("tides_spin")
    extras.add_force(force)
    particles["Earth"].params["k2"] = K2_EARTH / 2.0
    particles["Earth"].params["I"] = C_EARTH * M_EARTH * R_EARTH**2
    spin_period_years = EARTH_SPIN_PERIOD_HOURS / 24.0 / 365.25
    spin = 2.0 * np.pi / spin_period_years
    obl = np.radians(EARTH_OBLIQUITY_DEG)
    particles["Earth"].params["Omega"] = [0.0, spin * np.sin(obl), spin * np.cos(obl)]
    particles["Earth"].params["tau"] = tau_seconds / YEAR_SECONDS
    angular_momentum = np.asarray(sim.angular_momentum()) + np.asarray(extras.spin_angular_momentum())
    rotation = rebound.Rotation.to_new_axes(newz=angular_momentum)
    extras.rotate_simulation(rotation)
    extras.initialize_spin_ode(force)
    return extras


def create_simulation(config: SimulationConfig, beta: float, eccentricity: float = 1e-7):
    rebound, reboundx = _require_rebound()
    model = get_moon_model(config.moon_type)
    axes = satellite_semimajor_axes(beta, config.n_moons, model, mode=config.spacing_mode)
    sim = rebound.Simulation()
    sim.units = ("AU", "yr", "Msun")
    sim.integrator = "ias15"
    sim.collision = "direct"
    sim.collision_resolve = "halt"
    sim.add(m=M_SUN, hash="Sun")
    sim.add(m=M_EARTH, a=A_EARTH, e=E_EARTH, r=R_EARTH,
            primary=sim.particles[0], hash="Earth")
    for i, axis in enumerate(axes):
        sim.add(m=model.mass_solar, a=axis, e=eccentricity, r=model.radius_au,
                f=golden_phase(i), primary=sim.particles["Earth"], hash=f"M{i+1}")
    sim.move_to_com()
    extras = None
    if config.tau_seconds > 0:
        extras = _setup_tides(sim, config.tau_seconds, rebound, reboundx)
    r_roche = roche_radius(model)
    p_ref = np.sqrt((1.8 * r_roche) ** 3 / M_EARTH)
    sim.dt = p_ref / 20.0
    sim.ri_ias15.min_dt = p_ref / 50.0
    return sim, extras, axes, p_ref


def run_simulation(config: SimulationConfig, beta: float, eccentricity: float = 1e-7) -> dict[str, float | str]:
    rebound, _ = _require_rebound()
    start = time.perf_counter()
    model = get_moon_model(config.moon_type)
    sim, extras, initial_axes, p_ref = create_simulation(config, beta, eccentricity)
    # Keep the REBOUNDx Extras object alive for the full integration.
    particles = sim.particles
    p_inner = np.sqrt(initial_axes[0] ** 3 / M_EARTH)
    tmax = config.max_orbits * p_inner
    r_roche = roche_radius(model)
    r_hill = hill_radius()
    emax = np.zeros(config.n_moons)
    emin = np.ones(config.n_moons)
    check_far = 1000.0 * p_ref
    check_near = 40.0 * p_ref
    check_step = check_far
    fate = "stable"

    try:
        while sim.t < tmax:
            sim.integrate(min(sim.t + check_step, tmax))
            near = False
            for i in range(config.n_moons):
                a, e, distance, _ = _relative_elements(
                    rebound, particles["Earth"], particles[2 + i], M_EARTH + model.mass_solar
                )
                emax[i] = max(emax[i], e)
                emin[i] = min(emin[i], e)
                if a < 0:
                    fate = "a<0"
                    break
                if e > 1:
                    fate = "e>1"
                    break
                if distance < r_roche:
                    fate = "d<R_roche"
                    break
                boundary_value = a if config.outer_boundary == "semimajor-axis" else a * (1.0 + e)
                if boundary_value > 0.4 * r_hill:
                    fate = "a>0.4R_H" if config.outer_boundary == "semimajor-axis" else "apo>0.4R_H"
                    break
                if distance < 3.0 * r_roche or boundary_value > 0.35 * r_hill:
                    near = True
            if fate != "stable":
                break
            check_step = check_near if near else check_far
        if fate == "stable" and sim.t >= tmax:
            fate = "tmax"
    except rebound.Collision:
        fate = "collision"

    final_axes = np.zeros(config.n_moons)
    for i in range(config.n_moons):
        try:
            final_axes[i], _, _, _ = _relative_elements(
                rebound, particles["Earth"], particles[2 + i], M_EARTH + model.mass_solar
            )
        except Exception:
            final_axes[i] = np.nan

    row: dict[str, float | str] = {"dH": beta, "e0": eccentricity}
    for i in range(config.n_moons):
        row[f"eM{i+1}_max"] = emax[i]
        row[f"eM{i+1}_min"] = emin[i]
        row[f"dA{i+1}_over_A0"] = (final_axes[i] - initial_axes[i]) / initial_axes[i]
    row["t_end(PM)"] = sim.t / p_inner
    row["wall_time"] = time.perf_counter() - start
    row["fate"] = fate
    return row


def run_simulation_worker(payload: tuple[dict[str, Any], float, float]):
    config_dict, beta, eccentricity = payload
    return run_simulation(SimulationConfig(**config_dict), beta, eccentricity)


def config_payload(config: SimulationConfig) -> dict[str, Any]:
    return asdict(config)


def run_time_series(config: SimulationConfig, beta: float, years: float,
                    output_every_orbits: float = 5.0) -> dict[str, np.ndarray]:
    rebound, _ = _require_rebound()
    model = get_moon_model(config.moon_type)
    sim, extras, axes, _ = create_simulation(config, beta, 1e-7)
    # Keep the REBOUNDx Extras object alive for the full integration.
    particles = sim.particles
    p_inner = np.sqrt(axes[0] ** 3 / M_EARTH)
    times = np.arange(0.0, years + 0.5 * output_every_orbits * p_inner,
                      output_every_orbits * p_inner)
    positions = np.zeros((len(times), config.n_moons, 2))
    semimajor = np.zeros((len(times), config.n_moons))
    eccentricity = np.zeros((len(times), config.n_moons))
    periods = np.zeros((len(times), config.n_moons))
    r_hill = hill_radius()
    for j, target in enumerate(times):
        sim.integrate(target)
        for i in range(config.n_moons):
            moon = particles[2 + i]
            earth = particles["Earth"]
            positions[j, i] = [(moon.x - earth.x) / r_hill, (moon.y - earth.y) / r_hill]
            a, e, _, p = _relative_elements(rebound, earth, moon, M_EARTH + model.mass_solar)
            semimajor[j, i] = a / r_hill
            eccentricity[j, i] = e
            periods[j, i] = p
    return {
        "time_years": times,
        "positions_rh": positions,
        "semimajor_rh": semimajor,
        "eccentricity": eccentricity,
        "periods_years": periods,
        "beta": np.asarray(beta),
        "tau_seconds": np.asarray(config.tau_seconds),
        "moon_type": np.asarray(model.name),
    }
