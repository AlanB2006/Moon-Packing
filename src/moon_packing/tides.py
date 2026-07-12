from __future__ import annotations

from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from .constants import (
    A_EARTH, C_EARTH, EARTH_OBLIQUITY_DEG, EARTH_SPIN_PERIOD_HOURS,
    E_EARTH, K2_EARTH, M_EARTH, M_SUN, R_EARTH, YEAR_SECONDS,
    get_moon_model,
)
from .dynamics import hill_radius, roche_radius
from .simulation import _require_rebound

G_AU = 4.0 * np.pi**2


def _f(index: int, e: float) -> float:
    if index == 1:
        return 1 + 31/2*e**2 + 255/8*e**4 + 185/16*e**6 + 25/64*e**8
    if index == 2:
        return 1 + 15/2*e**2 + 45/8*e**4 + 5/16*e**6
    if index == 3:
        return 1 + 15/4*e**2 + 15/8*e**4 + 5/64*e**6
    if index == 4:
        return 1 + 3/2*e**2 + 1/8*e**4
    if index == 5:
        return 1 + 3*e**2 + 3/8*e**4
    raise ValueError(index)


def secular_derivatives_archived(_t, state, satellite_mass, tau_years):
    """Barnes-style CTL equations used by the supplied secular script.

    The state is moon a/e, planet a/e, and planetary spin. Keeping this
    implementation separate and explicit makes the Figure 1 provenance clear.
    """
    a_m, e_m, a_p, e_p, omega = state
    psi = np.radians(EARTH_OBLIQUITY_DEG)
    b_m = np.sqrt(max(1.0 - e_m**2, 1e-15))
    b_p = np.sqrt(max(1.0 - e_p**2, 1e-15))
    n_m = np.sqrt(G_AU * (M_EARTH + satellite_mass) / a_m**3)
    n_p = np.sqrt(G_AU * (M_SUN + M_EARTH) / a_p**3)

    z_pm = (3*G_AU**2*K2_EARTH*satellite_mass**2*(M_EARTH+satellite_mass)
            * R_EARTH**5/a_m**9*tau_years)
    z_ps = (3*G_AU**2*K2_EARTH*M_SUN**2*(M_EARTH+M_SUN)
            * R_EARTH**5/a_p**9*tau_years)

    da_m = (2*a_m**2/(G_AU*M_EARTH*satellite_mass) * z_pm
            * (np.cos(psi)*_f(2, e_m)/b_m**12*omega/n_m - _f(1, e_m)/b_m**15))
    de_m = (11*a_m*e_m/(2*G_AU*M_EARTH*satellite_mass) * z_pm
            * (np.cos(psi)*_f(4, e_m)/b_m**10*omega/n_m - _f(3, e_m)/b_m**13))
    da_p = (2*a_p**2/(G_AU*M_SUN*M_EARTH) * z_ps
            * (np.cos(psi)*_f(2, e_p)/b_p**12*omega/n_p - _f(1, e_p)/b_p**15))
    de_p = (11*a_p*e_p/(2*G_AU*M_SUN*M_EARTH) * z_ps
            * (np.cos(psi)*_f(4, e_p)/b_p**10*omega/n_p - _f(3, e_p)/b_p**13))

    inertia = C_EARTH*M_EARTH*R_EARTH**2
    d_omega = z_pm/(2*inertia*n_m) * (
        2*np.cos(psi)*_f(2, e_m)/b_m**12
        - (1+np.cos(psi)**2)*_f(5, e_m)/b_m**9*omega/n_m
    )
    d_omega += z_ps/(2*inertia*n_p) * (
        2*np.cos(psi)*_f(2, e_p)/b_p**12
        - (1+np.cos(psi)**2)*_f(5, e_p)/b_p**9*omega/n_p
    )
    return [da_m, de_m, da_p, de_p, d_omega]


def generate_secular_validation(years: float = 30000, tau_seconds: float = 698,
                                n_output: int = 301) -> np.ndarray:
    model = get_moon_model("Luna")
    a0 = 2.0 * roche_radius(model)
    omega0 = 2*np.pi/(EARTH_SPIN_PERIOD_HOURS/24/365.25)
    times = np.linspace(0.0, years, n_output)
    solution = solve_ivp(
        secular_derivatives_archived,
        (0.0, years),
        [a0, 1e-4, A_EARTH, E_EARTH, omega0],
        t_eval=times,
        args=(model.mass_solar, tau_seconds/YEAR_SECONDS),
        rtol=1e-11,
        atol=1e-13,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return np.column_stack([solution.t, solution.y[0]/hill_radius(), solution.y[1], solution.y[4]])


def generate_reboundx_validation(years: float = 30000, tau_seconds: float = 698,
                                 n_output: int = 301) -> np.ndarray:
    rebound, reboundx = _require_rebound()
    model = get_moon_model("Luna")
    sim = rebound.Simulation()
    sim.units = ("AU", "yr", "Msun")
    sim.integrator = "ias15"
    sim.add(m=M_SUN, hash="Sun")
    sim.add(m=M_EARTH, a=A_EARTH, e=E_EARTH, r=R_EARTH,
            primary=sim.particles[0], hash="Earth")
    sim.add(m=model.mass_solar, a=2*roche_radius(model), e=1e-4,
            r=model.radius_au, M=np.pi/3,
            primary=sim.particles["Earth"], hash="Moon")
    sim.move_to_com()
    extras = reboundx.Extras(sim)
    force = extras.load_force("tides_spin")
    extras.add_force(force)
    earth = sim.particles["Earth"]
    earth.params["k2"] = K2_EARTH/2
    earth.params["I"] = C_EARTH*M_EARTH*R_EARTH**2
    omega0 = 2*np.pi/(EARTH_SPIN_PERIOD_HOURS/24/365.25)
    obl = np.radians(EARTH_OBLIQUITY_DEG)
    earth.params["Omega"] = [0, omega0*np.sin(obl), omega0*np.cos(obl)]
    earth.params["tau"] = tau_seconds/YEAR_SECONDS
    total_l = np.asarray(sim.angular_momentum()) + np.asarray(extras.spin_angular_momentum())
    extras.rotate_simulation(rebound.Rotation.to_new_axes(newz=total_l))
    extras.initialize_spin_ode(force)

    times = np.linspace(0.0, years, n_output)
    output = np.zeros((n_output, 4))
    for i, target in enumerate(times):
        sim.integrate(target)
        orbit = sim.particles["Moon"].orbit(primary=sim.particles["Earth"])
        omega = np.linalg.norm(np.asarray(sim.particles["Earth"].params["Omega"]))
        output[i] = [target, orbit.a/hill_radius(), orbit.e, omega]
    return output


def save_validation(output_dir: str | Path, years: float = 30000, tau_seconds: float = 698):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    secular = generate_secular_validation(years, tau_seconds)
    np.savetxt(output_dir/"secular_tides.csv", secular, delimiter=",",
               header="time_years,a_over_RH,eccentricity,omega_per_year", comments="")
    reboundx_data = generate_reboundx_validation(years, tau_seconds)
    np.savetxt(output_dir/"reboundx_tides.csv", reboundx_data, delimiter=",",
               header="time_years,a_over_RH,eccentricity,omega_per_year", comments="")
