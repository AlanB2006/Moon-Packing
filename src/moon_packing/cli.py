from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys

import numpy as np

from .constants import get_moon_model
from .io import write_results
from .plotting import reproduce_figures
from .simulation import (
    SimulationConfig, config_payload, run_simulation_worker, run_time_series,
)
from .tides import save_validation


def _float_range(start: float, stop: float, step: float) -> np.ndarray:
    if step <= 0:
        raise ValueError("step must be positive")
    n = int(np.floor((stop - start) / step + 0.5)) + 1
    return np.round(start + np.arange(max(n, 0)) * step, 12)


def _common_simulation_arguments(parser: argparse.ArgumentParser):
    parser.add_argument("--moon-type", required=True, choices=["Luna", "Pluto", "Ceres", "luna", "pluto", "ceres"])
    parser.add_argument("--n-moons", type=int, required=True)
    parser.add_argument("--tau", type=float, default=0.0, help="Tidal time lag in seconds")
    parser.add_argument("--spacing-mode", choices=["archived", "recursive"], default="archived")
    parser.add_argument("--outer-boundary", choices=["semimajor-axis", "apocenter"], default="semimajor-axis")


def run_grid_main(argv=None):
    parser = argparse.ArgumentParser(description="Run a beta or beta/eccentricity grid of packed-moon simulations.")
    _common_simulation_arguments(parser)
    parser.add_argument("--beta-start", type=float, required=True)
    parser.add_argument("--beta-stop", type=float, required=True)
    parser.add_argument("--beta-step", type=float, required=True)
    parser.add_argument("--e-start", type=float, default=1e-7)
    parser.add_argument("--e-stop", type=float)
    parser.add_argument("--e-step", type=float)
    parser.add_argument("--max-orbits", type=float, default=1e7)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    get_moon_model(args.moon_type)
    betas = _float_range(args.beta_start, args.beta_stop, args.beta_step)
    if args.e_stop is None:
        eccentricities = np.asarray([args.e_start])
    else:
        if args.e_step is None:
            parser.error("--e-step is required when --e-stop is given")
        eccentricities = _float_range(args.e_start, args.e_stop, args.e_step)
    config = SimulationConfig(args.moon_type, args.n_moons, args.tau, args.max_orbits,
                              args.spacing_mode, args.outer_boundary)
    payloads = [(config_payload(config), float(beta), float(e))
                for e in eccentricities for beta in betas]
    print(f"Running {len(payloads)} integrations with {args.workers} worker(s).")
    rows = []
    if args.workers == 1:
        for i, payload in enumerate(payloads, 1):
            rows.append(run_simulation_worker(payload))
            print(f"\r{i}/{len(payloads)}", end="", flush=True)
    else:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_simulation_worker, payload) for payload in payloads]
            for i, future in enumerate(as_completed(futures), 1):
                rows.append(future.result())
                print(f"\r{i}/{len(payloads)}", end="", flush=True)
    print()
    write_results(args.output, rows, args.n_moons)
    print(args.output)


def run_timeseries_main(argv=None):
    parser = argparse.ArgumentParser(description="Run and save one packed-moon time series.")
    _common_simulation_arguments(parser)
    parser.add_argument("--beta", type=float, required=True)
    parser.add_argument("--years", type=float, default=25.0)
    parser.add_argument("--output-every-orbits", type=float, default=5.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    config = SimulationConfig(args.moon_type, args.n_moons, args.tau, spacing_mode=args.spacing_mode,
                              outer_boundary=args.outer_boundary)
    data = run_time_series(config, args.beta, args.years, args.output_every_orbits)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output, **data)
    print(args.output)


def figures_main(argv=None):
    parser = argparse.ArgumentParser(description="Reproduce manuscript figures.")
    parser.add_argument("--figures", nargs="+", type=int, default=[2, 3, 4, 8, 9, 10])
    parser.add_argument("--root", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    outputs = reproduce_figures(args.figures, args.root, args.output_dir)
    for output in outputs:
        print(output)
    return 0 if outputs else 1


def validate_tides_main(argv=None):
    parser = argparse.ArgumentParser(description="Generate REBOUNDx and secular tidal-validation tables.")
    parser.add_argument("--years", type=float, default=30000.0)
    parser.add_argument("--tau", type=float, default=698.0)
    parser.add_argument("--output-dir", type=Path, default=Path("data/generated/tides"))
    args = parser.parse_args(argv)
    save_validation(args.output_dir, args.years, args.tau)
    print(args.output_dir)
