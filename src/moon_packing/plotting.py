from __future__ import annotations

from pathlib import Path
import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from .constants import get_moon_model
from .dynamics import beta_max, resonance_beta
from .io import read_results

TAUS = [0, 100, 698]
MOON_COLORS = ["red", "blue", "saddlebrown", "purple", "orange"]
PAPER_CONFIG = {"Luna": {"figure": 2, "sweep": "2_Luna_{tau}s_final.txt", "map": "output_L_2moons_{flag}_9.0_{tau}.txt", "xmax": 9.0}, "Pluto": {"figure": 3, "sweep": "3_Pluto_{tau}s_final.txt", "map": "output_P_3moons_{flag}_11.0_{tau}.txt", "xmax": 11.0}, "Ceres": {"figure": 4, "sweep": "5_Ceres_{tau}s_final.txt", "map": "output_C_5moons_{flag}_14.0_{tau}.txt", "xmax": 14.0}}
RESONANCE_TICKS = ((7 / 2, r"$7\!:\!2$"), (4, r"$4\!:\!1$"), (9 / 2, r"$9\!:\!2$"), (5, r"$5\!:\!1$"), (6, r"$6\!:\!1$"), (7, r"$7\!:\!1$"), (8, r"$8\!:\!1$"), (9, r"$9\!:\!1$"))


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _sweep_path(root: Path, moon_name: str, tau: int) -> Path:
    template = PAPER_CONFIG[moon_name]["sweep"]
    return root / "data" / "sweeps" / moon_name / template.format(tau=tau)


def _map_path(root: Path, moon_name: str, tau: int) -> Path:
    template = PAPER_CONFIG[moon_name]["map"]
    flag = "T0" if tau == 0 else "T1"
    return root / "data" / "eccentricity_maps" / moon_name / template.format(flag=flag, tau=tau)


def _save(fig, output: Path) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    return output


def _resonance_tick_pairs(model, n_moons: int, xmin: float, xmax: float) -> list[tuple[float, str]]:
    tick_pairs = [(resonance_beta(period_ratio, model), label) for period_ratio, label in RESONANCE_TICKS]
    tick_pairs.append((beta_max(n_moons, model), r"$\beta_{\max}$"))
    return sorted((location, label) for location, label in tick_pairs if xmin <= location <= xmax)


def _add_resonance_top_axis(ax, model, n_moons: int):
    xmin, xmax = ax.get_xlim()
    tick_pairs = _resonance_tick_pairs(model, n_moons, xmin, xmax)
    top_ax = ax.secondary_xaxis("top")
    if tick_pairs:
        locations, labels = zip(*tick_pairs)
        top_ax.set_xticks(locations, labels)
    else:
        top_ax.set_xticks([])
    top_ax.tick_params(axis="x", which="major", direction="out", length=6, width=1, labelsize=8, pad=3)
    plt.setp(top_ax.get_xticklabels(), rotation=45, ha="left", rotation_mode="anchor")
    return top_ax


def plot_lifetime_eccentricity(moon_type: str, root: str | Path | None = None, output: str | Path | None = None) -> Path:
    root = Path(root) if root else _repo_root()
    model = get_moon_model(moon_type)
    cfg = PAPER_CONFIG[model.name]
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 6.2), sharex="col")
    panel = list("abcdef")

    for col, tau in enumerate(TAUS):
        table = read_results(_sweep_path(root, model.name, tau))
        beta = table["dH"]
        lifetime = table["t_end(PM)"]
        axes[0, col].scatter(beta, lifetime, s=2, c="black", rasterized=True)
        axes[0, col].set_yscale("log")
        axes[0, col].set_ylim(1e3, 1.6e7)
        axes[0, col].set_title(rf"$\tau={tau}\,$s")

        for i in range(1, table.n_moons + 1):
            axes[1, col].scatter(beta, table[f"eM{i}_max"], s=2, color=MOON_COLORS[i - 1], rasterized=True)

        axes[1, col].set_yscale("log")
        axes[1, col].set_ylim(1e-3, 1.5)
        axes[1, col].set_xlabel(r"$\beta\ (R_{H,m})$")

        for row in range(2):
            ax = axes[row, col]
            ax.set_xlim(2 * np.sqrt(3), cfg["xmax"])
            ax.minorticks_on()
            ax.tick_params(which="both", direction="out")
            ax.text(0.04, 0.88, panel[row * 3 + col], transform=ax.transAxes, fontsize=13, fontweight="bold")
            ax.axvline(beta_max(model.manuscript_count, model), ls="--", lw=1, color="0.25")

        if col == 0:
            axes[0, col].set_ylabel(r"Lifetime ($P_1$)")
            axes[1, col].set_ylabel(r"$e_{\max}$")
        else:
            axes[0, col].tick_params(labelleft=False)
            axes[1, col].tick_params(labelleft=False)

    handles = [plt.Line2D([], [], marker="s", ls="", color=MOON_COLORS[i], label=str(i + 1)) for i in range(model.manuscript_count)]
    fig.legend(handles=handles, title="Moon index", ncol=model.manuscript_count, loc="upper left", bbox_to_anchor=(0.08, 1.01), frameon=False)
    fig.subplots_adjust(wspace=0.08, hspace=0.12, top=0.88)
    output = Path(output) if output else root / "figures" / "reproduced" / f"figure{cfg['figure']:02d}_{model.name.lower()}_lifetime.pdf"
    return _save(fig, output)


def _heat_grid(table, moon_index: int, xmax: float):
    beta = table["dH"]
    e0 = table["e0"]
    values = table.eccentricity_range(moon_index).copy()
    unstable = np.isin(table.fate, ["a<0", "a>0.4R_H", "apo>0.4R_H", "e>1", "collision", "d<R_roche"])
    values[unstable] = 1.5
    xi = np.arange(3.5, xmax + 0.005, 0.01)
    yi = np.arange(0.0, 0.3001, 0.01)
    xg, yg = np.meshgrid(xi, yi)
    zi = griddata((beta, e0), values, (xg, yg), method="nearest")
    zi = np.clip(zi, 1e-2, None)
    return xi, yi, np.log10(zi)


def plot_eccentricity_map(moon_type: str, root: str | Path | None = None, output: str | Path | None = None) -> Path:
    root = Path(root) if root else _repo_root()
    model = get_moon_model(moon_type)
    cfg = PAPER_CONFIG[model.name]
    fig, axes = plt.subplots(3, 2, figsize=(10.5, 11.5), sharex=True, sharey=True)
    image = None
    panel = list("abcdef")
    beta_limit = beta_max(model.manuscript_count, model)

    for row, tau in enumerate(TAUS):
        table = read_results(_map_path(root, model.name, tau))
        moon_indices = [1, table.n_moons]

        for col, moon_index in enumerate(moon_indices):
            xi, yi, log_delta_e = _heat_grid(table, moon_index, cfg["xmax"])
            ax = axes[row, col]
            image = ax.pcolormesh(xi, yi, log_delta_e, cmap="gnuplot", vmin=-2, vmax=0, shading="auto", rasterized=True)
            ax.text(0.03, 0.89, panel[row * 2 + col], transform=ax.transAxes, fontsize=13, fontweight="bold")
            ax.axvline(beta_limit, color="black", ls="--", lw=1)

            for period_ratio, _ in RESONANCE_TICKS:
                location = resonance_beta(period_ratio, model)
                if 3.5 < location < cfg["xmax"]:
                    ax.axvline(location, color="black", ls=":", lw=0.55, alpha=0.35)

            ax.set_xlim(3.5, cfg["xmax"])
            ax.set_ylim(0.0, 0.3)
            ax.minorticks_on()
            ax.tick_params(which="both", direction="out")

            if row == 0:
                _add_resonance_top_axis(ax, model, model.manuscript_count)
                ax.set_title(f"Moon {moon_index}", pad=52)

            if row == 2:
                ax.set_xlabel(r"$\beta\ (R_{H,m})$")

            if col == 0:
                ax.set_ylabel(r"$e_0$")

        axes[row, 1].text(0.97, 0.92, rf"$\tau={tau}\,$s", transform=axes[row, 1].transAxes, ha="right", va="top", fontsize=10, bbox={"facecolor": "white", "alpha": 0.75, "edgecolor": "none", "pad": 1.5})

    if image is not None:
        cbar = fig.colorbar(image, ax=axes.ravel().tolist(), pad=0.02, fraction=0.035)
        cbar.set_label(r"$\log_{10}\Delta e$")

    fig.subplots_adjust(wspace=0.08, hspace=0.08, right=0.88, top=0.86)
    fig_no = {"Luna": 8, "Pluto": 9, "Ceres": 10}[model.name]
    output = Path(output) if output else root / "figures" / "reproduced" / f"figure{fig_no:02d}_{model.name.lower()}_eccentricity_map.pdf"
    return _save(fig, output)


def plot_tides_validation(root: str | Path | None = None, output: str | Path | None = None) -> Path:
    root = Path(root) if root else _repo_root()
    data_dir = root / "data" / "generated" / "tides"
    reboundx_path = data_dir / "reboundx_tides.csv"
    secular_path = data_dir / "secular_tides.csv"

    if not reboundx_path.exists() or not secular_path.exists():
        raise FileNotFoundError("Run moon-packing-validate-tides before reproducing Figure 1.")

    reboundx_data = np.genfromtxt(reboundx_path, delimiter=",", names=True)
    secular_data = np.genfromtxt(secular_path, delimiter=",", names=True)
    fig, axes = plt.subplots(3, 1, figsize=(4.7, 8.2), sharex=True)
    fields = [("a_over_RH", r"$a/R_H$"), ("eccentricity", r"$e_{sat}$"), ("omega_per_year", r"$\Omega$ (yr$^{-1}$)")]

    for ax, (field, label) in zip(axes, fields):
        ax.semilogx(reboundx_data["time_years"], reboundx_data[field], ".", ms=2, label="REBOUNDx")
        ax.semilogx(secular_data["time_years"], secular_data[field], ".", ms=2, label="Secular")
        ax.set_ylabel(label)
        ax.minorticks_on()
        ax.tick_params(which="both", direction="out")

    axes[0].legend(frameon=False)
    axes[-1].set_xlabel("Time (yr)")
    fig.subplots_adjust(hspace=0.12)
    output = Path(output) if output else root / "figures" / "reproduced" / "figure01_tides_validation.pdf"
    return _save(fig, output)


def plot_time_series(moon_type: str, beta_values: list[float], root: str | Path | None = None, output: str | Path | None = None) -> Path:
    root = Path(root) if root else _repo_root()
    model = get_moon_model(moon_type)
    paths = [root / "data" / "generated" / "timeseries" / f"{model.name.lower()}_beta{beta:g}_tau698.npz" for beta in beta_values]
    missing = [str(path) for path in paths if not path.exists()]

    if missing:
        raise FileNotFoundError("Generate the requested time series first:\n" + "\n".join(missing))

    fig, axes = plt.subplots(4, 3, figsize=(11.5, 12.5), sharex="col")

    for col, (path, beta) in enumerate(zip(paths, beta_values)):
        data = np.load(path)
        t = data["time_years"]
        pos = data["positions_rh"]
        semi = data["semimajor_rh"]
        ecc = data["eccentricity"]
        periods = data["periods_years"]

        for i in range(model.manuscript_count):
            axes[0, col].scatter(pos[:, i, 0], pos[:, i, 1], s=0.4, color=MOON_COLORS[i], rasterized=True)
            axes[1, col].scatter(t, semi[:, i], s=0.4, color=MOON_COLORS[i], rasterized=True)
            axes[2, col].scatter(t, ecc[:, i], s=0.4, color=MOON_COLORS[i], rasterized=True)
            if i > 0:
                axes[3, col].scatter(t, periods[:, i] / periods[:, i - 1], s=0.4, color=MOON_COLORS[i], rasterized=True)

        axes[0, col].set_aspect("equal", adjustable="box")
        axes[0, col].set_xlim(-0.44, 0.44)
        axes[0, col].set_ylim(-0.44, 0.44)
        axes[0, col].set_title(rf"$\beta={beta:g}$")
        axes[2, col].set_yscale("log")
        axes[3, col].set_xlabel("Time (yr)")

        for row in range(4):
            axes[row, col].minorticks_on()
            axes[row, col].tick_params(which="both", direction="out")

    labels = [r"$Y/R_H$", r"$a/R_H$", r"$e$", r"$P_{i+1}/P_i$"]

    for row, label in enumerate(labels):
        axes[row, 0].set_ylabel(label)

    for col in range(3):
        axes[0, col].set_xlabel(r"$X/R_H$")

    fig.subplots_adjust(wspace=0.1, hspace=0.12)
    fig_no = {"Luna": 5, "Pluto": 6, "Ceres": 7}[model.name]
    output = Path(output) if output else root / "figures" / "reproduced" / f"figure{fig_no:02d}_{model.name.lower()}_timeseries.pdf"
    return _save(fig, output)


def reproduce_figures(figures: list[int], root: str | Path | None = None, output_dir: str | Path | None = None) -> list[Path]:
    root = Path(root) if root else _repo_root()
    output_dir = Path(output_dir) if output_dir else root / "figures" / "reproduced"
    outputs: list[Path] = []

    for figure in figures:
        try:
            if figure == 1:
                outputs.append(plot_tides_validation(root, output_dir / "figure01_tides_validation.pdf"))
            elif figure in {2, 3, 4}:
                name = {2: "Luna", 3: "Pluto", 4: "Ceres"}[figure]
                outputs.append(plot_lifetime_eccentricity(name, root, output_dir / f"figure{figure:02d}_{name.lower()}_lifetime.pdf"))
            elif figure in {5, 6, 7}:
                name = {5: "Luna", 6: "Pluto", 7: "Ceres"}[figure]
                betas = {"Luna": [4.4, 8.2, 8.55], "Pluto": [4.3, 9.0, 9.7], "Ceres": [6.0, 13.2, 13.55]}[name]
                outputs.append(plot_time_series(name, betas, root, output_dir / f"figure{figure:02d}_{name.lower()}_timeseries.pdf"))
            elif figure in {8, 9, 10}:
                name = {8: "Luna", 9: "Pluto", 10: "Ceres"}[figure]
                outputs.append(plot_eccentricity_map(name, root, output_dir / f"figure{figure:02d}_{name.lower()}_eccentricity_map.pdf"))
            else:
                warnings.warn(f"No recipe is defined for figure {figure}")
        except FileNotFoundError as exc:
            warnings.warn(str(exc))

    return outputs
