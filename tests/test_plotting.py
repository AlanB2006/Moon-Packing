from pathlib import Path

import matplotlib.pyplot as plt

from moon_packing.constants import get_moon_model
from moon_packing.plotting import _add_resonance_top_axis, _resonance_tick_pairs, plot_eccentricity_map, plot_lifetime_eccentricity

ROOT = Path(__file__).resolve().parents[1]


def test_lifetime_plot(tmp_path):
    out = plot_lifetime_eccentricity("Luna", root=ROOT, output=tmp_path / "life.pdf")
    assert out.exists() and out.stat().st_size > 1000


def test_heatmap_plot(tmp_path):
    out = plot_eccentricity_map("Luna", root=ROOT, output=tmp_path / "map.pdf")
    assert out.exists() and out.stat().st_size > 1000


def test_resonance_tick_pairs_are_sorted_and_visible():
    model = get_moon_model("Luna")
    tick_pairs = _resonance_tick_pairs(model, model.manuscript_count, 3.5, 9.0)
    locations = [location for location, _ in tick_pairs]
    assert locations == sorted(locations)
    assert all(3.5 <= location <= 9.0 for location in locations)
    assert r"$\beta_{\max}$" in [label for _, label in tick_pairs]


def test_top_axis_contains_resonance_labels():
    model = get_moon_model("Luna")
    fig, ax = plt.subplots()
    ax.set_xlim(3.5, 9.0)
    top_ax = _add_resonance_top_axis(ax, model, model.manuscript_count)
    labels = [label.get_text() for label in top_ax.get_xticklabels()]
    plt.close(fig)
    assert r"$7\!:\!2$" in labels
    assert r"$9\!:\!1$" in labels
    assert r"$\beta_{\max}$" in labels
