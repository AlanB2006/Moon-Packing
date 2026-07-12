from pathlib import Path

from moon_packing.plotting import plot_eccentricity_map, plot_lifetime_eccentricity

ROOT = Path(__file__).resolve().parents[1]


def test_lifetime_plot(tmp_path):
    out = plot_lifetime_eccentricity("Luna", root=ROOT, output=tmp_path / "life.pdf")
    assert out.exists() and out.stat().st_size > 1000


def test_heatmap_plot(tmp_path):
    out = plot_eccentricity_map("Luna", root=ROOT, output=tmp_path / "map.pdf")
    assert out.exists() and out.stat().st_size > 1000
