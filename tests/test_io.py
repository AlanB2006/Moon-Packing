from pathlib import Path

from moon_packing.io import read_results

ROOT = Path(__file__).resolve().parents[1]


def test_read_luna_sweep():
    table = read_results(ROOT / "data/sweeps/Luna/2_Luna_0s_final.txt")
    assert len(table) > 100
    assert table.n_moons == 2
    assert table["dH"].shape == table.fate.shape
    assert "tmax" in set(table.fate)


def test_read_ceres_map():
    table = read_results(ROOT / "data/eccentricity_maps/Ceres/output_C_5moons_T1_14.0_698.txt")
    assert table.n_moons == 5
    assert len(table) > 30000
