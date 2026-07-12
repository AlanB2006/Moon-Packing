from moon_packing.cli import _float_range


def test_float_range_inclusive():
    values = _float_range(0.0, 0.3, 0.1)
    assert list(values) == [0.0, 0.1, 0.2, 0.3]
