import numpy as np
import pytest

from moon_packing.constants import get_moon_model
from moon_packing.dynamics import beta_max, hill_radius, resonance_beta, roche_radius, satellite_semimajor_axes


def test_basic_scales():
    luna = get_moon_model("L")
    assert 0 < roche_radius(luna) < hill_radius()
    assert beta_max(2, luna) > 2*np.sqrt(3)


def test_archived_axes_increase():
    model = get_moon_model("Pluto")
    axes = satellite_semimajor_axes(7.0, 3, model, mode="archived")
    assert np.all(np.diff(axes) > 0)


def test_recursive_is_available_and_distinct():
    model = get_moon_model("Ceres")
    archived = satellite_semimajor_axes(8.0, 5, model, mode="archived")
    recursive = satellite_semimajor_axes(8.0, 5, model, mode="recursive")
    assert archived.shape == recursive.shape
    assert not np.allclose(archived[2:], recursive[2:])


def test_resonance_beta():
    model = get_moon_model("Luna")
    assert resonance_beta(2.0, model) > 0
    with pytest.raises(ValueError):
        resonance_beta(1.0, model)
