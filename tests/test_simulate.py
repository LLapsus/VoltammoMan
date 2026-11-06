"""The simulator against textbook results (Randles-Sevcik, Nicholson's psi - dEp table)."""

import numpy as np
import pytest

from voltammoman.analysis import analyze_curve
from voltammoman.simulate import F, R, potential_program, randles_sevcik


def test_potential_program_like_nova():
    E = potential_program(0, 0.9, -0.3)
    assert len(E) == 984
    assert E[0] == pytest.approx(10 / 4096)
    assert E.max() == pytest.approx(0.90087890625)
    assert E.min() == pytest.approx(-0.30029296875)
    assert E[-1] == 0
    assert len(potential_program(0, 0.3, -0.5)) == 656


def test_reversible_wave(reversible):
    E, I = reversible
    v = analyze_curve(E, I).values
    assert 54e-3 < v["dEp"] < 62e-3                         # ~57-59 mV at 25 °C, 2.44 mV steps
    assert v["IpA"] == pytest.approx(randles_sevcik(), rel=0.04)
    assert v["E_half"] == pytest.approx(0.25, abs=0.005)


@pytest.mark.parametrize("psi,dEp_mV", [(1.0, 84), (0.25, 141), (0.1, 212)])
def test_quasi_reversible_matches_nicholson(sim, psi, dEp_mV):
    # Nicholson (1965): psi = k0 / sqrt(pi D n v F / R T)
    k0 = psi * np.sqrt(np.pi * 7e-6 * 0.1 * F / (R * 298.15))
    E, I = sim(k0=k0)
    v = analyze_curve(E, I).values
    assert v["dEp"] * 1e3 == pytest.approx(dEp_mV, abs=8)


def test_noise_is_reproducible(sim):
    a = sim(k0=1e-2, noise=1e-7, seed=3)[1]
    b = sim(k0=1e-2, noise=1e-7, seed=3)[1]
    np.testing.assert_array_equal(a, b)


def test_vertices_must_straddle_start():
    with pytest.raises(ValueError):
        potential_program(0, 0.5, 0.2)
