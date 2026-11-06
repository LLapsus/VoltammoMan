"""Peak evaluation on simulated curves with known answers."""

import numpy as np
import pytest

from conftest import BACKGROUND, RED_FIRST
from voltammoman.analysis import analyze_curve
from voltammoman.cycle import select_cycle, split_branches
from voltammoman.simulate import potential_program


@pytest.mark.parametrize("k0", [1.0, 1e-2, 1e-3])
@pytest.mark.parametrize("system", ["ox_first", "red_first"])
def test_background_does_not_change_the_peaks(sim, k0, system):
    kw = {} if system == "ox_first" else RED_FIRST
    clean = analyze_curve(*sim(k0=k0, **kw)).values
    noisy = analyze_curve(*sim(k0=k0, **kw, **BACKGROUND, noise=5e-8, seed=1)).values
    assert noisy["flags"] == ""
    assert abs(noisy["EpA"] - clean["EpA"]) < 8e-3
    assert abs(noisy["EpC"] - clean["EpC"]) < 8e-3
    assert noisy["IpA"] == pytest.approx(clean["IpA"], rel=0.06)
    assert noisy["IpC"] == pytest.approx(clean["IpC"], rel=0.06)


@pytest.mark.parametrize("k0", [1.0, 1e-2, 1e-3])
def test_nicholson_ratio_is_one_for_a_stable_couple(sim, k0):
    v = analyze_curve(*sim(k0=k0, **BACKGROUND)).values
    assert v["ratio_nicholson"] == pytest.approx(1.0, abs=0.12)


def test_first_peak(sim):
    assert analyze_curve(*sim(k0=1e-2)).values["first_peak"] == "anodic"
    assert analyze_curve(*sim(k0=1e-2, **RED_FIRST)).values["first_peak"] == "cathodic"


@pytest.mark.parametrize("method", ["flattest", "linear", "tangent"])
def test_all_methods_find_reversible_peaks(reversible, method):
    v = analyze_curve(*reversible, baseline=method).values
    assert 50e-3 < v["dEp"] < 65e-3
    assert v["IpA"] > 0 > v["IpC"]


def test_no_peak_is_flagged(sim):
    v = analyze_curve(*sim(k0=1e-2, c_red=0.0, c_ox=0.0, **BACKGROUND, noise=3e-8, seed=2)).values
    flags = v["flags"].split(";")
    assert "A:weak_peak" in flags
    assert any(f.startswith("C:") for f in flags)   # weak_peak or no_baseline, depending on the noise
    assert np.isnan(v["ratio_nicholson"])


def test_reversed_start_direction(sim, reversible):
    # the same system scanned negative first: 0 V -> -0.3 V -> 0.9 V -> 0 V
    v1 = analyze_curve(*reversible).values
    v2 = analyze_curve(*sim(k0=1.0, E=potential_program(0, -0.3, 0.9))).values
    assert v2["EpA"] == pytest.approx(v1["EpA"], abs=5e-3)
    assert v2["EpC"] == pytest.approx(v1["EpC"], abs=5e-3)
    assert v2["IpA"] == pytest.approx(v1["IpA"], rel=0.05)


def test_cycle_starting_at_vertex(reversible):
    E, I = reversible
    k = int(np.argmin(E))                     # rotate the cycle to start at the lower vertex
    E2, I2 = np.r_[E[k:], E[:k]], np.r_[I[k:], I[:k]]
    an, ca = split_branches(E2, I2)
    step = 10 / 4096
    assert an.E[0] == E.min() and abs(ca.E[-1] - E.min()) <= step
    v = analyze_curve(E2, I2).values
    assert v["dEp"] == pytest.approx(analyze_curve(E, I).values["dEp"], abs=5e-3)


def test_multi_cycle_selection(sim):
    E, I1 = sim(k0=1e-2)
    I2 = I1 * 2
    E3, I3 = np.r_[E, E, E], np.r_[I1, I1, I2]
    Es, Is = select_cycle(E3, I3, cycle=-1)
    np.testing.assert_array_equal(Is, I2)
    Es, Is = select_cycle(E3, I3, cycle=0)
    np.testing.assert_array_equal(Is, I1)
    scan = np.repeat([1, 2, 3], len(E))
    Es, Is = select_cycle(E3, I3, scan=scan, cycle=1)
    np.testing.assert_array_equal(Is, I1)
