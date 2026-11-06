import pytest

from voltammoman.simulate import potential_program, simulate_cv

BACKGROUND = dict(Cdl=3e-4, slope=5e-6, offset=2e-6)
RED_FIRST = dict(E=potential_program(0, 0.3, -0.5), E0=-0.17, c_ox=1e-6, c_red=0.0)


@pytest.fixture(scope="session")
def reversible():
    """Fast electron transfer, oxidation first, no background."""
    return simulate_cv(k0=1.0)


@pytest.fixture(scope="session")
def sim():
    """Cached simulations keyed by their parameters."""
    cache = {}

    def get(**kw):
        key = tuple(sorted((k, v if not hasattr(v, "tobytes") else v.tobytes()) for k, v in kw.items()))
        if key not in cache:
            cache[key] = simulate_cv(**kw)
        return cache[key]

    return get
