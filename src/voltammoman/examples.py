"""A set of synthetic example voltammograms (NOVA export format) for trying the tool out."""

from __future__ import annotations

from pathlib import Path

from .io import write_curve
from .simulate import potential_program, simulate_cv

# Two marker-like systems: oxidation first (ferrocyanide-like) and reduction first
# (hexaammineruthenium-like), each from fast to slow electron transfer, with a charging
# current, a sloped background and a little noise, as on real carbon electrodes.
_BACKGROUND = dict(Cdl=3e-4, slope=4e-6, noise=3e-8)
EXAMPLES = {
    **{f"sim-ox-k0_{k}": dict(E0=0.22, k0=float(k), c_red=1e-6, c_ox=0.0,
                              E=potential_program(0, 0.9, -0.3), **_BACKGROUND)
       for k in ("1e-1", "1e-2", "3e-3", "1e-3", "3e-4")},
    **{f"sim-red-k0_{k}": dict(E0=-0.17, k0=float(k), c_red=0.0, c_ox=1e-6,
                               E=potential_program(0, 0.3, -0.5), **_BACKGROUND)
       for k in ("1e-1", "1e-2", "1e-3")},
    "sim-no-peak": dict(E0=0.22, k0=1e-2, c_red=0.0, c_ox=0.0, E=potential_program(0, 0.9, -0.3), **_BACKGROUND),
}


def write_examples(directory) -> list[Path]:
    """Write all example curves into ``directory`` (created if needed); returns their paths."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for seed, (name, kw) in enumerate(EXAMPLES.items()):
        E, I = simulate_cv(seed=seed, **kw)
        paths.append(write_curve(directory / f"{name}.txt", E, I))
    return paths
