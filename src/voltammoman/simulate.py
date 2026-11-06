"""A small simulator of cyclic voltammograms, used for tests and examples.

Model: one-electron (or n-electron) couple O + n e- <-> R at a planar electrode with
semi-infinite linear diffusion (equal diffusion coefficients), Butler-Volmer kinetics, a
double-layer capacitance, optional uncompensated resistance, a linear background, a constant
current offset and Gaussian noise. Diffusion is solved by an implicit (backward Euler)
finite-difference scheme.

The default potential program reproduces NOVA exports: steps of 2.44140625 mV, start at 0 V,
first vertex, second vertex and back to 0 V; the first exported point is one step after the
start. Several cycles are simulated and only the last one is returned, as an exported
steady-state cycle.

Units: potentials V, currents A, concentrations mol/cm^3 (1 mM = 1e-6 mol/cm^3), lengths cm,
time s, capacitance F/cm^2, resistance ohm.
"""

from __future__ import annotations

import numpy as np
from scipy.linalg import solve_banded

F = 96485.33212   # C/mol
R = 8.314462618   # J/(mol K)
NOVA_STEP = 10 / 4096  # 2.44140625 mV


def potential_program(E_start=0.0, E_vertex1=0.9, E_vertex2=-0.3, step=NOVA_STEP):
    """Staircase potential program E_start -> E_vertex1 -> E_vertex2 -> E_start.

    Vertices are rounded outwards to whole steps from E_start. As in NOVA exports, the first
    point is one step after E_start and the last point is E_start.
    """
    def to_steps(x):
        return int(np.ceil(abs(x - E_start) / step - 1e-9)) * int(np.sign(x - E_start))

    k1, k2 = to_steps(E_vertex1), to_steps(E_vertex2)
    if k1 == 0 or k2 == 0 or np.sign(k1) == np.sign(k2):
        raise ValueError("the vertices must lie on opposite sides of E_start")
    d1, d2 = int(np.sign(k1)), int(np.sign(k2))
    seq = list(range(d1, k1 + d1, d1))           # start -> vertex 1
    seq += list(range(k1 - d1, k2 - d1, -d1))    # vertex 1 -> vertex 2
    seq += list(range(k2 - d2, -d2, -d2))        # vertex 2 -> start
    return E_start + np.array(seq) * step


def simulate_cv(E=None, scan_rate=0.1, E0=0.25, k0=1e-2, alpha=0.5, n=1, D=7e-6,
                c_ox=0.0, c_red=1e-6, area=0.1, T=298.15, Cdl=0.0, Ru=0.0, slope=0.0,
                offset=0.0, noise=0.0, cycles=2, substeps=8, seed=None):
    """Simulate a cyclic voltammogram; returns ``(E, I)`` of the last cycle.

    Parameters
    ----------
    E : potential program (default :func:`potential_program` with its defaults)
    scan_rate : V/s
    E0 : formal potential, V
    k0 : standard heterogeneous rate constant, cm/s (~1 reversible, ~1e-3 quasi-reversible)
    alpha : transfer coefficient
    n : number of electrons
    D : diffusion coefficient, cm^2/s (both species)
    c_ox, c_red : bulk concentrations of O and R, mol/cm^3
    area : electrode area, cm^2
    T : temperature, K
    Cdl : double-layer capacitance, F/cm^2 (charging current Cdl * area * dE/dt)
    Ru : uncompensated resistance, ohm (applied with the current of the previous step)
    slope : linear background, A/V (current slope * E added everywhere)
    offset : constant current offset, A
    noise : standard deviation of Gaussian noise, A
    cycles : number of simulated cycles; only the last one is returned
    substeps : time steps of the diffusion solver per potential step
    seed : random seed for the noise
    """
    E_prog = potential_program() if E is None else np.asarray(E, dtype=float)
    E_full = np.r_[E_prog[-1], E_prog]  # point 0 = start potential (= last point of the cycle)
    step_E = np.abs(np.diff(E_full))
    dt_point = np.median(step_E) / scan_rate
    dt = dt_point / substeps
    f = n * F / (R * T)
    c_tot = c_ox + c_red

    total_time = cycles * len(E_prog) * dt_point
    dx = np.sqrt(D * dt)  # lambda = 1
    L = 6 * np.sqrt(D * total_time)
    N = int(np.ceil(L / dx)) + 2
    lam = D * dt / dx**2
    c = np.full(N, c_ox)  # concentration of O; R = c_tot - c everywhere (equal D)

    # banded matrix (upper, diagonal, lower) for backward Euler
    ab = np.zeros((3, N))
    ab[0, 1:] = -lam
    ab[1, :] = 1 + 2 * lam
    ab[2, :-1] = -lam
    ab[1, -1], ab[2, -2] = 1.0, 0.0   # Dirichlet bulk boundary
    ab[0, 1] = -2 * lam               # surface node: half cell, flux from node 1

    out = np.empty(len(E_prog))
    i_prev = 0.0
    for cyc in range(cycles):
        for p in range(len(E_prog)):
            e_prev, e_now = (E_full[p], E_full[p + 1])
            dEdt = (e_now - e_prev) / dt_point
            i_f = 0.0
            for q in range(substeps):
                e = e_prev + (e_now - e_prev) * (q + 1) / substeps - i_prev * Ru
                eta = np.clip(f * (e - E0), -200, 200)
                kf, kb = k0 * np.exp(-alpha * eta), k0 * np.exp((1 - alpha) * eta)
                g = 2 * dt / dx
                ab[1, 0] = 1 + 2 * lam + g * (kf + kb)
                rhs = c.copy()
                rhs[0] += g * kb * c_tot
                rhs[-1] = c_ox
                c = solve_banded((1, 1), ab, rhs)
                rate = kf * c[0] - kb * (c_tot - c[0])  # net reduction rate, mol/(cm^2 s)
                i_f = -n * F * area * rate             # anodic current positive
            i_c = Cdl * area * dEdt
            i_prev = i_f + i_c
            out[p] = i_prev
    I = out + slope * E_prog + offset
    if noise:
        I = I + np.random.default_rng(seed).normal(0, noise, len(I))
    return E_prog.copy(), I


def randles_sevcik(n=1, area=0.1, D=7e-6, c=1e-6, scan_rate=0.1, T=298.15):
    """Peak current of a reversible wave, A (0.4463 n F A c sqrt(n F v D / R T))."""
    return 0.4463 * n * F * area * c * np.sqrt(n * F * scan_rate * D / (R * T))
