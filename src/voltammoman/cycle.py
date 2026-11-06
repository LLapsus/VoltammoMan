"""Splitting a recorded voltammogram into cycles and into anodic / cathodic branches."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Branch:
    """One sweep direction of a cycle, in the order in which it was recorded.

    ``E`` increases along an anodic branch and decreases along a cathodic one. ``t`` is the
    position of every point in the selected cycle (0, 1, 2, ...), used to tell which peak
    was recorded first. When the cycle starts inside the potential window, the branch that
    was interrupted at the start is joined from the end and the beginning of the cycle,
    so its ``t`` jumps.
    """

    E: np.ndarray
    I: np.ndarray
    t: np.ndarray
    direction: int  # +1 anodic (increasing E), -1 cathodic

    def __len__(self):
        return len(self.E)


def _directions(E):
    """Sign of every potential step; zero steps take the sign of the previous step."""
    d = np.sign(np.diff(E))
    for k in range(len(d)):
        if d[k] == 0:
            d[k] = d[k - 1] if k > 0 else 0
    if d[0] == 0:  # a flat start takes the first non-zero direction
        nz = np.flatnonzero(d)
        d[: nz[0] if len(nz) else len(d)] = d[nz[0]] if len(nz) else 1
    return d


def turning_points(E) -> np.ndarray:
    """Indices where the scan direction changes (the vertices)."""
    d = _directions(np.asarray(E, dtype=float))
    return np.flatnonzero(d[1:] != d[:-1]) + 1


def select_cycle(E, I, scan=None, cycle=-1):
    """Return ``(E, I)`` of one cycle of a possibly multi-cycle recording.

    With a ``scan`` column (scan number from NOVA), the cycle is chosen by it. Otherwise the
    recording is cut each time the potential passes the starting potential again in the
    starting direction. ``cycle`` is a zero-based index; negative values count from the end,
    and the default ``-1`` is the last *complete* cycle (with two vertices).
    """
    E, I = np.asarray(E, dtype=float), np.asarray(I, dtype=float)
    if scan is not None:
        scan = np.asarray(scan)
        ids = list(dict.fromkeys(scan[~np.isnan(scan)]))
        if len(ids) > 1:
            m = scan == ids[cycle]
            return E[m], I[m]
        return E, I
    if len(turning_points(E)) <= 2:
        return E, I
    d0 = _directions(E)[0]
    rel = d0 * (E - E[0])
    cuts = [0] + [k for k in range(1, len(E)) if rel[k - 1] < 0 <= rel[k]] + [len(E)]
    pieces = [(a, b) for a, b in zip(cuts[:-1], cuts[1:]) if b - a > 2]
    complete = [p for p in pieces if len(turning_points(E[p[0]:p[1]])) >= 2]
    chosen = (complete or pieces)[cycle]
    return E[chosen[0]:chosen[1]], I[chosen[0]:chosen[1]]


def split_branches(E, I, step_tol=1.5):
    """Split one cycle into its anodic and cathodic :class:`Branch`.

    Handles cycles that start at a vertex (two segments) or inside the window (three segments,
    e.g. 0 V -> E_max -> E_min -> 0 V as in NOVA exports). In the latter case the cycle is a
    steady-state one: its end continues into its beginning, so the interrupted branch is joined
    from the last and the first segment.
    """
    E, I = np.asarray(E, dtype=float), np.asarray(I, dtype=float)
    t = np.arange(len(E))
    tp = turning_points(E)
    bounds = np.r_[0, tp, len(E) - 1]
    segs = []  # (start, end) inclusive, direction
    for a, b in zip(bounds[:-1], bounds[1:]):
        if b > a:
            segs.append((a, b, int(np.sign(E[b] - E[a]))))
    if not any(s[2] > 0 for s in segs) or not any(s[2] < 0 for s in segs):
        raise ValueError("the curve has no complete forward and reverse sweep")
    step = np.median(np.abs(np.diff(E))) or 1e-3
    closed = abs(E[-1] - E[0]) <= step_tol * step

    def branch(direction):
        same = [s for s in segs if s[2] == direction]
        if len(same) >= 2 and closed and same[0][0] == 0 and same[-1][1] == len(E) - 1:
            (a1, b1, _), (a2, b2, _) = same[-1], same[0]
            start = a2 + 1 if abs(E[a2] - E[b1]) < step / 2 else a2  # skip a duplicated joint point
            idx = np.r_[np.arange(a1, b1 + 1), np.arange(start, b2 + 1)]
        else:
            a, b, _ = max(same, key=lambda s: s[1] - s[0])
            idx = np.arange(a, b + 1)
        return Branch(E[idx], I[idx], t[idx], direction)

    return branch(+1), branch(-1)
