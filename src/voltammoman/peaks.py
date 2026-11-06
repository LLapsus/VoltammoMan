"""Peak detection and baseline correction on one branch of a voltammogram.

Three baseline methods are available; all fit a straight line to the part of the branch
recorded *before* the peak (the capacitive / background region) and extrapolate it under the
peak:

``flattest`` (default)
    least-squares line through the most linear window of width ``width`` V in the background
    region (at least ``gap`` V before the peak); among windows that are about equally linear,
    the one closest to the peak is used. Most accurate on simulated curves and closest to
    manually evaluated NOVA values.
``linear``
    least-squares line through a fixed window ``base_from``..``base_to`` V before the peak
    (simple and predictable; too high when the faradaic wave starts early, e.g. slow kinetics).
``tangent``
    tangent to the (smoothed) curve at the point of the background region where it rises least
    steeply towards the peak, i.e. at the foot of the wave. This mimics a hand-drawn baseline.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.signal import find_peaks, savgol_filter

from .cycle import Branch

BASELINES = ("flattest", "linear", "tangent")


@dataclass
class PeakSettings:
    """Tunable parameters of the peak search (potentials in V)."""

    edge: float = 0.03       # peaks closer than this to a vertex are rejected
    search: float = 0.02     # the peak is refined within this distance of the raw extremum
    gap: float = 0.10        # tangent/flattest: background region ends this far before the peak
    width: float = 0.10      # flattest: width of the sliding window
    smooth: int = 11         # tangent: Savitzky-Golay window length (points) for the derivative
    base_from: float = 0.35  # linear: window starts this far before the peak ...
    base_to: float = 0.15    # ... and ends this far before it
    broad: float = 0.20      # flag a peak whose half-height point is further than this from Ep
    min_rel: float = 0.05    # flag a peak lower than this fraction of the branch current range ...
    min_snr: float = 10.0    # ... or lower than this multiple of the noise level
    min_points: int = 8      # minimal number of points of a baseline fit


@dataclass
class Peak:
    """Result for one branch. ``Ip`` is signed (anodic > 0, cathodic < 0) and baseline-corrected."""

    Ep: float = np.nan
    Ip: float = np.nan
    I_raw: float = np.nan             # current at the peak, relative to zero
    half_width: float = np.nan        # |Ep - Ep/2|, Ep/2 = potential at half height before the peak
    baseline: tuple | None = None     # (slope, intercept) of the baseline, I = slope * E + intercept
    base_range: tuple | None = None   # potential range used for the baseline
    index: int | None = None          # index into the branch
    t: int | None = None              # position of the peak in the cycle
    flags: list = field(default_factory=list)


def _line_fit(E, I):
    return tuple(np.polyfit(E, I, 1))


def _smooth(I, window):
    """Savitzky-Golay smoothing (quadratic) with an odd window no longer than the data."""
    n = len(I)
    win = window if window % 2 else window + 1
    win = min(win, n if n % 2 else n - 1)
    return savgol_filter(I, win, 2) if win >= 5 else I


def _baseline(E, I, k, sign, method, s: PeakSettings):
    """Baseline (slope, intercept) and its potential range for a peak at branch index ``k``."""
    progress = sign * (E - E[0])          # distance travelled since the start of the branch
    before = sign * (E[k] - E)            # distance before the peak
    if method == "linear":
        m = (before >= s.base_to) & (before <= s.base_from) & (progress >= s.edge)
        if m.sum() < s.min_points:
            return None, None
        return _line_fit(E[m], I[m]), (E[m].min(), E[m].max())

    region = np.flatnonzero((before >= s.gap) & (progress >= s.edge))
    if len(region) < s.min_points:
        return None, None
    if method == "flattest":
        # Windows are compared on the smoothed current, so that noise does not decide which
        # window is the most linear; the line itself is fitted to the raw current.
        Is = _smooth(I, s.smooth)
        Er = E[region]
        span = abs(Er[-1] - Er[0])
        if span <= s.width:
            sel = region
        else:
            wins = []  # (distance of the window from the peak, rms of a line fit, indices)
            for lo in np.arange(0, span - s.width + 1e-12, s.width / 4):
                w = region[(sign * (Er - Er[0]) >= lo) & (sign * (Er - Er[0]) <= lo + s.width)]
                if len(w) < s.min_points:
                    continue
                coef = np.polyfit(E[w], Is[w], 1)
                rms = np.sqrt(np.mean((Is[w] - np.polyval(coef, E[w])) ** 2))
                wins.append((before[w].min(), rms, w))
            if not wins:
                return None, None
            # Among the windows that are about as linear as the best one (within the noise left
            # after smoothing), take the one closest to the peak.
            noise = np.std(I[region] - Is[region])
            tol = 1.5 * min(w[1] for w in wins) + 0.3 * noise
            sel = min((w for w in wins if w[1] <= tol), key=lambda w: w[0])[2]
        return _line_fit(E[sel], I[sel]), (E[sel].min(), E[sel].max())

    if method == "tangent":
        Is = _smooth(I, s.smooth)
        dIdE = np.gradient(Is, E)
        # Along the sweep, d(sign*I)/d(progress) = dI/dE for both directions, so the point where
        # the curve rises least steeply towards the peak is the minimum of dI/dE.
        j = region[np.argmin(dIdE[region])]
        slope = dIdE[j]
        return (slope, Is[j] - slope * E[j]), (E[j], E[j])

    raise ValueError(f"unknown baseline method {method!r}; use one of {BASELINES}")


def find_peak(branch: Branch, method: str = "flattest", settings: PeakSettings | None = None) -> Peak:
    """Locate the peak on one branch and measure it against a baseline.

    The most prominent local extremum (maximum on an anodic, minimum on a cathodic branch) away
    from the vertices is taken as the peak. A baseline is constructed (see the module docstring),
    the peak is moved to the extremum of the baseline-corrected current within ``search`` V and
    the baseline is rebuilt for the final position. Problems are reported in ``flags``:
    ``no_peak``, ``peak_at_vertex``, ``no_baseline``, ``nonpositive_height``, ``weak_peak``
    (lower than ``min_rel`` of the current range of the branch or ``min_snr`` times the noise),
    ``broad_peak``.
    """
    s = settings or PeakSettings()
    E, I, sign = branch.E, branch.I, branch.direction
    out = Peak()
    progress, remaining = sign * (E - E[0]), sign * (E[-1] - E)
    inner = np.flatnonzero((progress >= s.edge) & (remaining >= s.edge))
    if len(inner) < 3:
        out.flags.append("no_peak")
        return out
    cand, props = find_peaks(sign * I[inner], prominence=0)
    if len(cand) == 0:
        out.flags.append("no_peak")
        return out
    k = inner[cand[np.argmax(props["prominences"])]]

    coef, rng = _baseline(E, I, k, sign, method, s)
    if coef is not None:
        near = np.flatnonzero(np.abs(E - E[k]) <= s.search)
        k = near[np.argmax(sign * (I[near] - np.polyval(coef, E[near])))]
        coef, rng = _baseline(E, I, k, sign, method, s)

    out.Ep, out.I_raw, out.index, out.t = E[k], I[k], int(k), int(branch.t[k])
    if progress[k] < s.edge or remaining[k] < s.edge:
        out.flags.append("peak_at_vertex")
        out.Ep = np.nan
        return out
    if coef is None:
        out.flags.append("no_baseline")
        return out
    out.baseline, out.base_range = coef, rng
    corr = I - np.polyval(coef, E)
    Ip = corr[k]
    if sign * Ip <= 0:
        out.flags.append("nonpositive_height")
        return out
    out.Ip = Ip
    noise = np.std(I - _smooth(I, s.smooth))
    if abs(Ip) < max(s.min_rel * (I.max() - I.min()), s.min_snr * noise):
        out.flags.append("weak_peak")
    # half-height point on the rising side of the peak
    below = np.flatnonzero(sign * corr[:k] <= sign * Ip / 2)
    if len(below):
        j = below[-1]
        e_half = np.interp(sign * Ip / 2, [sign * corr[j], sign * corr[j + 1]], [E[j], E[j + 1]])
        out.half_width = abs(out.Ep - e_half)
    if not len(below) or out.half_width > s.broad:
        out.flags.append("broad_peak")
    return out
