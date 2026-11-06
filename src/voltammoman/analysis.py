"""Evaluation of whole voltammograms: both peaks, their separation and reversibility."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .cycle import Branch, select_cycle, split_branches
from .peaks import Peak, PeakSettings, find_peak

# Columns of the result table, in order (see README for their meaning).
RESULT_COLUMNS = [
    "file", "n_points", "E_min", "E_max", "first_peak",
    "EpA", "EpC", "dEp", "E_half", "IpA", "IpC", "ratio", "ratio_nicholson",
    "wA", "wC", "IpA_raw", "IpC_raw", "baseline", "flags", "error",
]


@dataclass
class CurveResult:
    """All values for one curve (``values``) plus the intermediate objects needed for plotting."""

    values: dict
    E: np.ndarray
    I: np.ndarray
    anodic: Branch | None = None
    cathodic: Branch | None = None
    peak_a: Peak | None = None
    peak_c: Peak | None = None
    reference: dict = field(default_factory=dict)


def nicholson_ratio(forward: Branch, reverse: Branch, pf: Peak, pr: Peak, ref_offset: float = 0.05):
    """Reverse-to-forward peak-current ratio by Nicholson's empirical formula.

    R. S. Nicholson, Anal. Chem. 38 (1966) 1406::

        i_rev / i_fwd = (i_rev)0 / i_fwd + 0.485 (i_sp)0 / i_fwd + 0.086

    ``i_fwd`` is the baseline-corrected forward peak current, ``(i_rev)0`` the reverse peak current
    and ``(i_sp)0`` the current at the switching potential, both measured from a reference level
    instead of a baseline. Nicholson uses zero current as the reference; here each branch is
    referred to its own current at ``ref_offset`` V after the forward branch starts, i.e. well
    before the forward wave, which removes the charging current and any constant offset.

    Returns ``(ratio, reference)`` where ``reference`` holds the levels used (for plotting).
    """
    s = forward.direction
    invalid = {"no_peak", "peak_at_vertex", "nonpositive_height", "no_baseline", "weak_peak"}
    if not (np.isfinite(pf.Ip) and np.isfinite(pr.I_raw)) or invalid & set(pr.flags):
        return np.nan, {}
    e_ref = forward.E[0] + s * ref_offset
    order_f, order_r = np.argsort(forward.E), np.argsort(reverse.E)
    i_ref_f = np.interp(e_ref, forward.E[order_f], forward.I[order_f])
    i_ref_r = np.interp(e_ref, reverse.E[order_r], reverse.I[order_r])
    i_fwd = s * pf.Ip
    i_sp = s * (forward.I[-1] - i_ref_f)
    i_rev = -s * (pr.I_raw - i_ref_r)
    ratio = i_rev / i_fwd + 0.485 * i_sp / i_fwd + 0.086
    return ratio, dict(E_ref=e_ref, I_ref_forward=i_ref_f, I_ref_reverse=i_ref_r,
                       E_switch=forward.E[-1], I_switch=forward.I[-1])


def analyze_curve(E, I, baseline: str = "flattest", settings: PeakSettings | None = None,
                  scan=None, cycle: int = -1, name: str | None = None) -> CurveResult:
    """Evaluate one voltammogram.

    Returns a :class:`CurveResult` whose ``values`` contain (potentials in V, currents in A):

    - ``EpA``, ``EpC``: anodic / cathodic peak potential; ``dEp = EpA - EpC``;
      ``E_half = (EpA + EpC) / 2``,
    - ``IpA``, ``IpC``: baseline-corrected peak currents (``IpC < 0``),
    - ``ratio = IpA / |IpC|`` from the baseline-corrected currents,
    - ``ratio_nicholson``: the same ratio estimated by Nicholson's formula (more reliable for
      the peak recorded second, whose baseline is not a straight line),
    - ``wA``, ``wC``: half-peak widths |Ep - Ep/2|,
    - ``first_peak``: which peak was recorded first in the cycle ("anodic"/"cathodic"),
    - ``flags``: problems found, prefixed by ``A:``/``C:`` for the anodic/cathodic peak.
    """
    E, I = select_cycle(E, I, scan=scan, cycle=cycle)
    v = dict.fromkeys(RESULT_COLUMNS, np.nan)
    v.update(file=name, n_points=len(E), E_min=float(np.min(E)), E_max=float(np.max(E)),
             baseline=baseline, flags="", error="")
    res = CurveResult(values=v, E=E, I=I)
    an, ca = split_branches(E, I)
    pa, pc = find_peak(an, baseline, settings), find_peak(ca, baseline, settings)
    res.anodic, res.cathodic, res.peak_a, res.peak_c = an, ca, pa, pc

    v.update(EpA=pa.Ep, EpC=pc.Ep, IpA=pa.Ip, IpC=pc.Ip, wA=pa.half_width, wC=pc.half_width,
             IpA_raw=pa.I_raw, IpC_raw=pc.I_raw)
    v["dEp"] = pa.Ep - pc.Ep
    v["E_half"] = (pa.Ep + pc.Ep) / 2
    v["ratio"] = pa.Ip / -pc.Ip
    if pa.t is not None and pc.t is not None:
        anodic_first = pa.t < pc.t
        v["first_peak"] = "anodic" if anodic_first else "cathodic"
        fwd, rev, pf, pr = (an, ca, pa, pc) if anodic_first else (ca, an, pc, pa)
        r, res.reference = nicholson_ratio(fwd, rev, pf, pr)
        v["ratio_nicholson"] = 1 / r if anodic_first else r  # always IpA / |IpC|
    else:
        v["first_peak"] = ""
    v["flags"] = ";".join([f"A:{f}" for f in pa.flags] + [f"C:{f}" for f in pc.flags])
    return res
