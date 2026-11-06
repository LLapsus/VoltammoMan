"""Control plots of an evaluated voltammogram."""

from __future__ import annotations

import numpy as np

from .analysis import CurveResult


def plot_result(res: CurveResult, ax=None, title: str | None = None):
    """Plot the curve with both peaks, their baselines and heights; returns the figure.

    Currents are shown in µA. Dashed lines are the baselines, vertical segments the peak heights,
    the grey cross the reference levels used by Nicholson's formula.
    """
    if ax is None:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4.5))
    else:
        fig = ax.figure
    ax.plot(res.E, res.I * 1e6, color="0.2", lw=1, label="curve")
    for peak, branch, color, lab in ((res.peak_a, res.anodic, "C3", "anodic"),
                                     (res.peak_c, res.cathodic, "C0", "cathodic")):
        if peak is None or branch is None or not np.isfinite(peak.Ep):
            continue
        if peak.baseline is not None:
            lo, hi = sorted((peak.base_range[0], peak.base_range[1]))
            e = np.linspace(min(lo, peak.Ep), max(hi, peak.Ep), 50)
            ax.plot(e, np.polyval(peak.baseline, e) * 1e6, "--", color=color, lw=1)
            ax.plot([lo, hi], np.polyval(peak.baseline, [lo, hi]) * 1e6, color=color, lw=3, alpha=.35)
            base = np.polyval(peak.baseline, peak.Ep)
            ax.plot([peak.Ep, peak.Ep], [base * 1e6, peak.I_raw * 1e6], color=color, lw=1.5)
        ax.plot(peak.Ep, peak.I_raw * 1e6, "o", color=color, ms=5,
                label=f"{lab}: Ep = {peak.Ep * 1e3:.0f} mV" +
                      (f", Ip = {peak.Ip * 1e6:.2f} µA" if np.isfinite(peak.Ip) else ""))
    if res.reference:
        r = res.reference
        ax.plot([r["E_ref"]] * 2, [r["I_ref_forward"] * 1e6, r["I_ref_reverse"] * 1e6], "x", color="0.5")
        ax.plot(r["E_switch"], r["I_switch"] * 1e6, "s", color="0.5", ms=4)
    v = res.values
    sub = []
    if np.isfinite(v.get("dEp", np.nan)):
        sub.append(f"ΔEp = {v['dEp'] * 1e3:.0f} mV")
    if np.isfinite(v.get("ratio_nicholson", np.nan)):
        sub.append(f"IpA/IpC (Nicholson) = {v['ratio_nicholson']:.2f}")
    if v.get("flags"):
        sub.append(f"flags: {v['flags']}")
    ax.set_title((title or v.get("file") or "") + ("\n" + ", ".join(sub) if sub else ""), fontsize=9)
    ax.set_xlabel("E [V]")
    ax.set_ylabel("I [µA]")
    ax.axhline(0, color="0.8", lw=.5, zorder=0)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    return fig


def save_plots(results, directory, dpi: int = 110):
    """Save one PNG per :class:`CurveResult` into ``directory``; returns the list of paths."""
    from pathlib import Path

    from matplotlib.figure import Figure  # no pyplot: works headless and leaves notebooks alone

    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    out = []
    for res in results:
        fig = Figure(figsize=(7, 4.5))
        plot_result(res, ax=fig.subplots())
        p = directory / (Path(res.values["file"] or "curve").stem + ".png")
        fig.savefig(p, dpi=dpi)
        out.append(p)
    return out
