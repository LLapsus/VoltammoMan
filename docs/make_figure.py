"""Generate docs/parameters.png: what VoltammoMan measures and which settings it uses.

Run from the repository root:  python docs/make_figure.py
"""

from pathlib import Path

import numpy as np
from matplotlib.figure import Figure

from voltammoman.analysis import analyze_curve
from voltammoman.peaks import PeakSettings
from voltammoman.simulate import potential_program, simulate_cv

OUT = Path(__file__).with_name("parameters.png")
UA = 1e6
RED, BLUE, GREY = "#c0392b", "#2471a3", "0.45"


def arrow(ax, x0, y0, x1, y1, color, text=None, tx=None, ty=None, **kw):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.2, shrinkA=0, shrinkB=0))
    if text:
        ax.text(tx if tx is not None else (x0 + x1) / 2, ty if ty is not None else (y0 + y1) / 2,
                text, color=color, fontsize=9, **kw)


def main():
    E, I = simulate_cv(E=potential_program(0, 0.9, -0.3), E0=0.22, k0=3e-3, c_red=1e-6,
                       Cdl=3e-4, slope=4e-6, noise=3e-8, seed=3)
    s = PeakSettings()
    res = analyze_curve(E, I, settings=s)
    v, pa, pc = res.values, res.peak_a, res.peak_c

    fig = Figure(figsize=(12, 5.0))
    ax1, ax2 = fig.subplots(1, 2)

    # ---------------- left: results -------------------------------------------------------
    ax = ax1
    ax.plot(res.E, res.I * UA, color="0.15", lw=1.1)
    ax.axhline(0, color="0.85", lw=.6, zorder=0)
    for p, br, c, sgn, lab in ((pa, res.anodic, RED, 1, "A"), (pc, res.cathodic, BLUE, -1, "C")):
        lo, hi = sorted(p.base_range)
        e = np.linspace(min(lo, p.Ep) - .03, max(hi, p.Ep) + .03, 50)
        ax.plot(e, np.polyval(p.baseline, e) * UA, "--", color=c, lw=1)
        ax.plot([lo, hi], np.polyval(p.baseline, [lo, hi]) * UA, color=c, lw=4, alpha=.3,
                solid_capstyle="butt")
        base = np.polyval(p.baseline, p.Ep)
        arrow(ax, p.Ep, base * UA, p.Ep, p.I_raw * UA, c,
              f"Ip{lab} = {p.Ip * UA:.2f} µA", tx=p.Ep + sgn * .025, ty=(base + p.I_raw) / 2 * UA,
              va="center", ha="left" if sgn > 0 else "right")
        ax.plot(p.Ep, p.I_raw * UA, "o", color=c, ms=6, zorder=5)
        ax.axvline(p.Ep, color=c, ls=":", lw=.8)
        # half-peak width: from the half-height point before the peak to Ep
        e_half = p.Ep - sgn * p.half_width
        i_half = np.polyval(p.baseline, e_half) + p.Ip / 2
        ax.plot([e_half, p.Ep], [i_half * UA] * 2, color=c, lw=1)
        ax.plot(e_half, i_half * UA, "|", color=c, ms=8)
        ax.text(min(e_half, p.Ep) + p.half_width / 2, i_half * UA + sgn * 1.2, f"w{lab}",
                color=c, fontsize=9, ha="center", va="bottom" if sgn > 0 else "top")
    ax.text(pa.Ep + .01, -26, f"EpA = {pa.Ep * 1e3:.0f} mV",
            color=RED, fontsize=9, rotation=90, va="bottom")
    ax.text(pc.Ep - .01, 22, f"EpC = {pc.Ep * 1e3:.0f} mV", color=BLUE, fontsize=9,
            rotation=90, va="top", ha="right")
    y = 30
    arrow(ax, pc.Ep, y, pa.Ep, y, "0.2", f"ΔEp = {v['dEp'] * 1e3:.0f} mV",
          ty=y + 1.2, ha="center", va="bottom")
    ax.plot(v["E_half"], y, "v", color="0.2", ms=6)
    ax.text(v["E_half"] + .01, y - 1.5, f"E_half = {v['E_half'] * 1e3:.0f} mV", fontsize=8.5,
            color="0.2", va="top")
    if res.reference:
        r = res.reference
        ax.plot([r["E_ref"]] * 2, [r["I_ref_forward"] * UA, r["I_ref_reverse"] * UA], "x",
                color=GREY, ms=7)
        ax.plot(r["E_switch"], r["I_switch"] * UA, "s", color=GREY, ms=5)
        ax.text(r["E_ref"] + .02, r["I_ref_reverse"] * UA - 1, "Nicholson\nreferences",
                color=GREY, fontsize=8, va="top")
    ax.set_title("What is measured (one row of the result table)\n"
                 f"ratio = IpA/|IpC| = {v['ratio']:.2f},   ratio_nicholson = "
                 f"{v['ratio_nicholson']:.2f}", fontsize=10)
    ax.set_xlabel("E [V]")
    ax.set_ylabel("I [µA]")
    ax.set_ylim(-28, 36)

    # ---------------- right: settings (anodic branch) -------------------------------------
    ax = ax2
    br = res.anodic
    ax.plot(br.E, br.I * UA, color="0.15", lw=1.1, label="anodic branch")
    Emin, Emax = br.E.min(), br.E.max()
    ymin, ymax = -12, 34
    ax.set_ylim(ymin, ymax)
    # edge zones at both vertices
    for x0, x1 in ((Emin, Emin + s.edge), (Emax - s.edge, Emax)):
        ax.axvspan(x0, x1, color="0.6", alpha=.25, lw=0)
    ax.text(Emin + s.edge / 2, ymax - 1, "edge", ha="center", va="top", fontsize=8.5, rotation=90)
    ax.text(Emax - s.edge / 2, ymax - 1, "edge", ha="center", va="top", fontsize=8.5, rotation=90)
    # linear window
    ax.axvspan(pa.Ep - s.base_from, pa.Ep - s.base_to, ymin=.0, ymax=.08, color="#7d3c98",
               alpha=.25, lw=0)
    ax.text(pa.Ep - (s.base_from + s.base_to) / 2, ymin + 1.8, "linear: fixed window",
            color="#7d3c98", fontsize=8.5, ha="center", va="center")
    arrow(ax, pa.Ep - s.base_from, ymin + 9, pa.Ep, ymin + 9, "#7d3c98", "base_from",
          tx=pa.Ep - s.base_from + .01, ty=ymin + 9.5, va="bottom")
    arrow(ax, pa.Ep - s.base_to, ymin + 5.5, pa.Ep, ymin + 5.5, "#7d3c98", "base_to",
          tx=pa.Ep - s.base_to + .01, ty=ymin + 6, va="bottom")
    # flattest: searched region and chosen window
    y = 25
    ax.annotate("", xy=(Emin + s.edge, y), xytext=(pa.Ep - s.gap, y),
                arrowprops=dict(arrowstyle="|-|", color="#1e8449", lw=1.2, mutation_scale=4))
    ax.text((Emin + s.edge + pa.Ep - s.gap) / 2, y + .8,
            "flattest: windows of width `width` are tried here", color="#1e8449", fontsize=8.5,
            ha="center", va="bottom")
    arrow(ax, pa.Ep - s.gap, y - 4, pa.Ep, y - 4, "#1e8449", "gap", ty=y - 3.4, ha="center",
          va="bottom")
    lo, hi = sorted(pa.base_range)
    ax.plot([lo, hi], np.polyval(pa.baseline, [lo, hi]) * UA, color=RED, lw=4, alpha=.3,
            solid_capstyle="butt")
    e = np.linspace(lo - .03, pa.Ep + .05, 50)
    ax.plot(e, np.polyval(pa.baseline, e) * UA, "--", color=RED, lw=1, label="baseline")
    yb = np.polyval(pa.baseline, hi) * UA + 4
    arrow(ax, lo, yb, hi, yb, RED, "chosen window (`width`)", ty=yb + .6, ha="center", va="bottom")
    # search
    ax.axvspan(pa.Ep - s.search, pa.Ep + s.search, color=RED, alpha=.12, lw=0)
    ax.text(pa.Ep + s.search + .01, pa.I_raw * UA, "±search: peak refined\nhere", color=RED,
            fontsize=8.5, va="center")
    ax.plot(pa.Ep, pa.I_raw * UA, "o", color=RED, ms=6, zorder=5)
    ax.set_title("What can be set (shown on the anodic branch; defaults in V)\n"
                 f"edge {s.edge}, search {s.search}, gap {s.gap}, width {s.width}, "
                 f"base_from {s.base_from}, base_to {s.base_to}", fontsize=10)
    ax.set_xlabel("E [V]")
    ax.set_ylabel("I [µA]")
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(.08, .78))

    fig.tight_layout()
    fig.savefig(OUT, dpi=110)
    print(OUT, res.values)


if __name__ == "__main__":
    main()
