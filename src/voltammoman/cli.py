"""Command-line interface: ``vom analyze``, ``vom simulate``, ``vom examples``, ``vom gui``."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .peaks import BASELINES, PeakSettings


def _add_settings(p):
    g = p.add_argument_group("peak search settings (potentials in V)")
    d = PeakSettings()
    g.add_argument("--edge", type=float, default=d.edge, help="reject peaks closer than this to a vertex")
    g.add_argument("--gap", type=float, default=d.gap,
                   help="flattest/tangent: the baseline region ends this far before the peak")
    g.add_argument("--width", type=float, default=d.width, help="flattest: width of the baseline window")
    g.add_argument("--base-from", type=float, default=d.base_from, help="linear: window start before the peak")
    g.add_argument("--base-to", type=float, default=d.base_to, help="linear: window end before the peak")
    g.add_argument("--smooth", type=int, default=d.smooth, help="smoothing window (points)")


def _settings(a) -> PeakSettings:
    return PeakSettings(edge=a.edge, gap=a.gap, width=a.width, base_from=a.base_from,
                        base_to=a.base_to, smooth=a.smooth)


def cmd_analyze(a):
    import numpy as np
    import pandas as pd

    from .analysis import RESULT_COLUMNS
    from .batch import analyze_file, find_files, save_table

    files = find_files(a.paths or ["."], a.pattern, a.recursive)
    if not files:
        print(f"no files matching {a.pattern!r} found", file=sys.stderr)
        return 1
    settings = _settings(a)
    rows, results = [], []
    for k, f in enumerate(files, 1):
        try:
            res = analyze_file(f, a.baseline, settings, a.cycle)
            rows.append(res.values)
            results.append(res)
        except Exception as exc:  # noqa: BLE001
            row = dict.fromkeys(RESULT_COLUMNS, np.nan)
            row.update(file=f.name, baseline=a.baseline, flags="", error=f"{type(exc).__name__}: {exc}")
            rows.append(row)
        if not a.quiet:
            print(f"\r[{k}/{len(files)}] {f.name:<60.60}", end="", file=sys.stderr)
    if not a.quiet:
        print(file=sys.stderr)
    df = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    out = save_table(df, a.output)
    n_err, n_flag = (df["error"] != "").sum(), (df["flags"].fillna("") != "").sum()
    print(f"{len(df)} curves evaluated -> {out}  ({n_flag} with flags, {n_err} errors)")
    if a.plots:
        from .plot import save_plots
        paths = save_plots(results, a.plots)
        print(f"{len(paths)} control plots -> {a.plots}/")
    return 0


def cmd_simulate(a):
    from .io import write_curve
    from .simulate import potential_program, simulate_cv

    E = potential_program(a.start, a.vertex1, a.vertex2)
    E, I = simulate_cv(E, scan_rate=a.scan_rate, E0=a.E0, k0=a.k0, alpha=a.alpha, n=a.n, D=a.D,
                       c_ox=a.c_ox * 1e-6, c_red=a.c_red * 1e-6, area=a.area, Cdl=a.Cdl * 1e-6,
                       Ru=a.Ru, slope=a.slope * 1e-6, offset=a.offset * 1e-6, noise=a.noise * 1e-6,
                       cycles=a.cycles, seed=a.seed)
    print(write_curve(a.output, E, I))
    return 0


def cmd_examples(a):
    from .examples import write_examples
    paths = write_examples(a.directory)
    print(f"{len(paths)} example curves -> {a.directory}/")
    return 0


def cmd_gui(a):
    from .gui import launch
    launch(port=a.port, share=a.share, inbrowser=not a.no_browser)
    return 0


def build_parser():
    p = argparse.ArgumentParser(prog="vom", description="VoltammoMan: batch evaluation of cyclic voltammograms.")
    p.add_argument("--version", action="version", version=f"VoltammoMan {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    an = sub.add_parser("analyze", help="evaluate curve files and write a table",
                        description="Evaluate all curve files (default: *.txt in the current directory).")
    an.add_argument("paths", nargs="*", help="files and/or directories (default: current directory)")
    an.add_argument("-o", "--output", default="voltammoman_results.xlsx",
                    help="output table, .xlsx or .csv (default: %(default)s)")
    an.add_argument("--plots", nargs="?", const="voltammoman_plots", default=None, metavar="DIR",
                    help="also save one control plot per curve (default dir: %(const)s)")
    an.add_argument("-b", "--baseline", choices=BASELINES, default="flattest", help="baseline method")
    an.add_argument("-p", "--pattern", default="*.txt", help="file pattern in directories (default: %(default)s)")
    an.add_argument("-r", "--recursive", action="store_true", help="search directories recursively")
    an.add_argument("--cycle", type=int, default=-1, help="cycle to evaluate in multi-cycle files (default: last)")
    an.add_argument("-q", "--quiet", action="store_true", help="no progress output")
    _add_settings(an)
    an.set_defaults(func=cmd_analyze)

    si = sub.add_parser("simulate", help="simulate one curve in the NOVA export format")
    si.add_argument("output", help="output .txt file")
    si.add_argument("--start", type=float, default=0.0, help="start potential [V]")
    si.add_argument("--vertex1", type=float, default=0.9, help="first vertex [V]")
    si.add_argument("--vertex2", type=float, default=-0.3, help="second vertex [V]")
    si.add_argument("--scan-rate", type=float, default=0.1, help="[V/s]")
    si.add_argument("--E0", type=float, default=0.22, help="formal potential [V]")
    si.add_argument("--k0", type=float, default=1e-2, help="rate constant [cm/s]")
    si.add_argument("--alpha", type=float, default=0.5, help="transfer coefficient")
    si.add_argument("--n", type=int, default=1, help="number of electrons")
    si.add_argument("--D", type=float, default=7e-6, help="diffusion coefficient [cm^2/s]")
    si.add_argument("--c-ox", type=float, default=0.0, help="bulk concentration of the oxidised form [mM]")
    si.add_argument("--c-red", type=float, default=1.0, help="bulk concentration of the reduced form [mM]")
    si.add_argument("--area", type=float, default=0.1, help="electrode area [cm^2]")
    si.add_argument("--Cdl", type=float, default=0.0, help="double-layer capacitance [µF/cm^2]")
    si.add_argument("--Ru", type=float, default=0.0, help="uncompensated resistance [ohm]")
    si.add_argument("--slope", type=float, default=0.0, help="linear background [µA/V]")
    si.add_argument("--offset", type=float, default=0.0, help="current offset [µA]")
    si.add_argument("--noise", type=float, default=0.0, help="noise standard deviation [µA]")
    si.add_argument("--cycles", type=int, default=2, help="simulated cycles (the last one is written)")
    si.add_argument("--seed", type=int, default=None, help="random seed for the noise")
    si.set_defaults(func=cmd_simulate)

    ex = sub.add_parser("examples", help="write a set of synthetic example curves")
    ex.add_argument("directory", nargs="?", default="voltammoman_examples")
    ex.set_defaults(func=cmd_examples)

    gu = sub.add_parser("gui", help="start the graphical interface in the web browser (needs gradio)")
    gu.add_argument("--port", type=int, default=None, help="port (default: first free from 7860)")
    gu.add_argument("--share", action="store_true", help="create a public gradio.live link (data leave the computer!)")
    gu.add_argument("--no-browser", action="store_true", help="do not open the browser automatically")
    gu.set_defaults(func=cmd_gui)
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
