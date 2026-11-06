"""VoltammoMan: batch evaluation of cyclic voltammograms.

Typical use::

    from voltammoman import analyze_files, save_table
    df = analyze_files("path/to/curves")          # one row per file
    save_table(df, "results.xlsx")

or from the command line: ``vom analyze path/to/curves -o results.xlsx --plots``.
"""

__version__ = "0.1.0"

from .analysis import CurveResult, analyze_curve, nicholson_ratio  # noqa: E402
from .batch import analyze_file, analyze_files, find_files, save_table  # noqa: E402
from .io import read_curve, write_curve  # noqa: E402
from .peaks import BASELINES, PeakSettings, find_peak  # noqa: E402
from .simulate import potential_program, randles_sevcik, simulate_cv  # noqa: E402

__all__ = [
    "BASELINES", "CurveResult", "PeakSettings", "analyze_curve", "analyze_file", "analyze_files",
    "find_files", "find_peak", "nicholson_ratio", "potential_program", "randles_sevcik",
    "read_curve", "save_table", "simulate_cv", "write_curve",
]
