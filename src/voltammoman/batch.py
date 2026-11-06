"""Batch evaluation of many voltammogram files."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .analysis import RESULT_COLUMNS, CurveResult, analyze_curve
from .io import read_curve
from .peaks import PeakSettings

# Human-readable description of every result column (also written to the XLSX output).
COLUMN_DESCRIPTIONS = {
    "file": "file name",
    "n_points": "number of points of the evaluated cycle",
    "E_min": "lower vertex potential [V]",
    "E_max": "upper vertex potential [V]",
    "first_peak": "which peak was recorded first in the cycle (forward peak)",
    "EpA": "anodic peak potential [V]",
    "EpC": "cathodic peak potential [V]",
    "dEp": "peak separation EpA - EpC [V]",
    "E_half": "(EpA + EpC) / 2 [V]",
    "IpA": "anodic peak current above the baseline [A]",
    "IpC": "cathodic peak current below the baseline [A] (negative)",
    "ratio": "IpA / |IpC| from the baseline-corrected currents",
    "ratio_nicholson": "IpA / |IpC| estimated by Nicholson's formula (recommended for reversibility)",
    "wA": "anodic half-peak width |EpA - EpA/2| [V]",
    "wC": "cathodic half-peak width |EpC - EpC/2| [V]",
    "IpA_raw": "current at the anodic peak, without baseline [A]",
    "IpC_raw": "current at the cathodic peak, without baseline [A]",
    "baseline": "baseline method",
    "flags": "problems found (A: anodic, C: cathodic peak)",
    "error": "error message if the file could not be evaluated",
}


def find_files(paths, pattern: str = "*.txt", recursive: bool = False) -> list[Path]:
    """Expand files and directories into a sorted list of curve files."""
    if isinstance(paths, (str, Path)):
        paths = [paths]
    files = []
    for p in map(Path, paths):
        if p.is_dir():
            files += sorted(p.rglob(pattern) if recursive else p.glob(pattern))
        elif p.exists():
            files.append(p)
        else:
            raise FileNotFoundError(p)
    return files


def analyze_file(path, baseline: str = "flattest", settings: PeakSettings | None = None,
                 cycle: int = -1) -> CurveResult:
    """Read and evaluate one file."""
    cv = read_curve(path)
    return analyze_curve(cv["E"].to_numpy(), cv["I"].to_numpy(), baseline=baseline, settings=settings,
                         scan=cv["scan"].to_numpy() if "scan" in cv else None, cycle=cycle,
                         name=Path(path).name)


def analyze_files(paths, pattern: str = "*.txt", recursive: bool = False, baseline: str = "flattest",
                  settings: PeakSettings | None = None, cycle: int = -1, progress=None) -> pd.DataFrame:
    """Evaluate all curve files in ``paths`` (files and/or directories); one row per file.

    Files that cannot be read or evaluated get a row with the ``error`` column filled in, so a
    single bad file never stops the batch. ``progress`` is an optional callable taking
    ``(index, total, file)``.
    """
    files = find_files(paths, pattern, recursive)
    rows = []
    for k, f in enumerate(files):
        if progress:
            progress(k, len(files), f)
        try:
            rows.append(analyze_file(f, baseline, settings, cycle).values)
        except Exception as exc:  # noqa: BLE001 - report and continue
            row = dict.fromkeys(RESULT_COLUMNS, np.nan)
            row.update(file=f.name, baseline=baseline, flags="", error=f"{type(exc).__name__}: {exc}")
            rows.append(row)
    return pd.DataFrame(rows, columns=RESULT_COLUMNS)


def save_table(df: pd.DataFrame, path) -> Path:
    """Save results as ``.csv`` or ``.xlsx`` (the latter with a sheet describing the columns)."""
    path = Path(path)
    if path.suffix.lower() == ".xlsx":
        with pd.ExcelWriter(path) as xw:
            df.to_excel(xw, sheet_name="results", index=False)
            pd.DataFrame({"column": list(COLUMN_DESCRIPTIONS), "meaning": list(COLUMN_DESCRIPTIONS.values())}
                         ).to_excel(xw, sheet_name="columns", index=False)
    else:
        df.to_csv(path, index=False)
    return path
