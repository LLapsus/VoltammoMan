"""Reading and writing voltammograms in the NOVA (Metrohm Autolab) ASCII export format.

A typical NOVA export looks like this (``;`` separator, decimal point, CRLF line endings)::

    Potential applied (V);Smoothed WE(1).Current (A)
    0.00244140625;-3.44506730248918E-06
    0.0048828125;-3.43248639787946E-06
    ...

:func:`read_curve` also accepts other common variants: ``,``/tab separators, a decimal comma,
extra columns (time, scan number, ...) and other column names containing "potential" and
"current".
"""

from __future__ import annotations

import csv
import io
import math
from pathlib import Path

import numpy as np
import pandas as pd

NOVA_HEADER = ("Potential applied (V)", "Smoothed WE(1).Current (A)")


# --------------------------------------------------------------------------
# Reading
# --------------------------------------------------------------------------

def _pick_column(columns, keyword, prefer=()):
    """Name of the first column containing ``keyword``, preferring those containing ``prefer``."""
    cands = [c for c in columns if keyword in c.lower()]
    for p in prefer:
        best = [c for c in cands if p in c.lower()]
        if best:
            return best[0]
    return cands[0] if cands else None


def read_curve(path) -> pd.DataFrame:
    """Read one exported voltammogram.

    Returns a DataFrame with columns ``E`` (potential, V) and ``I`` (current, A), plus ``scan``
    when the file contains a scan-number column. Rows with non-numeric values are dropped.
    """
    text = Path(path).read_text(encoding="utf-8-sig", errors="replace")
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError(f"{path}: empty file")
    try:
        sep = csv.Sniffer().sniff(lines[0], delimiters=";,\t").delimiter
    except csv.Error:
        sep = ";" if ";" in lines[0] else ("\t" if "\t" in lines[0] else ",")
    body = "\n".join(lines[1:])
    # A decimal comma is only possible when the separator is not a comma.
    decimal = "," if sep != "," and "," in body and "." not in body else "."
    df = pd.read_csv(io.StringIO("\n".join(lines)), sep=sep, decimal=decimal)
    df.columns = [str(c).strip() for c in df.columns]

    e_col = _pick_column(df.columns, "potential", prefer=("applied", "we(1)"))
    i_col = _pick_column(df.columns, "current", prefer=("we(1).current", "smoothed"))
    if e_col is None or i_col is None:  # no recognisable header: take the first two columns
        if df.shape[1] < 2:
            raise ValueError(f"{path}: expected at least two columns (potential, current)")
        e_col, i_col = df.columns[:2]
    out = pd.DataFrame({"E": pd.to_numeric(df[e_col], errors="coerce"),
                        "I": pd.to_numeric(df[i_col], errors="coerce")})
    scan_col = _pick_column(df.columns, "scan")
    if scan_col is not None:
        out["scan"] = pd.to_numeric(df[scan_col], errors="coerce")
    out = out.dropna(subset=["E", "I"]).reset_index(drop=True)
    if len(out) < 10:
        raise ValueError(f"{path}: fewer than 10 data points")
    return out


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------

def format_potential(x: float) -> str:
    """Potential as NOVA writes it: the shortest exact decimal, ``0`` for zero (e.g. ``0.8984375``)."""
    x = float(x)
    if x == 0:
        return "0"
    s = repr(x)
    return s[:-2] if s.endswith(".0") else s


def format_current(x: float) -> str:
    """Current as NOVA writes it: 15 significant digits, no trailing zeros (``-3.44506730248918E-06``)."""
    x = float(x)
    if x == 0 or not math.isfinite(x):
        return "0" if x == 0 else repr(x)
    mant, exp = f"{x:.14E}".split("E")
    mant = mant.rstrip("0").rstrip(".")
    return f"{mant}E{exp}"


def write_curve(path, E, I, header=NOVA_HEADER) -> Path:
    """Write a voltammogram in the NOVA export format (``;`` separator, CRLF line endings)."""
    path = Path(path)
    rows = [";".join(header)]
    rows += [f"{format_potential(e)};{format_current(i)}" for e, i in zip(np.asarray(E), np.asarray(I))]
    path.write_bytes(("\r\n".join(rows) + "\r\n").encode("ascii"))
    return path
