"""Load MFiX monitor CSV output into pandas.

MFiX writes one text file per monitor region (keyword ``monitor_name(N)``).
The format is a short header block of lines starting with ``#`` followed by
whitespace- or comma-separated data rows whose first column is time [s] and
remaining columns are the monitored quantities in SI units.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


_HEADER_PREFIXES = ("#", "!")


def _scan_monitor_layout(path: Path) -> tuple[int, list[str]]:
    """Return (n_skip, column_names) for an MFiX monitor CSV/DAT file.

    MFiX 26.1.2 writes monitor files in one of two layouts:

    (A) A short block of ``#``-prefixed metadata lines, then a *quoted*,
        comma-separated header line (e.g. ``"Time","p_g"``), then data rows.

    (B) A short block of ``#``/``!`` metadata lines where the last metadata
        line itself holds the column names (whitespace-separated), followed
        directly by data rows.

    When MFiX is launched multiple times with the same ``run_name``, each
    launch APPENDS a fresh header block to the existing monitor file. This
    function skips to the row immediately after the *last* header line, so
    only the most recent run's rows are returned.
    """
    with path.open("r") as fh:
        lines = fh.readlines()

    # Find the index of the last header row (either '"Time"'-style quoted
    # header or last '#'-prefixed metadata line that looks like column names).
    last_header_idx = -1
    names: list[str] = []
    for i, raw in enumerate(lines):
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped.startswith(_HEADER_PREFIXES):
            tokens = stripped.lstrip("#! ").split()
            if tokens and not tokens[0].replace(".", "", 1).replace(
                "-", "", 1
            ).replace("e", "", 1).replace("E", "", 1).isdigit():
                last_header_idx = i
                names = tokens
            else:
                last_header_idx = i
                names = []
            continue
        if '"' in stripped:
            last_header_idx = i
            names = [tok.strip().strip('"') for tok in stripped.split(",") if tok.strip()]
            continue
    n_skip = last_header_idx + 1
    return n_skip, names


def read_monitor(path: str | Path) -> pd.DataFrame:
    """Read a single MFiX monitor file into a DataFrame indexed by time [s].

    Parameters
    ----------
    path : str | Path
        Path to a monitor output file (``*.csv``, ``*.dat``, or no extension).

    Returns
    -------
    pandas.DataFrame
        Index ``time`` [s]; one column per monitored field.
    """
    p = Path(path)
    n_skip, names = _scan_monitor_layout(p)
    kwargs: dict = {"skiprows": n_skip, "sep": r"[,\s]+", "engine": "python"}
    if names:
        kwargs["header"] = None
        kwargs["names"] = names
    df = pd.read_csv(p, **kwargs)
    time_col = df.columns[0]
    df = df.rename(columns={time_col: "time"}).set_index("time")
    df.attrs["source"] = str(p)
    return df


def read_monitors(paths: Iterable[str | Path]) -> dict[str, pd.DataFrame]:
    """Read a collection of monitor files, keyed by file stem."""
    return {Path(p).stem: read_monitor(p) for p in paths}
