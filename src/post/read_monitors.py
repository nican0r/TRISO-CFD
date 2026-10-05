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


def _read_header(path: Path) -> tuple[int, list[str]]:
    """Return (n_skip, column_names).

    The last ``#`` line before the first data row is treated as the column
    header. If no ``#`` line looks like a header, columns are left to pandas
    (numeric 0, 1, ...).
    """
    names: list[str] = []
    n_skip = 0
    with path.open("r") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                n_skip += 1
                continue
            if stripped.startswith(_HEADER_PREFIXES):
                tokens = stripped.lstrip("#! ").split()
                if tokens and not tokens[0].replace(".", "", 1).replace("-", "", 1).replace("e", "", 1).replace("E", "", 1).isdigit():
                    names = tokens
                n_skip += 1
                continue
            break
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
    n_skip, names = _read_header(p)
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
