"""Load MFiX VTK output with PyVista.

MFiX writes one ``.vtu`` per output interval plus a ``.pvd`` collection file
that pairs each ``.vtu`` with its physical time [s]. ``vtk_data='C'`` keeps the
fields on cell centres, so ``mesh.cell_data`` holds ``EP_G``, ``P_G``,
``Gas_Velocity``, etc.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path

import pyvista as pv


@dataclass(frozen=True)
class Frame:
    """A single VTU snapshot at a known simulation time."""

    time: float  # [s]
    path: Path
    mesh: pv.UnstructuredGrid


def read_vtu(path: str | Path) -> pv.UnstructuredGrid:
    """Read one ``.vtu`` file into a PyVista UnstructuredGrid."""
    return pv.read(str(path))


def read_pvd(pvd_path: str | Path) -> list[tuple[float, Path]]:
    """Parse a ``.pvd`` collection and return ``[(time_s, vtu_path), ...]``.

    Entries are returned in the order they appear in the file.
    """
    p = Path(pvd_path)
    tree = ET.parse(p)
    root = tree.getroot()
    out: list[tuple[float, Path]] = []
    for ds in root.iter("DataSet"):
        t = float(ds.attrib["timestep"])
        rel = ds.attrib["file"]
        out.append((t, (p.parent / rel).resolve()))
    return out


def load_pvd(pvd_path: str | Path) -> list[Frame]:
    """Load every frame in a ``.pvd`` into memory as ``Frame`` records.

    All meshes are read eagerly. For very large runs iterate ``read_pvd`` and
    call ``read_vtu`` per frame instead of using this helper.
    """
    return [Frame(t, p, read_vtu(p)) for t, p in read_pvd(pvd_path)]


def cell_field(mesh: pv.UnstructuredGrid, *candidates: str):
    """Return the first cell-data array whose name matches one of ``candidates``.

    Comparison is case-insensitive. Raises ``KeyError`` if none match.
    """
    lut = {k.lower(): k for k in mesh.cell_data.keys()}
    for c in candidates:
        if c.lower() in lut:
            return mesh.cell_data[lut[c.lower()]]
    raise KeyError(
        f"None of {candidates!r} found in cell_data; available: {list(mesh.cell_data.keys())}"
    )
