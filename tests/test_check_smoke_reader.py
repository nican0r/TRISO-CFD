"""Tests for the raw-appended VTK XML reader (.claude/skills/mfix-smoke-validate/check_smoke.py).

MFiX 26.1.2 writes header_type="UInt64" (8-byte block sizes); 23.2 omits
header_type, so the blocks carry the VTK default UInt32 (4-byte) sizes.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / ".claude/skills/mfix-smoke-validate"))
from check_smoke import read_vtk_xml  # noqa: E402

XYZ = np.array([[0.0, 0.01, 0.0], [0.001, 0.02, -0.002]], ">f4")
DIA = np.array([8.9e-4, 8.8e-4], ">f4")


def _vtp(tmp_path, header_attr, size_dtype):
    blocks = [XYZ.tobytes(), DIA.tobytes()]
    offsets, off = [], 0
    for b in blocks:
        offsets.append(off)
        off += np.dtype(size_dtype).itemsize + len(b)
    head = (
        f'<VTKFile type="PolyData" version="0.1" byte_order="BigEndian"{header_attr}>\n'
        f'<PolyData><Piece NumberOfPoints="  2 ">\n'
        f'<Points><DataArray type="Float32" Name="coordinates" NumberOfComponents="3" '
        f'format="appended" offset="  {offsets[0]} " /></Points>\n'
        f'<PointData><DataArray type="Float32" Name="Diameter" format="appended" '
        f'offset="  {offsets[1]} " /></PointData>\n'
        f'</Piece></PolyData>\n<AppendedData encoding="raw">\n_'
    ).encode()
    body = b"".join(np.array([len(b)], size_dtype).tobytes() + b for b in blocks)
    path = tmp_path / "p.vtp"
    path.write_bytes(head + body + b"\n</AppendedData></VTKFile>\n")
    return path


@pytest.mark.parametrize("header_attr, size_dtype", [
    (' header_type="UInt64"', ">u8"),   # MFiX 26.1.2
    ("", ">u4"),                         # MFiX 23.2: VTK default UInt32
    (' header_type="UInt32"', ">u4"),
])
def test_reads_both_block_size_headers(tmp_path, header_attr, size_dtype):
    out = read_vtk_xml(_vtp(tmp_path, header_attr, size_dtype))
    np.testing.assert_allclose(out["coordinates"], XYZ.astype(float))
    np.testing.assert_allclose(out["Diameter"], DIA.astype(float))
