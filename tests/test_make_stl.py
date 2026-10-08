"""Tests for src.make_stl.build_spouted_bed_stl."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest

from src.make_stl import build_spouted_bed_stl, build_spouted_bed_stl_implicit


@pytest.fixture()
def stl_path(tmp_path) -> Path:
    out = tmp_path / "test_bed.stl"
    build_spouted_bed_stl(
        R_c_m=0.038,
        R_i_m=0.00475,
        cone_half_angle_deg=30.0,
        H_dom_m=0.40,
        theta_segments=60,
        out_path=out,
    )
    return out


def _parse_vertices(text: str) -> np.ndarray:
    pat = re.compile(r"vertex\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)")
    return np.array(
        [[float(a), float(b), float(c)] for a, b, c in pat.findall(text)]
    )


def test_facet_count_matches_profile(stl_path):
    # Profile has 2 segments (cone + cylinder); each ring of theta_segments=60
    # quads produces 2 triangles -> 2 * 60 * 2 = 240 side-wall triangles.
    text = stl_path.read_text()
    nfacets = len(re.findall(r"\bfacet normal\b", text))
    assert nfacets == 240


def test_bounding_box(stl_path):
    verts = _parse_vertices(stl_path.read_text())
    xmin, ymin, zmin = verts.min(axis=0)
    xmax, ymax, zmax = verts.max(axis=0)
    # x,z should span [-R_c, +R_c] exactly on the cylinder course
    assert xmax == pytest.approx(0.038, rel=1e-6)
    assert xmin == pytest.approx(-0.038, rel=1e-6)
    assert zmax == pytest.approx(0.038, rel=1e-6)
    assert zmin == pytest.approx(-0.038, rel=1e-6)
    # axial span 0 .. H_dom
    assert ymin == pytest.approx(0.0, abs=1e-12)
    assert ymax == pytest.approx(0.40, rel=1e-6)


def test_orifice_rim_radius(stl_path):
    """At y=0 the inner rim of the cone must sit at r = R_i (orifice edge)."""
    verts = _parse_vertices(stl_path.read_text())
    r = np.sqrt(verts[:, 0] ** 2 + verts[:, 2] ** 2)
    at_floor = np.abs(verts[:, 1]) < 1e-9
    r_floor = r[at_floor]
    # Allow tiny tessellation chord error around the circle.
    assert np.min(r_floor) == pytest.approx(0.00475, rel=1e-3)


def test_cone_slope_matches_half_angle(stl_path):
    """Points between y=0 and y=cone_top should lie on r = R_i + y*tan(30 deg)."""
    verts = _parse_vertices(stl_path.read_text())
    r = np.sqrt(verts[:, 0] ** 2 + verts[:, 2] ** 2)
    y = verts[:, 1]
    tan30 = np.tan(np.radians(30.0))
    cone_top = (0.038 - 0.00475) / tan30
    on_cone = (y > 1e-9) & (y < cone_top - 1e-6)
    if np.any(on_cone):
        r_predicted = 0.00475 + y[on_cone] * tan30
        # Tessellation chord ~ R(y) * (1 - cos(pi/60)) ~ 0.14 % of R(y);
        # allow 1% rel to cover both ends.
        assert np.max(np.abs(r[on_cone] - r_predicted) / r_predicted) < 0.02


def test_cylinder_radius_above_cone(stl_path):
    """Points above the cone top should sit on r = R_c (within tessellation)."""
    verts = _parse_vertices(stl_path.read_text())
    r = np.sqrt(verts[:, 0] ** 2 + verts[:, 2] ** 2)
    y = verts[:, 1]
    tan30 = np.tan(np.radians(30.0))
    cone_top = (0.038 - 0.00475) / tan30
    on_cyl = y > cone_top + 1e-6
    # Chord error on the cylinder: R_c * (1 - cos(pi/60)) ~ 5.2e-5 m.
    assert np.all(np.abs(r[on_cyl] - 0.038) < 1e-4)


def test_header_and_footer(stl_path):
    text = stl_path.read_text()
    assert text.startswith("solid spouted_bed_V_3d\n")
    assert text.rstrip().endswith("endsolid spouted_bed_V_3d")


def test_invalid_cone_rejected(tmp_path):
    with pytest.raises(ValueError):
        build_spouted_bed_stl(
            R_c_m=0.038,
            R_i_m=0.00475,
            cone_half_angle_deg=1.0,  # cone too narrow -> h_cone > H_dom
            H_dom_m=0.05,
            theta_segments=20,
            out_path=tmp_path / "bad.stl",
        )


# ---------------------------------------------------------------------------
# Tests for the Day-10 fix 3 implicit-boolean / marching-cubes STL generator.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def implicit_stl_path(tmp_path_factory) -> Path:
    """Build a coarse implicit-SDF STL once per test session (slow)."""
    out = tmp_path_factory.mktemp("implicit") / "spouted_implicit.stl"
    build_spouted_bed_stl_implicit(
        R_c_m=0.025,
        R_i_m=0.002,
        cone_half_angle_deg=30.0,
        H_dom_m=0.20,
        L_stub_m=0.010,
        dx_sample_m=0.001,   # 1 mm sampling (coarse to keep the test fast)
        out_path=out,
    )
    return out


def test_implicit_header_and_footer(implicit_stl_path):
    text = implicit_stl_path.read_text()
    assert text.startswith("solid spouted_bed_V_3d_implicit\n")
    assert text.rstrip().endswith("endsolid spouted_bed_V_3d_implicit")


def test_implicit_facet_count_lower_bound(implicit_stl_path):
    text = implicit_stl_path.read_text()
    nfacets = len(re.findall(r"\bfacet normal\b", text))
    # Loose lower bound; at 1 mm sampling we expect ~5k-20k triangles.
    assert nfacets >= 1000, f"expected >=1000 facets, got {nfacets}"


def test_implicit_bounding_box(implicit_stl_path):
    verts = _parse_vertices(implicit_stl_path.read_text())
    # y-span should cover [-L_stub, H_dom] with small overrun due to the
    # 2*dx margin and marching-cubes grid alignment.
    y = verts[:, 1]
    assert y.min() <= -0.010 + 0.002   # within 2 mm of -L_stub
    assert y.max() >= 0.20 - 0.002     # within 2 mm of H_dom


def test_implicit_no_nan_inf(implicit_stl_path):
    verts = _parse_vertices(implicit_stl_path.read_text())
    assert np.all(np.isfinite(verts))


def test_implicit_orifice_rim_radius(implicit_stl_path):
    """At y=0 the minimum vertex radius should be close to R_i."""
    verts = _parse_vertices(implicit_stl_path.read_text())
    r = np.sqrt(verts[:, 0] ** 2 + verts[:, 2] ** 2)
    # Find vertices near y=0 (the orifice plane).
    near_floor = np.abs(verts[:, 1]) < 1.5e-3  # within 1.5*dx of y=0
    r_near = r[near_floor]
    assert r_near.size > 0
    # Within 2 * dx_sample = 2 mm of R_i = 2 mm.
    assert np.min(r_near) == pytest.approx(0.002, abs=2e-3)
