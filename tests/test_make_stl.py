"""Tests for src.make_stl.build_spouted_bed_stl_implicit."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
import yaml

from src.make_stl import build_spouted_bed_stl_implicit

# ---------------------------------------------------------------------------
# Production geometry values loaded once from params/geometry.yaml.
# All downstream fixtures and tests that touch the implicit-SDF path must use
# these so a change to the YAML immediately propagates to the test suite.
# ---------------------------------------------------------------------------
_REPO_ROOT = Path(__file__).resolve().parent.parent
_GEO = yaml.safe_load((_REPO_ROOT / "params" / "geometry.yaml").read_text())
_V3D = _GEO["spouted_bed_V_3d"]

# Production L_stub_m (R1 value: 0.020 m).  Used by the implicit-STL fixture
# and by test_production_L_stub_minimum to catch regression to the old 0.010 m
# value that crashed the first PIC attempt.
_PROD_L_STUB_M: float = float(_V3D["L_stub_m"])

# Minimum acceptable stub length [m].  Chosen midway between the broken
# pre-R1 value (0.010 m) and the current working value (0.020 m) so that
# drift in either dangerous direction triggers a test failure.
_L_STUB_MIN_M: float = 0.015


def _parse_vertices(text: str) -> np.ndarray:
    pat = re.compile(r"vertex\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)\s+([-+eE.0-9]+)")
    return np.array(
        [[float(a), float(b), float(c)] for a, b, c in pat.findall(text)]
    )


# ---------------------------------------------------------------------------
# Tests for the Day-10 fix 3 implicit-boolean / marching-cubes STL generator.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def implicit_stl_path(tmp_path_factory) -> Path:
    """Build a coarse implicit-SDF STL once per test session (slow).

    L_stub_m is taken from the production params/geometry.yaml
    (spouted_bed_V_3d.L_stub_m) so this fixture automatically tracks any
    future YAML change.  The test-only dx_sample_m is kept at 1 mm (coarser
    than the production 0.5 mm) to keep the test fast; it is not taken from
    the YAML because it only controls STL resolution and does not affect the
    physics that the tests below verify.
    """
    out = tmp_path_factory.mktemp("implicit") / "spouted_implicit.stl"
    build_spouted_bed_stl_implicit(
        R_c_m=0.025,
        R_i_m=0.002,
        cone_half_angle_deg=30.0,
        H_dom_m=0.20,
        L_stub_m=_PROD_L_STUB_M,       # production value from geometry.yaml
        dx_sample_m=0.001,              # 1 mm sampling (coarse to keep the test fast)
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
    # L_stub is the production value loaded from params/geometry.yaml so this
    # assertion automatically tightens or widens with any future YAML change.
    y = verts[:, 1]
    dx_test = 0.001   # test-fixture dx_sample_m (1 mm)
    assert y.min() <= -_PROD_L_STUB_M + 2 * dx_test   # within 2 mm of -L_stub
    assert y.max() >= 0.20 - 2 * dx_test               # within 2 mm of H_dom


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


def test_implicit_no_degenerate_facets(implicit_stl_path):
    """No written facet may have |cross(V2-V1, V3-V1)| below MFiX's safe limit.

    Regression guard for the Day-10 AWS Batch SIGFPE: PyVista's marching cubes
    on an SDF whose zero isosurface touches grid points emits "sliver" facets
    with coincident or collinear vertices.  Those slivers pass MFiX's STL
    reader filter (get_stl_data.f:973-992) when their written normal has unit
    length, then crash IS_POINT_INSIDE_FACET (get_stl_data.f:2386) with
    division by zero at Inv_denom = 1 / ||cross||^2.  The STL writer must
    either drop the sliver entirely or write [0, 0, 0] for the normal so
    MFiX's NORM > TOL_STL = 1e-6 branch catches it.
    """
    text = implicit_stl_path.read_text()
    facet_re = re.compile(
        r"facet normal\s+(\S+)\s+(\S+)\s+(\S+).*?"
        r"vertex\s+(\S+)\s+(\S+)\s+(\S+).*?"
        r"vertex\s+(\S+)\s+(\S+)\s+(\S+).*?"
        r"vertex\s+(\S+)\s+(\S+)\s+(\S+)",
        re.DOTALL,
    )
    survivors = 0
    for m in facet_re.finditer(text):
        f = [float(x) for x in m.groups()]
        normal = np.array(f[:3])
        v1 = np.array(f[3:6])
        v2 = np.array(f[6:9])
        v3 = np.array(f[9:12])
        normal_mag = np.linalg.norm(normal)
        cross_mag = np.linalg.norm(np.cross(v2 - v1, v3 - v1))
        # MFiX reads the written normal; if NORM > TOL_STL=1e-6 the facet
        # is accepted.  For the facet to then be safe in IS_POINT_INSIDE_FACET,
        # the vertex cross product must be well above float precision.
        if normal_mag > 1.0e-6 and cross_mag < 1.0e-12:
            survivors += 1
    assert survivors == 0, (
        f"{survivors} degenerate facets slipped past the STL writer; MFiX "
        f"would SIGFPE on them in cartesian_grid/get_stl_data.f:2386."
    )


def _facet_centroids_normals(text: str) -> tuple[np.ndarray, np.ndarray]:
    verts = _parse_vertices(text).reshape(-1, 3, 3)
    e1 = verts[:, 1] - verts[:, 0]
    e2 = verts[:, 2] - verts[:, 0]
    n = np.cross(e1, e2)
    n /= np.linalg.norm(n, axis=1)[:, None]
    return verts.mean(axis=1), n


def test_implicit_no_membrane_at_orifice_plane(implicit_stl_path):
    """No horizontal facet may span the orifice bore at the stub/cone seam.

    Regression guard for the Day-10 R5 finding (notes/day10_r5_orifice_clog.md):
    the stub and cone SDF primitives were both capped at exactly y = 0, the
    sampling grid lands on y = 0.0, and the union SDF was identically zero over
    the orifice disc, so marching cubes emitted a membrane (144 facets,
    r = 0.24-1.84 mm in the R4 production STL).  MFiX treated it as a wall and
    the inlet gas squeezed through a thin ring gap.
    """
    c, n = _facet_centroids_normals(implicit_stl_path.read_text())
    r = np.hypot(c[:, 0], c[:, 2])
    dx_test = 0.001
    interior = (c[:, 1] > -_PROD_L_STUB_M + 2 * dx_test) & (c[:, 1] < 0.20 - 2 * dx_test)
    horizontal = np.abs(n[:, 1]) > 0.7
    bad = interior & horizontal & (r < 0.002 - dx_test)
    assert bad.sum() == 0, (
        f"{bad.sum()} horizontal facets inside the orifice bore "
        f"(y = {np.round(c[bad, 1] * 1e3, 2)} mm): spurious STL membrane"
    )


def test_implicit_no_membrane_at_cone_top(implicit_stl_path):
    """Same seam guard at the cone/freeboard junction (y = h_cone)."""
    c, n = _facet_centroids_normals(implicit_stl_path.read_text())
    r = np.hypot(c[:, 0], c[:, 2])
    dx_test = 0.001
    interior = (c[:, 1] > -_PROD_L_STUB_M + 2 * dx_test) & (c[:, 1] < 0.20 - 2 * dx_test)
    horizontal = np.abs(n[:, 1]) > 0.7
    bad = interior & horizontal & (r < 0.025 - 2 * dx_test)
    assert bad.sum() == 0, (
        f"{bad.sum()} horizontal facets inside the vessel "
        f"(y = {np.round(c[bad, 1] * 1e3, 2)} mm): spurious STL membrane"
    )


# ---------------------------------------------------------------------------
# Day-10 R2 Option 2: tapered-stub STL
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def tapered_stl_path(tmp_path_factory) -> Path:
    """Build an implicit-SDF STL with Option-2 taper (R_stub -> R_i taper).

    Uses the same R_stub = 4 mm / R_i = 2 mm / L_stub = 20 mm numbers as the
    production geometry, with taper_length_m = 2 mm.  The 1 mm sampling
    (coarser than production 0.5 mm) keeps the test fast.
    """
    out = tmp_path_factory.mktemp("tapered") / "tapered.stl"
    build_spouted_bed_stl_implicit(
        R_c_m=0.025,
        R_i_m=0.002,
        cone_half_angle_deg=30.0,
        H_dom_m=0.20,
        L_stub_m=_PROD_L_STUB_M,
        dx_sample_m=0.001,
        out_path=out,
        R_stub_m=0.004,
        taper_length_m=0.002,
    )
    return out


def test_tapered_stub_radius_profile(tapered_stl_path):
    """Verify the stub radius narrows linearly through the taper.

    At y = -1.5 mm (deep in bulk stub, below taper): max r ~ R_stub = 4 mm.
    At y = -0.5 mm (75% through the 2 mm taper): max r ~ R_stub + 0.75 *
        (R_i - R_stub) = 4 + 0.75*(2-4) = 2.5 mm (tolerance for 1 mm sampling).
    At y = 0: max r ~ R_i = 2 mm.
    """
    verts = _parse_vertices(tapered_stl_path.read_text())
    r = np.sqrt(verts[:, 0] ** 2 + verts[:, 2] ** 2)
    y = verts[:, 1]
    dx = 0.001  # test-fixture sampling grid

    # Bulk stub (y around -1.5 mm, below the 2 mm taper that starts at y=-2 mm):
    bulk_mask = np.abs(y - (-0.0015)) < 0.5 * dx
    if bulk_mask.any():
        # Within one sampling-cell of 4 mm.
        assert r[bulk_mask].max() == pytest.approx(0.004, abs=1.5 * dx)

    # Mid-taper (y around -0.5 mm, 75% through the 2 mm taper): r_bound =
    # R_stub + (R_i - R_stub)*(y + taper)/taper = 4 + (-2)*(1.5)/2 = 2.5 mm.
    mid_mask = np.abs(y - (-0.0005)) < 0.5 * dx
    if mid_mask.any():
        r_max_mid = r[mid_mask].max()
        # Expected ~2.5 mm; chord/sampling error ~1 mm.
        assert r_max_mid == pytest.approx(0.0025, abs=1.5 * dx)
        # And it must be strictly less than the bulk stub radius.
        assert r_max_mid < 0.004

    # Orifice plane (y = 0): r_bound = R_i = 2 mm.
    orifice_mask = np.abs(y) < 0.5 * dx
    if orifice_mask.any():
        assert r[orifice_mask].max() == pytest.approx(0.002, abs=1.5 * dx)


def test_tapered_rejects_taper_longer_than_stub(tmp_path):
    """taper_length_m > L_stub_m must raise (taper would spill out of stub)."""
    with pytest.raises(ValueError):
        build_spouted_bed_stl_implicit(
            R_c_m=0.025,
            R_i_m=0.002,
            cone_half_angle_deg=30.0,
            H_dom_m=0.20,
            L_stub_m=0.020,
            dx_sample_m=0.002,
            out_path=tmp_path / "bad.stl",
            R_stub_m=0.004,
            taper_length_m=0.025,  # > L_stub_m
        )


# ---------------------------------------------------------------------------
# Regression guard: production L_stub_m must stay at or above the minimum
# that allows stable PIC spouting.  The pre-R1 value of 0.010 m crashed the
# first PIC launch (inlet jet too short to develop before hitting the cone).
# This test is completely independent of the implicit-STL fixture so it runs
# even when the slow marching-cubes tests are skipped.
# ---------------------------------------------------------------------------


def test_production_L_stub_minimum():
    """params/geometry.yaml spouted_bed_V_3d.L_stub_m must be >= 0.015 m.

    0.015 m is the midpoint between:
      - 0.010 m: the pre-R1 (broken) value that caused the first PIC run to
        crash (inlet jet collapsed before reaching the cone orifice plane).
      - 0.020 m: the current R1 production value.

    A change to any value below 0.015 m in geometry.yaml now fails this test,
    making the regression immediately visible instead of silently passing the
    suite while the simulation crashes at runtime.
    """
    assert _PROD_L_STUB_M >= _L_STUB_MIN_M, (
        f"params/geometry.yaml spouted_bed_V_3d.L_stub_m = {_PROD_L_STUB_M!r} m "
        f"is below the minimum safe value of {_L_STUB_MIN_M!r} m.  "
        f"The pre-R1 value of 0.010 m caused the first PIC launch to crash; "
        f"restoring a short stub will break the simulation.  "
        f"Update spouted_bed_V_3d.L_stub_m to >= {_L_STUB_MIN_M!r} m."
    )
