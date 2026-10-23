"""Where the point-source gas leaves its cells (C26 grid-layout A/B, 2026-10-20).

Reads a C26 run dir (manifest.json + BACKGROUND VTK).  For each BACKGROUND
frame in the window, from the staggered gas velocities (VTK Gas_Velocity is
the raw U_G, V_G, W_G on each cell's east, north and top faces; MFiX 26.1.2
cartesian_grid/vtk_out.f l.197):

  up      net gas flux out through the north faces of the point-source cells
          (manifest ps_cell_*), / Q = m_dot / rho_g.  1 - up leaves through
          their side faces (mass balance of the source cells; exact in a
          gas-only run, approximate with solids, which displace gas).
  closure signed V_G * A over every first-row north face, / Q (1 when steady).
  cross   net flux across the x = 0 and z = 0 faces inside the 2 x 2 source
          block (axis-on-face layout only), / Q: 0 for a symmetric source.
  centroid  upward-gas-flux centroid (x, z) in rows 2-9 (y = 4.2-18.9 mm),
          over cells whose north face lies wholly inside the cone.

North-face areas are the cell square clipped by the cone circle
r = R_t + y tan(half-angle) at the face height, using straight chords
between the edge intersections, as MFiX builds cut faces.  This is an
approximation of MFiX's STL cut faces; the closure checks it.

    .venv/bin/python notes/c26_source_split.py <run_dir> [t0 t1]
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / '.claude/skills/mfix-smoke-validate'))
from check_smoke import cell_centers, read_pvd, read_vtk_xml  # noqa: E402

T_WINDOW = (0.01, 0.05)   # [s] gas-only tests: steady from 0.005 s (2026-10-20)
CENTROID_ROWS = range(2, 10)


def clipped_square_area(x0, x1, z0, z1, R):
    """Area of [x0,x1] x [z0,z1] inside the circle r < R, with straight chords
    between the points where the square's edges cross the circle."""
    corners = [(x0, z0), (x1, z0), (x1, z1), (x0, z1)]
    poly = []
    for (ax, az), (bx, bz) in zip(corners, corners[1:] + corners[:1]):
        if math.hypot(ax, az) < R:
            poly.append((ax, az))
        dx, dz = bx - ax, bz - az
        a, b, c = dx * dx + dz * dz, 2 * (ax * dx + az * dz), ax * ax + az * az - R * R
        disc = b * b - 4 * a * c
        if disc > 0:
            for s in sorted(((-b - math.sqrt(disc)) / (2 * a), (-b + math.sqrt(disc)) / (2 * a))):
                if 0 < s < 1:
                    poly.append((ax + s * dx, az + s * dz))
    if len(poly) < 3:
        return 0.0
    p = np.array(poly)
    return 0.5 * abs(np.dot(p[:, 0], np.roll(p[:, 1], -1)) - np.dot(p[:, 1], np.roll(p[:, 0], -1)))


def main(run, t0, t1):
    run = Path(run)
    m = json.loads((run / 'manifest.json').read_text())
    dx, tan = m['dx_m'], math.tan(math.radians(m['cone_half_angle_deg']))
    Q = m['ps_massflow_g_kg_s'] / m['rho_g_kg_m3']
    x_min, z_min, y_min = m['x_min_m'], m['z_min_m'], m['y_min_m']
    j_ps = int(round((m['ps_cell_y_s_m'] - y_min) / dx))
    ps_i = range(int(round((m['ps_cell_x_w_m'] - x_min) / dx)), int(round((m['ps_cell_x_e_m'] - x_min) / dx)))
    ps_k = range(int(round((m['ps_cell_z_b_m'] - z_min) / dx)), int(round((m['ps_cell_z_t_m'] - z_min) / dx)))
    face_on_axis = m.get('axis_on_cell_face', False)

    frames = [(t, f) for t, f in read_pvd(run / 'BACKGROUND.pvd') if t0 <= t <= t1]
    if not frames:
        sys.exit(f'no BACKGROUND frames in [{t0}, {t1}] s')
    area_cache = {}

    def north_area(i, j, k):
        if (i, j, k) not in area_cache:
            y_n = y_min + (j + 1) * dx
            R = m['R_t_m'] + max(y_n, 0.0) * tan
            xa, za = x_min + i * dx, z_min + k * dx
            area_cache[(i, j, k)] = clipped_square_area(xa, xa + dx, za, za + dx, R)
        return area_cache[(i, j, k)]

    ups, clos, cx, cz, cents = [], [], [], [], {r: [] for r in CENTROID_ROWS}
    for _, f in frames:
        g = read_vtk_xml(run / f)
        c = cell_centers(g)
        ijk = np.floor((c - [x_min, y_min, z_min]) / dx).astype(int)
        U, V, W = g['Gas_Velocity'].T
        cell = {tuple(t): n for n, t in enumerate(ijk)}
        up = sum(V[cell[(i, j_ps, k)]] * north_area(i, j_ps, k)
                 for i in ps_i for k in ps_k if (i, j_ps, k) in cell)
        row = np.flatnonzero(ijk[:, 1] == j_ps)
        clo = sum(V[n] * north_area(*ijk[n]) for n in row)
        ups.append(up / Q)
        clos.append(clo / Q)
        if face_on_axis:
            # faces at x = 0 (east faces of the west source cells) and z = 0
            iw, kb = ps_i[0], ps_k[0]
            ax = sum(U[cell[(iw, j_ps, k)]] for k in ps_k if (iw, j_ps, k) in cell) * dx * dx
            az = sum(W[cell[(i, j_ps, kb)]] for i in ps_i if (i, j_ps, kb) in cell) * dx * dx
            cx.append(ax / Q)
            cz.append(az / Q)
        for r in CENTROID_ROWS:
            j = j_ps + r
            sel = np.flatnonzero((ijk[:, 1] == j) & (V > 0))
            y_n = y_min + (j + 1) * dx
            R = m['R_t_m'] + y_n * tan
            xs = x_min + (ijk[sel, 0] + 0.5) * dx
            zs = z_min + (ijk[sel, 2] + 0.5) * dx
            inside = np.hypot(np.abs(xs) + dx / 2, np.abs(zs) + dx / 2) < R
            q = V[sel][inside]
            if q.sum() > 0:
                cents[r].append(((q * xs[inside]).sum() / q.sum(), (q * zs[inside]).sum() / q.sum()))

    ups = np.array(ups)
    print(f'{run.name}: {len(frames)} BACKGROUND frames, t = {frames[0][0]:.3f}-{frames[-1][0]:.3f} s; '
          f'PS {len(ps_i)}x{len(ps_k)} cell(s) in row j = {j_ps}; Q = {Q:.4e} m^3/s')
    print(f'  up share (north faces of the source cells) {ups.mean():.3f} '
          f'(frames {ups.min():.3f}-{ups.max():.3f}); sideways {1 - ups.mean():.3f}')
    print(f'  closure (whole first-row top plane)        {np.mean(clos):.3f}')
    if face_on_axis:
        print(f'  net flux across x = 0 / Q {np.mean(cx):+.3f}, across z = 0 / Q {np.mean(cz):+.3f}')
    centroids = {}
    for r in CENTROID_ROWS:
        if cents[r]:
            a = np.array(cents[r]).mean(0) * 1e3
            y_c = m['ps_cell_y_s_m'] + (r + 0.5) * dx
            centroids[y_c] = (a[0] * 1e-3, a[1] * 1e-3)
            print(f'  upward-flux centroid y = {y_c*1e3:5.1f} mm: x {a[0]:+.2f} mm, z {a[1]:+.2f} mm')
    return {'n_frames': len(frames), 'up_share': float(ups.mean()), 'closure': float(np.mean(clos)),
            'cross_x': float(np.mean(cx)) if cx else None, 'cross_z': float(np.mean(cz)) if cz else None,
            'centroids_m': centroids}


if __name__ == '__main__':
    t0, t1 = (float(sys.argv[2]), float(sys.argv[3])) if len(sys.argv) > 3 else T_WINDOW
    main(sys.argv[1], t0, t1)
