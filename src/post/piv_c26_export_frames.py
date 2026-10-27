"""Export the C26 synthetic-PIV camera frames as image files for PIVlab, with
the reference results of the in-repo engines on the same frames.

The synthetic PIV (src/post/piv_c26.py) renders its camera frames in memory.
To check the PIV step with an independent program (PIVlab, which C26 used),
this writes exactly those frames to disk:

    <out-dir>/<run>/spout/frame_NNNN.tif   spout-map view (C26 Fig. 4 region)
    <out-dir>/<run>/wall/frame_NNNN.tif    wall-box view (C26 Fig. 1c box)
    <out-dir>/<run>/export.json            pixel size, frame times, frame
                                           interval, view geometry, the
                                           project's PIV grid and wall box
                                           in PIVlab image coordinates
    <out-dir>/<run>/reference_native.json / reference_map_native.csv
                                           wall velocity, visible-parcel wall
                                           velocity and spout map of the
                                           in-repo engine over the SAME frames
    <out-dir>/<run>/pivlab_layout_native/{spout,wall}/pair_NNNN.txt
                                           the engine's per-pair fields in the
                                           PIVlab ASCII layout (x, y, u, v in
                                           px and px/frame, y down), to
                                           self-test src.post.piv_c26_pivlab_compare
    <out-dir>/shift_check/                 with --shift-check: one synthetic
                                           pair of spheres moved by a known
                                           shift, to check PIVlab's sign and
                                           scale before the run frames

Images: 16-bit grayscale TIFF, intensity 0-1 -> 0-65535, top row first
(PIVlab / MATLAB image convention: x right, y DOWN, 1-based pixel indices).
The project renders row 0 at the bottom, so a project (row, col) grid point
is PIVlab pixel (x, y) = (col + 1, ny - row); project dr (up) = -v_PIVlab,
dc = u_PIVlab.  Compare PIVlab's export with the references using
src.post.piv_c26_pivlab_compare.  Raw frames are gitignored.

    MPLBACKEND=Agg .venv/bin/python -m src.post.piv_c26_export_frames \\
        cases/C26_3d/c02_apex_v232/aws_out [--max-frames 61] [--stride 1] \\
        [--shift-check]

All values SI; image quantities in pixels [px].
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np
import tifffile

from src.post import piv_c26 as P

OUT = P.REPO / "results" / "pivlab_export"
VIEWS = ("spout", "wall")


def to_tiff(img, path):
    """Project image (row 0 at the bottom, intensity 0-1) -> 16-bit TIFF, top
    row first."""
    a = np.clip(np.flipud(np.asarray(img, float)), 0.0, 1.0)
    tifffile.imwrite(path, np.rint(a * 65535).astype(np.uint16), photometric="minisblack")


def pivlab_xy(rows, cols, ny):
    """Project grid (rows, cols) [px, row 0 at the bottom, pixel centre at the
    integer] -> PIVlab image coordinates (x, y) [px, 1-based, y down]."""
    return np.asarray(cols, float) + 1.0, ny - np.asarray(rows, float)


def write_pivlab_ascii(path, rows, cols, dr, dc, ny, empty=None):
    """Write a field in PIVlab's 'Export > ASCII' layout: header lines, then
    one row per vector 'x, y, u, v' [px, px, px/frame, px/frame], x right,
    y down; empty windows (if given) are written as NaN, as PIVlab's
    rejected vectors are before interpolation."""
    x, y = pivlab_xy(rows, cols, ny)
    u, v = np.asarray(dc, float), -np.asarray(dr, float)
    if empty is not None:
        u, v = np.where(empty, np.nan, u), np.where(empty, np.nan, v)
    data = np.column_stack([x.ravel(), y.ravel(), u.ravel(), v.ravel()])
    with Path(path).open("w") as f:
        f.write("# PIVlab-layout export of the project engine (src.post.piv_c26_export_frames)\n")
        f.write("x [px], y [px], u [px/frame], v [px/frame]\n")
        np.savetxt(f, data, delimiter=", ", fmt="%.4f")


def view_record(view, step, wins, in_roi=None, box=None):
    ny, nx = view.shape
    x_pl, y_pl = pivlab_xy(view.rows, view.cols, ny)
    rec = {"shape_px": [int(ny), int(nx)], "x0_m": view.x0, "y0_m": view.y0,
           "grid_rows_px": [int(r) for r in view.rows[:, 0]], "grid_cols_px": [int(c) for c in view.cols[0, :]],
           "pivlab_grid_x_px": [float(v) for v in x_pl[0, :]], "pivlab_grid_y_px": [float(v) for v in y_pl[:, 0]],
           "grid_x_m": [float(v) for v in view.x], "grid_y_m": [float(v) for v in view.y],
           "windows_px": wins, "step_px": step}
    if box is not None:
        xc, yc, half, px = box
        # box in PIVlab pixel coordinates (x right, y down, 1-based)
        c0, c1 = (np.array([xc - half, xc + half]) - view.x0) / px
        r0, r1 = (np.array([yc - half, yc + half]) - view.y0) / px
        rec["box_m"] = {"xc": xc, "yc": yc, "half": half}
        rec["pivlab_box_px"] = {"x_min": float(c0 + 0.5), "x_max": float(c1 + 0.5),
                                "y_min": float(ny - r1 + 0.5), "y_max": float(ny - r0 + 0.5)}
        rec["n_grid_points_in_box"] = int(in_roi.sum())
    return rec


def export_run(run, out, stride=1, max_frames=None, write_layout=True, params=None, pixel_m=None,
               name_suffix=""):
    p = params or P.load_params()
    if pixel_m is not None:                      # camera resolution check (notes/c26_pivlab_rerun.md)
        p = copy.deepcopy(p)
        p["camera"]["pixel_m"] = float(pixel_m)
    run = Path(run)
    name = (run.parent.name if run.name == "aws_out" else run.name) + name_suffix
    out_run = Path(out) / name
    for v in VIEWS:
        (out_run / v).mkdir(parents=True, exist_ok=True)
    g = P.run_geometry(run, p, stride, max_frames)
    frames = g["frames"]
    px, dt = g["px"], float(np.median(np.diff([t for t, _ in frames])))
    written = set()
    refs = {}

    for engine in ("native",):
        layout = out_run / f"pivlab_layout_{engine}"
        if write_layout:
            for v in VIEWS:
                (layout / v).mkdir(parents=True, exist_ok=True)

        def on_pair(info, layout=layout):
            i = info["i"]
            for v in VIEWS:
                img_a, img_b, dr, dc, emp = info[v]
                view = g[v]
                for j, img in ((i, img_a), (i + 1, img_b)):
                    if (v, j) not in written:
                        to_tiff(img, out_run / v / f"frame_{j:04d}.tif")
                        written.add((v, j))
                if write_layout:
                    write_pivlab_ascii(layout / v / f"pair_{i:04d}.txt", view.rows, view.cols, dr, dc,
                                       view.shape[0], empty=None)

        res = P.analyse_run(run, p, stride, max_frames, on_pair=on_pair)
        row = P.summary_row(res)
        ref = {"engine": engine, "n_pairs": res["n_pairs"],
               "wall_piv_m_s": float(np.mean(res["wall_piv_m_s"])),
               "wall_piv_sd_m_s": float(np.std(res["wall_piv_m_s"])),
               "wall_piv_per_pair_m_s": [float(v) for v in res["wall_piv_m_s"]],
               "wall_visible_parcels_m_s": row["wall_visible_parcels_m_s"],
               "map_mean_m_s": row["map_mean_m_s"], "map_max_m_s": row["map_max_m_s"],
               "empty_window_frac": row["empty_window_frac"],
               "metrics_vs_exp": res.get("metrics"), "metrics_c26_model_vs_exp": res.get("metrics_c26_model"),
               "caption": P.caption(res)}
        (out_run / f"reference_{engine}.json").write_text(json.dumps(ref, indent=1))
        np.savetxt(out_run / f"reference_map_{engine}.csv", np.flipud(res["map"]), delimiter=",", fmt="%.5f",
                   header=f"{engine} engine spout map [m/s] on the Fig. 4 grid, top row first. " + P.caption(res))
        np.savetxt(out_run / f"reference_mappiv_{engine}.csv", np.flipud(res["map_piv_grid"]), delimiter=",",
                   fmt="%.5f", header=f"{engine} engine spout map [m/s] on the project PIV grid, top row first")
        refs[engine] = ref
        print(f"{engine:8s} wall {ref['wall_piv_m_s']:.4e} m/s, map mean {ref['map_mean_m_s']:.4f} m/s, "
              f"metrics vs exp {ref['metrics_vs_exp']}")

    meta = {
        "run_dir": str(run), "run": name, "case": str(g["m"]["case_id"]), "flow_slpm": g["flow"],
        "pixel_m": px, "frame_dt_s": dt, "n_frames": len(frames), "n_pairs": len(frames) - 1,
        "frame_times_s": [float(t) for t, _ in frames], "frame_files": [str(f) for _, f in frames],
        "stride": stride, "max_frames": max_frames,
        "image": {"format": "uint16 TIFF, intensity 0-1 -> 0-65535, top row first (y down)",
                  "convention": "PIVlab pixel (x, y) = (col + 1, ny - row); dr_project = -v_pivlab, dc = u_pivlab"},
        "piv": {k: p["piv"][k] for k in ("windows_px", "step_px", "uod_threshold", "uod_eps_px", "min_window_std")},
        "empty_windows": p["spout_map"]["empty_windows"],
        "map": {"y_map_m": g["y_map"], "bed_surface_m": g["y_s"], "H_m": g["H"], "R_c_m": g["R_c"],
                "panel_flow_slpm": g["pflow"], "fig4_grid": [int(g["nrow"]), int(g["ncol"])]},
        "views": {"spout": view_record(g["spout"], g["step"], g["wins"]),
                  "wall": view_record(g["wall"], g["step"], g["wins"], g["in_roi"],
                                      (g["xc"], g["yc"], g["half"], px))},
        "engines": list(refs),
    }
    (out_run / "export.json").write_text(json.dumps(meta, indent=1))
    print(f"wrote {len(written)} frames ({len(frames)} per view), export.json and references to {out_run}")
    return out_run, meta, refs


def export_shift_check(out, shift_px=(1.5, -2.5), n=600, size=200, seed=0):
    """A pair of rendered sphere images moved by shift_px = (dx, dy) in
    project axes (x right, y up) [px], as in tests/test_piv_c26.py, written
    as TIFFs with the expected PIVlab result."""
    px = 1.0e-4
    rng = np.random.default_rng(seed)
    xyz = np.column_stack([rng.uniform(0, size * px, n), rng.uniform(0, size * px, n), rng.uniform(0, 0.01, n)])
    d = np.full(n, 0.6e-3)
    a, _ = P.render(xyz, d, (0.0, 0.0), (size, size), px)
    b, _ = P.render(xyz + np.array([shift_px[0] * px, shift_px[1] * px, 0.0]), d, (0.0, 0.0), (size, size), px)
    out = Path(out) / "shift_check"
    out.mkdir(parents=True, exist_ok=True)
    to_tiff(a, out / "frame_0000.tif")
    to_tiff(b, out / "frame_0001.tif")
    meta = {"shift_project_px": {"dx_right": shift_px[0], "dy_up": shift_px[1]},
            "expected_pivlab_px_per_frame": {"u": shift_px[0], "v": -shift_px[1]},
            "note": "PIVlab's y axis points down: v must come out as -dy_up. Mean over all vectors "
                    "away from the border should match within ~0.05 px."}
    (out / "export.json").write_text(json.dumps(meta, indent=1))
    print(f"shift check written to {out}: expected PIVlab u = {shift_px[0]} px, v = {-shift_px[1]} px")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run", type=Path, nargs="?", help="run dir with manifest.json and PARTICLES_PIV.pvd")
    ap.add_argument("--out-dir", type=Path, default=OUT)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--max-frames", type=int, default=None)
    ap.add_argument("--no-layout", action="store_true", help="skip the PIVlab-layout ASCII of the engine")
    ap.add_argument("--shift-check", action="store_true", help="also write the known-shift pair")
    ap.add_argument("--pixel-m", type=float, default=None, help="camera pixel size [m] (default params camera.pixel_m)")
    ap.add_argument("--name-suffix", default="", help="appended to the export dir name, e.g. _px50")
    a = ap.parse_args(argv)
    if a.shift_check:
        export_shift_check(a.out_dir)
    if a.run is not None:
        export_run(a.run, a.out_dir, a.stride, a.max_frames, not a.no_layout, pixel_m=a.pixel_m, name_suffix=a.name_suffix)
    elif not a.shift_check:
        ap.error("give a run dir and/or --shift-check")
    return 0


if __name__ == "__main__":
    sys.exit(main())
