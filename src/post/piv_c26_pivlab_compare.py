"""Compare a PIVlab analysis of the exported C26 frames with the in-repo
engine (native) on the same frames.

Input: the export dir written by src.post.piv_c26_export_frames (export.json,
reference_<engine>.json, reference_map_<engine>.csv) and a directory of
PIVlab 'File > Export > ASCII' files, one per frame pair, in frame order
(sorted by name).  The PIVlab session must be UNCALIBRATED: x, y in pixels
(1-based, y down) and u, v in px/frame; the first four numeric columns are
read, extra columns (vorticity, ...) are ignored, header lines skipped, NaN
allowed.  The conversion to m/s uses pixel_m and frame_dt_s from export.json.

For --view wall: per pair, the mean of +v (PIVlab y down = C26's downward
positive) over the vectors inside the 4 mm box, then the time mean; compared
with the engine's wall_piv_m_s.
For --view spout: per pair, the vectors as exported (run PIVlab's
post-processing with 'interpolate missing data' first, so every grid point
carries a vector, as the project's empty_windows = interpolate); time mean of
the vectors over the pairs with any finite vector; |mean vector| resampled
onto the digitized Fig. 4 grid; CC / SSI / MSE against the C26 experiment
map and Pearson CC against the engine's map.

Self-test (checks the parser and the coordinate conventions, must reproduce
the engine to ~1e-4): point --pivlab at the exporter's own
pivlab_layout_native/<view>/.

    .venv/bin/python -m src.post.piv_c26_pivlab_compare \\
        results/pivlab_export/c02_momentum24 --view wall \\
        --pivlab ~/PIVlab_out/c02_wall/

All values SI.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from src.post import piv_c26 as P

NUM = re.compile(r"[-+]?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?|nan|NaN|Inf|-Inf")


def read_pivlab_ascii(path):
    """(x, y, u, v) 1-D arrays from one PIVlab ASCII export; the first four
    numeric fields of each data row, any of comma / semicolon / tab / space
    separators, non-numeric header lines skipped."""
    rows = []
    with Path(path).open() as f:
        for line in f:
            parts = [t for t in re.split(r"[,;\t ]+", line.strip()) if t]
            if len(parts) < 4:
                continue
            try:
                vals = [float(t) for t in parts[:4]]
            except ValueError:
                continue
            rows.append(vals)
    if not rows:
        raise ValueError(f"{path}: no numeric rows")
    a = np.array(rows)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def to_grid(x, y, u, v):
    """Scatter rows -> (xs, ys, U, V) on the regular PIVlab grid (ys increasing
    = downward in the image)."""
    xs, ys = np.unique(x), np.unique(y)
    U = np.full((ys.size, xs.size), np.nan)
    V = np.full((ys.size, xs.size), np.nan)
    i = np.searchsorted(ys, y)
    j = np.searchsorted(xs, x)
    U[i, j], V[i, j] = u, v
    return xs, ys, U, V


def natural_key(p):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", p.name)]


def pivlab_files(d):
    files = sorted([f for f in Path(d).iterdir() if f.suffix.lower() in (".txt", ".dat", ".csv", ".asc")],
                   key=natural_key)
    if not files:
        raise ValueError(f"no PIVlab ASCII files in {d}")
    return files


def load_series(files):
    """Stack of per-pair fields on a common grid: xs, ys, U (n, ny, nx), V."""
    first = to_grid(*read_pivlab_ascii(files[0]))
    xs, ys = first[0], first[1]
    U, V = [first[2]], [first[3]]
    for f in files[1:]:
        gx, gy, u, v = to_grid(*read_pivlab_ascii(f))
        if gx.shape != xs.shape or gy.shape != ys.shape or not (np.allclose(gx, xs) and np.allclose(gy, ys)):
            raise ValueError(f"{f}: grid differs from {files[0]}")
        U.append(u); V.append(v)
    return xs, ys, np.array(U), np.array(V)


def compare_wall(meta, xs, ys, U, V, refs):
    box = meta["views"]["wall"]["pivlab_box_px"]
    px, dt = meta["pixel_m"], meta["frame_dt_s"]
    inb = ((xs[None, :] >= box["x_min"]) & (xs[None, :] <= box["x_max"])
           & (ys[:, None] >= box["y_min"]) & (ys[:, None] <= box["y_max"]))
    if not inb.any():
        raise ValueError("no PIVlab vectors inside the wall box; check the view / units")
    per_pair = np.array([np.nanmean(V[k][inb]) for k in range(V.shape[0])]) * px / dt   # downward +
    out = {"view": "wall", "n_pairs_pivlab": int(V.shape[0]), "n_vectors_in_box": int(inb.sum()),
           "n_vectors_in_box_project": meta["views"]["wall"]["n_grid_points_in_box"],
           "wall_pivlab_m_s": float(np.nanmean(per_pair)), "wall_pivlab_sd_m_s": float(np.nanstd(per_pair)),
           "nan_frac_in_box": float(np.isnan(V[:, inb]).mean())}
    for e, r in refs.items():
        out[f"wall_{e}_m_s"] = r["wall_piv_m_s"]
        out[f"wall_rel_diff_vs_{e}"] = out["wall_pivlab_m_s"] / r["wall_piv_m_s"] - 1.0
        n = min(len(per_pair), len(r["wall_piv_per_pair_m_s"]))
        out[f"wall_per_pair_cc_vs_{e}"] = P.pearson_cc(np.nan_to_num(per_pair[:n]),
                                                       np.array(r["wall_piv_per_pair_m_s"][:n]))
    return out, per_pair


def compare_spout(meta, xs, ys, U, V, refs, params):
    px, dt = meta["pixel_m"], meta["frame_dt_s"]
    view, mp = meta["views"]["spout"], meta["map"]
    ny = view["shape_px"][0]
    any_vec = np.isfinite(U).reshape(U.shape[0], -1).any(axis=1) & np.isfinite(V).reshape(V.shape[0], -1).any(axis=1)
    Um, Vm = np.nanmean(U[any_vec], axis=0), np.nanmean(V[any_vec], axis=0)
    # PIVlab grid -> project metres; rows must increase upward for the resampler
    rows = ny - ys                      # project rows, decreasing along ys
    x_m = view["x0_m"] + (xs - 1.0 + 0.5) * px
    y_m = view["y0_m"] + (rows + 0.5) * px
    order = np.argsort(y_m)
    mag = np.hypot(Um, -Vm)[order] * px / dt                 # |mean vector| [m/s]
    mag = np.where(np.isfinite(mag), mag, 0.0)
    grid = P.fig4_grid_map(mag, SimpleNamespace(x=x_m, y=y_m[order]), mp["y_map_m"], mp["H_m"], mp["R_c_m"],
                           mp["fig4_grid"][0], mp["fig4_grid"][1])
    out = {"view": "spout", "n_pairs_pivlab": int(U.shape[0]), "n_pairs_with_vectors": int(any_vec.sum()),
           "grid_pivlab": f"{ys.size}x{xs.size}", "grid_project": f"{len(view['grid_rows_px'])}x{len(view['grid_cols_px'])}",
           "nan_frac_pivlab": float(np.isnan(V).mean()),
           "map_mean_pivlab_m_s": float(grid.mean()), "map_max_pivlab_m_s": float(grid.max())}
    exp_path = P.fig4_path("exp", mp["panel_flow_slpm"])
    if exp_path.exists() and float(mp["panel_flow_slpm"]) == float(meta["flow_slpm"]):
        exp = P.load_fig4("exp", mp["panel_flow_slpm"])
        m = P.map_metrics(exp, grid, params["metrics"])
        out.update({f"{k}_pivlab_vs_exp": v for k, v in m.items()})
    for e, r in refs.items():
        out[f"map_mean_{e}_m_s"] = r["map_mean_m_s"]
        rmap = np.flipud(np.loadtxt(r["_map_path"], delimiter=",", comments="#"))
        out[f"map_cc_pivlab_vs_{e}"] = P.pearson_cc(grid, rmap)
        out[f"map_rms_diff_vs_{e}_m_s"] = float(np.sqrt(np.mean((grid - rmap) ** 2)))
        if r.get("metrics_vs_exp"):
            for k, v in r["metrics_vs_exp"].items():
                out[f"{k}_{e}_vs_exp"] = v
    return out, grid


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("export_dir", type=Path, help="results/pivlab_export/<run>")
    ap.add_argument("--pivlab", type=Path, required=True, help="dir of PIVlab ASCII exports, one per pair")
    ap.add_argument("--view", choices=("wall", "spout"), required=True)
    ap.add_argument("--tag", default="pivlab", help="output name suffix")
    a = ap.parse_args(argv)
    params = P.load_params()
    meta = json.loads((a.export_dir / "export.json").read_text())
    refs = {}
    for e in meta["engines"]:
        f = a.export_dir / f"reference_{e}.json"
        if f.exists():
            refs[e] = json.loads(f.read_text())
            refs[e]["_map_path"] = a.export_dir / f"reference_map_{e}.csv"
    files = pivlab_files(a.pivlab)
    xs, ys, U, V = load_series(files)
    if U.shape[0] != meta["n_pairs"]:
        print(f"warning: {U.shape[0]} PIVlab pairs vs {meta['n_pairs']} exported pairs", file=sys.stderr)
    if a.view == "wall":
        out, series = compare_wall(meta, xs, ys, U, V, refs)
    else:
        out, series = compare_spout(meta, xs, ys, U, V, refs, params)
        np.savetxt(a.export_dir / f"map_{a.tag}_{a.view}.csv", np.flipud(series), delimiter=",", fmt="%.5f",
                   header=f"{a.tag} spout map [m/s] on the Fig. 4 grid, top row first, from {a.pivlab}")
    out["pivlab_dir"] = str(a.pivlab)
    path = a.export_dir / f"compare_{a.tag}_{a.view}.csv"
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out))
        w.writeheader(); w.writerow(out)
    for k, v in out.items():
        print(f"  {k}: {v:.4g}" if isinstance(v, float) else f"  {k}: {v}")
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
