"""Digitize the time-averaged velocity-magnitude maps of C26 Fig. 4.

C26 Fig. 4 (Conry et al. 2026, p.6) shows side-view maps of particle
velocity magnitude above the bed at 10 / 12 / 14 slpm: (a) experiment PIV,
(b) PIV of C26's simulation.  Table 3 scores (b) against (a).  This script
reads the figure's JPEG straight from the C26 PDF (checking its SHA-256), maps
each pixel to a velocity through the jet colorbar (nearest colour, linear
scale 4e-4 to 0.4 m/s), and writes each map on the comparison grid
(params/piv_c26.yaml fig4.grid_cols cells across the vessel):

    data/digitized/c26_fig4_{exp,model}_{10,12,14}slpm.csv

Rows are stored top row first, as in the figure.  Each model map is resampled
onto the grid of the experiment panel at the same flow.  It then prints CC,
SSI and MSE of model vs experiment through src/post/piv_c26.py's metrics, next
to C26 Table 3.

    .venv/bin/python -m src.post.digitize_c26_fig4
"""

from __future__ import annotations

import argparse
import hashlib
import io
import re
import sys
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.spatial import cKDTree

from src.make_case_c26 import load_sweep
from src.post.piv_c26 import DIGITIZED, REPO, fig4_path, load_params, map_metrics


def extract_jpeg(pdf_path, width, height, sha256):
    """Raw bytes of the DCT-encoded image of the given size in the PDF."""
    raw = Path(pdf_path).read_bytes()
    for m in re.finditer(rb"<<([^<>]*?/Filter\s*/DCTDecode[^<>]*?)>>\s*stream\r?\n", raw):
        d = m.group(1)
        w, h, n = (re.search(rb"/" + k + rb"\s+(\d+)", d) for k in (b"Width", b"Height", b"Length"))
        if w and h and n and int(w.group(1)) == width and int(h.group(1)) == height:
            data = raw[m.end():m.end() + int(n.group(1))]
            got = hashlib.sha256(data).hexdigest()
            if got != sha256:
                raise ValueError(f"Fig. 4 image hash {got} != expected {sha256}: different PDF?")
            return data
    raise ValueError(f"no {width}x{height} JPEG in {pdf_path}")


def colorbar_lut(img, rows, cols):
    """(colours, values-fraction) of the colorbar, tick marks median-filtered out."""
    lut = img[rows[0]:rows[1] + 1, cols[0]:cols[1] + 1].mean(axis=0)
    return ndimage.median_filter(lut, size=(9, 1), mode="nearest")


def invert_colormap(rgb, lut, values, max_dist):
    """Velocity per pixel from the nearest colorbar colour; NaN off the colormap."""
    dist, i = cKDTree(lut).query(rgb.reshape(-1, 3))
    out = values[i].reshape(rgb.shape[:2])
    out[dist.reshape(rgb.shape[:2]) > max_dist] = np.nan
    return out


def fill_nan(a):
    """Replace NaNs by the nearest valid value."""
    bad = np.isnan(a)
    if not bad.any():
        return a
    idx = ndimage.distance_transform_edt(bad, return_distances=False, return_indices=True)
    return a[tuple(idx)]


def resample(a, shape):
    """Box-filter resample to (rows, cols)."""
    im = Image.fromarray(a.astype(np.float32))   # mode F
    return np.asarray(im.resize((shape[1], shape[0]), Image.Resampling.BOX), dtype=float)


def digitize(p=None):
    p = p or load_params()
    f4 = p["fig4"]
    jpg = extract_jpeg(REPO / f4["pdf"], f4["image_width_px"], f4["image_height_px"], f4["image_sha256"])
    img = np.asarray(Image.open(io.BytesIO(jpg)).convert("RGB"), dtype=float)
    lut = colorbar_lut(img, f4["colorbar_rows"], f4["colorbar_cols"])
    values = np.linspace(f4["v_min_m_s"], f4["v_max_m_s"], len(lut))
    maps, info = {}, {}
    for flow in sorted(f4["panels_exp"]):
        for kind in ("exp", "model"):
            c0, c1, r0, r1 = f4[f"panels_{kind}"][flow]
            v = invert_colormap(img[r0:r1 + 1, c0:c1 + 1], lut, values, f4["max_color_distance"])
            info[kind, flow] = (v.shape, float(np.isnan(v).mean()))
            maps[kind, flow] = fill_nan(v)
        ncol = int(f4["grid_cols"])
        nrow = int(round(ncol * maps["exp", flow].shape[0] / maps["exp", flow].shape[1]))
        for kind in ("exp", "model"):
            maps[kind, flow] = resample(maps[kind, flow], (nrow, ncol))
    return maps, info


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out-dir", type=Path, default=DIGITIZED)
    a = ap.parse_args(argv)
    p = load_params()
    maps, info = digitize(p)
    t3 = load_sweep()["piv_targets_table3"]
    case_of_flow = {10: "c01", 12: "c02", 14: "c03"}
    a.out_dir.mkdir(parents=True, exist_ok=True)
    for (kind, flow), v in sorted(maps.items()):
        src_shape, nan_frac = info[kind, flow]
        label = "experiment PIV (Fig. 4a)" if kind == "exp" else "C26 model PIV (Fig. 4b)"
        path = a.out_dir / fig4_path(kind, flow).name
        hdr = (f"C26 Fig. 4 {label}, {flow} slpm: time-averaged particle velocity magnitude [m/s]\n"
               f"digitized by src/post/digitize_c26_fig4.py (params/piv_c26.yaml fig4) from the "
               f"{src_shape[1]}x{src_shape[0]} px panel; {100 * nan_frac:.2f} % off-colormap pixels filled\n"
               f"grid {v.shape[1]} cols (left wall to right wall) x {v.shape[0]} rows (top row first, "
               f"bottom row at the masked bed surface); colorbar linear 4e-4..0.4 m/s")
        np.savetxt(path, v, delimiter=",", fmt="%.4f", header=hdr)
        print(f"wrote {path.relative_to(REPO)}  {v.shape}  mean {v.mean():.4f} m/s  "
              f"off-colormap {100 * nan_frac:.2f} %")
    print("\nC26 model vs experiment on the digitized grid (Table 3 in brackets):")
    for flow, cid in case_of_flow.items():
        met = map_metrics(maps["exp", flow], maps["model", flow], p["metrics"])
        print(f"  {flow} slpm: " + "  ".join(f"{k} {met[k]:.4f} [{t3[cid][k]}]" for k in met))
    return 0


if __name__ == "__main__":
    sys.exit(main())
