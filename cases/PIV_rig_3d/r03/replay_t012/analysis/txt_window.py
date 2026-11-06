"""Replay .TXT rows only (after the last time reset): residual kinds and levels by 0.5 ms."""
import re, sys, collections
import numpy as np
row = re.compile(r"^\s*([\d.]+)\s+([\d.E+-]+)\s+(\S)\s+(\d+)\s+([A-Za-z]+\d*)=([\d.E+-]+)")
rows = []
for line in open(sys.argv[1], errors="replace"):
    m = row.match(line)
    if m: rows.append((float(m[1]), float(m[2]), int(m[4]), m[5], float(m[6])))
t = np.array([r[0] for r in rows]); cut = np.where(np.diff(t) < -1e-6)[0]
start = cut[-1] + 1 if len(cut) else 0
rows = rows[start:]
print("replay rows", len(rows), "t", rows[0][0], "->", rows[-1][0])
bins = collections.defaultdict(list)
for r in rows: bins[int(r[0] / 0.0005)].append(r)
print("  t_ms  steps  dt_med   NITmed NITmax  worst-residual share  resid_med")
for b in sorted(bins):
    rr = bins[b]; kinds = collections.Counter(r[3] for r in rr)
    print(f"{b*0.5:7.1f} {len(rr):5d} {np.median([r[1] for r in rr]):9.2e} {np.median([r[2] for r in rr]):5.0f} {max(r[2] for r in rr):5d}  {dict(kinds)}  {np.median([r[4] for r in rr]):.2e}")
