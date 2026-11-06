"""dt history of a rig run from its MFiX .TXT file, in time bins."""
import re, sys
import numpy as np
path, step = sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
tmax = float(sys.argv[3]) if len(sys.argv) > 3 else None
tmin = float(sys.argv[4]) if len(sys.argv) > 4 else 0.0
row = re.compile(r"^\s*([\d.]+)\s+([\d.E+-]+)\s+(\S)\s+(\d+)\s+(\S+)\s+([\d.]+)s\s+([\d.]+)s\s+(\d+)\s+(\d+)")
t, dt, cv, nit, pip, res = [], [], [], [], [], []
for line in open(path, errors="replace"):
    m = row.match(line)
    if not m: continue
    t.append(float(m[1])); dt.append(float(m[2])); cv.append(m[3]); nit.append(int(m[4])); res.append(m[5].split("=")[0]); pip.append(int(m[9]))
t, dt, nit, pip = map(np.array, (t, dt, nit, pip)); cv = np.array(cv); res = np.array(res)
print(f"rows {len(t)}  t_end {t[-1]:.4f}  CV flags {dict(zip(*np.unique(cv, return_counts=True)))}")
print("   t0     t1     steps  dt_med    dt_min    dt_max   NIT_med NIT_p90 nonconv  PIP  worst-res")
edges = np.arange(tmin, (tmax or t[-1]) + step, step)
for a, b in zip(edges[:-1], edges[1:]):
    s = (t >= a) & (t < b)
    if not s.any(): continue
    u, c = np.unique(res[s], return_counts=True)
    print(f"{a:6.3f} {b:6.3f} {s.sum():7d} {np.median(dt[s]):9.2e} {dt[s].min():9.2e} {dt[s].max():9.2e} {np.median(nit[s]):6.0f} {np.percentile(nit[s],90):6.0f} {(cv[s]!='C').sum():6d} {pip[s][-1]:6d} {u[c.argmax()]}")
