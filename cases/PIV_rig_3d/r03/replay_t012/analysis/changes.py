"""Largest frame-to-frame field changes in the replay, and where. Also the fastest parcels' histories."""
import re, sys, importlib.util
import numpy as np
spec = importlib.util.spec_from_file_location("cs", sys.argv[2] + "/.claude/skills/mfix-smoke-validate/check_smoke.py")
cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)
d = sys.argv[1]; t0, t1 = float(sys.argv[3]), float(sys.argv[4])
fr = [(float(a), b) for a, b in re.findall(r'timestep="\s*([\d.E+-]+)".*?file="([^"]+)"', open(f"{d}/BACKGROUND.pvd").read())]
fr = [f for f in fr if t0 <= f[0] <= t1]
prev = None
print(" t_ms | max dEP_G  r,y(mm) EP | max dVg  r,y(mm) EP V | max dP(excl. uniform) r,y | median-P-shift")
for t, f in fr:
    b = cs.read_vtk_xml(f"{d}/{f}"); xyz = cs.cell_centers(b); r = np.hypot(xyz[:,0], xyz[:,2]); y = xyz[:,1]
    ep, v, p = b["EP_G"], b["Gas_Velocity"], b["P_G"]
    if prev is not None:
        de = np.abs(ep - prev[0]); dv = np.linalg.norm(v - prev[1], axis=1); dp = p - prev[2]
        m = np.median(dp); dpe = np.abs(dp - m)
        i, j, k = de.argmax(), dv.argmax(), dpe.argmax()
        print(f"{t*1e3:6.2f} | {de[i]:.3f} {r[i]*1e3:5.1f},{y[i]*1e3:5.1f} {ep[i]:.2f} | {dv[j]:6.1f} {r[j]*1e3:5.1f},{y[j]*1e3:5.1f} {ep[j]:.2f} {np.linalg.norm(v[j]):5.1f} | {dpe[k]:6.0f} {r[k]*1e3:5.1f},{y[k]*1e3:5.1f} | {m:7.1f}")
    prev = (ep, v, p)
