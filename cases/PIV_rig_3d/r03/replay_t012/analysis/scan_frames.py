"""Per-frame scan of the r03 replay: where the pressure/velocity extremes sit."""
import re, sys, importlib.util
import numpy as np
spec = importlib.util.spec_from_file_location("cs", sys.argv[2] + "/.claude/skills/mfix-smoke-validate/check_smoke.py")
cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)
d = sys.argv[1]; t0, t1 = float(sys.argv[3]), float(sys.argv[4])
dx = 2.4286e-3; P0 = 101325.0
def frames(pvd):
    return [(float(a), b) for a, b in re.findall(r'timestep="\s*([\d.E+-]+)".*?file="([^"]+)"', open(f"{d}/{pvd}").read())]
bgf, ptf = frames("BACKGROUND.pvd"), frames("PARTICLES.pvd")
print("  t_ms | Pmin(Pa) r,y(mm) EPG | Pmax(Pa) r,y(mm) EPG | maxV(no src) r,y EPG | n_clip | slab_P  | vp_max r,y(mm) | axis v@18 EP@18")
for t, f in bgf:
    if not (t0 <= t <= t1): continue
    b = cs.read_vtk_xml(f"{d}/{f}")
    xyz = cs.cell_centers(b); r = np.hypot(xyz[:, 0], xyz[:, 2]); y = xyz[:, 1]
    p = b["P_G"] - P0; ep = b["EP_G"]; v = np.linalg.norm(b["Gas_Velocity"], axis=1)
    src = (r < 0.5*dx) & (y < 0)
    q = np.where(~src)[0]
    i, j, k = q[p[q].argmin()], q[p[q].argmax()], q[v[q].argmax()]
    slab = (y >= 33.24e-3) & (y < 33.24e-3 + dx)
    ax18 = (r < 0.5*dx) & (y > 16.6e-3) & (y < 19.05e-3)
    tp, fp = min(ptf, key=lambda a: abs(a[0] - t))
    pt = cs.read_vtk_xml(f"{d}/{fp}"); vp = np.linalg.norm(pt["Velocity"], axis=1); ip = vp.argmax()
    c = pt["coordinates"]; rp = np.hypot(c[ip, 0], c[ip, 2])
    print(f"{t*1e3:7.2f} | {p[i]:7.0f} {r[i]*1e3:5.1f},{y[i]*1e3:5.1f} {ep[i]:.2f} | {p[j]:7.0f} {r[j]*1e3:5.1f},{y[j]*1e3:5.1f} {ep[j]:.2f} | {v[k]:6.1f} {r[k]*1e3:5.1f},{y[k]*1e3:5.1f} {ep[k]:.2f} | {(ep<=0.4201).sum():5d} | {p[slab].mean():7.1f} | {vp[ip]:6.2f} {rp*1e3:5.1f},{c[ip,1]*1e3:6.1f} | {v[ax18].mean():6.1f} {ep[ax18].mean():.2f}")
