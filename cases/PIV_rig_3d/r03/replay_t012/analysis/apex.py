"""Test: do parcels in the cells beside the jet base line up with the failed iterations?"""
import re, sys, importlib.util, collections
import numpy as np
spec = importlib.util.spec_from_file_location("cs", sys.argv[2] + "/.claude/skills/mfix-smoke-validate/check_smoke.py")
cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)
d = sys.argv[1]; dx = 2.4286e-3
def pvd(n): return [(float(a), b) for a, b in re.findall(r'timestep="\s*([\d.E+-]+)".*?file="([^"]+)"', open(f"{d}/{n}").read())]
bg = [f for f in pvd("BACKGROUND.pvd") if f[0] >= 0.1199]; pt = pvd("PARTICLES.pvd")
fails = collections.Counter()
for line in open(f"{d}/run.log", errors="replace"):
    if ":-(" in line:
        m = re.match(r"\s*t=\s*([\d.]+)", line); fails[round(float(m[1]) * 1e4)] += 1   # 0.1 ms resolution
print("t_ms | fails in (t-0.25,t] | parcels y<dx: n, r_mm, v_mps(ids) | apex non-axis cells y<dx: min EP_G, r_mm, |v_g| | axis y=1.2 EP_G, v")
ptc = {}
for t, f in bg:
    tp, fp = min(pt, key=lambda a: abs(a[0] - t))
    P = cs.read_vtk_xml(f"{d}/{fp}"); c = P["coordinates"]; v = np.linalg.norm(P["Velocity"], axis=1); ids = P["Particle_ID"]
    rp = np.hypot(c[:, 0], c[:, 2]); s = np.where(c[:, 1] < dx)[0]
    B = cs.read_vtk_xml(f"{d}/{f}"); x = cs.cell_centers(B); r = np.hypot(x[:, 0], x[:, 2]); y = x[:, 1]; e = B["EP_G"]; V = np.linalg.norm(B["Gas_Velocity"], axis=1)
    low = (y >= 0) & (y < dx) & (r > 0.5 * dx); i = np.where(low)[0]; k = i[e[i].argmin()]; ax = np.where((y >= 0) & (y < dx) & (r <= 0.5 * dx))[0][0]
    tk = round(t * 1e4); nf = sum(fails[q] for q in range(tk - 2, tk + 1))
    desc = " ".join(f"{rp[q]*1e3:.1f}/{c[q,1]*1e3:.1f}/{v[q]:.2f}#{int(ids[q])}" for q in s[:5])
    print(f"{t*1e3:7.2f} | {nf:3d} | {len(s):2d} {desc:44s} | {e[k]:.2f} {r[k]*1e3:4.1f} {V[k]:5.1f} | {e[ax]:.2f} {V[ax]:5.1f}")
print("low-cell list (r_mm, y_mm, EP_G@last frame):", [(round(float(r[q])*1e3,2), round(float(y[q])*1e3,2), round(float(e[q]),2)) for q in np.where((y>=0)&(y<dx))[0]])
