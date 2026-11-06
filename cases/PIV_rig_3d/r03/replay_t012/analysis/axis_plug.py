"""Axis EP_G by frame against failed iterations: does the stall line up with the plug clearing?"""
import re, sys, importlib.util, collections
import numpy as np
spec = importlib.util.spec_from_file_location("cs", sys.argv[2] + "/.claude/skills/mfix-smoke-validate/check_smoke.py")
cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)
d = sys.argv[1]; dx = 2.4286e-3
bg = [(float(a), b) for a, b in re.findall(r'timestep="\s*([\d.E+-]+)".*?file="([^"]+)"', open(f"{d}/BACKGROUND.pvd").read())]
bg = [f for f in bg if f[0] >= 0.1199]
fails = collections.Counter()
for line in open(f"{d}/run.log", errors="replace"):
    if ":-(" in line:
        fails[round(float(re.match(r"\s*t=\s*([\d.]+)", line)[1]) * 1e4)] += 1
ys = [8.5, 10.9, 13.4, 15.8, 18.2, 20.6, 23.1, 25.5, 27.9, 30.4]
print("  t_ms fails | axis EP_G at y=" + " ".join(f"{y:5.1f}" for y in ys) + " | peak P, y | total open-area jet flow: v_y*EP at y=25.5")
for t, f in bg:
    b = cs.read_vtk_xml(f"{d}/{f}"); x = cs.cell_centers(b); r = np.hypot(x[:,0], x[:,2]); y = x[:,1]; e = b["EP_G"]; p = b["P_G"] - 101325
    ax = np.where(r < 0.5*dx)[0]; o = ax[np.argsort(y[ax])]
    vals = [e[o[np.argmin(np.abs(y[o]*1e3 - q))]] for q in ys]
    src = o[y[o] < 0]; rest = np.setdiff1d(o, src); k = rest[p[rest].argmax()]
    q = o[np.argmin(np.abs(y[o]*1e3 - 25.5))]
    tk = round(t*1e4); nf = sum(fails[w] for w in range(tk-2, tk+1))
    flag = "*" if nf else " "
    print(f"{t*1e3:7.2f} {nf:3d}{flag} |                " + " ".join(f"{v:5.2f}" for v in vals) + f" | {p[k]:6.0f} {y[k]*1e3:5.1f} | {b['Gas_Velocity'][q,1]:6.1f}")
