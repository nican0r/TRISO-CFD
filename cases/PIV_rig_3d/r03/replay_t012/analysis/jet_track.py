"""Jet-tip history along the axis: tip height (contiguous EP_G > 0.8 from the source), peak pressure, plug above the tip."""
import re, sys, importlib.util
import numpy as np
spec = importlib.util.spec_from_file_location("cs", ".claude/skills/mfix-smoke-validate/check_smoke.py")
cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)
d, t0, t1, every = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
dx = 2.4286e-3; P0 = 101325.0
fr = [(float(a), b) for a, b in re.findall(r'timestep="\s*([\d.E+-]+)".*?file="([^"]+)"', open(f"{d}/BACKGROUND.pvd").read())]
fr = [f for f in fr if t0 <= f[0] <= t1][::every]
print("  t_ms | tip_y mm | v_tip m/s | Pmax Pa @y | plug above tip: EP min, cells<0.6 | v_axis@6mm | n_clip")
for t, f in fr:
    b = cs.read_vtk_xml(f"{d}/{f}"); xyz = cs.cell_centers(b); r = np.hypot(xyz[:, 0], xyz[:, 2]); y = xyz[:, 1]
    ax = np.where((r < 0.5 * dx) & (y > 0) & (y < 0.08))[0]; ax = ax[np.argsort(y[ax])]
    ep = b["EP_G"][ax]; vy = b["Gas_Velocity"][ax, 1]; p = b["P_G"][ax] - P0
    n = 0
    while n < len(ax) and ep[n] > 0.8: n += 1
    tip = y[ax][n - 1] * 1e3 if n else 0.0
    above = ep[n:n + 8]
    print(f"{t*1e3:7.2f} | {tip:6.1f} | {vy[n-1] if n else 0:6.1f} | {p.max():6.0f} @{y[ax][p.argmax()]*1e3:5.1f} | {above.min() if len(above) else 1:.2f} {int((above < 0.6).sum()):2d} | {vy[2]:6.1f} | {(b['EP_G'] <= 0.4201).sum():4d}")
