"""Bed/jet state of one rig frame pair: axis profile plus particle summary."""
import sys, importlib.util
import numpy as np
spec = importlib.util.spec_from_file_location("cs", ".claude/skills/mfix-smoke-validate/check_smoke.py")
cs = importlib.util.module_from_spec(spec); spec.loader.exec_module(cs)
d, idx = sys.argv[1], sys.argv[2]
dx = 2.4286e-3
bg = cs.read_vtk_xml(f"{d}/BACKGROUND_{idx}.vtu")
pt = cs.read_vtk_xml(f"{d}/PARTICLES_{idx}.vtp")
xyz = cs.cell_centers(bg)
vg = bg["Gas_Velocity"]; ep = bg["EP_G"]; p = bg["P_G"]
r = np.hypot(xyz[:, 0], xyz[:, 2])
axis = r < 0.5 * dx
o = np.argsort(xyz[axis, 1])
print("axis cells: y_mm  EP_G  v_y  P_g-101325")
for y, e, v, q in zip(xyz[axis, 1][o], ep[axis][o], vg[axis, 1][o], p[axis][o]):
    if y < 0.075:
        print(f"  {y*1e3:6.1f} {e:6.3f} {v:7.2f} {q-101325:8.0f}")
print("P_g extremes:", f"min {p.min()-101325:.0f} at {xyz[p.argmin()]*1e3}", f"max {p.max()-101325:.0f} at {xyz[p.argmax()]*1e3}")
sp = np.linalg.norm(vg, axis=1)
print("max |v_g|", f"{sp.max():.1f} at {xyz[sp.argmax()]*1e3}, EP_G there {ep[sp.argmax()]:.3f}; EP_G min {ep.min():.3f}")
pc = pt["coordinates"]
print("parcels", len(pc), "y range mm", f"{pc[:,1].min()*1e3:.2f} {pc[:,1].max()*1e3:.2f}")
print("keys", list(pt.keys())[:12])
