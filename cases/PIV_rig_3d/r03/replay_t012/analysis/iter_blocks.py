"""Per-step iteration blocks from a full_log run.log: dt, NIT, first/last residuals, failed or not."""
import re, sys
import numpy as np
head = re.compile(r"^\s*Time =\s*([\d.]+)\s+Dt =\s*([\d.E+-]+)")
row = re.compile(r"^\s*(\d+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+([\d.E+-]+)\s+(\S+)")
steps, cur = [], None
for line in open(sys.argv[1], errors="replace"):
    m = head.match(line)
    if m:
        cur = dict(t=float(m[1]), dt=float(m[2]), p=[], u=[], v=[], fail=False); steps.append(cur); continue
    if cur is None: continue
    m = row.match(line)
    if m and int(m[1]) == len(cur["p"]) + 1:
        cur["p"].append(float(m[2])); cur["u"].append(float(m[3])); cur["v"].append(float(m[4])); continue
    if ":-(" in line: cur["fail"] = True
t0, t1 = float(sys.argv[2]), float(sys.argv[3])
S = [s for s in steps if t0 <= s["t"] <= t1 and s["p"]]
print(f"steps {len(S)} failed {sum(s['fail'] for s in S)}")
mode = sys.argv[4] if len(sys.argv) > 4 else "bins"
if mode == "bins":
    print(" dt bin        steps fail  NIT med/p90/max   P0_1 med   P0_1 p90   share NIT>=10")
    edges = [0, 1.5e-7, 3e-7, 6e-7, 1.2e-6, 2.4e-6, 4.8e-6, 7e-6, 1.1e-5]
    for a, b in zip(edges[:-1], edges[1:]):
        q = [s for s in S if a <= s["dt"] < b]
        if not q: continue
        n = np.array([len(s["p"]) for s in q]); p1 = np.array([s["p"][0] for s in q])
        print(f" {a:8.1e}-{b:8.1e} {len(q):5d} {sum(s['fail'] for s in q):4d}  {np.median(n):4.0f} {np.percentile(n,90):4.0f} {n.max():4d}   {np.median(p1):8.1e}  {np.percentile(p1,90):8.1e}   {(n>=10).mean():.2f}")
else:
    for s in S:
        if s["fail"] or len(s["p"]) >= int(mode):
            print(f"t={s['t']:.5f} dt={s['dt']:.2e} nit={len(s['p']):2d} {'FAIL' if s['fail'] else 'ok  '} P0:", " ".join(f"{x:.0e}" for x in s["p"][:25]))
