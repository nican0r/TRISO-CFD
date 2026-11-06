import re, sys
p = sys.argv[1]
s = open(p, newline="").read()
def sub(pat, rep):
    global s
    s, n = re.subn(pat, rep, s, count=1, flags=re.M | re.I)
    if n != 1: raise SystemExit(f"no match: {pat}")
sub(r"^(\s*run_type\s*=\s*)'NEW'", r"\1'RESTART_1'")
sub(r"^(\s*tstop\s*=\s*).*$", r"\g<1>0.135")
sub(r"^(\s*res_dt\s*=\s*.*)$", "\\1\n  full_log = .True.")
sub(r"^(\s*vtk_dt\(1\)\s*=\s*).*$", r"\g<1>0.00025")
sub(r"^(\s*vtk_dt\(2\)\s*=\s*).*$", r"\g<1>0.00025")
for i in range(1, 5):
    sub(rf"^(\s*monitor_dt\({i}\)\s*=\s*).*$", r"\g<1>0.0002")
open(p, "w", newline="").write(s)
