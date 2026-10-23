"""CGP drag check: which diameter does MFiX use for coarse-grained parcel drag,
and how much would the other choice change the bed?

Part 1 searches the installed MFiX source (the solver named in a C26 case's
run.sh) for the lines that set CGP drag, parcel size and volume, and contact.
Part 2 evaluates the TGS law (src/correlations.py) with the case deck's gas
properties and statistical weight, at the real and at the parcel diameter.

The MFiX 23.2 source (C26's version) is not available here, so Part 1 checks
26.1.2 only; the notes ("CGP drag check") weigh the 23.2 evidence.

Usage:
    .venv/bin/python notes/c26_cgp_drag_check.py [--case c01_massonly_axisface]
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path

import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from src.correlations import (  # noqa: E402
    G_STD, terminal_velocity, tgs_drag_coefficient, umf_tgs,
)

CASES = REPO / 'cases' / 'C26_3d'

# (file under model/, pattern, must be present, what it means)
SOURCE_CHECKS = [
    ('des/drag_gp_des.f', r'IF\(CGDEM\)\s*DPM\s*=\s*2\.0d0\*DES_CGP_RPR\(NP\)', True,
     'drag law diameter = real particle diameter'),
    ('des/drag_gp_des.f', r'CALL DRAG_TGS\(DgA,EPgtmp,Mu,ROPg,VREL,DPM,DPA,phis\)', True,
     'TGS called with that diameter (DPM; DPA = DPM for one phase)'),
    ('des/drag_gp_des.f', r'^\s*F_gstmp\s*=\s*DgA\s*\*\s*PVOL\(NP\)\s*$', True,
     'drag coefficient x parcel volume (model A)'),
    ('des/make_arrays_des.f', r'^\s*PVOL\(L\)\s*=\s*\(4\.0D0/3\.0D0\)\*PI\*DES_RADIUS\(L\)\*\*3\s*$', True,
     'parcel volume from the parcel radius; PMASS = PVOL * rho'),
    ('des/generate_particles_mod.f',
     r'DES_CGP_RPR\(PIP\)\s*=\s*DES_RADIUS\(PIP\)/DES_CGP_STW\(PIP\)\*\*\(1\.0d0/3\.0d0\)', True,
     'real radius = parcel radius / W^(1/3) (IC seeding, constant W)'),
    ('check_data/check_solids_common_discrete.f',
     r'D_p0\(M\)\s*=\s*D_p0\(M\)\s*\*\s*CGP_STAT_WT\(M\)\*\*\(1\.d0/3\.d0\)', True,
     'parcel diameter = d W^(1/3) (C26 Eq. 1)'),
    ('drag_gs.f', r'RE = DPA\*VREL\*ROPg/Mug', True, 'TGS Re (mirrored in correlations.py)'),
    ('drag_gs.f', r'F_STOKES = 18D0\*Mug\*EPg/DPM\*\*2', True, 'TGS Stokes scale (mirrored)'),
]
CONTACT_FILES = ['des/calc_force_dem.f', 'des/calc_collision_wall_mod.f',
                 'des/cfwallcontact.f', 'des/cffctowall.f']
CGP_TOKEN = r'^[^!]*(\bCGDEM\b|DES_CGP_|CGP_STAT_WT)'   # code, not comments


def mfix_src_for(case_dir):
    """MFiX source tree of the solver that run.sh launches."""
    m = re.search(r'"([^"]+)/bin/mfixsolver\w*"', (case_dir / 'run.sh').read_text())
    if not m:
        sys.exit(f'no mfixsolver path in {case_dir / "run.sh"}')
    src = Path(m.group(1)) / 'share' / 'mfix' / 'src'
    if not (src / 'model').is_dir():
        sys.exit(f'no MFiX source tree at {src}')
    return src


def find(path, pattern):
    rx = re.compile(pattern, re.I)
    return [(i, ln.rstrip()) for i, ln in enumerate(path.read_text(errors='replace').splitlines(), 1)
            if rx.search(ln)]


def deck_value(text, key):
    m = re.search(rf'^\s*{re.escape(key)}\s*=\s*(\S+)', text, re.M | re.I)
    if not m:
        sys.exit(f'{key} not in deck')
    return m.group(1).strip("'\"")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--case', default='c01_massonly_axisface')
    a = ap.parse_args(argv)
    case_dir = CASES / a.case
    src = mfix_src_for(case_dir)
    model = src / 'model'
    version = (src / 'VERSION').read_text().strip() if (src / 'VERSION').exists() else '?'

    print(f'== Part 1: MFiX {version} source ({src})')
    ok = True
    for rel, pat, want, meaning in SOURCE_CHECKS:
        hits = find(model / rel, pat)
        st = 'PASS' if bool(hits) == want else 'FAIL'
        ok &= st == 'PASS'
        where = ', '.join(f'{rel}:{i}' for i, _ in hits) or rel
        print(f'  [{st}] {meaning}\n         {where}')
    for rel in CONTACT_FILES:
        hits = find(model / rel, CGP_TOKEN)
        st = 'PASS' if not hits else 'FAIL'
        ok &= st == 'PASS'
        print(f'  [{st}] no CGP term in contact routine {rel}'
              + ''.join(f'\n         {rel}:{i}: {t.strip()}' for i, t in hits))

    deck_path = next(case_dir.glob('*.mfx'))
    deck = deck_path.read_text()
    manifest = json.loads((case_dir / 'manifest.json').read_text())
    ysz = yaml.safe_load((REPO / 'params' / 'particles.yaml').read_text())['case_C26_ysz']
    rho_g = float(deck_value(deck, 'ro_g0'))
    mu_g = float(deck_value(deck, 'mu_g0'))
    w = float(deck_value(deck, 'cgp_stat_wt(1)'))
    d_real = float(deck_value(deck, 'd_p0(1)'))
    rho_p = float(deck_value(deck, 'ro_s0(1)'))
    clip = float(deck_value(deck, 'des_epg_clip'))
    mdot = float(deck_value(deck, 'ps_massflow_g(1)'))
    drag = deck_value(deck, 'drag_type')
    if abs(d_real - ysz['d_p_mean_m']) > 1e-9 or rho_p != ysz['rho_p_kg_per_m3']:
        sys.exit('deck particle size/density differ from params/particles.yaml case_C26_ysz')
    d_parcel = d_real * w ** (1.0 / 3.0)
    r_top = min(manifest['R_t_m'] + manifest['H_static_m']
                * math.tan(math.radians(manifest['cone_half_angle_deg'])), manifest['R_c_m'])
    u_top = mdot / rho_g / (math.pi * r_top ** 2)

    print(f'\n== Part 2: TGS law, deck {deck_path.name} (drag_type {drag})')
    print(f'  gas rho {rho_g:.4f} kg/m3, mu {mu_g:.4e} Pa s; rho_p {rho_p:.0f} kg/m3; W {w:g}')
    print(f'  d real {d_real * 1e6:.1f} um, d parcel {d_parcel * 1e6:.1f} um '
          f'(ratio W^(1/3) = {d_parcel / d_real:.3f})')
    print(f'  {"eps":>6} {"U_mf real":>10} {"U_mf parcel":>12} {"ratio":>6} '
          f'{"DgA real/parcel at U_mf,real":>29}')
    eps_list = sorted(set([*ysz['eps_mf_range'], clip]))
    for eps in eps_list:
        u_r = umf_tgs(d_real, rho_p, rho_g, mu_g, eps)
        u_p = umf_tgs(d_parcel, rho_p, rho_g, mu_g, eps)
        dga = (tgs_drag_coefficient(eps, u_r / eps, d_real, rho_g, mu_g)
               / tgs_drag_coefficient(eps, u_r / eps, d_parcel, rho_g, mu_g))
        print(f'  {eps:6.2f} {u_r:8.3f} m/s {u_p:10.3f} m/s {u_p / u_r:6.2f} {dga:29.2f}')
    ut_r = terminal_velocity(d_real, rho_p, rho_g, mu_g, G_STD)
    ut_p = terminal_velocity(d_parcel, rho_p, rho_g, mu_g, G_STD)
    print(f'  single-particle terminal velocity: real {ut_r:.2f} m/s, parcel {ut_p:.2f} m/s '
          f'(ratio {ut_p / ut_r:.2f})')
    print(f'  bed-top superficial gas velocity at {a.case}: {u_top:.3f} m/s '
          f'(r_top {r_top * 1e3:.1f} mm, ps_massflow_g {mdot:.4e} kg/s)')
    print(f'\nOVERALL source checks: {"PASS" if ok else "FAIL"}')
    return 0 if ok else 1


if __name__ == '__main__':
    sys.exit(main())
