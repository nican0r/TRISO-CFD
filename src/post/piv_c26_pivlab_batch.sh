#!/usr/bin/env bash
# PIVlab re-run of the C26 synthetic PIV for one exported run: every variant
# over both views, then the comparison with the in-repo engines and C26.
#
#   src/post/piv_c26_pivlab_batch.sh <run>            e.g. c02_apex_v232
#   src/post/piv_c26_pivlab_batch.sh <run> <variant>  one variant only
#
# Needs results/pivlab_export/<run>/ from src.post.piv_c26_export_frames and
# MATLAB R2026b with PIVlab 3.14 in ~/Documents/MATLAB/PIVlab (notes/
# c26_pivlab_rerun.md).  Scores: <export>/<variant>/score.json
# (src/post/piv_c26_pivlab_score.m).  Variants (src/post/piv_c26_pivlab_run.m):
#   pivlab_default : C26 / CES25 settings + PIVlab GUI defaults (CLAHE on)
#   pivlab_noclahe : as default, CLAHE off (diagnostic: PIVlab returns a fixed
#                    spurious vector for the constant background CLAHE makes
#                    of the particle-free regions of the rendered frames)
#   pivlab_ces25   : CES25 binarization (dots, > 0.9, 1 px dilation) + defaults
#   pivlab_ces25_noclahe : CES25 binarization, CLAHE off
#   pivlab_c26img  : the frames as C26 Fig. 4d shows them going into PIVlab:
#                    8-bit white parcel-sized disks on black (CES25
#                    binarization of the render, no dilation) + defaults
#   pivlab_c26img_noclahe : same frames, CLAHE off
set -euo pipefail
REPO="$(cd "$(dirname "$0")/../.." && pwd)"
RUN="${1:?run name, e.g. c02_apex_v232}"
ONLY="${2:-}"
MATLAB="${MATLAB:-/Applications/MATLAB_R2026b.app/bin/matlab}"
EXP="$REPO/results/pivlab_export/$RUN"
[ -f "$EXP/export.json" ] || { echo "no export at $EXP" >&2; exit 1; }

declare -A OPTS=(
  [pivlab_default]=""
  [pivlab_noclahe]=", clahe=0"
  [pivlab_ces25]=", binarize=1"
  [pivlab_ces25_noclahe]=", binarize=1, clahe=0"
  [pivlab_c26img]=", binarize=1, binarize_thresh=0, binarize_dilate_px=0, image_class='uint8'"
  [pivlab_c26img_noclahe]=", binarize=1, binarize_thresh=0, binarize_dilate_px=0, image_class='uint8', clahe=0"
)
VARIANTS=(pivlab_default pivlab_noclahe pivlab_ces25 pivlab_ces25_noclahe pivlab_c26img pivlab_c26img_noclahe)
[ -n "$ONLY" ] && VARIANTS=("$ONLY")

CMD="addpath('$REPO/src/post');"
for v in "${VARIANTS[@]}"; do
  for view in wall spout; do
    CMD+=" piv_c26_pivlab_run('$EXP','$view','$EXP/$v/$view'${OPTS[$v]}, quiet=true);"
    CMD+=" disp('$RUN $v $view done');"
  done
done
echo "== MATLAB: $RUN ${VARIANTS[*]} ($(date))"
"$MATLAB" -batch "$CMD" 2>&1 | grep -v "Trial License\|^\s*$\|^Try \|recommendations\|^WARNING: .*interpass"

cd "$REPO"
for v in "${VARIANTS[@]}"; do
  for view in wall spout; do
    echo "== compare: $RUN $v $view"
    MPLBACKEND=Agg .venv/bin/python -m src.post.piv_c26_pivlab_compare "$EXP" --view "$view" \
        --pivlab "$EXP/$v/$view" --tag "$v" | grep -E "wall_pivlab_m_s|wall_rel_diff|per_pair_cc|nan_frac|map_mean_pivlab|CC_pivlab|SSI_pivlab|MSE_m2_s2_pivlab|map_cc_pivlab_vs|^wrote"
  done
done
# C26-style scoring in MATLAB with MATLAB defaults (the reference numbers;
# the Python compare above is only a cross-check against the repo engines)
CMD="addpath('$REPO/src/post');"
for v in "${VARIANTS[@]}"; do CMD+=" piv_c26_pivlab_score('$EXP','$v');"; done
echo "== MATLAB scoring: $RUN"
"$MATLAB" -batch "$CMD" 2>&1 | grep -v "Trial License\|^\s*$"
echo "== done $RUN ($(date))"
