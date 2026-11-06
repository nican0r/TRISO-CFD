#!/usr/bin/env bash
# A/B on r03's stalled state: run the baseline deck from t = 0.12 s to 0.1294 s (inside the stall), then
# continue that state to 0.1300 s with leq_it(1) at the default 20 and at 200.
set -uo pipefail
R="$(cd "$(dirname "$0")/.." && pwd)"
run() { (cd "$1" && date +%s > .t_start && docker run --rm --platform=linux/amd64 -e LOCAL_SMOKE=1 -e MFX_FILE="RIG_r03.mfx" -e NP=4 -e NODESI=2 -e NODESJ=1 -e NODESK=2 -v "$PWD":/work triso-cfd/mfix:23.2 > /dev/null 2> docker.err; echo "$1 exit $?"; date +%s > .t_end); }
run "$R/replay_ab_state"
echo "state run fails: $(/usr/bin/grep -c ':-(' "$R/replay_ab_state/run.log")"
for arm in cap20 cap200; do
  d="$R/replay_ab_$arm"; mkdir "$d"
  for f in "$R"/replay_ab_state/*; do case "$f" in *BACKGROUND_*|*PARTICLES_*|*run.log|*.LOG|*.OUT|*.pvd|*boundary.vtk|*FACETS*|*MESH.vtu|*VTU_FRAME*|*chain.sh|*docker.err) ;; *) cp "$f" "$d/";; esac; done
  sed -i '' -E 's/^([[:space:]]*tstop[[:space:]]*=[[:space:]]*).*$/\10.1300/' "$d/RIG_r03.mfx"
  [ "$arm" = cap200 ] && sed -i '' -E 's/^([[:space:]]*full_log = \.True\.)$/\1\'$'\n''  leq_it(1) = 200/' "$d/RIG_r03.mfx"
  /usr/bin/grep -n -E "tstop|leq_it|run_type|full_log" "$d/RIG_r03.mfx"
  run "$d"
done
