#!/usr/bin/env bash
# Capture a stalled state of r03 and A/B the pressure-solve sweep cap on it.
# 1) baseline numerics from t = 0.12 s with res_dt = 0.1 ms; snapshot each completed .RES pair; stop the
#    container at the first failed iteration.  2) continue the last snapshot for 0.8 ms with leq_it(1) = 20 and 200.
set -uo pipefail
R="$(cd "$(dirname "$0")/.." && pwd)"; S="$R/replay_cap_state"; mkdir -p "$S/snap"
img=(docker run --rm --platform=linux/amd64 -e LOCAL_SMOKE=1 -e MFX_FILE="RIG_r03.mfx" -e NP=4 -e NODESI=2 -e NODESJ=1 -e NODESK=2)
(cd "$S" && "${img[@]}" --name r03_cap_state -v "$PWD":/work triso-cfd/mfix:23.2 > /dev/null 2> docker.err) &
last=""
while sleep 1; do
  docker ps --format '{{.Names}}' | /usr/bin/grep -q '^r03_cap_state$' || { [ -n "$last" ] && break; [ -f "$S/run.log" ] && ! pgrep -f "name r03_cap_state" > /dev/null && break; continue; }
  if [ -f "$S/run.log" ] && /usr/bin/grep -q ':-(' "$S/run.log"; then docker kill r03_cap_state > /dev/null; echo "first failure seen; container stopped"; break; fi
  m=$(stat -f %m "$S/RIG_R03.RES" 2>/dev/null || echo 0); now=$(date +%s)
  if [ "$m" != "$last" ] && [ $((now - m)) -ge 1 ] && /usr/bin/grep -q "Time =" "$S/run.log" 2>/dev/null; then
    cp "$S/RIG_R03.RES" "$S/snap/RIG_R03.RES.tmp" && cp "$S/RIG_R03_DES.RES" "$S/snap/RIG_R03_DES.RES.tmp" \
      && mv "$S/snap/RIG_R03.RES.tmp" "$S/snap/RIG_R03.RES" && mv "$S/snap/RIG_R03_DES.RES.tmp" "$S/snap/RIG_R03_DES.RES" \
      && /usr/bin/grep "Time =" "$S/run.log" | tail -1 > "$S/snap/at.txt"; last="$m"
  fi
done
wait
echo "capture: fails $(/usr/bin/grep -c ':-(' "$S/run.log"); last step $(/usr/bin/grep 'Time =' "$S/run.log" | tail -1); snapshot near $(cat "$S/snap/at.txt" 2>/dev/null)"
/usr/bin/grep -q ':-(' "$S/run.log" || { echo "no stall captured"; exit 0; }
t=$(/usr/bin/grep ':-(' "$S/run.log" | head -1 | sed -E 's/^[[:space:]]*t=[[:space:]]*([0-9.]+).*/\1/')
tstop=$(awk -v t="$t" 'BEGIN{printf "%.4f", t + 0.0008}')
for arm in cap20 cap200; do
  d="$R/replay_cap_$arm"; mkdir "$d"
  for f in "$S"/*; do case "$f" in *BACKGROUND_*|*PARTICLES_*|*run.log|*.LOG|*.OUT|*.pvd|*boundary.vtk|*FACETS*|*MESH.vtu|*VTU_FRAME*|*chain.sh|*docker.err|*/snap) ;; *) cp "$f" "$d/";; esac; done
  cp "$S/snap/RIG_R03.RES" "$S/snap/RIG_R03_DES.RES" "$d/"
  sed -i '' -E "s/^([[:space:]]*tstop[[:space:]]*=[[:space:]]*).*$/\1$tstop/" "$d/RIG_r03.mfx"
  sed -i '' -E 's/^([[:space:]]*res_dt[[:space:]]*=[[:space:]]*).*$/\10.01/' "$d/RIG_r03.mfx"
  [ "$arm" = cap200 ] && sed -i '' -E 's/^([[:space:]]*full_log = \.True\.)$/\1\'$'\n''  leq_it(1) = 200/' "$d/RIG_r03.mfx"
  /usr/bin/grep -n -E "tstop|leq_it|res_dt" "$d/RIG_r03.mfx"
  (cd "$d" && "${img[@]}" -v "$PWD":/work triso-cfd/mfix:23.2 > /dev/null 2> docker.err; echo "$arm exit $?")
done
