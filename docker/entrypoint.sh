#!/usr/bin/env bash
# AWS Batch entrypoint for MFiX 26.1.2.
#
# Required env vars (set via the Batch job's containerOverrides):
#   CASE_S3        s3://bucket/cases/<case>/        (prefix with case files)
#   RESULTS_S3     s3://bucket/results/<case>/      (prefix to upload outputs to)
#   MFX_FILE       <name>.mfx                       (deck filename inside case dir)
#   NP             <int>                            (MPI rank count; default 1)
#
# Optional:
#   CASE_S3_ARRAY  s3://bucket/sweep/<prefix>/      (array mode; takes precedence
#                                                   over CASE_S3. The entrypoint
#                                                   reads AWS_BATCH_JOB_ARRAY_INDEX
#                                                   and syncs
#                                                   <CASE_S3_ARRAY>case_<NN>/)
#   RESULTS_S3_ARRAY  s3://bucket/sweep/<prefix>/results/
#   EXTRA_MFIX_ARGS   extra argv passed to mfixsolver (e.g. "-r")
#   NODESI/NODESJ/NODESK  explicit DMP partitioning (must multiply to NP).
#                         If all are unset, entrypoint auto-picks: NODESJ=NP,
#                         NODESI=NODESK=1 (axial decomposition — biased toward
#                         the typical tallest axis in spouted-bed geometries);
#                         otherwise any unset one defaults to 1.
#
# Local smoke-test mode: if LOCAL_SMOKE=1 is set, the entrypoint skips all S3
# interaction and simply runs the solver on /work (which should be bind-mounted
# from the host at `docker run -v`).

set -euo pipefail

log() { printf '[entrypoint %s] %s\n' "$(date -u +%H:%M:%S)" "$*"; }

# -- DMP decomposition --------------------------------------------------------

inject_dmp_keywords() {
    # mfixsolver_dmp (pymfix.py wrapper) does not accept command-line keyword
    # overrides. Instead, patch NODESI/NODESJ/NODESK directly into the deck
    # (idempotent: remove any existing nodesi/j/k lines first).
    local mfx=$1 np=$2
    if [[ "$np" -le 1 ]]; then return 0; fi
    local ni=${NODESI:-}
    local nj=${NODESJ:-}
    local nk=${NODESK:-}
    if [[ -z "$ni" && -z "$nj" && -z "$nk" ]]; then
        ni=1; nj=$np; nk=1
    else
        ni=${ni:-1}; nj=${nj:-1}; nk=${nk:-1}
    fi
    local prod=$(( ni * nj * nk ))
    if [[ "$prod" -ne "$np" ]]; then
        log "ERROR: NODESI*NODESJ*NODESK = $ni*$nj*$nk = $prod but NP = $np"
        exit 3
    fi
    log "DMP decomposition: nodesi=$ni nodesj=$nj nodesk=$nk"
    # Strip any pre-existing nodes* lines (case-insensitive, leading whitespace OK).
    sed -i -E '/^[[:space:]]*nodes[ijk][[:space:]]*=.*$/Id' "$mfx"
    {
        echo ""
        echo "### DMP partitioning (injected by entrypoint.sh)"
        echo "  nodesi = $ni"
        echo "  nodesj = $nj"
        echo "  nodesk = $nk"
    } >> "$mfx"
}

run_mfix() {
    local mfx=$1 np=$2
    inject_dmp_keywords "$mfx" "$np"
    if [[ "$np" -le 1 ]]; then
        log "exec mfixsolver -f $mfx"
        mfixsolver -f "$mfx" ${EXTRA_MFIX_ARGS:-} 2>&1 | tee run.log
        return ${PIPESTATUS[0]}
    fi
    log "exec mpirun -np $np mfixsolver_dmp -f $mfx"
    mpirun --allow-run-as-root -np "$np" \
        --bind-to none \
        --mca btl_vader_single_copy_mechanism none \
        --mca btl self,vader,tcp \
        mfixsolver_dmp -f "$mfx" ${EXTRA_MFIX_ARGS:-} 2>&1 | tee run.log
    return ${PIPESTATUS[0]}
}

# -- Local smoke-test path ----------------------------------------------------

if [[ "${LOCAL_SMOKE:-0}" == "1" ]]; then
    log "LOCAL_SMOKE=1 — skipping S3 sync, running on /work bind-mount"
    cd /work
    NP="${NP:-1}"
    MFX_FILE="${MFX_FILE:?MFX_FILE must be set}"
    set +e
    run_mfix "$MFX_FILE" "$NP"
    exit $?
fi

# -- AWS Batch path -----------------------------------------------------------

NP="${NP:-1}"
MFX_FILE="${MFX_FILE:?MFX_FILE must be set}"

if [[ -n "${CASE_S3_ARRAY:-}" ]]; then
    IDX_RAW="${AWS_BATCH_JOB_ARRAY_INDEX:?CASE_S3_ARRAY set but no array index}"
    IDX=$(printf '%02d' "$IDX_RAW")
    CASE_S3="${CASE_S3_ARRAY%/}/case_${IDX}/"
    RESULTS_S3="${RESULTS_S3_ARRAY:?RESULTS_S3_ARRAY required with CASE_S3_ARRAY}"
    RESULTS_S3="${RESULTS_S3%/}/case_${IDX}/"
    log "array index $IDX_RAW -> $CASE_S3"
else
    CASE_S3="${CASE_S3:?CASE_S3 or CASE_S3_ARRAY must be set}"
    RESULTS_S3="${RESULTS_S3:?RESULTS_S3 must be set}"
fi

WORKDIR=/work/run
mkdir -p "$WORKDIR"
cd "$WORKDIR"

log "sync in: $CASE_S3 -> $WORKDIR"
aws s3 sync "$CASE_S3" "$WORKDIR/" --no-progress

if [[ ! -f "$MFX_FILE" ]]; then
    log "ERROR: $MFX_FILE not found in $WORKDIR after sync"
    ls -la "$WORKDIR" || true
    exit 2
fi

log "launching: NP=$NP, MFX_FILE=$MFX_FILE"
START=$(date +%s)
set +e
run_mfix "$MFX_FILE" "$NP"
RC=$?
set -e
END=$(date +%s)
log "mfix exit=$RC wall=$((END - START))s"

log "sync out: $WORKDIR -> $RESULTS_S3"
aws s3 sync "$WORKDIR/" "$RESULTS_S3" --no-progress \
    --exclude '*.pyc' --exclude '__pycache__/*'

exit $RC
