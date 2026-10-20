# Running MFiX on AWS Batch

This project runs MFiX 26.1.2 on AWS Batch (EC2, on-demand) for workloads that
are too long or too parallel for the laptop (Day 11-12 CFD-DEM, Day 13 16-run
parametric sweep). The laptop is still used to render cases and post-process
results; Batch only runs the solver.

## What the stack looks like

```
  laptop                      AWS
  ------                      ---
  render case ---s3 sync---> s3://triso-cfd-mfix/cases/<case>/           (single)
                             s3://triso-cfd-mfix/sweep/<name>/case_NN/   (array)
                                            |
  submit-job --------------> Batch job queue (triso-cfd-q)
                                            |
                                      Compute env (EC2 on-demand
                                      c7a.4xlarge/8xlarge,
                                      BEST_FIT_PROGRESSIVE)
                                            |
                                      Container (ECR image
                                      triso-cfd/mfix:26.1.2)
                                            |
                             s3://triso-cfd-mfix/results/<case>/              (single)
                             s3://triso-cfd-mfix/sweep/<name>/results/case_NN/ (array)
  poll + s3 sync <----------            CloudWatch logs
```

Pieces, with their file / resource name:

| Resource | Name | Where it is defined |
| --- | --- | --- |
| Docker image | `triso-cfd/mfix:26.1.2` | `docker/Dockerfile`, `docker/entrypoint.sh` |
| ECR repo | `triso-cfd/mfix` | provisioned by `src/aws/provision.py` |
| S3 bucket | `triso-cfd-mfix` | provisioned by `src/aws/provision.py` |
| IAM roles | service, instance, job, execution | provisioned by `src/aws/provision.py` |
| Compute env | `triso-cfd-ce` | provisioned by `src/aws/provision.py` |
| Job queue | `triso-cfd-q` | provisioned by `src/aws/provision.py` |
| Job definition | `triso-cfd-mfix-26-1-2` | provisioned by `src/aws/provision.py` |
| Submit script | `.venv/bin/python -m src.aws.submit` | `src/aws/submit.py` |
| Poll script | `.venv/bin/python -m src.aws.poll` | `src/aws/poll.py` |
| Shared config | `src/aws/config.py` | all names, region, instance types |
| Run manifest | `results/run_log.csv` (one 'submitted' row per job, carries `aws_job_id`) | appended by `src/aws/submit.py:append_run_log_row` (submit.py:61) |

Note on tracking: there is no stand-alone AWS job manifest file. Each case
directory carries its own `cases/<case>/manifest.json` (written by the day-N
driver) with physical parameters, and `results/run_log.csv` is the project-
wide tracked manifest (CLAUDE.md rule 15: "Submit jobs via the tracked
manifest system"). `submit.py` appends a `status=submitted` row *before*
dispatch, so the aws_job_id and S3 prefix are logged even if Batch rejects
the submission.

Everything is on-demand. No spot, no EFA, no multi-node parallel, no
ParallelCluster — MFiX's DMP decomposition caps around 16 ranks at ~10^5 cells
so a single c7a instance is the right granularity, and on-demand means no run
gets reclaimed mid-flight (Day 11-12's 20-60 h CFD-DEM run is the reason).

## One-time setup

1. Install Docker (already present on this laptop) and configure the AWS CLI
   with credentials that can create S3 buckets, IAM roles, ECR repos, and
   Batch resources.
2. Provision the stack and build + push the image:

   ```bash
   .venv/bin/python -m src.aws.provision --build
   ```

   Idempotent — safe to re-run. `--build` is only needed the first time or
   when `docker/Dockerfile` / `docker/entrypoint.sh` changes.

3. Verify:

   ```bash
   .venv/bin/python -m src.aws.provision --status
   ```

## Per-run workflow

### Single run (Day 10, Day 11-12)

```bash
# 1. render the case locally (same driver that worked on the laptop)
.venv/bin/python -m src.run_V_3d_day10

# 2. dry-run plan (no S3, no submit — useful for sanity-check)
.venv/bin/python -m src.aws.submit \
    --case-dir cases/V_3d/day10 \
    --mfx V_3d_day10.mfx \
    --job-name d10-pic \
    --vcpu 16 --memory 30000 \
    --dry-run

# 3. real submit (syncs case to S3, submits Batch job, appends run_log row)
.venv/bin/python -m src.aws.submit \
    --case-dir cases/V_3d/day10 \
    --mfx V_3d_day10.mfx \
    --job-name d10-pic \
    --vcpu 16 --memory 30000
# prints: [submit] job id: <job-id>

# 4. watch until done, then pull results into cases/.../aws_out
.venv/bin/python -m src.aws.poll <job-id> --watch --pull cases/V_3d/day10/aws_out
```

Day 11-12 is identical, just with a larger `--vcpu` / `--memory` and the
Day-11-12 driver.

### Array run (Day 13)

```bash
# 1. render all 16 sweep leaves under cases/sweep13/case_00/ ... case_15/
.venv/bin/python -m src.run_V_3d_sweep13   # (day 13 driver)

# 2. submit as a Batch array
.venv/bin/python -m src.aws.submit \
    --array-root cases/sweep13 \
    --mfx V_3d_sweep.mfx \
    --job-name d13-sweep \
    --array-size 16 \
    --vcpu 16 --memory 30000

# 3. poll the parent job id; children are <id>:0, <id>:1, ...
.venv/bin/python -m src.aws.poll <parent-id> --watch
```

The entrypoint reads `AWS_BATCH_JOB_ARRAY_INDEX`, formats it as `case_NN`,
and syncs only its own leaf from S3. The compute env's `MAX_VCPUS = 64`
(in `src/aws/config.py`) lets ~4 children run in parallel; raise it to 256
to run all 16 at once (same total cost, just finishes sooner).

## Entrypoint contract

The container's `/opt/entrypoint.sh` reads these env vars (passed via Batch's
`containerOverrides`):

| Var | Required | Meaning |
| --- | --- | --- |
| `CASE_S3` | yes (single) | `s3://<bucket>/cases/<case>/` to sync in |
| `RESULTS_S3` | yes (single) | `s3://<bucket>/results/<case>/` to sync out |
| `CASE_S3_ARRAY` | yes (array) | `s3://<bucket>/sweep/<name>/` root for array children |
| `RESULTS_S3_ARRAY` | yes (array) | corresponding results root |
| `MFX_FILE` | yes | deck filename inside the case dir |
| `NP` | no | MPI rank count (default 1; `src/aws/submit.py` sets it = `--vcpu`) |
| `EXTRA_MFIX_ARGS` | no | extra argv for `mfixsolver` (e.g. `-r` to restart) |
| `NODESI` / `NODESJ` / `NODESK` | no | explicit DMP partitioning; product must equal `NP`. If all unset, `docker/entrypoint.sh:inject_dmp_keywords` auto-picks `NODESI=1, NODESJ=NP, NODESK=1` (axial decomposition for spouted-bed geometries) and patches the keywords into the deck before launch — `mfixsolver_dmp` does not accept CLI keyword overrides. |
| `LOCAL_SMOKE` | no | if `1`, skip S3 and run on `/work` bind-mount |

## Smoke-testing the image

Before the first cloud submit — or after any Dockerfile change — run the
container locally against the tut_tfm2d case with a shortened `tstop`:

```bash
mkdir /tmp/mfix-smoke && cp cases/tut_tfm2d/example_2d_fluidized_bed.mfx /tmp/mfix-smoke/
sed -i.bak 's/tstop               = 5.0/tstop               = 0.02/' \
    /tmp/mfix-smoke/example_2d_fluidized_bed.mfx
docker run --rm --platform=linux/amd64 \
    -e LOCAL_SMOKE=1 \
    -e MFX_FILE=example_2d_fluidized_bed.mfx \
    -e NP=1 \
    -v /tmp/mfix-smoke:/work \
    triso-cfd/mfix:26.1.2
ls /tmp/mfix-smoke   # expect BACKGROUND_0000.vtu, *.RES, *.SP1-8, *.LOG
rm -rf /tmp/mfix-smoke
```

The initial setup was validated this way: pre-processing emitted
`NUMBER OF FLUID CELLS = 14400`, the solver reached `t = 0.02 s` in ~18 s wall,
and all expected output files appeared.

## Job sizing from the case manifest (C26)

`src/aws/submit.py` reads `aws_vcpu`, `aws_memory_mb`, `nranks` and
`nodes_i/j/k` from the case's `manifest.json` for any flag left unset.
Explicit flags win. Manifests without `aws_vcpu`, which covers every V_3d case,
keep the old defaults (16 vCPU / 30000 MB).

C26 (`src/run_C26_3d.py`) writes 4 vCPU, 6000 MB and 4 ranks (2×1×2,
gridmap i = k = [13,14]). The account quota is 32 vCPU, so 8 C26 jobs run at
once, packed 4 per c7a.4xlarge. The basis is in `src/aws/config.py`
(`C26_MEMORY_MB`) and `notes/c26_conry2026_case.md` ("AWS run cost").

```bash
.venv/bin/python -m src.aws.submit --case-dir cases/C26_3d/c02 \
    --mfx C26_c02.mfx --job-name c26-c02 --dry-run   # drop --dry-run to submit
```

## Costs (on-demand, us-east-1 c7a)

Single-node, no EFA, no placement group.

| Day | Instance | Hours | $/hr | Est. $ |
| --- | --- | --- | --- | --- |
| 10 (single PIC) | c7a.4xlarge (16 vCPU) | 1-2 | $0.821 | $0.82-1.64 |
| 11-12 (CFD-DEM) | c7a.8xlarge (32 vCPU, 16 ranks) | 20-60 | $1.642 | $33-100 |
| 13 (16-run sweep) | c7a.4xlarge × 16 | ~1 each | $0.821 | $13-26 |
| **Total** | | | | **~$47-128** |

Add ~$0.10/month S3 storage; CloudWatch log retention is free in the default
30-day window.

## Common operations

```bash
# list recent jobs on the queue
aws batch list-jobs --job-queue triso-cfd-q --max-results 20

# tail live logs for a running job (fetch logStreamName via describe-jobs)
aws logs tail /aws/batch/job --follow \
    --log-stream-name-prefix <logStreamName>

# cancel a running job
aws batch terminate-job --job-id <id> --reason "manual cancel"

# tear down (if ever needed) — reverse of provision, delete in this order:
#   job queue -> compute env -> job definition (deregister) ->
#   S3 bucket (empty first) -> ECR repo -> IAM roles
```

## Not implemented

- No multi-node MPI. If a run ever outgrows a single c7a.8xlarge (~16 effective
  ranks), the next step is `type: MULTI_NODE` with EFA — not addressed here.
- No `retryStrategy.attempts > 1`. `src/aws/provision.py:provision_job_definition`
  registers the job def with `retryStrategy={"attempts": 1}` (= no retry).
  On-demand doesn't need it; a Batch retry would just mask a real MFiX failure.
- No cost dashboard. The run_log.csv row carries `aws_job_id`; wall-clock and
  $/run can be computed from Batch's `startedAt` / `stoppedAt` fields on demand.
