# Running MFiX on AWS Batch

This project runs MFiX 23.2 on AWS Batch (EC2, on-demand; 26.1.2 until
2026-10-24, see "MFiX version" below) for workloads that
are too long or too parallel for the laptop (on this branch, the Case C26 3 s
CFD-CGDEM runs c01–c09; the Day 10–13 Case V/D runs and their drivers are on
branch `case-v-d-day10-13`). The laptop is still used to render cases and
post-process results; Batch only runs the solver.

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
                                      triso-cfd/mfix:23.2)
                                            |
                             s3://triso-cfd-mfix/results/<case>/              (single)
                             s3://triso-cfd-mfix/sweep/<name>/results/case_NN/ (array)
  poll + s3 sync <----------            CloudWatch logs
```

Pieces, with their file / resource name:

| Resource | Name | Where it is defined |
| --- | --- | --- |
| Docker image | `triso-cfd/mfix:23.2` (`:26.1.2` still in ECR) | `docker/Dockerfile`, `docker/entrypoint.sh` |
| ECR repo | `triso-cfd/mfix` | provisioned by `src/aws/provision.py` |
| S3 bucket | `triso-cfd-mfix` | provisioned by `src/aws/provision.py` |
| IAM roles | service, instance, job, execution | provisioned by `src/aws/provision.py` |
| Compute env | `triso-cfd-ce` | provisioned by `src/aws/provision.py` |
| Job queue | `triso-cfd-q` | provisioned by `src/aws/provision.py` |
| Job definition | `triso-cfd-mfix-23-2` (named from `IMAGE_TAG`; `triso-cfd-mfix-26-1-2` still registered) | `src/aws/provision.py` (`--image-only` under `aws login`) |
| Submit script | `.venv/bin/python -m src.aws.submit` | `src/aws/submit.py` |
| Poll script | `.venv/bin/python -m src.aws.poll` | `src/aws/poll.py` |
| CLI wrappers | `aws batch submit-job` / `describe-jobs` via subprocess | `src/aws/cli.py` (used by submit, poll and `notes/c26_watch_post.py`) |
| Shared config | `src/aws/config.py` | all names, region, instance types |
| Run manifest | `results/run_log.csv` (one 'submitted' row per job, carries `aws_job_id`) | appended by `src/aws/submit.py::append_run_log_row` |

Note on tracking: there is no stand-alone AWS job manifest file. Each case
directory carries its own `cases/<case>/manifest.json` (written by the case
driver, e.g. `src/run_C26_3d.py`) with physical parameters, and `results/run_log.csv` is the project-
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
   Batch resources. Since 2026-10-22 the credentials come from `aws login`;
   see "Credentials" below for what that means for the Python tooling.
2. Provision the stack and build + push the image:

   ```bash
   .venv/bin/python -m src.aws.provision --build
   ```

   This boto3 path does not work under `aws login` (see "Credentials").
   On the existing stack, a new image or MFiX release only needs the
   AWS-CLI path, which builds, pushes and registers the job definition:

   ```bash
   .venv/bin/python -m src.aws.provision --image-only --build
   ```

   Idempotent — safe to re-run. `--build` is only needed the first time or
   when `docker/Dockerfile` / `docker/entrypoint.sh` changes.

3. Verify:

   ```bash
   .venv/bin/python -m src.aws.provision --status
   ```

## Credentials: `aws login` and the CLI wrappers (2026-10-22)

Since 2026-10-22 the laptop authenticates with the AWS CLI v2 `aws login`
provider (CLI 2.38.4, universal build, reinstalled after the macOS 27 upgrade
broke the x86-only install and the old access key stopped working). boto3 in
`.venv` cannot read that provider without the optional `botocore[crt]`
package. User decision (2026-10-22): do not add the dependency; shell out to
the CLI instead.

- `src/aws/cli.py` runs `aws <service> <op> --region <REGION> --output json`
  and parses stdout. `batch_submit_job(payload)` passes the same API payload
  boto3's `submit_job(**payload)` took, as `--cli-input-json`;
  `batch_describe_jobs(ids)` returns the `jobs` list. A non-zero exit raises
  `RuntimeError` with the CLI's stderr.
- `src/aws/submit.py`, `src/aws/poll.py` and `notes/c26_watch_post.py` use
  the wrappers and no longer import boto3. `tests/test_aws_cli.py` mocks
  `subprocess.run` and checks the command shape, the JSON handling, the
  error path and that submit / poll stay boto3-free.
- `src/aws/provision.py` (one-off stack setup) still uses boto3, except
  `--image-only` (2026-10-24), which pushes an image and registers its job
  definition through `src/aws/cli.py`. The full provisioning path will
  fail under `aws login` until it is converted the same way;
  `requirements.txt` keeps boto3 / botocore pinned for it. Do not install
  `awscrt` / `botocore[crt]` into `.venv` to work around this.
- Sessions expire after some hours. Every wrapper then fails with

  ```
  aws: [ERROR]: Your session has expired. Please reauthenticate using 'aws login'.
  ```

  Run `aws login` (it opens a browser) and retry the command. Jobs already
  on Batch are unaffected; only the laptop-side calls stop.
- The first `aws s3 sync` of a submit has failed once with exit 1 and passed
  on retry (2026-10-22); treat one sync failure as transient before
  debugging.
- Verified 2026-10-22: nine C26 `c0X_momentum24` jobs submitted through the
  CLI path; `poll` and the watcher read their status.

## MFiX version (23.2 since 2026-10-24)

The image runs MFiX 23.2, the release C26 ran (notes/c26_conry2026_case.md),
so a C26 deck can be run on C26's own solver. **23.2 is the project default
for all new runs** (user decision 2026-10-24): the default image tag and job
definition (`src/aws/config.py`) and the default deck release
(`params/numerics.yaml mfix_version`). 26.1.2 is the release of every run
before 2026-10-24 and the laptop's native DMP install.

- **What the image holds** (`docker/Dockerfile`): NETL's Linux install command
  (`mfix mfix-gui mfix-src mfix-doc mfix-solver`, all `=23.2`) plus gcc /
  gfortran 13, openmpi 4.1, cmake and make. Stock 23.2 ships only a serial
  `mfixsolver` (no MPI, no OpenMP; `openmpi` first appears as a dependency
  in 23.3.2), so the Dockerfile builds a batch DMP solver from the 23.2
  source with MFiX's own tool (`python -m mfixgui.build_mfixsolver --batch
  --dmp`) and installs it as `mfixsolver_dmp`. The entrypoint contract is
  unchanged: NP = 1 runs the stock serial solver, NP > 1 runs `mpirun
  mfixsolver_dmp`.
- **conda, not mamba**: mamba 2.x in `miniforge3:latest` fetches conda-forge
  as shards and only for the named packages, so the MFiX channel's
  conda-forge dependencies (flask, matplotlib) "do not exist". The
  Dockerfile uses `conda create`, as NETL's command does.
- **The channel link expires.** It is the per-download URL from NETL's
  registration e-mail; the 26.1.2 one returns HTTP 401 now. Rebuilding the
  image later needs a fresh link; the pushed `triso-cfd/mfix:23.2` in ECR
  does not.
- **Decks.** 23.2 aborts on unknown keywords (Error 2000) and four keywords
  in every 26.1.2 C26 deck do not exist in 23.2. `src.run_C26_3d` writes
  23.2 decks by default (see notes/c26_conry2026_case.md "MFiX 23.2 decks
  and AWS image"); `--mfix-version 26.1.2` re-renders the old release. The
  manifest records `mfix_version`; `src/aws/submit.py` refuses a deck
  whose version differs from `IMAGE_TAG`, also under `--dry-run`, and the
  driver refuses to render into a run dir holding a case for another
  release. Manifests without the field are 26.1.2 decks, so every existing
  case dir (C26 before 2026-10-24, V_3d) needs `TRISO_IMAGE_TAG=26.1.2` to
  be resubmitted. Only the C26 generator has a 23.2 mode.
- **Old release on demand.** `TRISO_IMAGE_TAG=26.1.2` selects the 26.1.2
  image and job definition (`triso-cfd-mfix-26-1-2`), both still on AWS
  (for reruns of pre-2026-10-24 cases only). The
  `results/run_log.csv` row of every submit carries `mfix=<tag>` in
  `key_params`.
- **Local runs of a 23.2 deck.** The laptop has no 23.2 DMP solver (its
  `mfix-23.2` conda env is the stock serial one), so the case's `run.sh`
  runs the image with `LOCAL_SMOKE=1` on the case dir: the same binary as
  on AWS, amd64 under emulation, so slower than the native 26.1.2. Another
  machine needs the image first (`docker pull` from ECR, or
  `src.aws.provision --image-only --build` with a valid channel link).

## Per-run workflow

### Single run (C26)

```bash
# 1. render the case locally (same driver that worked on the laptop); the
#    deck is for 23.2, the image's release, by default
.venv/bin/python -m src.run_C26_3d --case c04 --run-dir-suffix _massonly_v232

# 2. dry-run plan (no S3, no submit; sizing comes from the manifest)
.venv/bin/python -m src.aws.submit \
    --case-dir cases/C26_3d/c04_massonly_v232 \
    --mfx C26_c04_massonly_v232.mfx \
    --job-name c26-c04-massonly-v232 \
    --dry-run

# 3. real submit (syncs case to S3, submits Batch job, appends run_log row):
#    the same command without --dry-run
# prints: [submit] job id: <job-id>

# 4. watch until done, then pull results into cases/.../aws_out
.venv/bin/python -m src.aws.poll <job-id> --watch --pull cases/C26_3d/c04_massonly_v232/aws_out
```

### Array run

```bash
# 1. render one complete case dir per index: <sweep_root>/case_00/, case_01/, ...
#    (the Day-13 driver src/run_V_3d_day13.py, on branch case-v-d-day10-13,
#    writes cases/D_sweep13/ this way)

# 2. submit as a Batch array
.venv/bin/python -m src.aws.submit \
    --array-root cases/<sweep_root> \
    --mfx <deck>.mfx \
    --job-name <sweep-name> \
    --array-size N \
    --vcpu 16 --memory 30000

# 3. poll the parent job id; children are <id>:0, <id>:1, ...
.venv/bin/python -m src.aws.poll <parent-id> --watch
```

The entrypoint reads `AWS_BATCH_JOB_ARRAY_INDEX`, formats it as `case_NN`,
and syncs only its own leaf from S3. The compute env's `MAX_VCPUS = 64`
(in `src/aws/config.py`) caps how many children run at once (4 at 16 vCPU);
raising it runs more at once (same total cost, just finishes sooner).

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
| `NP` | no | MPI rank count (default 1; `src/aws/submit.py` sets it to `--nranks`, else the manifest's `nranks`, else `--vcpu`) |
| `EXTRA_MFIX_ARGS` | no | extra argv for the solver. 23.2: serial `mfixsolver` takes `-l`/`-c`/`-p`; the DMP build takes `-l LEVEL` and `KEYWORD=VALUE`; any other option stops it |
| `NODESI` / `NODESJ` / `NODESK` | no | explicit DMP partitioning; product must equal `NP`. If all unset, `docker/entrypoint.sh:inject_dmp_keywords` auto-picks `NODESI=1, NODESJ=NP, NODESK=1` (axial decomposition for spouted-bed geometries) and patches the keywords into the deck before launch (26.1.2's pymfix `mfixsolver_dmp` took no CLI keyword overrides). A deck that already sets the requested decomposition is left untouched, so a local `LOCAL_SMOKE` run does not rewrite the case deck. |
| `LOCAL_SMOKE` | no | if `1`, skip S3 and run on `/work` bind-mount |

## Smoke-testing the image

Before the first cloud submit — or after any Dockerfile change — run the
image locally on a C26 deck for the image's release, 4 ranks (the DMP
solver), about 10 steps. `--tstop` is a driver flag, so no sed edit of the
deck is needed:

```bash
.venv/bin/python -m src.run_C26_3d --case c02 --tstop 1e-4 \
    --run-dir-suffix _imgsmoke --out-root /tmp/mfix-smoke
/tmp/mfix-smoke/c02_imgsmoke/run.sh     # docker, LOCAL_SMOKE=1, NP=4
grep -E 'Successful seeding|Total number of particles in the system|Final simulation time' \
    /tmp/mfix-smoke/c02_imgsmoke/run.log
rm -rf /tmp/mfix-smoke
```

Validated this way on 2026-10-24 (23.2 image): "Successful seeding",
27,777-27,779 parcels (26.1.2, same deck: 27,778; the PSD draw varies by a
few parcels per run), `Final simulation time = 0.0001 s`, ~13 s wall
including emulated start-up. The old recipe (`cases/tut_tfm2d`) no longer
works: that deck is a 26.1.2 GUI file and 23.2 rejects three of its
keywords (`bc_default_walls`, `des_rigid_motion`, `mu_g_model`).

Two image faults this smoke caught on 2026-10-24:
- the DMP build needs conda-forge `binutils` (no unprefixed `ar` otherwise:
  "Error running link command" on `libmfixcore.a`);
- Open MPI 4.1's `mpirun` aborts with "plm_rsh_agent: ssh : rsh ... could
  not be found" without `OMPI_MCA_plm=isolated` (set in the Dockerfile).

## Job sizing from the case manifest (C26)

`src/aws/submit.py` reads `aws_vcpu`, `aws_memory_mb`, `nranks` and
`nodes_i/j/k` from the case's `manifest.json` for any flag left unset.
Explicit flags win. Manifests without `aws_vcpu`, which covers every V_3d case,
keep the old defaults (16 vCPU / 30000 MB).

C26 (`src/run_C26_3d.py`) writes 4 vCPU, 4 ranks (2×1×2; gridmap
i = k = [12,12] on the 24-cell layout that is the default since 2026-10-22,
[13,13] on the 26-cell axis-face layout, [13,14] on the old 27-cell layout)
and 6000 MB, scaled up by parcel count (c09: 8000 MB). The
account quota is 32 vCPU, so 8 C26 jobs run at once, packed 4 per
c7a.4xlarge. The basis is in `src/aws/config.py` (`C26_MEMORY_MB`) and
`notes/c26_conry2026_case.md` ("AWS run cost").

```bash
.venv/bin/python -m src.aws.submit --case-dir cases/C26_3d/c02_massonly_v232 \
    --mfx C26_c02_massonly_v232.mfx --job-name c26-c02-massonly-v232 --dry-run   # drop --dry-run to submit
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
