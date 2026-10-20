"""AWS Batch job orchestration for MFiX runs.

Modules:
  config      - project-wide names/constants (bucket, queue, job def, region)
  provision   - idempotent setup of S3 bucket, ECR repo, IAM roles, Batch
                compute env + queue + job definition
  submit      - s3 sync an already-rendered case dir, batch submit-job
                (single/array), append a status=submitted row to
                results/run_log.csv
  poll        - poll Batch job state, optionally pull results on success
                (does not edit results/run_log.csv)
"""
