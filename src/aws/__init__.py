"""AWS Batch job orchestration for MFiX runs.

Modules:
  config      - project-wide names/constants (bucket, queue, job def, region)
  provision   - idempotent setup of S3 bucket, ECR repo, IAM roles, Batch
                compute env + queue + job definition
  submit      - render case locally, s3 sync, batch submit-job (single/array),
                append row to results/run_log.csv with status=submitted
  poll        - poll Batch job state, pull results on success, update log row
"""
