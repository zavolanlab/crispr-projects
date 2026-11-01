#!/bin/bash
# Slurm job state for Snakemake's cluster-generic executor: prints running, success or failed.
id=${1%%;*}
state=$(sacct -j "$id" --format=State --noheader --parsable2 2>/dev/null | head -1)
[ -n "$state" ] || state=$(squeue -j "$id" -h -o %T 2>/dev/null)
case "$state" in
  COMPLETED) echo success ;;
  ""|PENDING|RUNNING|CONFIGURING|COMPLETING|REQUEUED|RESIZING|SUSPENDED|REQUEUE_HOLD|REQUEUE_FED) echo running ;;
  *) echo failed ;;
esac
