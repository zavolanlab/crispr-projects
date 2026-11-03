#!/bin/bash
# Slurm job state for Snakemake's cluster-generic executor: prints running, success or failed.
# sacct can fail briefly (slurmdbd); it is retried, and squeue is asked for jobs it does not know yet. A job that
# neither knows is reported as failed, so that Snakemake does not wait for it forever.
id=${1%%;*}
state=""
for attempt in 1 2 3; do
  state=$(sacct -j "$id" --format=State --noheader --parsable2 2>/dev/null | head -1)
  [ -n "$state" ] && break
  state=$(squeue -j "$id" -h -o %T 2>/dev/null)
  [ -n "$state" ] && break
  sleep 5
done
case "${state%% *}" in
  COMPLETED) echo success ;;
  PENDING|RUNNING|CONFIGURING|COMPLETING|REQUEUED|RESIZING|SUSPENDED|REQUEUE_HOLD|REQUEUE_FED) echo running ;;
  *) echo failed ;;
esac
