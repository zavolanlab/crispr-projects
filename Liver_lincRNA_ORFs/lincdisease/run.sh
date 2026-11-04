#!/bin/bash
# Run the lincdisease workflow (workflow/Snakefile) on sciCORE. Start it on a login node: downloads need the network,
# the other steps are submitted to Slurm. See workflow/README.md.
#
#   run.sh [--root DIR] COMMAND [ARGS]
#
#   plan  [TARGET...]   what would run (changes nothing); plan -v also says why
#   start [TARGET...]   run in the background (default target: all); safe to log out
#   status              is it running, recent log lines, Slurm jobs
#   stop                stop the workflow and cancel its Slurm jobs
#   redo STEP...        rerun these steps and everything downstream of them
#   unlock              after a crash or a killed login session ("Directory cannot be locked")
#   adopt [TARGET...]   accept existing outputs as up to date without running them (default: all), e.g. after an
#                       edit that cannot change them (comments, logging); say why in the commit message
#
# The paths of config/config.yaml are relative to the root directory: --root DIR, else the current directory.
# Targets: all (annotation tables and HTML table), govaere (Govaere 2020 CRAMs), or any output file / rule name.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$PWD
if [ "${1:-}" = --root ]; then ROOT=$(cd "${2:?--root needs a directory}" && pwd); shift 2; fi
ENV=${SNAKEMAKE_ENV:-$ROOT/conda_envs/snakemake}
[ -x "$ENV/bin/snakemake" ] || { echo "no Snakemake in $ENV: run from the project root, pass --root DIR, or set SNAKEMAKE_ENV"; exit 1; }
export PATH=$ENV/bin:$PATH
cfg() { "$ENV/bin/python" -c "import yaml,sys; print(yaml.safe_load(open(sys.argv[1]))[sys.argv[2]])" "$HERE/config/config.yaml" "$1"; }
RESULTS=$(cfg results_dir)
[ -e "$ROOT/$RESULTS" ] || { echo "$ROOT has no $RESULTS: run from the project root or pass --root DIR (see workflow/README.md)"; exit 1; }
OUT=$ROOT/$(cfg out_dir)
LOGS=$OUT/logs/workflow
PIDFILE=$LOGS/snakemake.pid
mkdir -p "$LOGS"
cd "$OUT"   # Slurm jobs start in the submission directory

# --config takes several values: it goes first, so that targets after the last option are not read as config
SMK=(snakemake --config "root=$ROOT" --snakefile "$HERE/workflow/Snakefile" --profile "$HERE/workflow/profile"
     --cluster-generic-status-cmd "$HERE/workflow/scripts/slurm_status.sh")
smk() { "${SMK[@]}" "$@"; }
# the pid file may outlive the login session and its pid be reused: check that the process is Snakemake
running() { [ -s "$PIDFILE" ] && { tr '\0' ' ' < "/proc/$(cat "$PIDFILE")/cmdline"; } 2>/dev/null | grep -q snakemake; }
start() {
  if running; then echo "already running (pid $(cat "$PIDFILE")); see run.sh status"; exit 1; fi
  log=$LOGS/$(date +%Y%m%d_%H%M%S).log
  ln -sfn "$log" "$LOGS/latest.log"
  setsid nohup bash -c '
    "$@"; rc=$?
    if [ $rc = 0 ]; then msg="lincdisease workflow finished"; else msg="lincdisease workflow FAILED (exit $rc)"; fi
    echo "$(date) $msg" >> '"$log"'
    command -v sandbox-notify > /dev/null && sandbox-notify "$msg" || true' _ "${SMK[@]}" "$@" \
    > "$log" 2>&1 < /dev/null &
  echo $! > "$PIDFILE"
  echo "started (pid $!); log: $log"
  echo "follow it with: run.sh status   or   tail -f $LOGS/latest.log"
}

cmd=${1:-plan}; shift || true
case "$cmd" in
  plan)
    if [ "${1:-}" = -v ]; then shift; smk -n "${@:-all}"; exit; fi
    smk -n "${@:-all}" --quiet rules
    echo "Details, with the reason for each job: run.sh plan -v [TARGET...]"
    ;;
  start) start "${@:-all}" ;;
  redo)
    [ $# -gt 0 ] || { echo "usage: run.sh redo STEP..."; exit 1; }
    start all --forcerun "$@"
    ;;
  status)
    if running; then echo "RUNNING (pid $(cat "$PIDFILE"))"; else echo "not running"; fi
    if [ -e "$LOGS/latest.log" ]; then
      echo "--- $(readlink "$LOGS/latest.log")"
      grep -E "steps \(|Error|error|FAILED|finished" "$LOGS/latest.log" | tail -15 || true
    fi
    echo "--- Slurm jobs"; squeue -h -o "%i %j %T %M" | grep -E " ld\." || echo "none"
    ;;
  stop)
    if running; then kill -INT -- "-$(cat "$PIDFILE")"; echo "stopping; Snakemake cancels its Slurm jobs"; else echo "not running"; fi
    ;;
  unlock) smk --unlock ;;
  adopt) smk --touch "${@:-all}" ;;
  *) sed -n 2,18p "$0"; exit 1 ;;
esac
