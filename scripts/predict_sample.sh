#!/usr/bin/env bash
set -euo pipefail
source "${GMRS_ROOT}/scripts/lib.sh"
bfile=""; covariates=""; model="$GMRS_ROOT/model/GMRS_model.rds"; out="$GMRS_ROOT/results/predictions/prediction.tsv"; keep_run=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bfile) [[ $# -ge 2 ]] || die "--bfile requires a value"; bfile="$2"; shift 2;;
    --covariates) [[ $# -ge 2 ]] || die "--covariates requires a value"; covariates="$2"; shift 2;;
    --model) [[ $# -ge 2 ]] || die "--model requires a value"; model="$2"; shift 2;;
    --out) [[ $# -ge 2 ]] || die "--out requires a value"; out="$2"; shift 2;;
    --keep-run) keep_run=1; shift;;
    *) die "unknown prediction argument: $1";;
  esac
done
[[ -n "$bfile" && -n "$covariates" ]] || die "--bfile and --covariates are required"
bfile="$(resolve_path "$bfile")"; covariates="$(resolve_path "$covariates")"; model="$(resolve_path "$model")"; out="$(resolve_path "$out")"
[[ -f "$model" ]] || die "model not found: $model"
[[ -f "$covariates" ]] || die "covariate file not found: $covariates"
run_id="predict-$(date +%Y%m%dT%H%M%S)-$$"; run_dir="$GMRS_ROOT/results/runs/$run_id"; score_root="$run_dir/scores"
mkdir -p "$score_root" "$(dirname "$out")"
cleanup(){ status=$?; if [[ $status -eq 0 && "$keep_run" != 1 ]]; then rm -rf -- "$run_dir"; elif [[ $status -ne 0 ]]; then log "prediction failed; diagnostics retained: $run_dir"; fi; return $status; }
trap cleanup EXIT
SCORE_ROOT="$score_root" "$GMRS_ROOT/scripts/score_cohort.sh" all "$bfile"
"${CONDA_BIN:-conda}" run --no-capture-output -n "${R_CONDA_ENV:-r_environment_name}" Rscript \
  "$GMRS_ROOT/scripts/predict_model.R" --model "$model" --score-root "$score_root" --covariates "$covariates" --out "$out"
log "prediction saved: $out"
if [[ "$keep_run" == 1 ]]; then log "QC and intermediate scores retained: $run_dir"; fi
