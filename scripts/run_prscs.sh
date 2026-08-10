#!/usr/bin/env bash
set -euo pipefail
source "${GMRS_ROOT}/scripts/lib.sh"

target="${1:-all}"
run_one() {
  local layer="$1" feature="$2" chr="$3" summary="$4" n="$5"
  summary="$(resolve_path "$summary")"
  [[ -f "$summary" ]] || die "summary file not found: $summary"
  [[ "$feature" != */* ]] || die "feature_id cannot contain '/': $feature"
  if [[ "$n" == AUTO ]]; then
    n="$(awk -v c="$SUMMARY_N_COLUMN" 'NF>=c {v=$c} END {print v}' "$summary")"
  fi
  [[ "$n" =~ ^[0-9]+([.][0-9]+)?$ ]] || die "invalid sample size for $feature: $n"
  local outdir ref_bim chr_args=()
  outdir="$(weight_dir "$layer")"; mkdir -p "$outdir/logs"
  if [[ "$chr" == ALL ]]; then
    ref_bim="${BIM_PREFIX_TEMPLATE/\{chr\}/}"
  else
    ref_bim="${BIM_PREFIX_TEMPLATE/\{chr\}/$chr}"
    chr_args=(--chrom="$chr")
  fi
  ref_bim="$(resolve_path "$ref_bim")"
  if [[ "$chr" == ALL ]]; then
    [[ -f "${ref_bim}1.bim" ]] || die "reference BIM not found: ${ref_bim}1.bim"
  else
    [[ -f "${ref_bim}.bim" ]] || die "reference BIM not found: ${ref_bim}.bim"
  fi
  local prefix="$outdir/$feature"
  if compgen -G "${prefix}_pst_eff_a1_b0.5_phiauto*.txt" >/dev/null && [[ "${FORCE:-0}" != 1 ]]; then
    log "skip existing weight: $layer/$feature"; return
  fi
  if [[ "${FORCE:-0}" == 1 ]]; then
    local stale=()
    mapfile -t stale < <(compgen -G "${prefix}_pst_eff_a1_b0.5_phiauto*.txt" || true)
    ((${#stale[@]} == 0)) || rm -f -- "${stale[@]}"
  fi
  log "PRS-CS: $layer/$feature (chr=$chr, N=$n)"
  MKL_NUM_THREADS="$THREADS" NUMEXPR_NUM_THREADS="$THREADS" OMP_NUM_THREADS="$THREADS" \
    "$PYTHON_BIN" "$PRSCS_PY" --ref_dir="$LD_REF_DIR" --bim_prefix="$ref_bim" \
      --sst_file="$summary" --n_gwas="$n" "${chr_args[@]}" --n_iter="$N_ITER" \
      --n_burnin="$N_BURNIN" --thin="$THIN" --out_dir="$prefix" \
      >"$outdir/logs/$feature.log" 2>&1
}

for layer in $(layers_for "$target"); do
  manifest="$(manifest_for "$layer")"
  while IFS=$'\t' read -r feature chr summary n _; do
    [[ "$feature" == feature_id ]] && continue
    run_one "$layer" "$feature" "$chr" "$summary" "$n"
  done < "$manifest"
done

