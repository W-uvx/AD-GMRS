#!/usr/bin/env bash
set -euo pipefail
source "${GMRS_ROOT}/scripts/lib.sh"

[[ -f "$PRSCS_PY" ]] || die "PRS-CS script not found: $PRSCS_PY"
[[ -d "$LD_REF_DIR" ]] || die "LD reference directory not found: $LD_REF_DIR"
command -v "$PLINK_BIN" >/dev/null 2>&1 || die "PLINK not found: $PLINK_BIN"
[[ -f "$(resolve_path "$TRAIN_BFILE").bed" && -f "$(resolve_path "$TRAIN_BFILE").bim" && -f "$(resolve_path "$TRAIN_BFILE").fam" ]] || die "Incomplete training PLINK prefix: $TRAIN_BFILE"
[[ "$BIM_PREFIX_TEMPLATE" == *'{chr}'* ]] || die "BIM_PREFIX_TEMPLATE must contain {chr}"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python not found: $PYTHON_BIN"
command -v sha256sum >/dev/null 2>&1 || die "sha256sum not found"
command -v sort >/dev/null 2>&1 || die "sort not found"
[[ -x "${CONDA_BIN:-}" ]] || command -v "${CONDA_BIN:-conda}" >/dev/null 2>&1 || die "Conda not found"
"${CONDA_BIN:-conda}" run -n "${R_CONDA_ENV:-r_environment_name}" Rscript --version >/dev/null 2>&1 || die "Rscript unavailable in Conda environment: ${R_CONDA_ENV:-r_environment_name}"
"${CONDA_BIN:-conda}" run --no-capture-output -n "${R_CONDA_ENV:-r_environment_name}" Rscript "${GMRS_ROOT}/scripts/validate_inputs.R"

for layer in gwas eqtl pqtl; do
  manifest="$(manifest_for "$layer")"
  [[ -f "$manifest" ]] || die "$layer manifest not found: $manifest"
  awk -F '\t' 'NR==1 {if ($1!="feature_id" || $2!="chromosome" || $3!="summary_file" || $4!="sample_size") exit 1}
    NR>1 {if (NF<4 || $1=="" || $2=="" || $3=="" || $4=="") exit 1}' "$manifest" || die "invalid manifest: $manifest"
  log "$layer manifest OK ($(($(wc -l < "$manifest") - 1)) features)"
done
log "validation completed"

