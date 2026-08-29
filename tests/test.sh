#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
example_dir="$root/tests/example"
output="$example_dir/prediction.tsv"
for suffix in bed bim fam; do
  [[ -s "$example_dir/example.$suffix" ]] || {
    echo "Missing example genotype file: $example_dir/example.$suffix" >&2
    exit 1
  }
done
[[ -s "$example_dir/covariates.tsv" ]] || {
  echo "Missing example covariates: $example_dir/covariates.tsv" >&2
  exit 1
}
[[ -s "$root/model/GMRS_model.rds" ]] || {
  echo "Missing bundled model: $root/model/GMRS_model.rds" >&2
  exit 1
}
rm -f -- "$output"
GMRS_CONFIG="$example_dir/config.sh" "$root/bin/gmrs" predict \
  --bfile "$example_dir/example" \
  --covariates "$example_dir/covariates.tsv" \
  --out "$output"
"$CONDA_BIN" run --no-capture-output -n "$R_CONDA_ENV" Rscript --vanilla -e '
  suppressPackageStartupMessages(library(data.table))
  args <- commandArgs(trailingOnly = TRUE)
  x <- fread(args[[1]])
  required <- c("FID", "IID", "PRS_raw", "TRS_raw", "ProRS_raw",
                "PRS", "TRS", "ProRS", "linear_predictor", "disease_probability")
  stopifnot(nrow(x) == 2L,
            identical(x$IID, c("EXAMPLE1", "EXAMPLE2")),
            all(required %in% names(x)),
            all(is.finite(x$disease_probability)),
            all(x$disease_probability >= 0 & x$disease_probability <= 1))
  print(x)
  cat("GMRS end-to-end example PASS\n")
' "$output"
