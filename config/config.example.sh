#!/usr/bin/env bash
# Copy this file to config/config.sh, then replace every /path/to/... placeholder.

PRSCS_PY="/path/to/software/PRScs.py"
PYTHON_BIN="/path/to/software/python"
LD_REF_DIR="/path/to/reference/ldblk_1kg_EUR"
BIM_PREFIX_TEMPLATE="/path/to/reference/genotypes/chr{chr}"
PLINK_BIN="/path/to/software/plink"
CONDA_BIN="/path/to/software/conda"
R_CONDA_ENV="r_environment_name"

# PLINK prefix only; do not append .bed/.bim/.fam.
TRAIN_BFILE="/path/to/training/genotypes/test_cohort"

GWAS_MANIFEST="manifests/gwas.tsv"
EQTL_MANIFEST="manifests/eqtl.tsv"
PQTL_MANIFEST="manifests/pqtl.tsv"

WEIGHT_ROOT="/path/to/output/weights"
GWAS_WEIGHT_SUBDIR="gwas"
EQTL_WEIGHT_SUBDIR="eqtl"
PQTL_WEIGHT_SUBDIR="pqtl"

N_ITER=1000
N_BURNIN=500
THIN=5
THREADS=10
SUMMARY_N_COLUMN=8

# Optional FID/IID keep file. Leave empty to use every genotype sample.
KEEP_FILE=""

# Optional evidence tables with feature_id, PP4 and Z columns.
EQTL_EVIDENCE=""
PQTL_EVIDENCE=""
INTACT_EPSILON=0.10
INTACT_GAMMA=1

POPULATION_PROFILE="/path/to/population/test_population_profile.tsv"
MIN_SNP_MATCH_RATE=0.70

