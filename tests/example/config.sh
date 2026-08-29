#!/usr/bin/env bash
PLINK_BIN="/storage/main/projects/bigcgpu-prj/wangjiahao/software/plink/plink"
CONDA_BIN="/storage/main/projects/bigcgpu-prj/wangjiahao/software/Anaconda3/Anaconda3/bin/conda"
R_CONDA_ENV="r4"
GWAS_MANIFEST="manifests/gwas.tsv"
EQTL_MANIFEST="manifests/eqtl.tsv"
PQTL_MANIFEST="manifests/pqtl.tsv"
WEIGHT_ROOT="../weight"
GWAS_WEIGHT_SUBDIR="PRS"
EQTL_WEIGHT_SUBDIR="WH_GTEx"
PQTL_WEIGHT_SUBDIR="WH_Iceland"
THREADS=2
KEEP_FILE=""
MIN_SNP_MATCH_RATE=0.70
