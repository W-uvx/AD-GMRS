# AD-GMRS

Atopic Dermatitis Genotype-based Molecular Risk Score for

AD-GMRS builds and applies a genotype-only disease-risk model from GWAS, eQTL and
pQTL summary statistics. It creates PRS-CS weights, calculates sample-level
scores with PLINK, aggregates transcriptomic and proteomic scores, trains an
elastic-net model, and predicts new samples from genotype files alone.

```text
GWAS/eQTL/pQTL summary statistics
              -> PRS-CS weights
              -> PLINK feature scores
              -> PRS + TRS + ProRS
              -> fitted risk model
new genotype  -> PRS + TRS + ProRS -> disease probability
```

## Requirements

- Bash 4+, Python 3, PRS-CS and PLINK 1.9
- A Conda environment containing R, `data.table` and `glmnet`

## 1. Configure paths

```bash
cp config/config.example.sh config/config.sh
```

Replace every `/path/to/...` value in `config/config.sh`:

| Variable | Required file or directory |
|---|---|
| `PRSCS_PY` | PRS-CS `PRScs.py` script |
| `PYTHON_BIN` | Python executable used to run PRS-CS |
| `LD_REF_DIR` | PRS-CS LD reference directory |
| `BIM_PREFIX_TEMPLATE` | Reference PLINK prefix containing `{chr}` |
| `PLINK_BIN` | PLINK executable |
| `CONDA_BIN` | Conda executable |
| `R_CONDA_ENV` | Conda environment name containing R packages |
| `TRAIN_BFILE` | Training PLINK prefix, without file extension |
| `WEIGHT_ROOT` | Output/root directory for posterior weights |
| `POPULATION_PROFILE` | Population Mean/SD table |

`TRAIN_BFILE=/path/to/training/genotypes/test_cohort` requires:

```text
/path/to/training/genotypes/test_cohort.bed
/path/to/training/genotypes/test_cohort.bim
/path/to/training/genotypes/test_cohort.fam
```

## 2. Configure manifests

Edit the three files under `manifests/`. They are tab-separated and have four
columns:

```text
feature_id  chromosome  summary_file  sample_size
```

Examples:

```text
test_gwas       ALL  /path/to/summary/test_gwas.tsv           100000
test_gene_1     1    /path/to/summary/eqtl/test_gene_1.tsv     AUTO
test_protein_1  1    /path/to/summary/pqtl/test_protein_1.tsv  AUTO
```

`feature_id` must exactly match the prefix used by its PRS-CS weight file.
`AUTO` reads the last value from `SUMMARY_N_COLUMN`; an explicit sample size is
safer when sample size is constant.

Summary statistics must be formatted for the installed PRS-CS version. A common
header is:

```text
SNP  A1  A2  BETA  P  N
```

PRS-CS posterior weights are headerless six-column files:

```text
chromosome  SNP  position  A1  A2  posterior_beta
```

## 3. Population standardisation file

`POPULATION_PROFILE` must be tab-separated:

```text
File            Mean   SD    Min   Max   N
test_gene_1     0.10   0.25  -1.0  1.2   10000
test_protein_1  0.05   0.18  -0.8  0.9   10000
```

Every selected eQTL and pQTL feature must have a finite positive SD. Features
with SD zero must be removed from the manifest. When a GWAS row is absent, its
Mean/SD is estimated once from the training cohort and stored in the model.

## 4. Optional evidence weights

Set `EQTL_EVIDENCE` and/or `PQTL_EVIDENCE` to tables containing:

```text
feature_id  PP4  Z
test_gene_1 0.8  4.2
```

## 5. Validate and generate weights

```bash
bin/AD-GMRS validate
bin/AD-GMRS weights all
```

Existing weights are reused when their filename starts with the matching
`feature_id` and ends in the PRS-CS posterior-weight suffix.

## 6. Build the training phenotype

The disease table requires `eid`, `Date` and `Category`. The baseline table
requires `eid` and one or more assessment-date columns. The genotype FAM must
contain IDs from the same cohort.

```bash
bin/AD-GMRS phenotype \
  --disease /path/to/phenotype/test_diagnoses.tsv \
  --baseline /path/to/phenotype/test_baseline_dates.csv \
  --fam /path/to/training/genotypes/test_cohort.fam \
  --codes TEST_CODE1,TEST_CODE2 \
  --outcome-name test_outcome \
  --out /path/to/output/test_phenotype.tsv \
  --mode prevalent
```

Output:

```text
FID  IID  test_outcome
```

`prevalent` defines a case when the first target diagnosis is on or before the
baseline date. `ever` defines a case when a target diagnosis occurs at any time.

## 7. Score the training cohort and fit the model

```bash
bin/AD-GMRS score all

bin/AD-GMRS train \
  --phenotype /path/to/output/test_phenotype.tsv \
  --outcome test_outcome \
  --population-profile /path/to/population/test_population_profile.tsv \
  --model /path/to/output/model/test_model.rds
```

Training writes:

```text
test_model.rds
test_model.score_distribution.tsv
test_model.oof.tsv
```

The phenotype may cover a subset of scored training subjects; training uses the
FID/IID intersection and reports its size. The model stores feature-level
Mean/SD, aggregation weights, composite-score Mean/SD, exact feature order,
weight-file SHA-256 checksums, fitted coefficients and package versions.

## 8. Predict from a new genotype

```bash
bin/AD-GMRS predict \
  --bfile /path/to/new_sample/test_sample \
  --model /path/to/output/model/test_model.rds \
  --out /path/to/output/predictions/test_sample.tsv \
  --keep-run
```

The only sample-specific input is the PLINK genotype prefix (`.bed/.bim/.fam`).
The output contains:

```text
FID IID PRS_raw TRS_raw ProRS_raw PRS TRS ProRS linear_predictor disease_probability
```

`PRS/TRS/ProRS` are population-standardised composite Z-scores. The prediction
run rejects missing features, mismatched SNP coverage, duplicate IDs and models
that require unavailable clinical covariates. It also rejects weights that are
not byte-for-byte equivalent to the canonical weights used for training.

`disease_probability` is the model probability under the outcome definition and
case prevalence of the training data. It is not automatically an absolute
population risk; external calibration is required when the deployment
population or sampling design differs.

