# BIRDbath pool-effect regression

This bundle contains the final regression endpoint, supporting R functions,
model documentation, a smoke test, and the generated results for the two
BIRDbath `P_reg` thresholds evaluated in the current high-coverage dataset.
Raw input data are not included.

## Model

The endpoint fits

```text
log2(theta_BIRD[v,p]) = alpha
                         + gamma * log2(theta_BIRDbath[v])
                         + beta_pool[p]
                         + error[v,p]
```

with the pool effects constrained to sum to zero. Thus, `beta_pool[p]` is the
estimated deviation of pool `p` from the average pool. Confidence intervals
and p-values use standard errors clustered by variant. See
`docs/final-pool-effect-model.md` for the complete methods paragraph.

## Contents

```text
R/                                  Model, input, and inference functions
scripts/estimate_pool_effects.R     Reusable analysis-ready endpoint
scripts/run_high_coverage_analysis.R
                                    Dataset-specific filtering wrapper
docs/final-pool-effect-model.md     Methods description
tests/final_pool_effect_model_smoke_test.R
                                    Dependency-light regression test
results/preg_gt_0.95/               Final P_reg > 0.95 results
results/preg_gt_0.3/                Sensitivity run at P_reg > 0.3
```

Each result directory contains:

- `pool-effects.tsv`: pool coefficients, clustered uncertainty, p-values,
  multiplicative interpretations, and sample counts.
- `model-summary.tsv`: overall model estimates and fit statistics.
- `analysis-summary.txt`: short human-readable model summary.
- `pool-effects.png`: coefficient plot with 95% confidence intervals.
- `bird-vs-birdbath-by-pool.png`: pool-faceted BIRD versus BIRDbath plots.

## Running the reusable endpoint

The input must already be ready for analysis. It must be a headerless,
tab-delimited file containing variant ID, BIRDbath theta, and then one BIRD
theta column per pool. The estimator does not recognize or filter on a
`P_reg` column.

```text
Rscript scripts/estimate_pool_effects.R READY_DATA.tsv OUTPUT_DIRECTORY
```

It can also be called from R with a ready data frame or matrix:

```r
source("scripts/estimate_pool_effects.R")
estimate_pool_effects(ready_data, output_dir)
```

The only non-base-R runtime dependency is `ggplot2`, used to produce the two
figures.

## Dataset-specific wrapper

`scripts/run_high_coverage_analysis.R` documents the preparation used for the
current high-coverage dataset: it reads the unfiltered source file, filters on
the final BIRDbath `P_reg` column, removes that column, and passes the prepared
data frame to the reusable endpoint. The source dataset itself is not included
in this bundle.

Within the original project layout, it can be run as:

```text
Rscript scripts/run_high_coverage_analysis.R 0.95 results/preg_gt_0.95
Rscript scripts/run_high_coverage_analysis.R 0.3 results/preg_gt_0.3
```

## Reproducibility check

From the bundle root:

```text
Rscript tests/final_pool_effect_model_smoke_test.R
```
