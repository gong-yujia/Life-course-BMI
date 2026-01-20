## FLOW-MR analysis (MVMR)

This directory contains the script used to perform multivariable Mendelian randomization based on the FLOW-MR framework.

### Input data preparation

Prior to running `run_FLOWMR.R`, input objects required by FLOW-MR were prepared, including:

- `Gamma_hat`: SNP–exposure association estimates across multiple exposures
- `Sd_hat`: corresponding standard errors
- `cor_mat`: correlation matrix between exposures

These objects were constructed from clumped genetic instruments using summary-level data.
Exposure–exposure correlation matrices (`cor_mat`) were estimated using the GRAPPLE framework, which accounts for sample overlap between GWAS summary statistics. Specifically, correlations
were derived from overlap coefficients estimated across summary-level data and were used to adjust for correlated estimation errors in the FLOW-MR model.

### FLOW-MR implementation

The script `run_FLOWMR.R` iterates over predefined time points and outcomes, and applies
the `BayesMediation` function implemented in the FLOWMR package to estimate direct and
indirect causal effects.

### Output

For each outcome and time point, the script generates a summary table containing posterior
estimates and uncertainty measures, saved as: flowmr_summary_<OUTCOME>_<TIME>.csv.

