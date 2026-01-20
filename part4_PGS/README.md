# prscs_array.slurm
This directory contains a SLURM script for estimating SNP effect size weights using PRS-CS. The script is designed to run PRS-CS chromosome-wise (chr1–22) for a given time point using GWAS summary statistics.

# Cox_analysis.R
This repository contains the code for conducting Cox regression analyses on polygenic risk scores (PGS) across multiple time points. The analysis evaluates the association between PGS and various outcomes, adjusting for key covariates.

## Key Features:
- Adjusts for age, genotyping batch, and principal components.
- Stratifies PRS into low, intermediate, and high risk groups.
- Outputs HR, CI, and p-values for each outcome.
- **Gender-specific PGS for Prostate Cancer and Breast Cancer**:
  - For **prostate cancer** and **breast cancer**, gender-specific PGS files are used, with separate calculations for males and females.
  - **Sex is not included as a covariate** in the model for these specific cancers, as the PGS already accounts for the gender difference.
