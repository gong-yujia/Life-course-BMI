# Part 2: Univariable Mendelian Randomization (UVMR)

This directory contains scripts for univariable Mendelian randomization (UVMR) analyses
based on harmonised exposure–outcome summary statistics. Multiple MR estimators are implemented
to assess the robustness of causal effect estimates.

All scripts assume that exposure–outcome harmonisation has been completed in advance
and that harmonised datasets are stored as `.rds` files.

---

## Directory structure

```text
part2_UVMR/
├── README.md
├── run_harmonised.R
├── UVMR_RAPS.R
├── UVMR_basic.R
├── UVMR_cml.R
├── UVMR_dIVW.R
├── UVMR_IVWrobust.R
└── run_Q_Egger.R

