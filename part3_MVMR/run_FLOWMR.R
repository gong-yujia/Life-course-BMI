#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(FLOWMR)
})

# ================== Path settings ==================
flowmr_input_dir  <- "path/to/flowMR_input"
flowmr_result_dir <- "path/to/flowMR_result"

# Create output root directory
dir.create(flowmr_result_dir, showWarnings = FALSE, recursive = TRUE)

msg <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n")

# ================== 1. Find all timepoint subdirectories ==================
time_dirs <- list.dirs(flowmr_input_dir, full.names = FALSE, recursive = FALSE)

if (length(time_dirs) == 0) {
  stop("No timepoint subdirectories found. Please check whether flowMR_input/ contains subfolders.")
}

# ================== 2. each timepoint ==================
for (time_tag in time_dirs) {

  msg("=== Processing timepoint: ", time_tag, " ===")

  this_input_dir  <- file.path(flowmr_input_dir, time_tag)

  this_result_dir <- file.path(flowmr_result_dir, time_tag)
  dir.create(this_result_dir, showWarnings = FALSE, recursive = TRUE)

  # Find all Gamma_Sd_*.rds under this timepoint directory
  rds_files <- list.files(
    this_input_dir,
    pattern = paste0("^Gamma_Sd_.*_", time_tag, "\\.rds$"),
    full.names = TRUE
  )

  if (length(rds_files) == 0) {
    msg("Timepoint ", time_tag, " has no Gamma_Sd_*_", time_tag, ".rds files. Skipping.")
    next
  }

  # ========== 3. each outcome ==========
  for (rds_in in rds_files) {

    base_name <- basename(rds_in)
    outcome_tag <- sub("^Gamma_Sd_", "", base_name)
    outcome_tag <- sub(paste0("_", time_tag, "\\.rds$"), "", outcome_tag)

    msg(">>> FLOW-MR: outcome = ", outcome_tag, " | timepoint = ", time_tag)

    # === Read object ===
    obj <- readRDS(rds_in)

    Gamma_hat <- obj$Gamma_hat   # 3 x K
    Sd_hat    <- obj$Sd_hat      # 3 x K
    cor_mat   <- obj$cor_mat     # 3 x 3
    SNPs      <- obj$snp         # SNP ID list

    # Print checks
    cat("Gamma_hat dim =", dim(Gamma_hat),
        "; Sd_hat dim =", dim(Sd_hat), "\n")
    cat("First few cols of Gamma_hat:\n")
    print(Gamma_hat[, 1:min(5, ncol(Gamma_hat)), drop = FALSE])
    cat("cor_mat:\n")
    print(cor_mat)

    # === Run FLOW-MR ===
    set.seed(25110101)
    fit <- FLOWMR::BayesMediation(
      Gamma_hat,
      Sd_hat,
      cor = cor_mat,
      inv = TRUE
    )

    print(fit$summary)

    # === Format summary table ===
    sum_df <- as.data.frame(fit$summary)
    sum_df$param <- rownames(fit$summary)
    sum_df <- sum_df[, c("param", setdiff(colnames(sum_df), "param"))]

    # === Write CSV: flowmr_summary_<OUTCOME>_<TIME>.csv ===
    out_csv <- file.path(
      this_result_dir,
      paste0("flowmr_summary_", outcome_tag, "_", time_tag, ".csv")
    )

    write.csv(sum_df, file = out_csv, row.names = FALSE)

    cat("Results saved to: ", out_csv, "\n\n")
  }

  msg("=== Finished timepoint: ", time_tag, " ===\n")
}

msg("All timepoints/outcomes completed.")

