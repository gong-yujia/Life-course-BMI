#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
  library(stringr)
})

# ===================== Paths =====================
base_dir   <- "."
output_dir <- file.path(base_dir, "MR_robust_results")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ===================== Locate outcome folders =====================
subfolders <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
subfolders <- subfolders[grepl("outcome\\.chr1_22\\.filtered_ready$", basename(subfolders))]

if (length(subfolders) == 0) {
  stop("No *_outcome.chr1_22.filtered_ready directories found. Please check base_dir.")
}

message("Number of outcome folders to process: ", length(subfolders))

# ===================== Main loop: per outcome =====================
for (folder in subfolders) {
  outcome_name <- str_extract(basename(folder), "^[^_]+") 
  message("\nProcessing outcome: ", outcome_name)

  # Harmonised RDS files under this outcome folder
  rds_files <- list.files(folder, pattern = "_harmonised\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) {
    rds_files <- list.files(folder, pattern = "_harmonised\\.rds$", full.names = TRUE, recursive = TRUE)
  }
  message("Found harmonised RDS files: ", length(rds_files))
  if (length(rds_files) == 0) next

  results_list <- list()

  for (rds_path in rds_files) {
    fname <- basename(rds_path)
    exposure_name <- str_remove(fname, "_harmonised\\.rds$")
    message("Analyzing: ", fname)

    dat <- tryCatch(
      readRDS(rds_path),
      error = function(e) {
        message("Failed to read RDS: ", e$message)
        NULL
      }
    )
    if (is.null(dat)) next

    # Compatibility: some objects may be wrapped as list(dat = <data>)
    if (!is.data.frame(dat) && !data.table::is.data.table(dat)) {
      if (is.list(dat) && !is.null(dat$dat)) dat <- dat$dat
    }
    if (is.null(dat) || nrow(dat) < 2) {
      message("Empty data or <2 rows, skipping.")
      next
    }

    # Required columns
    need_cols <- c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
    if (!all(need_cols %in% names(dat))) {
      message("Missing required columns: ",
              paste(setdiff(need_cols, names(dat)), collapse = ", "))
      next
    }

    # Cleaning: numeric conversion, remove non-finite values, deduplicate by SNP
    dt <- as.data.table(dat)[, ..need_cols]
    for (cc in setdiff(need_cols, "SNP")) {
      suppressWarnings(dt[, (cc) := as.numeric(get(cc))])
    }
    dt <- dt[
      is.finite(beta.exposure) & is.finite(se.exposure) &
      is.finite(beta.outcome)  & is.finite(se.outcome)
    ]
    dt <- unique(dt, by = "SNP")
    if (nrow(dt) < 2) {
      message("Less than 2 valid SNPs after cleaning, skipping.")
      next
    }

    # Build MRInput
    mr_input_obj <- mr_input(
      bx   = dt$beta.exposure,
      bxse = dt$se.exposure,
      by   = dt$beta.outcome,
      byse = dt$se.outcome,
      snps = dt$SNP
    )

    # Robust IVW
    res <- tryCatch(
      mr_ivw(mr_input_obj, model = "random", robust = TRUE),
      error = function(e) {
        message("mr_ivw failed: ", e$message)
        NULL
      }
    )
    if (is.null(res)) next

    nsnp_here <- nrow(dt)

    # Assemble result row
    res_row <- data.frame(
      exposure = exposure_name,
      nsnp     = nsnp_here,
      beta     = res@Estimate,
      se       = res@StdError,
      lci      = res@CILower,
      uci      = res@CIUpper,
      pval     = res@Pvalue,
      OR       = exp(res@Estimate),
      OR_lci   = exp(res@CILower),
      OR_uci   = exp(res@CIUpper),
      method   = "IVW (robust)",
      stringsAsFactors = FALSE
    )

    results_list[[length(results_list) + 1]] <- res_row
  }

  # Write results for this outcome
  if (length(results_list) > 0) {
    final_df <- rbindlist(results_list, fill = TRUE)
    out_path <- file.path(output_dir, paste0("MR_ivw_robust_", outcome_name, ".csv"))
    fwrite(final_df, out_path)
    message("Saved: ", out_path, " (", nrow(final_df), " rows)")
  } else {
    message("No valid results for outcome: ", outcome_name)
  }
}

message("\nAll tasks completed. Results saved in ", output_dir)
