#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(MendelianRandomization)
})

# ===================== Paths =====================
base_dir <- "."
out_dir  <- file.path(base_dir, "MR_divw_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===================== Locate outcome folders =====================
# Expected folder name pattern: *_outcome.chr1_22.filtered_ready
subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
subdirs <- subdirs[grepl("outcome\\.chr1_22\\.filtered_ready$", basename(subdirs))]

if (length(subdirs) == 0) {
  stop("No *_outcome.chr1_22.filtered_ready directories found. Please check base_dir.")
}
message("Number of outcome folders to process: ", length(subdirs))

# ===================== Main loop: per outcome =====================
for (folder in subdirs) {
  outcome_name <- str_extract(basename(folder), "^[^_]+") 
  message("\nProcessing outcome: ", outcome_name)

  rds_files <- list.files(folder, pattern = "_harmonised\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) {
    rds_files <- list.files(folder, pattern = "_harmonised\\.rds$", full.names = TRUE, recursive = TRUE)
  }
  message("Found harmonised RDS files: ", length(rds_files))
  if (length(rds_files) == 0) next

  res_list <- list()

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
    if (is.null(dat) || nrow(dat) == 0) {
      message("Empty data, skipping.")
      next
    }

    # Required columns
    need_cols <- c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
    if (!all(need_cols %in% names(dat))) {
      message("Missing required columns: ",
              paste(setdiff(need_cols, names(dat)), collapse = ", "))
      next
    }

    dt <- as.data.table(dat)

    # Use keep flags if present
    if ("mr_keep.exposure" %in% names(dt)) {
      dt <- dt[mr_keep.exposure == TRUE]
    } else if ("mr_keep" %in% names(dt)) {
      dt <- dt[mr_keep == TRUE]
    }
    if (nrow(dt) < 2) {
      message("Less than 2 SNPs after filtering, skipping.")
      next
    }

    # Convert to numeric; remove non-finite values; deduplicate by SNP
    for (cc in c("beta.exposure", "se.exposure", "beta.outcome", "se.outcome")) {
      suppressWarnings(dt[, (cc) := as.numeric(get(cc))])
    }
    dt <- dt[
      is.finite(beta.exposure) & is.finite(se.exposure) &
      is.finite(beta.outcome)  & is.finite(se.outcome)
    ]
    dt <- unique(dt, by = "SNP")
    if (nrow(dt) < 2) {
      message("Less than 2 SNPs after cleaning, skipping.")
      next
    }

    # Build MRInput
    mr_in <- mr_input(
      bx   = dt$beta.exposure,
      bxse = dt$se.exposure,
      by   = dt$beta.outcome,
      byse = dt$se.outcome,
      snps = dt$SNP
    )

    # Debiased IVW
    divw <- tryCatch(
      mr_divw(mr_in),
      error = function(e) {
        message("mr_divw failed: ", e$message)
        NULL
      }
    )
    if (is.null(divw)) next

    # Extract S4 slots
    nsnp <- length(dt$SNP)
    beta <- divw@Estimate
    se   <- divw@StdError
    lci  <- divw@CILower
    uci  <- divw@CIUpper
    pval <- divw@Pvalue

    res_row <- data.frame(
      exposure = exposure_name,
      nsnp     = nsnp,
      beta     = beta,
      se       = se,
      lci      = lci,
      uci      = uci,
      pval     = pval,
      OR       = exp(beta),
      OR_lci   = exp(lci),
      OR_uci   = exp(uci),
      method   = "Debiased IVW",
      stringsAsFactors = FALSE
    )
    res_list[[length(res_list) + 1]] <- res_row
  }

  # Write results for this outcome
  if (length(res_list) > 0) {
    final_df <- rbindlist(res_list, fill = TRUE)
    out_path <- file.path(out_dir, paste0("MR_divw_", outcome_name, ".csv"))
    fwrite(final_df, out_path)
    message("Saved: ", out_path, " (", nrow(final_df), " rows)")
  } else {
    message("No valid results for outcome: ", outcome_name)
  }
}

message("\nAll tasks completed. Results saved in ", out_dir)
