#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(MRcML)
})

# ===================== Paths =====================
base_dir <- "."
out_dir  <- file.path(base_dir, "MR_cml_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===================== Locate outcome folders =====================
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

    # Required columns (from TwoSampleMR harmonise_data output)
    need_base <- c("SNP", "beta.exposure", "se.exposure", "beta.outcome", "se.outcome")
    if (!all(need_base %in% names(dat))) {
      message("Missing required columns: ",
              paste(setdiff(need_base, names(dat)), collapse = ", "))
      next
    }

    dt <- as.data.table(dat)

    # Use keep flags if present
    if ("mr_keep.exposure" %in% names(dt)) {
      dt <- dt[mr_keep.exposure == TRUE]
    } else if ("mr_keep" %in% names(dt)) {
      dt <- dt[mr_keep == TRUE]
    }

    if (nrow(dt) < 3) {
      message("Less than 3 SNPs after filtering, skipping.")
      next
    }

    # Convert to numeric where applicable; remove non-finite values; deduplicate SNPs
    for (cc in c(
      "beta.exposure", "se.exposure", "beta.outcome", "se.outcome",
      "samplesize.exposure", "ncase.exposure", "ncontrol.exposure"
    )) {
      if (cc %in% names(dt)) suppressWarnings(dt[, (cc) := as.numeric(get(cc))])
    }

    dt <- dt[
      is.finite(beta.exposure) & is.finite(se.exposure) &
      is.finite(beta.outcome)  & is.finite(se.outcome)
    ]
    dt <- unique(dt, by = "SNP")

    if (nrow(dt) < 3) {
      message("Less than 3 SNPs after cleaning, skipping.")
      next
    }

    n <- NA_real_
    if ("samplesize.exposure" %in% names(dt)) {
      n <- suppressWarnings(median(dt$samplesize.exposure, na.rm = TRUE))
    }
    if (!is.finite(n) || is.na(n) || n <= 0) {
      if (all(c("ncase.exposure", "ncontrol.exposure") %in% names(dt))) {
        n <- suppressWarnings(median(dt$ncase.exposure + dt$ncontrol.exposure, na.rm = TRUE))
      }
    }
    if (!is.finite(n) || is.na(n) || n <= 0) {
      message("Cannot determine sample size n, skipping.")
      next
    }

    bx <- dt$beta.exposure
    by <- dt$beta.outcome
    sx <- dt$se.exposure
    sy <- dt$se.outcome
    nsnp <- nrow(dt)

    # Run MR-cML
    fit <- tryCatch(
      MRcML::mr_cML(
        b_exp = bx, b_out = by,
        se_exp = sx, se_out = sy,
        n = n,
        random_start = 100,
        random_seed = 123
      ),
      error = function(e) {
        message("mr_cML failed: ", e$message)
        NULL
      }
    )
    if (is.null(fit)) next

    # Extract BIC-based estimate
    theta <- as.numeric(fit$BIC_theta)
    se    <- as.numeric(fit$BIC_se)
    pval  <- as.numeric(fit$BIC_p)

    if (!is.finite(theta) || !is.finite(se) || se <= 0) {
      message("Invalid BIC estimate (theta/se), skipping.")
      next
    }

    lci <- theta - 1.96 * se
    uci <- theta + 1.96 * se

    ninv_bic <- tryCatch(length(fit$BIC_invalid), error = function(e) NA_integer_)

    df <- data.frame(
      exposure = exposure_name,
      nsnp     = nsnp,
      beta     = theta,
      se       = se,
      lci      = lci,
      uci      = uci,
      pval     = pval,
      OR       = exp(theta),
      OR_lci   = exp(lci),
      OR_uci   = exp(uci),
      method   = "cML (BIC)",
      n_invalid_BIC = ninv_bic,
      stringsAsFactors = FALSE
    )

    res_list[[length(res_list) + 1]] <- df
  }

  # Write results for this outcome
  if (length(res_list) > 0) {
    final_df <- rbindlist(res_list, fill = TRUE)
    out_path <- file.path(out_dir, paste0("MR_cml_BIC_", outcome_name, ".csv"))
    fwrite(final_df, out_path)
    message("Saved: ", out_path, " (", nrow(final_df), " rows)")
  } else {
    message("No valid results for outcome: ", outcome_name)
  }
}

message("\nAll tasks completed. Results saved in ", out_dir)
