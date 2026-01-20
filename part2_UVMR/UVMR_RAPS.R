#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(mr.raps)
})

# ========= Path setup =========
base_dir   <- "."  
out_dir    <- file.path(base_dir, "MR_raps_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ========= Find all outcome folders =========
subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
subdirs <- subdirs[grepl("outcome\\.chr1_22\\.filtered_ready$", basename(subdirs))]

if (length(subdirs) == 0) {
  stop("No *_outcome.chr1_22.filtered_ready directories found. Please check base_dir.")
}

message("Number of outcome folders to process: ", length(subdirs))

# ========= Main loop: Process each outcome =========
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

    dat <- tryCatch(readRDS(rds_path), error = function(e) { message("Failed to read: ", e$message); NULL })
    if (is.null(dat)) next

    # Some objects might be wrapped in list$dat
    if (!is.data.frame(dat) && !data.table::is.data.table(dat)) {
      if (is.list(dat) && !is.null(dat$dat)) dat <- dat$dat
    }
    if (is.null(dat) || nrow(dat) == 0) { message("Data is empty, skipping"); next }

    # Keep standard columns for harmonised data from TwoSampleMR
    need_cols <- c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome",
                   "mr_keep.exposure","mr_keep")
    miss <- setdiff(c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome"), names(dat))
    if (length(miss) > 0) {
      message("Missing required columns: ", paste(miss, collapse=", "), ", skipping")
      next
    }

    dt <- as.data.table(dat)

    # Handle 'keep' flags
    if ("mr_keep.exposure" %in% names(dt)) {
      dt <- dt[mr_keep.exposure == TRUE]
    } else if ("mr_keep" %in% names(dt)) {
      dt <- dt[mr_keep == TRUE]
    }
    if (nrow(dt) < 2) { message("After filtering, less than 2 SNPs, skipping"); next }

    # Convert to numeric, remove non-finite values, deduplicate SNPs
    for (cc in c("beta.exposure","se.exposure","beta.outcome","se.outcome")) {
      suppressWarnings(dt[, (cc) := as.numeric(get(cc))])
    }
    dt <- dt[is.finite(beta.exposure) & is.finite(se.exposure) &
             is.finite(beta.outcome)  & is.finite(se.outcome)]
    dt <- unique(dt, by = "SNP")
    if (nrow(dt) < 2) { message("After cleaning, less than 2 SNPs, skipping"); next }

    bx <- dt$beta.exposure
    by <- dt$beta.outcome
    sx <- dt$se.exposure
    sy <- dt$se.outcome

    # Try overdispersed robust, fall back to simple if it fails
    model_used <- NA_character_
    fit <- tryCatch(
      {
        model_used <- "RAPS robust (overdispersed, Huber)"
        mr.raps.overdispersed.robust(bx, by, sx, sy, diagnostics = FALSE, loss.function = "huber")
      },
      warning = function(w) {
        if (grepl("overdispersion parameter is very small", conditionMessage(w))) {
          message("Overdispersed parameter too small, falling back to simple")
          return(NULL)
        } else {
          warning(w); return(NULL)
        }
      },
      error = function(e) { message("Robust model failed: ", e$message); NULL }
    )

    if (is.null(fit)) {
      fit <- tryCatch(
        { model_used <- "RAPS simple"; mr.raps.simple(bx, by, sx, sy) },
        error = function(e) { message("Simple model also failed: ", e$message); NULL }
      )
    }
    if (is.null(fit)) next

    beta <- as.numeric(fit$beta.hat)
    se   <- as.numeric(fit$beta.se)
    z    <- beta / se
    pval <- 2 * pnorm(-abs(z))
    lci  <- beta - 1.96 * se
    uci  <- beta + 1.96 * se

    res_row <- data.frame(
      exposure = exposure_name,
      nsnp     = nrow(dt),
      beta     = beta,
      se       = se,
      lci      = lci,
      uci      = uci,
      pval     = pval,
      OR       = exp(beta),
      OR_lci   = exp(lci),
      OR_uci   = exp(uci),
      method   = model_used,
      stringsAsFactors = FALSE
    )
    res_list[[length(res_list) + 1]] <- res_row
  }

  # Save results for the current outcome
  if (length(res_list) > 0) {
    final_df <- rbindlist(res_list, fill = TRUE)
    out_path <- file.path(out_dir, paste0("MR_raps_", outcome_name, ".csv"))
    fwrite(final_df, out_path)
    message("Saved: ", out_path, " (", nrow(final_df), " rows)")
  } else {
    message("No valid results for: ", outcome_name)
  }
}

message("\nAll tasks completed. Results saved in ", out_dir)
