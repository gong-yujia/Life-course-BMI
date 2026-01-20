#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(stringr)
})

# ----------- Set up paths -------------
base_dir   <- "."  
out_dir    <- file.path(base_dir, "MR_basic_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------- Get all outcome subdirectories-------------
subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
subdirs <- subdirs[grepl("outcome\\.chr1_22\\.filtered_ready$", basename(subdirs))]

if (length(subdirs) == 0) {
  stop("No *_outcome.chr1_22.filtered_ready directories found. Please check the base_dir.")
}

message("Number of outcome folders to process: ", length(subdirs))

# ----------- Main loop: process each outcome -----------------
for (folder in subdirs) {
  outcome_name <- str_extract(basename(folder), "^[^_]+")  
  message("\nProcessing outcome: ", outcome_name)

  rds_files <- list.files(folder, pattern = "_harmonised\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) {
    # If no files, check subdirectories recursively
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

    # Basic cleaning
    need_cols <- c("SNP","beta.exposure","se.exposure","beta.outcome","se.outcome")
    if (!all(need_cols %in% names(dat))) {
      message("Missing required columns: ", paste(setdiff(need_cols, names(dat)), collapse=", "))
      next
    }

    dt <- as.data.table(dat)
    for (cc in c("beta.exposure","se.exposure","beta.outcome","se.outcome")) {
      suppressWarnings(dt[, (cc) := as.numeric(get(cc))])
    }
    dt <- dt[is.finite(beta.exposure) & is.finite(se.exposure) &
             is.finite(beta.outcome)  & is.finite(se.outcome)]
    dt <- unique(dt, by = "SNP")

    if (nrow(dt) < 2) { message("Less than 2 valid SNPs, skipping"); next }

    # TwoSampleMR::mr requires harmonised standard columns; no renaming done here
    nsnp_here <- length(unique(dt$SNP))

    # Run basic MR methods
    mrdf <- tryCatch(
      {
        mr(dt,
           method_list = c(
             "mr_ivw"
           )
        )
      },
      error = function(e) { message("mr() error: ", e$message); NULL }
    )
    if (is.null(mrdf) || nrow(mrdf) == 0) { next }

    mrdf <- as.data.table(mrdf)
    mrdf[, `:=`(
      exposure_file = exposure_name,
      outcome_dir   = outcome_name,
      nsnp          = nsnp_here
    )]

    res_list[[length(res_list) + 1]] <- mrdf
  }

  # Save results for the current outcome
  if (length(res_list) > 0) {
    final_df <- rbindlist(res_list, fill = TRUE)
    out_path <- file.path(out_dir, paste0("MR_basic_", outcome_name, ".csv"))
    fwrite(final_df, out_path)
    message("Saved: ", out_path, " (", nrow(final_df), " rows)")
  } else {
    message("No valid results for: ", outcome_name)
  }
}

message("\nAll tasks completed. Results saved in ", out_dir)
