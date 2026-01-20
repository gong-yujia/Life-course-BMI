#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
  library(tools)
})

# --------------------------- 0. Global settings ---------------------------
options(datatable.fread.datatable = FALSE)   
na_vec <- c(".", "NA", "")
msg <- function(...) message(format(Sys.time(), "[%H:%M:%S] "), ...)

# Keep only standard single nucleotide alleles
keep_snps <- function(dt, ea = "EA", nea = "NEA") {
  dt[
    nchar(get(ea))  == 1 & get(ea)  %in% c("A","T","C","G") &
      nchar(get(nea)) == 1 & get(nea) %in% c("A","T","C","G")
  ]
}

# --------------------------- 1. Parameter paths ---------------------------
# Replace with relative paths or placeholders
plink_path  <- "path/to/plink"  
bfile_path  <- "path/to/reference_LD/"  
p_threshold <- 1e-5
clump_kb    <- 10000
clump_r2    <- 0.001

# Current directory for exposure
exposure_root <- "."   
output_dir    <- file.path(exposure_root, "clump_result", "1e-5")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

outcome_dir   <- "path/to/outcome/" 
outcome_files <- list.files(outcome_dir, pattern = "outcome\\.chr1_22\\.filtered\\.tsv$", full.names = TRUE)
if (length(outcome_files) == 0) stop("No outcome files found")

# --------------------------- 2. Exposure files ---------------------------
# Replace with actual exposure data file path
exposure_files <- list.files(exposure_root, pattern = "^cleaned_.*\\.tsv$", full.names = TRUE)
if (!length(exposure_files)) stop("No cleaned_*.tsv files found")

# --------------------------- 3. Clumping ---------------------------
iv_summary <- data.frame(file = character(), n_IV = integer(), mean_F = numeric(), stringsAsFactors = FALSE)

for (file in exposure_files) {
  msg("Processing exposure: ", basename(file))
  
  df <- tryCatch({
    dat <- fread(file, na.strings = na_vec)
    setDT(dat)
    setnames(dat, names(dat), toupper(trimws(names(dat))))
    num_cols <- intersect(c("BETA","SE","P","N","EAF"), names(dat))
    if (length(num_cols)) dat[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
    dat
  }, error = function(e) { msg("Cannot read: ", e$message); NULL })
  if (is.null(df)) next
  
  required_cols <- c("RSID","BETA","SE","P","N","EA","NEA","EAF")
  if (!all(required_cols %in% names(df))) {
    msg("Missing columns, skipping: ", basename(file))
    iv_summary <- rbind(iv_summary, data.frame(file = basename(file), n_IV = 0, mean_F = NA))
    next
  }
  
  # P value threshold + F-statistic > 10 + standard alleles
  df <- df[df$P < p_threshold, ]
  df[, F_STAT := (BETA/SE)^2]
  df <- df[F_STAT > 10]
  df <- keep_snps(df)
  if (!nrow(df)) {
    msg("No valid SNPs after filtering")
    iv_summary <- rbind(iv_summary, data.frame(file = basename(file), n_IV = 0, mean_F = NA))
    next
  }
  
  # Format data for clumping
  dat <- format_data(as.data.frame(df), type = "exposure",
                     snp_col = "RSID", beta_col = "BETA", se_col = "SE",
                     eaf_col = "EAF", effect_allele_col = "EA", other_allele_col = "NEA",
                     pval_col = "P", samplesize_col = "N")
  setDT(dat)
  dat[, `:=`(rsid = SNP, pval = pval.exposure, id = id.exposure)]
  
  # Clumping (clump_p = p_threshold)
  clumped <- tryCatch(
    ld_clump_local(dat, bfile = bfile_path, plink_bin = plink_path,
                   clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = p_threshold),
    error = function(e) { msg("Clumping failed: ", e$message); NULL }
  )
  
  raw_name <- file_path_sans_ext(basename(file))
  prefix   <- sub("^cleaned_", "", raw_name)
  out_path <- file.path(output_dir, paste0("clumped_", prefix, ".txt"))
  
  if (!is.null(clumped) && nrow(clumped)) {
    setDT(clumped)
    req  <- c("SNP","P","A1","A2","BP","CHR")
    have <- intersect(req, names(clumped))
    if (!length(have)) {
      msg("Clumped output missing critical columns (SNP/P/BP/CHR); skipping save")
      iv_summary <- rbind(iv_summary, data.frame(file = basename(file), n_IV = 0, mean_F = NA))
      next
    }
    clumped <- clumped[, ..have]
    clumped <- clumped[!is.na(SNP)]
    clumped <- clumped[!duplicated(SNP)]
    if (!nrow(clumped)) {
      msg("No SNPs after cleaning clumped output")
      iv_summary <- rbind(iv_summary, data.frame(file = basename(file), n_IV = 0, mean_F = NA))
      next
    }
    
    fwrite(clumped, out_path, sep = "\t")
    msg("Saved → ", out_path)
    
    # mean_F
    iv_summary <- rbind(
      iv_summary,
      data.frame(
        file   = basename(file),
        n_IV   = nrow(clumped),
        mean_F = round(
          mean(
            merge(clumped, df[, .(RSID, F_STAT)], by.x = "SNP", by.y = "RSID",
                  all.x = FALSE, all.y = FALSE)$F_STAT,
            na.rm = TRUE
          ),
          2
        )
      )
    )
  } else {
    msg("No SNPs after clumping")
    iv_summary <- rbind(iv_summary, data.frame(file = basename(file), n_IV = 0, mean_F = NA))
  }
}

# Write IV summary
fwrite(iv_summary, file.path(output_dir, "clumping_IV_count_summary.txt"), sep = "\t")
msg("Clumping summary saved")

# ---------------- 4. Harmonise + Steiger ----------------
clumped_files <- list.files(output_dir, "^clumped_.*\\.txt$", full.names = TRUE)

steiger_summary <- data.table(
  outcome = character(), exposure = character(),
  n_before = integer(), n_after = integer(), n_removed = integer(),
  removed_snps = character()
)

for (outcome_file in outcome_files) {
  outcome_name <- file_path_sans_ext(basename(outcome_file))
  msg("Processing outcome: ", outcome_name)
  
  # outcome
  outcome_raw <- fread(outcome_file, na.strings = na_vec)
  setDT(outcome_raw)
  setnames(outcome_raw, names(outcome_raw), toupper(trimws(names(outcome_raw))))
  has_eaf <- "EAF" %in% names(outcome_raw)
  oc_num  <- intersect(c("BETA","SE","P","N","EAF"), names(outcome_raw))
  if (length(oc_num)) outcome_raw[, (oc_num) := lapply(.SD, as.numeric), .SDcols = oc_num]
  
  outcome_dat <- tryCatch(
    format_data(as.data.frame(outcome_raw), type = "outcome",
                snp_col = "RSID", beta_col = "BETA", se_col = "SE",
                effect_allele_col = "EA", other_allele_col = "NEA",
                pval_col = "P",
                eaf_col = if (has_eaf) "EAF" else NULL,
                samplesize_col = "N"),
    error = function(e) { msg("Outcome format_data failed: ", e$message); NULL }
  )
  if (is.null(outcome_dat)) next
  setDT(outcome_dat)  # <<< Key: convert to data.table for column selection
  
  outcome_dat <- outcome_dat[
    is.finite(beta.outcome) & is.finite(se.outcome) &
      !is.na(effect_allele.outcome) & !is.na(other_allele.outcome)
  ]
  if (!nrow(outcome_dat)) { msg("Outcome data is empty, skipping"); next }
  
  result_dir <- file.path(output_dir, paste0(outcome_name, "_ready"))
  dir.create(result_dir, showWarnings = FALSE)
  
  for (cf in clumped_files) {
    exp_name <- file_path_sans_ext(basename(cf))
    msg("Harmonising: ", exp_name)
    
    clumped_dt <- fread(cf, na.strings = na_vec)
    if (!("SNP" %in% names(clumped_dt))) { msg("Clumped file missing SNP column: ", cf); next }
    clumped_snps <- unique(na.omit(clumped_dt[["SNP"]]))
    if (!length(clumped_snps)) { msg("Clumped SNP list is empty"); next }
    
    cleaned_file <- file.path(exposure_root, paste0("cleaned_", sub("^clumped_", "", exp_name), ".tsv"))
    if (!file.exists(cleaned_file)) { msg("Exposure file not found: ", cleaned_file); next }
    
    # exposure
    exp_raw <- fread(cleaned_file, na.strings = na_vec)
    setDT(exp_raw)
    setnames(exp_raw, names(exp_raw), toupper(trimws(names(exp_raw))))
    ex_num <- intersect(c("BETA","SE","P","N","EAF"), names(exp_raw))
    if (length(ex_num)) exp_raw[, (ex_num) := lapply(.SD, as.numeric), .SDcols = ex_num]
    
    exp_dat <- tryCatch(
      format_data(as.data.frame(exp_raw), type = "exposure",
                  snp_col = "RSID", beta_col = "BETA", se_col = "SE",
                  eaf_col = "EAF", effect_allele_col = "EA", other_allele_col = "NEA",
                  pval_col = "P", samplesize_col = "N"),
      error = function(e) { msg("Exposure format_data failed: ", e$message); NULL }
    )
    if (is.null(exp_dat)) next
    setDT(exp_dat)  # <<< Key: convert to data.table
    
    exp_dat <- exp_dat[exp_dat$SNP %in% clumped_snps, ]
    exp_dat <- exp_dat[
      is.finite(beta.exposure) & is.finite(se.exposure) &
        !is.na(effect_allele.exposure) & !is.na(other_allele.exposure)
    ]
    if (nrow(exp_dat) == 0) { msg("No matching SNPs in exposure data"); next }
    
    # harmonise
    hdat <- tryCatch(
      harmonise_data(exp_dat, outcome_dat, action = 2),
      error = function(e) { msg("Harmonisation failed: ", e$message); NULL }
    )
    if (is.null(hdat)) next
    setDT(hdat)  # <<< Key: convert to data.table
    
    hdat <- hdat[
      is.finite(beta.exposure) & is.finite(se.exposure) &
        is.finite(beta.outcome)  & is.finite(se.outcome)
    ]
    if (!nrow(hdat)) { msg("Harmonised data is empty, skipping"); next }
    
    # Steiger (only when outcome contains EAF)
    if (has_eaf) {
      before_n   <- nrow(hdat)
      before_snp <- hdat$SNP
      
      hdat <- tryCatch(
        steiger_filtering(hdat),
        error = function(e) { msg("Steiger filtering failed: ", e$message); hdat }
      )
      
      after_n <- nrow(hdat)
      removed <- setdiff(before_snp, hdat$SNP)
      
      # Global summary
      steiger_summary <- rbind(
        steiger_summary,
        data.table(
          outcome = outcome_name,
          exposure = exp_name,
          n_before = before_n,
          n_after  = after_n,
          n_removed = before_n - after_n,
          removed_snps = if (length(removed)) paste(removed, collapse = ",") else ""
        ),
        fill = TRUE
      )
      if (length(removed)) {
        fwrite(data.table(SNP = removed),
               file.path(result_dir, paste0(exp_name, "_steiger_removed_snps.tsv")),
               sep = "\t")
      }
    } else {
      msg("Missing EAF → skipping Steiger filtering")
    }
    
    # Save
    saveRDS(hdat, file.path(result_dir, paste0(exp_name, "_harmonised.rds")))
    msg("Saved harmonised: ", exp_name)
  }
}

# Write Steiger global summary
if (nrow(steiger_summary)) {
  fwrite(steiger_summary, file.path(output_dir, "steiger_summary_all.tsv"), sep = "\t")
  msg("Steiger summary saved → steiger_summary_all.tsv")
} else {
  msg("No Steiger executed (possibly all outcomes lack EAF), skipping steiger_summary_all.tsv")
}

msg("All outcome analyses completed!")
