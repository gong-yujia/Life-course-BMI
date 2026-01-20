#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(stringr)
})

# ===================== Config =====================
min_snps_intercept <- 3  
root <- "."             

# I2GX (for reference / QC)
I2GX <- function(bx, bxse) {
  v <- stats::var(bx, na.rm = TRUE)
  e <- mean(bxse^2, na.rm = TRUE)
  v / (v + e)
}

# ===================== Collect all harmonised RDS files =====================
outcome_dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
outcome_dirs <- outcome_dirs[grepl("_outcome\\.chr1_22\\.filtered_ready$", basename(outcome_dirs))]

files <- unlist(lapply(outcome_dirs, function(d) {
  list.files(d, pattern = "_harmonised\\.rds$", recursive = TRUE, full.names = TRUE)
}))
files <- files[file.exists(files)]

if (!length(files)) stop("No *_harmonised.rds files found.")

# ===================== Output directory =====================
outdir <- file.path(root, "MR_Q_Egger")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

res_list <- list()

# ===================== Main loop =====================
for (f in files) {
  dat <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(dat)) next
  dat <- as.data.frame(dat)
                  
  if ("remove" %in% names(dat)) {
    dat <- dat[is.na(dat$remove) | dat$remove == FALSE, ]
  }

  if ("mr_keep" %in% names(dat)) {
    dat <- dat[dat$mr_keep %in% TRUE, ]
  } else if (all(c("mr_keep.exposure", "mr_keep.outcome") %in% names(dat))) {
    dat <- dat[dat$mr_keep.exposure %in% TRUE & dat$mr_keep.outcome %in% TRUE, ]
  }

  if (all(c("palindromic", "ambiguous") %in% names(dat))) {
    dat <- dat[!(dat$palindromic & dat$ambiguous), ]
  }

  need <- c("SNP", "beta.exposure", "se.exposure",
            "beta.outcome", "se.outcome", "exposure", "outcome")
  dat <- dat[stats::complete.cases(dat[, intersect(need, names(dat))]), ]

  nsnps <- length(unique(dat$SNP))
  expo  <- if ("exposure" %in% names(dat)) paste(unique(dat$exposure), collapse = ";") else NA_character_
  outc  <- if ("outcome"  %in% names(dat)) paste(unique(dat$outcome),  collapse = ";") else NA_character_

  Q_ivw <- Qp_ivw <- Q_egger <- Qp_egger <- NA_real_
  eg_int <- eg_se <- eg_p <- NA_real_
  I2 <- if (nsnps > 1) I2GX(dat$beta.exposure, dat$se.exposure) else NA_real_

  if (nsnps >= min_snps_intercept) {
    # Cochran's Q (IVW and Egger)
    het <- tryCatch(
      mr_heterogeneity(dat, method_list = c("mr_ivw", "mr_egger_regression")),
      error = function(e) NULL
    )
    if (!is.null(het)) {
      if (any(grepl("Inverse variance weighted", het$method, fixed = TRUE))) {
        idx <- grepl("Inverse variance weighted", het$method, fixed = TRUE)
        Q_ivw  <- het$Q[idx][1]
        Qp_ivw <- het$Q_pval[idx][1]
      }
      if (any(grepl("MR Egger", het$method))) {
        idx <- grepl("MR Egger", het$method)
        Q_egger  <- het$Q[idx][1]
        Qp_egger <- het$Q_pval[idx][1]
      }
    }

    # MR-Egger intercept test
    pleio <- tryCatch(mr_pleiotropy_test(dat), error = function(e) NULL)
    if (!is.null(pleio) && nrow(pleio) > 0) {
      eg_int <- pleio$egger_intercept[1]
      eg_se  <- pleio$se[1]
      eg_p   <- pleio$pval[1]
    }
  }

  res_list[[length(res_list) + 1]] <- data.frame(
    file = f,
    exposure = expo,
    outcome = outc,
    nsnps = nsnps,
    I2GX = I2,
    Q_ivw = Q_ivw,
    Qp_ivw = Qp_ivw,
    Q_egger = Q_egger,
    Qp_egger = Qp_egger,
    Egger_intercept = eg_int,
    Egger_intercept_se = eg_se,
    Egger_intercept_p = eg_p,
    stringsAsFactors = FALSE
  )
}

res <- rbindlist(res_list, fill = TRUE)

# ===== BH-FDR correction within outcome =====
if (nrow(res)) {
  setDT(res)

  res[, n_tests_outcome := sum(!is.na(Egger_intercept_p)), by = outcome]

  res[, Egger_intercept_FDR := {
    p <- Egger_intercept_p
    q <- rep(NA_real_, .N)
    if (sum(!is.na(p)) > 0) q[!is.na(p)] <- p.adjust(p[!is.na(p)], "BH")
    q
  }, by = outcome]

  res[, Qp_ivw_FDR := {
    p <- Qp_ivw
    q <- rep(NA_real_, .N)
    if (sum(!is.na(p)) > 0) q[!is.na(p)] <- p.adjust(p[!is.na(p)], "BH")
    q
  }, by = outcome]

  res[, Qp_egger_FDR := {
    p <- Qp_egger
    q <- rep(NA_real_, .N)
    if (sum(!is.na(p)) > 0) q[!is.na(p)] <- p.adjust(p[!is.na(p)], "BH")
    q
  }, by = outcome]
}

# ===== Save results =====
setcolorder(
  res,
  c("file", "outcome", "exposure", "nsnps", "I2GX",
    "Q_ivw", "Qp_ivw", "Qp_ivw_FDR",
    "Q_egger", "Qp_egger", "Qp_egger_FDR",
    "Egger_intercept", "Egger_intercept_se", "Egger_intercept_p",
    "Egger_intercept_FDR", "n_tests_outcome")
)

outfile <- file.path(outdir, "Q_and_Egger_intercept_summary.csv")
fwrite(res, outfile)
message("Output written to: ", outfile, " (rows: ", nrow(res), ")")
