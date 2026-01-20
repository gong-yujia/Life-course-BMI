library(data.table)
library(survival)
library(broom)

# ========== Preparation ==========
prs_base_dir   <- "/path/to/PRS_file"
pheno_base_dir <- "/path/to/phenotype/data"
covar_path     <- file.path(pheno_base_dir, "covariants_data", "covariants_all.csv")
lifestyle_path <- "/path/to/lifestyle/data"

timepoints <- c("7years", "8years")
outcomes   <- c("CAD", "HF", "T2DM")

covar     <- fread(covar_path)
lifestyle <- fread(lifestyle_path) 

# Output directory
outdir <- file.path(pheno_base_dir, "PRScs_result_re", "cox_result", "cox_result_interaction")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

all_results <- data.table()

# ========== 1. Main Loop ==========
for (tp in timepoints) {
  prs_path <- file.path(prs_base_dir, tp, "PRS_score_merged_phi1e2_allchr_zscore.profile")
  if (!file.exists(prs_path)) next
  prs <- fread(prs_path)
  setnames(prs, old = "FID", new = "eid", skip_absent = TRUE)

  for (outcome in outcomes) {
    message("Analyzing timepoint: ", tp, " | Outcome: ", outcome)
    pheno_path <- file.path(pheno_base_dir, paste0(outcome, "_for_cox.csv"))
    if (!file.exists(pheno_path)) next

    pheno <- fread(pheno_path)

    # Merging prs + pheno + covar + lifestyle
    dat <- Reduce(function(x, y) merge(x, y, by = "eid"), list(prs, pheno, covar))
    dat <- merge(dat, lifestyle, by = "eid")
    if (!("follow_up_time" %in% names(dat))) next
    dat <- dat[follow_up_time > 0]
    if (nrow(dat) == 0) next
    dat[, follow_up_time_years := follow_up_time / 365.25]

    # ===== PRS Grouping =====
    if (!("PRS_Z" %in% names(dat))) next
    qs <- quantile(dat$PRS_Z, probs = c(0, 0.2, 0.8, 1), na.rm = TRUE)
    if (length(unique(qs)) < 4) {
      pr <- frank(dat$PRS_Z, ties.method = "average", na.last = "keep") / sum(!is.na(dat$PRS_Z))
      dat[, PRS_group := cut(pr, breaks = c(0, 0.2, 0.8, 1),
                             labels = c("low", "intermediate", "high"),
                             include.lowest = TRUE)]
    } else {
      dat[, PRS_group := cut(PRS_Z, breaks = qs,
                             labels = c("low", "intermediate", "high"),
                             include.lowest = TRUE, right = TRUE)]
    }
    dat[, PRS_group := factor(PRS_group, levels = c("low", "intermediate", "high"))]

    # Factorizing: Gender / Education / Batch
    dat[, Sex_num            := factor(Sex_num, levels = c(1, 2), labels = c("Male", "Female"))]
    dat[, Education_group    := factor(Education_group, levels = c(1, 2), labels = c("Remainder", "Higher education"))]
    dat[, Genotype_batch_num := factor(Genotype_batch_num, levels = c(1, 2), labels = c("Axiom", "BiLEVE"))]

    # Lifestyle variables (alcohol removed)
    lifestyle_vars <- c("Physical_score", "Diet_score", "Sedentary_score", "Sleep_score", "Smoking_score")

    # Uniformly set 5 lifestyle variables as factors, reference level = healthy (1)
    for (vv in lifestyle_vars) {
      if (vv %in% names(dat)) {
        dat[[vv]] <- factor(dat[[vv]], levels = c(1, 0))  # 1 is the reference
      }
    }

    # Iterate over lifestyle variables for interaction, the other 4 are just covariates
    for (var in lifestyle_vars) {
      if (!var %in% names(dat)) next

      other_vars <- setdiff(lifestyle_vars, var)

      rhs <- paste(
        paste0("PRS_group * ", var),
        paste(other_vars, collapse = " + "),
        "Age + Sex_num + Genotype_batch_num + Education_group + Townsend_index",
        "PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10",
        sep = " + "
      )
      fml <- as.formula(paste0("Surv(follow_up_time_years, ", outcome, "_event) ~ ", rhs))

      fit <- try(coxph(fml, data = dat), silent = TRUE)
      if (inherits(fit, "try-error")) {
        message("⚠️ Skipping: Model did not converge or variables missing → ", tp, " / ", outcome, " / ", var)
        next
      }

      beta <- coef(fit)
      V    <- vcov(fit)
      cn   <- names(beta)

      # Capture the necessary coefficient names
      b_inter   <- "PRS_groupintermediate"
      b_high    <- "PRS_grouphigh"
      b_var0    <- grep(paste0("^", var, "[^:]*$"), cn, value = TRUE)  # Example: "Sleep_score0"
      b_ix_int  <- grep(paste0("^PRS_groupintermediate:", var, "[^:]*$|^", var, "[^:]*:PRS_groupintermediate$"), cn, value = TRUE)
      b_ix_high <- grep(paste0("^PRS_grouphigh:", var, "[^:]*$|^", var, "[^:]*:PRS_grouphigh$"), cn, value = TRUE)

      # If the lifestyle variable does not have two levels (e.g., all 1 or 0), skip
      if (length(b_var0) == 0) next

      # Construct 6 comparison groups (relative to baseline: low + healthy(1))
      groups <- list(
        "low + healthy(1)"          = c(),
        "intermediate + healthy(1)" = c(b_inter),
        "high + healthy(1)"         = c(b_high),
        "low + unhealthy(0)"        = c(b_var0),
        "intermediate + unhealthy(0)" = c(b_inter, b_var0, b_ix_int),
        "high + unhealthy(0)"         = c(b_high,  b_var0, b_ix_high)
      )

      # Calculate log(HR) and SE for each group (linear contrast)
      dt_out <- data.table()
      for (g in names(groups)) {
        terms <- groups[[g]]
        L <- setNames(rep(0, length(cn)), cn)
        if (length(terms)) L[terms] <- 1

        logHR <- sum(L * beta)
        SE    <- as.numeric(sqrt(t(L) %*% V %*% L))
        HR    <- exp(logHR)
        lower <- exp(logHR - 1.96 * SE)
        upper <- exp(logHR + 1.96 * SE)

        # For baseline (all 0), SE = 0, and the interval = 1
        if (all(L == 0)) {
          HR <- 1; lower <- 1; upper <- 1
        }

        dt_out <- rbind(dt_out, data.table(
          Timepoint = tp,
          Outcome   = outcome,
          lifestyle_var = var,
          Group6    = g,
          logHR     = logHR,
          HR        = HR,
          lower_CI  = lower,
          upper_CI  = upper
        ))
      }

      all_results <- rbindlist(list(all_results, dt_out), fill = TRUE)
    }
  }
}

fwrite(all_results, file = file.path(outdir, "PRS_Lifestyle_allOutcomes.csv"))
