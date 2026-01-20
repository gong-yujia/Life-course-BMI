library(data.table)
library(survival)
library(broom)

# Outcomes and Time Points
outcomes <- c("AD", "CAD", "HF", "MS", "T2DM", "colorectal", "PD")
time_points <- c("1y5m", "1year", "2years", "3months", "3years", "5years", "6months", "6weeks", "7years", "8years", "8months", "birth")

# Path setup
prs_base_path <- "/path/to/PRS_result"  
output_base_path <- "/path/to/cox_result"  

# Read covariates
covar <- fread("/path/to/covariants_all.csv")
setnames(covar, old = "eid", new = "eid", skip_absent = TRUE)

for(time_point in time_points) {
  cat("Processing time point:", time_point, "\n")
  
  prs_file <- file.path(prs_base_path, time_point, "PRS_score_merged_phi1e2_allchr_zscore.profile")
  if (!file.exists(prs_file)) {
    cat("PRS file not found for time point:", time_point, "- skipping\n")
    next
  }
  
  prs <- fread(prs_file)
  setnames(prs, old = "FID", new = "eid", skip_absent = TRUE)
  
  output_path <- file.path(output_base_path, time_point)
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
  }
  
  for(outcome in outcomes) {
    cat("  Processing outcome:", outcome, "\n")
    
    pheno_file <- paste0("/path/to/", outcome, "_for_cox.csv")
    if (!file.exists(pheno_file)) {
      cat("  Phenotype file not found for outcome:", outcome, "- skipping\n")
      next
    }
    
    pheno <- fread(pheno_file)
    setnames(pheno, old = "eid", new = "eid", skip_absent = TRUE)
    
    dat <- merge(prs, pheno, by = "eid")
    dat <- merge(dat, covar, by = "eid")
    
    dat$PRS_group <- cut(dat$PRS_Z,
                         breaks = quantile(dat$PRS_Z, probs = c(0, 0.2, 0.8, 1), na.rm = TRUE),
                         labels = c("low", "intermediate", "high"),
                         include.lowest = TRUE, right = TRUE)
    dat$PRS_group <- factor(dat$PRS_group, levels = c("low", "intermediate", "high"))
    dat$Sex_num <- factor(dat$Sex_num, levels = c(1, 2), labels = c("Male", "Female"))
    dat$Genotype_batch_num <- factor(dat$Genotype_batch_num, levels = c(1, 2), labels = c("Axiom", "BiLEVE"))
    dat <- dat[follow_up_time > 0]
    dat$follow_up_time_years <- dat$follow_up_time / 365.25

    event_var <- paste0(outcome, "_event")
    if (!(event_var %in% names(dat))) {
      cat("  Event variable", event_var, "not found in data, skipping\n")
      next
    }
    
    dat <- dat[complete.cases(dat[, c(event_var, "follow_up_time_years", "PRS_group", "Age", "Sex_num", 
                                      "Genotype_batch_num", 
                                      "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), with=FALSE]), ]
    
    cox_formula <- as.formula(
      paste0("Surv(follow_up_time_years, ", event_var, ") ~ PRS_group + Age + Sex_num + Genotype_batch_num + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
    )
    
    cox_model <- coxph(cox_formula, data = dat)
    
    full_result_file <- file.path(output_path, paste0("PRS_cox_", outcome, ".csv"))
    fwrite(broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE), full_result_file)
    
    result_cox <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)
    prs_hr <- result_cox[grep("^PRS_group", result_cox$term), c("term", "estimate", "conf.low", "conf.high", "p.value")]
    hr_result_file <- file.path(output_path, paste0("PRS_HR_", outcome, ".csv"))
    fwrite(prs_hr, hr_result_file)
    
    cat("    Completed:", outcome, "for time point:", time_point, "\n")
  }
  cat("Completed all outcomes for time point:", time_point, "\n")
}

cat("All time points and outcomes processed successfully!\n")
cat("Results saved to:", output_base_path, "\n")
