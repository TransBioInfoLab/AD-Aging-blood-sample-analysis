#######################################################################################################
# =================================================================================================== #
# Single CpG assocation test
# =================================================================================================== #
#######################################################################################################
# Load package
library(tidyverse)
library(survival)
library(janitor)
library(geepack)
library(gee)
library(lme4)
library(lmerTest)
# ===================================================================================================
# Single CpG test: cox regression
# ===================================================================================================
# Function: get_cox_coef
# ------------------------------------------------------------------------------
# Description:
#   Performs a Cox regression for a single CpG using survival data.
#
# Parameters:
#   cpg        : Numeric vector. Methylation values (or scaled values) for a single CpG.
#   pheno      : Data frame. Phenotype data containing survival information.
#   time_var   : Character. Name of the variable in 'pheno' representing time-to-event.
#   event_var  : Character. Name of the variable in 'pheno' representing the event indicator.
#   adjust_var : (Optional) Character vector. Additional covariates to adjust for in the regression.
#   ...        : Additional arguments passed to the coxph() function.
#
# Returns:
#   A data frame containing the regression coefficient(s) for the CpG along with a warning flag.
# ------------------------------------------------------------------------------
get_cox_coef <- function(cpg, pheno, time_var, event_var, adjust_var = NULL, ...) {
  
  # Construct the survival formula
  fo <- paste0("Surv(", time_var, ",", event_var, ") ~ cpg")
  if(!is.null(adjust_var)){
    fo <- paste0(fo, " + ", paste0(adjust_var, collapse = "+"))
  }
  formula <- as.formula(fo)
  
  data <- data.frame(cpg = cpg, pheno)
  
  # Check if any warnings occur during model fitting
  w <- check_warning(
    formula, 
    data,
    ...
  )
  
  # Fit the Cox regression model
  cox_mod <- coxph(
    formula,
    data = data,
    ...
  )
  
  # Extract summary coefficients
  cox_coef <- summary(cox_mod)$coefficients %>% data.frame()
  
  coef_df <- cox_coef[grepl("cpg", rownames(cox_coef)), ] %>% 
    clean_names()
  
  coef_df$warning <- w
  
  return(coef_df)
}
# ------------------------------------------------------------------------------
# Function: cox_coef
# ------------------------------------------------------------------------------
# Description:
#   Wrapper function to perform Cox regression across multiple CpGs in parallel.
#
# Parameters:
#   beta       : Matrix or data frame. CpG methylation values (rows = CpGs, columns = samples).
#   pheno      : Data frame. Phenotype data containing survival information.
#   time_var   : Character. Name of the time-to-event variable in 'pheno'.
#   event_var  : Character. Name of the event indicator variable in 'pheno'.
#   adjust_var : (Optional) Character vector. Covariates to adjust for in the regression.
#   parallel   : Logical. Whether to run the analysis in parallel (default: TRUE).
#   scale      : Logical. Whether to scale the CpG values (default: TRUE).
#   fdr_method : Character. Method for multiple testing correction (default: "fdr").
#   save       : Logical. Whether to save the results to a CSV file (default: FALSE).
#   dir.save   : (Optional) Character. Directory path where the CSV file will be saved.
#   cores      : Numeric. Number of CPU cores to use if running in parallel (default: 10).
#   prefix     : Character. Prefix for naming the output file (default: "Framingham").
#   ...        : Additional arguments passed to get_cox_coef().
#
# Returns:
#   A data frame of Cox regression results with FDR-adjusted p-values.
# ------------------------------------------------------------------------------
cox_coef <- function(beta, pheno, time_var, event_var, adjust_var = NULL, 
                     parallel = TRUE, 
                     scale = TRUE,
                     fdr_method = "fdr", 
                     save = FALSE, 
                     dir.save = NULL,
                     cores = 10,
                     prefix = "Framingham", ...) {
  
  if(parallel) doParallel::registerDoParallel(cores)
  
  results <- plyr::adply(
    beta,
    .margins = 1,
    .fun = function(cg) {
      if(scale) cg <- scale(cg)
      suppressWarnings({
        get_cox_coef(
          cpg = cg, 
          pheno = pheno,
          time_var = time_var,
          event_var = event_var,
          adjust_var = adjust_var,
          ...
        )
      })
    }, .id = "cpg", .parallel = parallel
  )
  
  colnames(results)[1] <- "cpg"
  # Add FDR-adjusted p-values to the results
  results <- add_fdr(results, method = fdr_method)
  
  if(save){
    write_csv(
      results,
      file.path(dir.save, paste0(prefix, "_single_cpg_cox_results.csv"))
    )
  }
  return(results)
}

# ------------------------------------------------------------------------------
# Function: check_warning
# ------------------------------------------------------------------------------
# Description:
#   Checks for warnings during model fitting with coxph().
#
# Parameters:
#   formula : Formula. The survival model formula.
#   data    : Data frame. Data used for model fitting.
#   ...     : Additional arguments passed to coxph().
#
# Returns:
#   Logical. TRUE if a warning occurred; FALSE otherwise.
# ------------------------------------------------------------------------------
check_warning <- function(formula, data, ...) {
  
  mod <- tryCatch({
    coxph(
      formula,
      data = data,
      ...
    )
  }, warning = function(w) return(NULL))
  
  if(length(class(mod) > 2)){
    w <- FALSE
  } else if (!is.null(mod) && class(mod) == "coxph"){
    w <- FALSE
  } else {
    w <- TRUE
  }
  
  return(w)
}

# ------------------------------------------------------------------------------
# Function: add_fdr
# ------------------------------------------------------------------------------
# Description:
#   Adds a column of FDR-adjusted p-values to the results data frame.
#
# Parameters:
#   results : Data frame. Results containing at least one column with p-values.
#   method  : Character. Multiple testing correction method (default: "fdr").
#
# Returns:
#   The input data frame with an additional column "fdr".
# ------------------------------------------------------------------------------
add_fdr <- function(results, method = "fdr") {
  
  pr <- results[[grep("^p", colnames(results), value = TRUE)]]
  results$fdr <- p.adjust(pr, method = method)
  
  return(results)
}

# ===================================================================================================
# Single CpG test: regression (linear model)
# ===================================================================================================
# Function: lm_test
# ------------------------------------------------------------------------------
# Description:
#   Performs linear regression on each CpG across samples.
#
# Parameters:
#   beta        : Matrix or data frame. CpG methylation values (rows = CpGs, columns = samples).
#   pheno       : Data frame. Phenotype data including the variable to be tested.
#   test_var    : Character. Name of the variable in 'pheno' to test against CpG methylation.
#   covariates  : (Optional) Character vector. Additional covariates to adjust for.
#   as_response : Logical. If TRUE, the CpG values are treated as the response variable (default: TRUE).
#   convert_to_M: Logical. If TRUE, transforms beta values to M-values using logit2() (default: TRUE).
#   parallel    : Logical. Whether to run the analysis in parallel (default: TRUE).
#   scale       : Logical. Whether to scale the CpG values before analysis (default: TRUE).
#   fdr_method  : Character. Method for multiple testing correction (default: "fdr").
#   save        : Logical. Whether to save the results to a CSV file (default: FALSE).
#   dir.save    : (Optional) Character. Directory path to save the output file.
#   cores       : Numeric. Number of CPU cores to use if running in parallel (default: 10).
#   prefix      : Character. Prefix for naming the output file (default: "Framingham").
#   ...         : Additional arguments passed to get_lm_coef().
#
# Returns:
#   A data frame with linear regression results and FDR-adjusted p-values.
# ------------------------------------------------------------------------------
lm_test <- function(beta, pheno, test_var, covariates = NULL, 
                    as_response = TRUE, convert_to_M = TRUE,
                    parallel = TRUE,
                    scale = TRUE,
                    fdr_method = "fdr", 
                    save = FALSE, 
                    dir.save = NULL,
                    cores = 10,
                    prefix = "Framingham", ...) {
  
  if(parallel) doParallel::registerDoParallel(cores)
  
  if(as_response) {
    res_var <- "cpg"
    test_var <- test_var
  } else {
    res_var <- test_var
    test_var <- "cpg"
  }
  
  if(convert_to_M) {
    mat <- minfi::logit2(beta)
  } else {
    mat <- beta
  }
  
  results <- plyr::adply(
    mat,
    .margins = 1,
    .fun = function(cg) {
      if(scale) cg <- scale(cg)
      data <- data.frame(cpg = cg, pheno)
      suppressWarnings({
        get_lm_coef(
          res_var = res_var, 
          test_var = test_var,
          data = data,
          covariates = covariates,
          ...
        )
      })
    }, .id = "cpg", .parallel = parallel
  )
  
  colnames(results)[1] <- "cpg"
  results <- add_fdr(results, method = fdr_method)
  
  if(save){
    write_csv(
      results,
      file.path(dir.save, paste0(prefix, "_single_cpg_test_results.csv"))
    )
  }
  return(results)
}

# ------------------------------------------------------------------------------
# Function: get_lm_coef
# ------------------------------------------------------------------------------
# Description:
#   Fits a linear regression model for a single CpG and extracts the coefficient for the test variable.
#
# Parameters:
#   res_var    : Character. Name of the response variable in the regression.
#   test_var   : Character. Name of the predictor variable to test.
#   data       : Data frame. Contains both the response and predictor variables, along with any covariates.
#   covariates : (Optional) Character vector. Additional covariates to include in the model.
#   ...        : Additional arguments passed to lm().
#
# Returns:
#   A data frame containing the regression coefficient(s) for the test variable.
# ------------------------------------------------------------------------------
get_lm_coef <- function(res_var, test_var, data, covariates = NULL, ...) {
  
  fo <- paste0(res_var, " ~ ", test_var)
  if(!is.null(covariates)){
    fo <- paste0(fo, " + ", paste0(covariates, collapse = "+"))
  }
  formula <- stats::as.formula(fo)
  
  tryCatch({
    lm_mod <- stats::lm(
      formula,
      data = data,
      ...
    )
    lm_coef <- data.frame(summary(lm_mod)$coefficients)
    coef_df <- lm_coef[grepl(test_var, rownames(lm_coef)), ]
    coef_df <- janitor::clean_names(coef_df)
  }, error = function(e){
    return(NULL)
  })
  
  coef_df
}

# ===================================================================================================
# Single CpG test: GEE
# ===================================================================================================
# Function: gee_test
# ------------------------------------------------------------------------------
# Description:
#   Performs a generalized estimating equations (GEE) analysis for a single CpG across samples.
#
# Parameters:
#   beta       : Matrix or data frame. CpG methylation values (rows = CpGs, columns = samples).
#   pheno      : Data frame. Phenotype data containing the variable to test.
#   test_var   : Character. Name of the predictor variable in 'pheno' to test against methylation.
#   group_var  : Character. Name of the grouping variable in 'pheno' for the GEE model.
#   covariates : (Optional) Character vector. Additional covariates to adjust for.
#   use_package: Character. Which GEE package to use ("geepack" or "gee", default "geepack").
#   as_response: Logical. If TRUE, treat methylation as the response variable (default: TRUE).
#   convert_to_M: Logical. If TRUE, converts beta values to M-values (default: TRUE).
#   parallel   : Logical. Whether to run the analysis in parallel (default: TRUE).
#   scale      : Logical. Whether to scale the CpG values (default: TRUE).
#   fdr_method : Character. Method for multiple testing correction (default: "fdr").
#   save       : Logical. Whether to save the results to a CSV file (default: FALSE).
#   dir.save   : (Optional) Character. Directory path to save the output file.
#   cores      : Numeric. Number of CPU cores to use if running in parallel (default: 10).
#   prefix     : Character. Prefix for naming the output file (default: "Framingham").
#   ...        : Additional arguments passed to get_gee_coef().
#
# Returns:
#   A data frame of GEE analysis results with FDR-adjusted p-values.
# ------------------------------------------------------------------------------
gee_test <- function(beta, pheno, test_var, 
                     group_var, covariates = NULL, 
                     use_package = "geepack",
                     as_response = TRUE, 
                     convert_to_M = TRUE,
                     parallel = TRUE,
                     scale = TRUE,
                     fdr_method = "fdr", 
                     save = FALSE, 
                     dir.save = NULL,
                     cores = 10,
                     prefix = "Framingham", ...) {
  
  if(parallel) doParallel::registerDoParallel(cores)
  
  if(as_response) {
    res_var <- "cpg"
    test_var <- test_var
  } else {
    res_var <- test_var
    test_var <- "cpg"
  }
  
  if(convert_to_M) {
    mat <- minfi::logit2(beta)
  } else {
    mat <- beta
  }
  
  results <- plyr::adply(
    mat,
    .margins = 1,
    .fun = function(cg) {
      if(scale) cg <- scale(cg)
      data <- data.frame(cpg = cg, pheno)
      suppressWarnings({
        get_gee_coef(
          res_var = res_var, 
          test_var = test_var,
          group_var = group_var,
          data = data,
          covariates = covariates,
          use_package = use_package,
          ...
        )
      })
    }, .id = "cpg", .parallel = parallel
  )
  
  colnames(results)[1] <- "cpg"
  results <- add_fdr(results, method = fdr_method)
  
  if(save){
    write_csv(
      results,
      file.path(dir.save, paste0(prefix, "_gee_single_cpg_test_results.csv"))
    )
  }
  return(results)
}

# ------------------------------------------------------------------------------
# Function: get_gee_coef
# ------------------------------------------------------------------------------
# Description:
#   Fits a GEE model for a single CpG and extracts the coefficient for the test variable.
#
# Parameters:
#   res_var    : Character. Name of the response variable.
#   test_var   : Character. Name of the predictor variable.
#   group_var  : Character. Name of the grouping variable.
#   data       : Data frame. Contains both the response and predictor variables, along with any covariates.
#   covariates : (Optional) Character vector. Additional covariates to adjust for.
#   use_package: Character. Package to use for GEE ("geepack" or "gee").
#   ...        : Additional arguments passed to the GEE function.
#
# Returns:
#   A data frame containing the GEE coefficient(s) for the test variable.
# ------------------------------------------------------------------------------
get_gee_coef <- function(res_var, test_var, group_var, data, covariates = NULL, use_package = "geepack", ...) {
  
  fo <- paste0(res_var, " ~ ", test_var)
  if(!is.null(covariates)){
    fo <- paste0(fo, " + ", paste0(covariates, collapse = "+"))
  }
  formula <- stats::as.formula(fo)
  data[[group_var]] <- factor(data[[group_var]])
  
  if(use_package == "geepack") {
    fn <- geepack::geeglm
  }
  if(use_package == "gee") {
    fn <- gee::gee
  }
  
  tryCatch({
    gee_mod <- fn(
      formula,
      data = data,
      id = data[[group_var]],
      ...
    )
    
    gee_coef <- data.frame(summary(gee_mod)$coefficients)
    coef_df <- gee_coef[grepl(test_var, rownames(gee_coef)), ]
    coef_df <- janitor::clean_names(coef_df)
  }, error = function(e){
    return(NULL)
  })
  
  coef_df
}

# ===================================================================================================
# Single CpG test: LME (Linear Mixed Effects)
# ===================================================================================================
# Function: lme_test
# ------------------------------------------------------------------------------
# Description:
#   Performs linear mixed-effects modeling for a single CpG across samples.
#
# Parameters:
#   beta          : Matrix or data frame. CpG methylation values (rows = CpGs, columns = samples).
#   pheno         : Data frame. Phenotype data containing the variable to test.
#   test_var      : Character. Name of the predictor variable in 'pheno' to test against methylation.
#   random_effect : (Optional) Character vector. Random effect variable(s) to include in the model.
#   covariates    : (Optional) Character vector. Additional covariates to adjust for.
#   as_response   : Logical. If TRUE, treat methylation as the response variable (default: TRUE).
#   convert_to_M  : Logical. If TRUE, convert beta values to M-values using logit2() (default: TRUE).
#   parallel      : Logical. Whether to run the analysis in parallel (default: TRUE).
#   scale         : Logical. Whether to scale the CpG values (default: TRUE).
#   fdr_method    : Character. Method for multiple testing correction (default: "fdr").
#   save          : Logical. Whether to save the results to a CSV file (default: FALSE).
#   dir.save      : (Optional) Character. Directory path to save the output file.
#   cores         : Numeric. Number of CPU cores to use if running in parallel (default: 10).
#   prefix        : Character. Prefix for naming the output file (default: "Framingham").
#   ...           : Additional arguments passed to get_lme_coef().
#
# Returns:
#   A data frame of LME analysis results with FDR-adjusted p-values.
# ------------------------------------------------------------------------------
lme_test <- function(beta, pheno, test_var, 
                     random_effect = NULL, covariates = NULL, 
                     as_response = TRUE, 
                     convert_to_M = TRUE,
                     parallel = TRUE,
                     scale = TRUE,
                     fdr_method = "fdr", 
                     save = FALSE, 
                     dir.save = NULL,
                     cores = 10,
                     prefix = "Framingham", ...) {
  
  if(parallel) doParallel::registerDoParallel(cores)
  
  if(as_response) {
    res_var <- "cpg"
    test_var <- test_var
  } else {
    res_var <- test_var
    test_var <- "cpg"
  }
  
  if(convert_to_M) {
    mat <- minfi::logit2(beta)
  } else {
    mat <- beta
  }
  
  results <- plyr::adply(
    mat,
    .margins = 1,
    .fun = function(cg) {
      if(scale) cg <- scale(cg)
      data <- data.frame(cpg = cg, pheno)
      suppressWarnings({
        get_lme_coef(
          res_var = res_var, 
          test_var = test_var,
          random_effect = random_effect,
          data = data,
          covariates = covariates,
          ...
        )
      })
    }, .id = "cpg", .parallel = parallel
  )
  
  colnames(results)[1] <- "cpg"
  results <- add_fdr(results, method = fdr_method)
  
  if(save){
    write_csv(
      results,
      file.path(dir.save, paste0(prefix, "_lme_single_cpg_test_results.csv"))
    )
  }
  return(results)
}

# ------------------------------------------------------------------------------
# Function: get_lme_coef
# ------------------------------------------------------------------------------
# Description:
#   Fits a linear mixed-effects model for a single CpG and extracts the coefficient for the test variable.
#
# Parameters:
#   res_var       : Character. Name of the response variable.
#   test_var      : Character. Name of the predictor variable to test.
#   data          : Data frame. Contains the response, predictor, covariates, and random effect variables.
#   random_effect : (Optional) Character vector. Variable(s) to be treated as random effects.
#   covariates    : (Optional) Character vector. Additional covariates to include in the model.
#   ...           : Additional arguments passed to lmerTest::lmer().
#
# Returns:
#   A data frame containing the LME coefficient(s) for the test variable.
# ------------------------------------------------------------------------------
get_lme_coef <- function(res_var, test_var, data, random_effect = NULL, covariates = NULL, ...) {
  
  fo <- paste0(res_var, " ~ ", test_var)
  if(!is.null(covariates)){
    fo <- paste0(fo, " + ", paste0(covariates, collapse = "+"))
  }
  if(!is.null(random_effect)) {
    fo <- paste0(fo, " + ", paste0(random_effect, collapse = "+"))
  }
  formula <- stats::as.formula(fo)
  
  tryCatch({
    lme_mod <- lmerTest::lmer(
      formula,
      data = data,
      ...
    )
    
    lme_coef <- data.frame(summary(lme_mod)$coefficients)
    coef_df <- lme_coef[grepl(test_var, rownames(lme_coef)), ]
    coef_df <- janitor::clean_names(coef_df)
  }, error = function(e){
    return(NULL)
  })
  
  coef_df
}

# ===================================================================================================
# Single CpG test: partial correlation
# ===================================================================================================
# Function: partial_corr
# ------------------------------------------------------------------------------
# Description:
#   Computes a partial correlation (default Spearman) between a specified dependent variable
#   and other variables while controlling for additional covariates.
#
# Parameters:
#   formula : Character or formula. A formula string specifying the variables (e.g., "y ~ x + z").
#   data    : Data frame. The dataset containing the variables.
#   dep_var : (Optional) Character. Name of the independent variable of interest.
#             If NULL, the second variable in the formula is used.
#   method  : Character. Correlation method to use ("spearman" by default).
#
# Returns:
#   A data frame summarizing the partial correlation coefficient, test statistic, p-value, and sample size for the specified variable.
# ------------------------------------------------------------------------------
partial_corr <- function(formula, data, dep_var = NULL, method = "spearman") {
  
  formula <- as.formula(formula)
  
  # Extract variable names from the formula
  vars <- all.vars(formula)
  
  res_var <- vars[1]
  if(is.null(dep_var)){
    dep_var <- vars[2] 
  }
  
  # Remove rows with missing values
  data_sel <- na.omit(data %>% dplyr::select(all_of(vars)))
  
  # Compute partial correlation using ppcor package
  par_corr <- ppcor::pcor(data_sel, method = method)
  
  par_summ <- data.frame(
    var = vars[-1],
    rho = par_corr$estimate[!vars %in% res_var, res_var],
    statistic = par_corr$statistic[!vars %in% res_var, res_var],
    p.value = par_corr$p.value[!vars %in% res_var, res_var],
    n = par_corr$n
  )
  
  par_ind <- par_summ %>% filter(var == dep_var)
  par_ind$var <- NULL
  
  return(par_ind)
}

# ===================================================================================================
# Adjust M values function
# ===================================================================================================
# Function: methyl_adj
# ------------------------------------------------------------------------------
# Description:
#   Adjusts methylation beta values by converting them to M-values, fitting a linear model
#   using specified adjustment variables, extracting the residuals, and optionally converting 
#   the residuals back to beta values.
#
# Parameters:
#   mat             : Matrix or data frame. Methylation beta values (rows = CpGs, columns = samples).
#   pheno           : Data frame. Phenotype or covariate data for adjustment.
#   adjust_var      : Character vector. Covariates to adjust for.
#   convert_to_M    : Logical. If TRUE, convert beta values to M-values using logit2() (default: TRUE).
#   return_to_beta  : Logical. If TRUE, convert the residual M-values back to beta values (default: TRUE).
#   parallel        : Logical. Whether to use parallel processing (default: TRUE).
#
# Returns:
#   A matrix of adjusted methylation values (either in M or beta scale).
# ------------------------------------------------------------------------------
methyl_adj <- function(mat, 
                       pheno, 
                       adjust_var, 
                       convert_to_M = TRUE, 
                       return_to_beta = TRUE,
                       parallel = TRUE) {
  
  if(convert_to_M){
    # Convert beta values to M-values
    mat <- minfi::logit2(mat)
  }
  mat <- as.matrix(mat)
  if(parallel) doParallel::registerDoParallel(10)
  
  resid_mat <- plyr::aaply(
    mat,
    .margins = 1,
    .fun = function(m) {
      dat <- data.frame(M = m, pheno)
      fo <- paste0("M ~ ", paste0(adjust_var, collapse = "+"))
      lm_mod <- lm(as.formula(fo), data = dat)
      r <- resid(lm_mod)
      return(r)
    }, .parallel = parallel
  )
  
  if(return_to_beta){
    # Convert adjusted M-values back to beta values
    resid_mat <- minfi::ilogit2(resid_mat)
  }
  
  return(resid_mat)
}

