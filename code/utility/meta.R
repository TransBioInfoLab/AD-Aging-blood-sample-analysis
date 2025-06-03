#######################################################################################################
# ===================================================================================================
# Function for meta analysis 
# ===================================================================================================
#######################################################################################################
library(metafor)
library(meta)
library(survival)
library(doParallel)
# ===================================================================================================
# meta analysis with metagen
# ===================================================================================================
# meta_wrapper: a wrapper function to perform meta-analysis across studies
# Parameters:
#   results_list: A list of data frames containing results from individual studies.
#   effect: The name of the column containing the effect estimates.
#   se: The name of the column containing the standard errors.
#   full_table: Logical. If TRUE, returns a full table with additional study-specific information.
#   sm: A string specifying the summary measure (default "HR" for hazard ratio, see metagen() in details).
#   test_var: The name of the identifier variable (default "cpg") to match across studies.
#   select_var: A vector of column names to retain in the full table output.

meta_wrapper <- function(results_list, effect, se, full_table = T, sm = "HR", test_var = "cpg",
                         select_var = c("seqnames", "start", "end", "width", "Islands.UCSC.Relation_to_Island",
                                        "UCSC_RefGene_Name", "UCSC_RefGene_Group", "GREAT_annotation")){
  
  # Extract all CpG identifiers across the results_list for the test variable.
  cpg <- unlist(purrr::map(results_list, ~.[[test_var]]))
  # Count how many times each CpG appears.
  cpg_tb <- table(cpg)
  # Select only those CpGs that appear in more than one study.
  cpg_select <- names(cpg_tb[cpg_tb > 1])
  
  # Register 30 parallel workers for the meta-analysis computation.
  doParallel::registerDoParallel(30)
  
  # For each common CpG, perform meta-analysis in parallel.
  results <- foreach::foreach(
    i = cpg_select,
    .errorhandling = "remove"
  ) %dopar% {
    # Extract results for the current CpG across studies.
    dat <- purrr::map(results_list, ~.[.[[test_var]] == i,])
    # Combine the data from different studies into one data frame.
    dat <- purrr::reduce(dat, rbind) 
    # Remove any NULL or empty elements.
    dat <- purrr::compact(dat)
    
    # Call the meta_fn function to perform the meta-analysis on the data for this CpG.
    meta_fn(dat = dat, 
            effect = effect,
            se = se,
            full_table = full_table,
            sm = sm,
            test_var = test_var,
            select_var = select_var)
  }
  
  # Stop the parallel cluster once processing is complete.
  doParallel::stopImplicitCluster()
  
  # Combine all individual meta-analysis results into one data.table.
  results <- data.table::rbindlist(results, fill = T)
  # Adjust the fixed-effect p-values for multiple testing (FDR correction) and insert before the column 'k'.
  results <- mutate(results, 
                    fdr = p.adjust(pVal.fixed, method = 'fdr'), .before = k)
  # Return the final meta-analysis results.
  results
}

# meta_fn: performs meta-analysis on a given dataset for a single CpG or region.
# Parameters:
#   dat: A data frame containing the results for a single CpG or region from multiple studies.
#   effect: The name of the column containing the effect estimates.
#   se: The name of the column containing the standard errors.
#   full_table: Logical. If TRUE, returns a full table with additional study-specific columns.
#   sm: A string specifying the summary measure (default "HR", see metagen() in details).
#   return_metagen: Logical. If TRUE, returns the meta-analysis object from metagen rather than the summary table.
#   test_var: The identifier variable (default "cpg") used to match the region.
#   select_var: A vector of column names to retain when reshaping the data for a full table.
#   ...: Additional arguments passed to the metagen function.

meta_fn <- function(dat, effect, se, full_table = T, sm = "HR", return_metagen = F, test_var = "cpg", 
                    select_var = c("seqnames", "start", "end", "width", "Islands.UCSC.Relation_to_Island",
                                   "UCSC_RefGene_Name", "UCSC_RefGene_Group", "GREAT_annotation"),
                    ...){
  
  # Perform meta-analysis using the metagen function from the meta package.
  # TE: effect estimate; seTE: standard error of the effect.
  f <- metagen(
    TE = dat[[effect]],
    seTE = dat[[se]],
    sm = sm,
    ...
  )
  
  # Create a data frame summarizing the meta-analysis results.
  mod_coef <- data.frame(
    cpg = unique(dat[[test_var]]), 
    estimate = f$TE.fixed,          
    se = f$seTE.fixed,              
    pVal.fixed = f$pval.fixed,     
    pVal.random = f$pval.random,   
    pVal.Q = f$pval.Q,             
    # Create a string indicating the direction of effect in each study using '+' or '-' signs.
    direction = paste0(ifelse(dat[[effect]] > 0, "+", "-"), collapse = ""),
    k = nrow(dat)                 # Number of studies/data points included in the analysis
  )
  # Rename the first column to match the test variable name.
  colnames(mod_coef)[1] <- test_var
  
  if(full_table) {
    # Remove extra columns that are not needed for the full table.
    dat$warning <- NULL
    dat$study <- NULL
    # Define the columns to be retained.
    col_names <- c(test_var, select_var)
    # Reshape the data to a wider format, so each study's information is in separate columns.
    dat_wider <- pivot_wider(
      dat, id_cols = all_of(col_names),
      names_from = study_id,
      names_glue = "{study_id}_{.value}",
      values_from = colnames(dat)[! colnames(dat) %in% c(col_names, "study_id")]
    )
    
    # Merge the meta-analysis summary with the reshaped data.
    mod_coef <- left_join(mod_coef, dat_wider, by = test_var)
  } 
  
  # Return the meta-analysis object if requested, else return the summarized coefficient table.
  if(return_metagen) {
    return(f)
  } else {
    return(mod_coef)
  }
  
}
