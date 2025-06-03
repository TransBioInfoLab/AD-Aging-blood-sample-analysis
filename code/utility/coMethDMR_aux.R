#######################################################################################################
# =================================================================================================== #
# Auxillar function for coMethDMR
# =================================================================================================== #
#######################################################################################################
library(coMethDMR)       
library(GenomicRanges)   
library(tidyverse)         
# Load custom utility functions for annotation and bacon processing
source("~/TBL Dropbox/Wei Zhang/AD-Aging-blood-sample-analysis/code/utility/annotation_and_bacon.R")
# ===================================================================================================
## Function to extract genomic ranges from an array based on a specified genome.
# Parameters:
#   array: An input array (e.g., methylation array) for which genomic annotations are required.
#   genome: A string specifying the genome build (e.g., "hg19", "hg38") used for annotation.
#   dir.data.aux: (Optional) A directory path containing auxiliary data used by the annotation function.

get_genomic_range <- function(array, genome, dir.data.aux = NULL) {
  
  # Get the annotated genomic ranges using a custom function get_anno_gr
  anno.gr <- get_anno_gr(array = array, genome = genome, dir.data.aux = dir.data.aux)
  
  # Convert the GRanges object to a data frame and filter out sex chromosomes and NA values.
  anno_df <- data.frame(anno.gr)  %>%
    dplyr::filter(!(seqnames %in% c("chrX", "chrY"))) %>%
    dplyr::filter(!is.na(seqnames))
  
  # Calculate the minimal start and maximal end positions for each chromosome,
  # and create a new column representing the range as a string.
  range_df <- anno_df %>%
    dplyr::select("seqnames", "start", "end") %>%
    dplyr::group_by(.data$seqnames) %>%
    dplyr::summarise(start = min(.data$start), end = max(.data$end)) %>%
    dplyr::mutate(range = paste0(.data$start, "-", .data$end))
  
  # Return a GRanges object constructed from the summarized range information.
  GenomicRanges::GRanges(
    seqnames = range_df$seqnames,
    ranges = range_df$range
  )
}
# ===================================================================================================
## Function to map the output region back to the input region based on overlaps.
# Parameters:
#   input_region: A character vector of genomic regions in a string format (e.g., "chr1:1000-5000").
#   output_region: A character vector of genomic regions from an output analysis in a similar string format.

map_input_output_region <- function(input_region, output_region) {
  # Split the input_region string into components (chromosome, start, end)
  input_df <- data.frame(
    str_split(input_region, ":|-", simplify = T)
  )
  colnames(input_df) <- c("seqnames", "start", "end")
  input_gr <- makeGRangesFromDataFrame(input_df)
  
  # Similarly, split the output_region string into components
  output_df <- data.frame(
    str_split(output_region, ":|-", simplify = T)
  )
  colnames(output_df) <- c("seqnames", "start", "end")
  output_gr <- makeGRangesFromDataFrame(output_df)
  
  # Find overlapping regions between input and output GRanges objects.
  overlap_region <- findOverlaps(output_gr, input_gr, type = "within")
  
  # Return a data frame mapping the input region corresponding to the output region overlap.
  data.frame(input_region = input_region[overlap_region@to], output_region = output_region)
}

map_overlap_region <- function(input_region, output_region) {
  # Split the input_region string into components (chromosome, start, end)
  input_df <- data.frame(
    str_split(input_region, ":|-", simplify = T)
  )
  colnames(input_df) <- c("seqnames", "start", "end")
  input_gr <- makeGRangesFromDataFrame(input_df)
  
  # Similarly, split the output_region string into components
  output_df <- data.frame(
    str_split(output_region, ":|-", simplify = T)
  )
  colnames(output_df) <- c("seqnames", "start", "end")
  output_gr <- makeGRangesFromDataFrame(output_df)
  
  # Find overlapping regions between input and output GRanges objects.
  overlap_region <- findOverlaps(input_gr, output_gr)
  
  data.frame(input_idx = overlap_region@from,
             output_idx = overlap_region@to,
             input_region = input_region[overlap_region@from],
             output_region = output_region[overlap_region@to])

}
# ===================================================================================================
## Function to compute residuals from a linear model using methylation data.
# Parameters:
#   mat: A matrix of methylation beta values, where rows are CpG sites and columns are samples.
#   pheno: A data frame containing phenotype information; rows correspond to samples.
#   betaToM: Logical. If TRUE, converts beta values to M-values; otherwise, uses beta values as is.
#   nCore: Integer specifying the number of cores for parallel processing.

getResid <- function(mat, pheno, betaToM = T, nCore = 1){
  # Optionally transform beta values to M-values.
  if(betaToM){
    M.val <- log2(mat) - log2(1-mat)
  } else M.val <- mat
  
  # Check that the sample size in the methylation matrix and phenotype data match.
  if(ncol(mat) != nrow(pheno)){
    message("Sample size of methylation matrix and pheno data is not the same.")
    # Align phenotype rows with the methylation matrix columns.
    pheno <- pheno[match(colnames(mat), rownames(pheno)),]
  } 
  
  M.val <- as.matrix(M.val)
  
  # Set up parallel processing if more than one core is specified.
  if(nCore != 1){
    doParallel::registerDoParallel(nCore)
    para <- T
  } else para <- F
  
  # Apply a linear model to each row of the M-value matrix.
  resid.mat <- plyr::aaply(
    .data = M.val,
    .margins = 1,
    .fun = function(one){
      dat <- data.frame (M = one, pheno)
      tryCatch({
        suppressMessages({
          # Fit a linear model of methylation (M) against all phenotype variables.
          res <- lm(M ~ ., data = dat, na.action = na.exclude) %>% resid()
          return(res)
        })
      }, error = function(e) {message(e); return(NULL)})
    }, .parallel = para
  )
  
  return(resid.mat)
}
# ===================================================================================================
## Function to identify significant regions based on the median methylation value.
# Parameters:
#   coMeth.list: A list where each element contains a vector of CpG IDs belonging to a co-methylated region.
#   M: A matrix of methylation values (either beta or M-values) with rows as CpGs and columns as samples.
#   pheno: A data frame of phenotype information corresponding to the samples.
#   test_var: A string indicating the column name in pheno (or the analysis) that is being tested.
#   adjust_var: (Optional) A vector of strings specifying additional variables for adjustment in the linear model.
#   nCore: Integer specifying the number of cores for parallel processing.

getSigRegion <- function(coMeth.list, M, pheno, test_var, adjust_var = NULL, nCore = 1){
  
  # Set up parallel processing if more than one core is specified.
  if(nCore != 1){
    doParallel::registerDoParallel(nCore)
    para <- T
  } else para <- F
  
  # Create a linear regression formula with the test variable; include adjustment variables if provided.
  fo <-  paste0("median.M ~ ", test_var)
  if(!is.null(adjust_var)){
    fo <- paste0(fo, " + ", paste0(adjust_var, collapse = "+"))
  }
  formula <- as.formula(fo)
  
  # Loop over each co-methylation region in the list.
  coMeth.df <- plyr::ldply(
    .data = coMeth.list,
    .fun = function(CpGs){
      # Extract methylation values for CpGs in the region.
      CpGinRegion <- M[rownames(M) %in% CpGs,]
      # Calculate the median methylation for the region for each sample.
      median.region <- matrixStats::colMedians(CpGinRegion)
      dat <- data.frame(median.M = median.region, pheno)
      tryCatch({
        suppressMessages({
          # Fit a linear model using the defined formula.
          mod <- lm(formula, data = dat)
          # Extract coefficients from the model summary.
          lm_coef <- summary(mod)$coefficients %>% data.frame()
          
          # Filter coefficients related to the test variable and clean column names.
          coef_df <- lm_coef[grepl(test_var, rownames(lm_coef)),] %>% 
            janitor::clean_names()
          
          # Add additional information: CpGs in the region and their count.
          coef_df$CpGsinRegion <- paste0(CpGs, collapse = ",")
          coef_df$n <- length(CpGs)
          return(coef_df)
        })
      }, error = function(e) {message(e); return(NULL)})
    }, .parallel = para, .id = "region"
  )
  
  # Adjust the p-values for multiple testing using the FDR method.
  coMeth.df <- coMeth.df %>% mutate(fdr = p.adjust(pr_t, method = "fdr"))
  
  return(coMeth.df)
}
