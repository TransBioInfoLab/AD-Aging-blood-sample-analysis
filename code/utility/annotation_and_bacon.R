#######################################################################################################
# =================================================================================================== #
# Function for annotation and inflation adjusted 
# =================================================================================================== #
#######################################################################################################
library(SummarizedExperiment)
library(tidyverse)
library(bacon)
library(GWASTools)
library(minfi)
library(rGREAT)
# ===================================================================================================
# Annotation Single CpGs
# ===================================================================================================
# Function: get_anno_gr
# ------------------------------------------------------------------------------
# Description:
#   Retrieves annotation information as a GRanges object for the specified methylation array and genome.
#
# Parameters:
#   array        : Character. An input array (e.g., methylation array) for which genomic annotations 
#                  are required. Supported options include "HM450", "EPICv1", "EPIC", and "EPICv2".
#   genome       : Character. A string specifying the genome build used for annotation, such as "hg19" 
#                  or "hg38".
#   dir.data.aux : (Optional) Character. A directory path containing auxiliary data used by the 
#                  annotation function, particularly for reading external manifest files when using 
#                  genome "hg38".
#
# Returns:
#   A GRanges object containing the genomic annotation data.

get_anno_gr <- function(array = "HM450", 
                        genome = "hg19",
                        dir.data.aux = NULL) {

  if(array == "HM450"){
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  } 
  if(array %in% c("EPICv1", "EPIC")){
    if(genome == "hg19") {
      library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
      anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    }
    if(genome == "hg38") {
      anno <- read_csv(
        file.path(dir.data.aux, "infinium-methylationepic-v-1-0-b5-manifest-file.csv"),
        show_col_types = F, 
        skip = 7
      )
    }
  }
  if(array == "EPICv2") {
    library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
    anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38)
  }
  
  if(array %in% c("EPICv1", "EPIC") & genome == "hg38") {
    anno.gr <- anno %>% makeGRangesFromDataFrame(
      seqnames.field = "CHR_hg38",
      start.field = "Start_hg38", end.field = "End_hg38", 
      strand.field = "Strand_hg38",
      keep.extra.columns = T,
      na.rm = T
    )
  } else {
    anno.gr <- anno %>% makeGRangesFromDataFrame(
      start.field = "pos", end.field = "pos", keep.extra.columns = T
    )
  }
  anno.gr
} 
# ------------------------------------------------------------------------------
# Function: annotate_results
# ------------------------------------------------------------------------------
# Description:
#   Annotates a CpG results data frame by merging genomic annotation (from get_anno_gr) with 
#   additional GREAT analysis results.
#
# Parameters:
#   result       : Data frame. A table containing CpG IDs and associated statistics to annotate.
#   array        : Character. The methylation array type (e.g., "HM450", "EPICv1", "EPIC").
#   genome       : Character. The genome version ("hg19" or "hg38").
#   dir.data.aux : Character. Directory path for auxiliary data files containing precomputed annotations.
#   save         : Logical. Whether to save the annotated result as a CSV file.
#   dir.save     : Character. Directory path where the annotated result file will be saved.
#   prefix       : Character. A prefix used for naming the output file.
#
# Returns:
#   A data frame with additional columns for genomic annotation and GREAT analysis.

annotate_results <- function(result, 
                             array = "HM450", 
                             genome = "hg19",
                             dir.data.aux = NULL, 
                             save = T,
                             dir.save = NULL,
                             prefix = "Framingham"){
  
  # load(file.path(dir.data.aux,"E073_15_coreMarks_segments.rda"))
  # load(file.path(dir.data.aux,"meta_analysis_cpgs.rda"))
  
  anno.gr <- get_anno_gr(array = array, genome = genome, dir.data.aux = dir.data.aux)
  anno_df <- as.data.frame(anno.gr)
  
  if(genome == "hg19") {

    if(array == "HM450"){
      load(file.path(dir.data.aux,"great_HM450_array_annotation.rda"))
    }
    if(array %in% c("EPICv1", "EPIC")){
      load(file.path(dir.data.aux,"great_EPIC_array_annotation.rda"))
    }
    result <- cbind(
      result,
      anno_df[result$cpg,c("seqnames", "start", "end", "width", "Relation_to_Island", "UCSC_RefGene_Name", "UCSC_RefGene_Group")] 
    )
    result <- dplyr::left_join(result, great, by = c("seqnames","start","end","cpg"))
  }
  
  # Add annotation

  if(genome == "hg38"){
    
    if(array %in% c("EPICv1", "EPIC")) {
      load(file.path(dir.data.aux,"great_EPIC_array_annotation.rda"))
      anno_df <-  anno_df %>% 
        mutate(cpg = Name) %>%
        dplyr::select(cpg, seqnames, start, end, width, 
                      UCSC_RefGene_Group, UCSC_RefGene_Name, Relation_to_UCSC_CpG_Island) %>%
        unique()
    } else {
      load(file.path(dir.data.aux,"great_EPICv2_array_annotation.rda"))
      anno_df <-  anno_df %>% 
        mutate(cpg = gsub("_.*", "",anno_df$Name)) %>%
        dplyr::select(cpg, seqnames, start, end, width, 
                      UCSC_RefGene_Group, UCSC_RefGene_Name, Relation_to_Island) %>%
        unique()
    }
    
    result <- left_join(
      result,
      anno_df
    )

    result <- dplyr::left_join(result, great %>% ungroup %>% dplyr::select(cpg, GREAT_annotation))
  }

  if(save){
    write_csv(
      result,
      file.path(dir.save, paste0(prefix, "_annotated_results.csv"))
    )
  }
  
  return(result)
}
# ===================================================================================================
# Annotation DMR
# ===================================================================================================
# Annotate the regions with E073 15-core marks segmentation states
# ------------------------------------------------------------------------------
# Description:
#   Annotates regions with ChromHMM segmentation states based on E073 15-core marks segmentation data.
#
# Parameters:
#   result       : Data frame. Regions to be annotated.
#   ChmmModels.gr: GRanges object. ChromHMM segmentation data containing state information.
#   region_var   : Character. The column name in 'result' representing region identifiers.
#
# Returns:
#   The original data frame with an added column "E073_15_coreMarks_segments_state" for segmentation states.

annotate_coreMarks_segments <- function(result, ChmmModels.gr, region_var) {
  # Create a GRanges object from the result data frame
  result.gr <- result %>% 
    makeGRangesFromDataFrame(start.field = "start", 
                             end.field = "end", 
                             seqnames.field = "seqnames")

  # Find overlaps between the result regions and the ChromHMM models
  hits <- findOverlaps(result.gr, ChmmModels.gr) %>% as.data.frame()
  hits$state <- ChmmModels.gr$state[hits$subjectHits]
  hits$region <- result[[region_var]][hits$queryHits]
  
  # Match the state information to each region based on the 'region' string
  result$E073_15_coreMarks_segments_state <- hits$state[match(result[[region_var]], hits$region)]
  
  return(result)
}
# ------------------------------------------------------------------------------
# Function: annotate_great
# ------------------------------------------------------------------------------
# Description:
#   Annotates regions with gene associations using GREAT analysis.
#
# Parameters:
#   result : Data frame. Regions to be annotated.
#   genome : Character. The genome version ("hg19" or "hg38") used in GREAT analysis.
#
# Returns:
#   A data frame with an additional column "GREAT_annotation" containing GREAT-derived gene information.
# ------------------------------------------------------------------------------
annotate_great <- function(result, genome) {
  # Create a GRanges object from the result data frame
  result.gr <- result %>% 
    makeGRangesFromDataFrame(start.field = "start", 
                             end.field = "end", 
                             seqnames.field = "seqnames")
  
  # Submit the GREAT job and retrieve gene associations
  job <- submitGreatJob(result.gr, species = genome)
  regionsToGenes_gr <- rGREAT::getRegionGeneAssociations(job)
  regionsToGenes <- as.data.frame(regionsToGenes_gr)
  
  # Create annotation strings for each region
  GREAT_annotation <- lapply(seq_len(length(regionsToGenes$annotated_genes)), function(i) {
    g <- ifelse(regionsToGenes$dist_to_TSS[[i]] > 0,
                paste0(regionsToGenes$annotated_genes[[i]], " (+", regionsToGenes$dist_to_TSS[[i]], ")"),
                paste0(regionsToGenes$annotated_genes[[i]], " (", regionsToGenes$dist_to_TSS[[i]], ")"))
    paste0(g, collapse = ";")
  })
  
  # Select key columns from GREAT output and combine with the annotation strings
  great <- dplyr::select(regionsToGenes, seqnames, start, end, width)
  great <- data.frame(great, GREAT_annotation = unlist(GREAT_annotation))
  
  # Merge the GREAT annotation with the original result data frame
  result <- dplyr::left_join(result, great, by = c("seqnames", "start", "end"))
  
  return(result)
}
# ------------------------------------------------------------------------------
# Function: annotate_region
# ------------------------------------------------------------------------------
# Description:
#   Annotates regions with CpG UCSC information by summarizing overlapping CpG details.
#
# Parameters:
#   result      : Data frame. Regions to be annotated.
#   array       : Character. The methylation array type (e.g., "HM450", "EPIC").
#   genome      : Character. The genome version ("hg19" or "hg38").
#   dir.data.aux: Character. Directory path for auxiliary annotation data.
#   cpg         : Vector. List of CpG identifiers to be considered.
#   region_var  : Character. (Optional) Column name in 'result' that contains region identifiers. Default is "region".
#   cores       : Numeric. Number of CPU cores to use for parallel processing. Default is 30.
#
# Returns:
#   A data frame with additional columns summarizing UCSC CpG annotation information.

annotate_region <- function(result, array, genome, dir.data.aux, cpg, region_var = "region", cores = 30) {

  anno.gr <- get_anno_gr(array = array, genome = genome, dir.data.aux = dir.data.aux)
  anno_df <- as.data.frame(anno.gr)
  
  if(genome == "hg38" & array %in% c("EPIC","EPICv1")) island_var <- "Relation_to_UCSC_CpG_Island"
  else island_var <- "Relation_to_Island"
  
  # Create a GRanges object from the result data frame
  result.gr <- result %>% 
    makeGRangesFromDataFrame(start.field = "start", 
                             end.field = "end", 
                             seqnames.field = "seqnames")
  
  doParallel::registerDoParallel(cores)
  
  anno <- plyr::ldply(
    1:length(result.gr),
    .fun = function (i) {
      rg <- result.gr[i]
      overlap <- findOverlaps(rg, anno.gr)
      hit <- subjectHits(overlap)
      anno_df_sub <- anno_df[hit,]
      cpgs <- intersect(anno_df_sub$Name, cpg)
      anno_df_sub <- anno_df_sub[match(cpgs, anno_df_sub$Name),]
      
      # cpg in region
      cpgs_in_region <- paste(cpgs, collapse = ",") 
      
      # UCSC_RefGene_Name
      UCSC_RefGene_Name <- paste0(unique(anno_df_sub[,"UCSC_RefGene_Name"]), collapse = ";")
      
      # UCSC_RefGene_Accession
      UCSC_RefGene_Accession <- paste0(unique(anno_df_sub[,"UCSC_RefGene_Accession"]), collapse = ";")
      
      # UCSC_RefGene_Group
      UCSC_RefGene_Group <- paste0(unique(anno_df_sub[,"UCSC_RefGene_Group"]), collapse = ";")
      
      # Relation_to_Island
      Relation_to_Island <- paste0(unique(anno_df_sub[,island_var]), collapse = ";")
      
      df <- data.frame(
        num_probes = length(cpgs),
        UCSC_RefGene_Name = UCSC_RefGene_Name,
        UCSC_RefGene_Accession = UCSC_RefGene_Accession,
        UCSC_RefGene_Group = UCSC_RefGene_Group,
        Relation_to_Island = Relation_to_Island,
        cpgs_in_region = cpgs_in_region
      )
      
      df[df == "NA"] <- NA
      
      df
    }, .parallel = T
  )

  cbind(result, anno)
}
# ------------------------------------------------------------------------------
# Function: annotate_enhancer
# ------------------------------------------------------------------------------
# Description:
#   Annotates regions with enhancer overlap information using an external enhancer dataset.
#
# Parameters:
#   result            : Data frame. Regions to be annotated.
#   nasser.enhancer.gr: GRanges object. Enhancer regions with associated cell type information.
#   cpg               : Vector. List of CpG identifiers used in the analysis.
#   genome            : Character. The genome version ("hg19" or "hg38").
#   array             : Character. The methylation array type.
#
# Returns:
#   A data frame with two additional columns:
#     - nasser_is_enhancer: Logical indicating enhancer overlap.
#     - nasser_is_enhancer_cell_types: Character string of associated cell types.


annotate_enhancer <- function(result, nasser.enhancer.gr, cpg, genome, array) {
  # Create a GRanges object from the result data frame
  result.gr <- result %>% 
    makeGRangesFromDataFrame(start.field = "start", 
                             end.field = "end", 
                             seqnames.field = "seqnames")
  
  # Find overlaps between the result regions and the enhancer regions
  hits <- findOverlaps(result.gr, nasser.enhancer.gr) %>% as.data.frame()
  
  # Initialize the enhancer annotation columns
  result$nasser_is_enhancer <- FALSE
  result$nasser_is_enhancer[unique(hits$queryHits)] <- TRUE
  
  result$nasser_is_enhancer_cell_types <- NA
  result$nasser_is_enhancer_cell_types[unique(hits$queryHits)] <- sapply(unique(hits$queryHits), function(x) {
    paste(unique(nasser.enhancer.gr$CellType[hits$subjectHits[hits$queryHits %in% x]]), collapse = ",")
  })
  
  return(result)
}

annotate_chmm <- function(result, dir.data.aux = dir.data.aux) {
  
  load(file.path(dir.data.aux,"E073_15_coreMarks_segments.rda"))
  
  message("Annotating E073_15_coreMarks_segments")
  
  result$region <- paste0(result$seqnames,":",result$start,"-", result$end)      
  result$start <- as.numeric(result$start)
  result$end <- as.numeric(result$end)
  result.gr <- result %>% makeGRangesFromDataFrame(
    start.field = "start",
    end.field = "end",
    seqnames.field = "seqnames"
  )
  hits <- findOverlaps(result.gr, ChmmModels.gr) %>% as.data.frame()
  hits$state <- ChmmModels.gr$state[hits$subjectHits]
  hits$region <- result$region[hits$queryHits]
  result$E073_15_coreMarks_segments_state <- hits$state[match(result$region,hits$region)]
  result$region <- NULL
  
  result
}
# ------------------------------------------------------------------------------
# Function: add_dmr_annotation
# ------------------------------------------------------------------------------
# Description:
#   Integrates multiple annotation methods (region, enhancer, coreMarks segmentation, GREAT) for DMRs.
#
# Parameters:
#   result       : Data frame. DMRs (differentially methylated regions) to be annotated.
#   cpg          : Vector. List of CpG identifiers used for regional annotation.
#   dir.data.aux : Character. Directory path for auxiliary data files (e.g., segmentation, enhancer data).
#   region_var   : (Optional) Character. Column name in 'result' containing region identifiers.
#                  If NULL, a region identifier is created using chromosome, start, and end.
#   array        : Character. The methylation array type (default "EPIC").
#   genome       : Character. Genome version ("hg19" or "hg38").
#   annotate     : Vector of characters. Specifies which annotation types to apply. Options include:
#                  "region", "enhancer", "E073_15_coreMarks_segments", "GREAT".
#
# Returns:
#   A data frame with multiple annotation columns added.

add_dmr_annotation <- function(result, 
                               cpg, 
                               dir.data.aux, 
                               region_var = NULL,
                               array = "EPIC", 
                               genome = "hg19",
                               annotate = c("region", "enhancer", "E073_15_coreMarks_segments", "GREAT")) {
  
  # Prepare the result data frame for genomic range creation

  if(is.null(region_var)) {
    result$seqnames <- paste0("chr", result$`#chrom`)
    result$region <- paste0(result$seqnames, ":", result$start, "-", result$end)
    region_var <- "region"
  } else {
    region <- str_split(result[[region_var]], ":|-", simplify = T)
    result$seqnames <- region[,1]
    result$start <- region[,2]
    result$end <- region[,3]
  }
  result$start <- as.numeric(result$start)
  result$end <- as.numeric(result$end)

  
  # Annotate with E073_15_coreMarks_segments if selected.
  if ("E073_15_coreMarks_segments" %in% annotate) {
    load(file.path(dir.data.aux, "E073_15_coreMarks_segments.rda"))
    message("Annotating E073_15_coreMarks_segments")
    result <- annotate_coreMarks_segments(result, ChmmModels.gr, region_var)
  }
  
  # Annotate with GREAT if selected.
  if ("GREAT" %in% annotate) {
    message("Annotating GREAT")
    result <- annotate_great(result, genome)
  }
  
  # Annotate with enhancer if selected.
  if ("enhancer" %in% annotate) {
    # Load enhancer data and filter
    data <- readr::read_tsv(file.path(dir.data.aux, "AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz"),
                            show_col_types = F)
    CellType.selected <- readxl::read_xlsx(file.path(dir.data.aux, "Nassser study selected biosamples.xlsx"), 
                                           col_names = FALSE) %>% dplyr::pull(1)
    data.filtered <- data %>% 
      dplyr::filter(CellType %in% CellType.selected) %>% 
      dplyr::filter(!isSelfPromoter) %>% 
      dplyr::filter(class != "promoter")
    nasser.enhancer.gr <- data.filtered %>% 
      makeGRangesFromDataFrame(start.field = "start",
                               end.field = "end",
                               seqnames.field = "chr",
                               keep.extra.columns = TRUE)
    
    message("Annotating enhancer")
    result <- annotate_enhancer(result, nasser.enhancer.gr, cpg, genome, array)
  }
  
  # Annotate with UCSC (island annotation) if selected.
  if ("region" %in% annotate) {
    message("Annotating region")
    result <- annotate_region(result, array, genome, dir.data.aux, cpg, region_var = region_var)
  }
  
  return(result)
}
# ===================================================================================================
# Bacon correction
# ===================================================================================================
# Function: bacon_adj
# ------------------------------------------------------------------------------
# Description:
#   Adjusts for bias and genomic inflation in EWAS data using the bacon method.
#   Supports correction using either z-scores or effect sizes with standard errors.
#
# Parameters:
#   data    : Data frame. The EWAS results containing statistics to be adjusted.
#   est_var : Character. Column name for effect size estimates.
#   z_var   : Character. Column name for z-scores.
#   std_var : Character. Column name for standard errors.
#   use_z   : Logical. If TRUE, bacon correction is performed using z-scores only; otherwise, effect sizes and SEs are used.
#   save    : Logical. Whether to save the bacon-corrected data and inflation statistics to files.
#   dir.save: Character. Directory path where the output files will be saved.
#   prefix  : Character. A prefix used for naming the output files.
#
# Returns:
#   A list containing:
#     - data.with.inflation: The corrected EWAS data frame.
#     - bacon.obj        : The bacon object from the initial analysis.
#     - inflation.stat   : A data frame summarizing inflation and bias statistics.

bacon_adj <- function(data, est_var, z_var, std_var, 
                      use_z = F, 
                      save = F, 
                      dir.save = NULL,
                      prefix = "Framingham"){
  
  ### 1. Compute genomic inflation factor before bacon adjustment
  data <- data %>% mutate(
    chisq = get(z_var)^2
  ) 
  
  # inflation factor - last term is median from chisq distrn with 1 df  
  inflationFactor <- median(data$chisq,na.rm = TRUE) / qchisq(0.5, 1)
  print("lambda")
  print(inflationFactor)

  ### 2. bacon analysis
  if(use_z){
    ### 2. bacon analysis
    z_scores <- data[[z_var]]
    bc <- bacon(
      teststatistics = z_scores,
      na.exclude = TRUE,
      verbose = F
    )
    # inflation factor
    print("lambda.bacon")
    print(inflation(bc))
    # bias
    print("estimate bias")
    print(bias(bc))
    print("estimates")
    print(bacon::estimates(bc))
    
    ### 3. Create final dataset
    data.with.inflation <- data %>% mutate(
      zScore.bacon = tstat(bc)[,1],
      pValue.bacon.z = pval(bc)[,1],
      fdr.bacon.z = p.adjust(pval(bc), method = "fdr"),
    ) %>% mutate(z.value = z_scores)
    
    print("o After bacon correction")
    print("Conventional lambda")
    lambda.con <- median((data.with.inflation$zScore.bacon) ^ 2,na.rm = TRUE)/qchisq(0.5, 1)
    print(lambda.con)
    
    # percent_null <- trunc ( bacon::estimates(bc)[1]*100, digits = 0)
    # percent_1  <- trunc ( bacon::estimates(bc)[2]*100, digits = 0 )
    # percent_2  <- 100 - percent_null - percent_1  
    bc2 <- bacon(
      teststatistics = data.with.inflation$zScore.bacon,
      na.exclude = TRUE,
      priors = list(
        sigma = list(alpha = 1.28,  beta = 0.36), 
        mu = list(lambda = c(0, 3, -3), tau = c(1000, 100, 100)), 
        epsilon = list(gamma = c(90, 5, 5)))
    )
  } else {
    est <- data[[est_var]]
    se <- data[[std_var]]
    
    
    bc <- bacon(
      teststatistics = NULL,
      effectsizes = est,
      standarderrors = se,
      na.exclude = TRUE,
      verbose = F
    )
    #  posteriors(bc)
    # inflation factor
    print("lambda.bacon")
    print(inflation(bc))
    # bias
    print("estimate bias")
    print(bias(bc))
    print("estimates")
    print(bacon::estimates(bc))
    
    ### 3. Create final dataset
    data.with.inflation <- data.frame(
      data,
      Estimate.bacon = bacon::es(bc),
      StdErr.bacon = bacon::se(bc),
      pValue.bacon = pval(bc),
      fdr.bacon = p.adjust(pval(bc), method = "fdr"),
      stringsAsFactors = FALSE
    )
    
    print("o After bacon correction")
    print("Conventional lambda")
    lambda.con <- median((data.with.inflation$Estimate.bacon/data.with.inflation$StdErr.bacon) ^ 2,na.rm = TRUE)/qchisq(0.5, 1)
    print(lambda.con)
    
    # percent_null <- trunc ( estimates(bc)[1]*100, digits = 0)
    # percent_1  <- trunc ( estimates(bc)[2]*100, digits = 0 )
    # percent_2  <- 100 - percent_null - percent_1  
    bc2 <- bacon(
      teststatistics = NULL,
      effectsizes =  data.with.inflation$Estimate.bacon,
      standarderrors = data.with.inflation$StdErr.bacon,
      na.exclude = TRUE,
      priors = list(
        sigma = list(alpha = 1.28,  beta = 0.36), 
        mu = list(lambda = c(0, 3, -3), tau = c(1000, 100, 100)), 
        epsilon = list(gamma = c(99, .5, .5)))
    )
  }
  
  print("inflation")
  print(inflation(bc2))
  print("estimates")
  print(bacon::estimates(bc2))
  data.with.inflation$chisq <- NULL
  
  inflation.stat <- data.frame(
    "Inflation.org" = inflationFactor,
    "Inflation.bacon" = inflation(bc),
    "Bias.bacon" = bias(bc),
    "Inflation.after.correction" = lambda.con,
    "Inflation.bacon.after.correction" = inflation(bc2),
    "Bias.bacon.after.correction" = bias(bc2)
  )
  
  
  if(save){
    readr::write_csv(
      data.with.inflation,
      file.path(dir.save, paste0(prefix, "_bacon_correction.csv"))
    )
    
    writexl::write_xlsx(
      inflation.stat,
      file.path(dir.save, paste0(prefix, "_inflation_stats.xlsx"))
    )
  }
  return(
    list(
      "data.with.inflation" = data.with.inflation,
      "bacon.obj" = bc,
      "inflation.stat" = inflation.stat
    )
  )
}



