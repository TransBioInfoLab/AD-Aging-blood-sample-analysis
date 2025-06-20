#######################################################################################################
# ===================================================================================================
# Function for pathway analysis 
# ===================================================================================================
#######################################################################################################
library(methylGSA)
library(missMethyl)
library(msigdbr)
# ===================================================================================================
# Wrapper function for missMethyl
# ===================================================================================================
missMethyl_wrapper <- function (dir.data.aux, use_msigdbr = T,  ...) {
  anno.gr <- get_anno_gr(array = "EPIC", genome = "hg19", dir.data.aux)
  anno_df <- DataFrame(anno.gr)
  
  if(use_msigdbr) {
    gs.list <- NULL
    gs.des.list <- NULL
    for (subcat in c("CP:BIOCARTA", "CP:KEGG", "CP:PID", "CP:REACTOME", "CP:WIKIPATHWAYS")) {
      gs <- get_msigdbr(subcategory = subcat, gene_id = "human_entrez_gene")
      gs.list <- c(gs.list, gs$GS.list)
      gs.des.list <- c(gs.des.list, gs$Des)
    }
    gsa <- gsameth(..., anno = anno_df, collection = gs.list)
    
    gsa <- gsa %>% 
      tibble::rownames_to_column("Pathway") %>%
      dplyr::mutate(Description = gs.des.list, .before = N)
  } else {
    
    gsa_ls <- plyr::llply(
      c("GO", "KEGG"),
      .fun = function(gs) {
        gometh(..., collection = gs, anno = anno_df, sig.genes = F)
      }
    )
    names(gsa_ls) <- c("GO", "KEGG")
    colnames(gsa_ls$GO)[2] <- "Description"
    gsa_ls$GO$Pathway_ID <- rownames(gsa_ls$GO)
    gsa_ls$KEGG$Pathway_ID <- rownames(gsa_ls$KEGG)
    gsa <- data.table::rbindlist(gsa_ls, idcol = "collection", fill = T) 
  }
  
  
  gsa
  
}
# ------------------------------------------------------------------------------
# Get msigdbr
get_msigdbr <- function(species = "Homo sapiens", category = "C2", subcategory = "CP", gene_id = "gene_symbol") {
  pathway_df <- msigdbr(species = species, collection = category, subcollection = subcategory)
  pathway <- unique(pathway_df$gs_name)

  pathway_des <- pathway_df %>% 
    dplyr::select(gs_name, gs_description) %>% 
    unique()
  pathway_des <- pathway_des[match(pathway, pathway_des$gs_name),] %>% pull(gs_description, name = gs_name)
  GS.list <- plyr::llply(
    pathway,
    .fun = function(p) {
      pathway_df %>% filter(gs_name %in% p) %>% pull(ncbi_gene)
    }
  )
  names(GS.list) <- pathway 
  
  list(GS.list = GS.list,
       Des = pathway_des)
}
# ===================================================================================================
# Wrapper function for methylGSA
# ===================================================================================================
# Function: methylGSA_wrapper
# ------------------------------------------------------------------------------
# Description:
#   A wrapper function to run methylGSA analysis on CpG p-values. For each gene set category,
#   it performs gene set analysis using the specified method (e.g., "GSEA" or "ORA"). If desired,
#   it also adds annotation for significant results using CpG annotations.
#
# Parameters:
#   cpg.pval         : Numeric vector or data frame. CpG-level p-values used for gene set analysis.
#   array            : Character. Methylation array type (e.g., "450K", "EPIC").
#   method           : Character. The analysis method to use (default is "GSEA"). Other options may include "ORA".
#   GS               : Vector. A list of gene set categories (GS) to be analyzed.
#   add_anno_for_sig : Logical. If TRUE (default), additional annotations are added for significantly enriched gene sets.
#   use_msigdbr      : Logical. If TRUE, gene sets are retrieved from the msigdbr package (default is FALSE).
#   fun              : Character. Name of the methylGSA function to be used (default is "methylRRA").
#   minsize          : Numeric. Minimum gene set size to consider (default is 100).
#   maxsize          : Numeric. Maximum gene set size to consider (default is 500).
#
# Returns:
#   A list where each element corresponds to one gene set category in GS. If add_anno_for_sig is TRUE,
#   each list element contains a sub-list with the raw methylGSA results (named "results") and the corresponding
#   annotations for significant gene sets.
# ------------------------------------------------------------------------------
methylGSA_wrapper <- function(cpg.pval, array, method = "GSEA", GS, add_anno_for_sig = T, 
                              use_msigdbr = F, fun = "methylRRA", minsize = 100, maxsize = 500){
  
  if(add_anno_for_sig) {
    CpGs <- getAnnot(array.type = array) %>% as.data.frame()
  }
  
  results <- plyr::llply(
    GS,
    .fun = function(gs.type) {
      if(use_msigdbr) {

        GS.list <- get_msigdbr(species = "Homo sapiens", category = "C2", subcategory = toupper(gs.type))$GS.list
        
      } else {
        GS.list <- NULL
      }
      
      fn <- eval(parse(text = fun))
      if(fun == "methylRRA") {
        res <- fn(
          cpg.pval = cpg.pval,
          method = method, 
          array.type = array,
          GS.type = gs.type,
          GS.list = GS.list,
          minsize = minsize,
          maxsize = maxsize
        )
      } else {
        res <- fn(
          cpg.pval = cpg.pval,
          array.type = array,
          GS.type = gs.type,
          GS.list = GS.list,
          minsize = minsize,
          maxsize = maxsize
        )
      }

  
      res_sig <- res %>% filter(.,padj < 0.05)
      
      # Add annotations
      if(add_anno_for_sig) {
        res_anno <- plyr::alply(
          res_sig,
          .margins = 1,
          .fun = function(pathway){

            if(method == "ORA") {
              pathway <- pathway %>% 
                separate_rows(overlap, sep = ",") 
              core <- "overlap"
            }
            if(method == "GSEA") {
              pathway <- pathway %>% 
                separate_rows(core_enrichment, sep = "/") 
              core <- "core_enrichment"
            }
           
            plyr::adply(
              pathway,
              .margins = 1,
              .fun = function(gene){
                cpgs <- CpGs %>%
                  filter(UCSC_RefGene_Name %in% gene[[core]]) %>%
                  dplyr::select(Name, UCSC_RefGene_Group)
                gene <- data.frame(gene, cpg = cpgs$Name)
              }, .id = "pathway"
            )
          }
        )
        names(res_anno) <- res_sig$Description
        return(
          c(list(results = res),
            res_anno)
        )
       
      }
      
      res
     
    }
  )
  
  names(results) <- GS
  results

}
# ------------------------------------------------------------------------------
# Get annotation (source function from methylGSA)
getAnnot = function(array.type, group = "all"){
  if(array.type=="450K"){
    FullAnnot = tryCatch({
      minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    },
    error = function(e){
      stop("IlluminaHumanMethylation450kanno.ilmn12.hg19 needs to
be installed and loaded before running methylglm/methylRRA")
    })
  }else{
    FullAnnot = tryCatch({
      minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    },
    error = function(e){
      stop("IlluminaHumanMethylationEPICanno.ilm10b4.hg19 needs to
be installed and loaded before running methylglm/methylRRA")
    })
  }

  FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group")]
  FullAnnot = FullAnnot[str_length(rownames(FullAnnot))==10,]
  FullAnnot = FullAnnot[!FullAnnot$UCSC_RefGene_Name=="",]
  ## get the first gene in each USCS_RefGene_Name
  temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Name,split=";"),
                '[', 1, FUN.VALUE=character(1))
  FullAnnot$UCSC_RefGene_Name = temp
  ## get the first gene group in each UCSC_RefGene_Group
  temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Group,split=";"),
                '[', 1, FUN.VALUE=character(1))
  FullAnnot$UCSC_RefGene_Group = temp

  if(group == "body"){
    FullAnnot =
      FullAnnot[FullAnnot$UCSC_RefGene_Group%in%c("Body", "1stExon"),]
  }

  if(group == "promoter1"){
    FullAnnot = FullAnnot[grepl("TSS",FullAnnot$UCSC_RefGene_Group),]
  }

  if(group == "promoter2"){
    FullAnnot =
      FullAnnot[FullAnnot$UCSC_RefGene_Group%in%c("TSS200", "TSS1500",
                                                  "1stExon", "5'UTR"),]
  }

  return(FullAnnot)
}