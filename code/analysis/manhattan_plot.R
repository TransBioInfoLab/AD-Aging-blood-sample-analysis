# Manhattan plot
dir.base <- "."
dir.results <- file.path(dir.base, "analysis_results")
dir.results.meta <- file.path(dir.results, "meta_analysis")
dir.data.aux <- file.path(dir.base, "../DATASETS/Aux_Sync/") 
dir.results.combp <- file.path(dir.results, "DMR")
dir.plot <- file.path(dir.results, "plots")

source(file.path(dir.base, "code/utility/plot.R")) 

results <- read_csv(
  file.path(dir.results.meta,
            "meta_analysis_FHS9_ADNI_results.csv")
)

cpgs.sig <- results %>% filter(pVal.fixed < 1e-05 & 
               direction %in% c("++", "--")) %>%
  pull(cpg) 

sig.dmr.df <- read_csv(
  file.path(dir.results.combp, "coMeth_combp_overlap_dmr.csv")
) 
sig.dmr_top20 <- sig.dmr.df %>% 
  slice_min(z_sidak_p, n = 20)
sig.dmr <- str_split(sig.dmr_top20$cpgs_in_region, ",") %>% unlist() ## 83

cpgs_df <- data.frame(
  "CpG" = results$cpg,
  "pVal.final" = results$pVal.fixed,
  "fdr" = results$fdr,
  "pos" = results$start,
  "chr" = as.numeric(gsub("chr", "", results$seqnames)),
  "GREAT" = gsub("\\(.*| ", "", results$GREAT_annotation) 
) 

annotated_cpg <- results %>% 
  filter(cpg %in% cpgs.sig) %>%
  slice_min(pVal.fixed, n = 20) %>%
  mutate(GREAT = gsub("\\;.*| ", "", GREAT_annotation),
         TSS = as.numeric(str_extract(GREAT, "(?<=\\().*(?=\\))"))) %>%
  filter(abs(TSS) < 2000 & abs(TSS) > 100) %>%
  pull(GREAT) %>%
  gsub("\\(.*| ", "",.)
#annotated_cpg[annotated_cpg == "RBM14-RBM4"]  <- "RBM14"
annotated_dmr <- results %>% 
  filter(cpg %in% sig.dmr) %>%
  mutate(GREAT = gsub("\\;.*| ", "", GREAT_annotation),
         TSS = as.numeric(str_extract(GREAT, "(?<=\\().*(?=\\))"))) %>%
  filter(abs(TSS) < 2000 & abs(TSS) > 100) %>%
  pull(GREAT) %>%
  gsub("\\(.*| ", "",.)

plot_manh(cpgs_df, annotated = unique(c(annotated_cpg, annotated_dmr)), colored = "")

ggplot2::ggsave(
  filename = file.path(dir.plot,"manhattan_plot.pdf"),
  width = 10,
  height = 6
)
