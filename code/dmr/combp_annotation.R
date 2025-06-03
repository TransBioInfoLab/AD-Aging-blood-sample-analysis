dir.base <- "."
dir.data <- file.path(dir.base, "datasets/combp")
dir.results <- file.path(dir.base, "analysis_results")
dir.results.combp <- file.path(dir.results, "DMR")
dir.results.meta <- file.path(dir.results, "meta_analysis")
dir.data.aux <- file.path(dir.base, "../DATASETS/Aux_Sync")
for(p in grep("dir.",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

# Auxilary function
source(file.path(dir.base, "code/utility/annotation_and_bacon.R"))

# Load data

results <- read_csv(
  file.path(dir.results.meta,
            "meta_analysis_FHS9_ADNI_results.csv")
)

combp_results <- readxl::read_xlsx(
  file.path(dir.data, "cnew.regions-p.bed.xlsx")
)

# Add annotation
results_anno <- add_dmr_annotation(
  result = combp_results,
  cpg = results$cpg,
  array = "EPIC",
  genome = "hg38",
  dir.data.aux = dir.data.aux
)

# Add direction
direction <- plyr::laply(
  str_split(results_anno$cpgs_in_region, ","),
  .fun = function(cpgs){
    est <- results %>% filter(cpg %in% cpgs) %>% pull(estimate)
    paste0(ifelse(est < 0, "-", "+"), collapse = "")
  }
)
results_anno$direction <- direction
# Calculate percentage of direction
pct.direction <- plyr::laply(
  str_split(direction, ""),
  .fun = function(d){
    sum(d == "-")/length(d)
  }
)
results_anno$pct_direction <- pct.direction


write_csv(
  results_anno,
  file.path(dir.results.combp, "combp_FHS9_ADNI_results_annotated.csv")
)


results_anno <- read_csv(
  file.path(dir.results.combp, "combp_FHS9_ADNI_results_annotated.csv")
)

results_sig <- results_anno %>% 
  filter(n_probes >= 3, z_sidak_p < 0.05, pct_direction %in% c(0,1), z_p < 1e-05)

write_csv(
  results_sig,
  file.path(dir.results.combp, "combp_filtered_results_annotated.csv")
)
