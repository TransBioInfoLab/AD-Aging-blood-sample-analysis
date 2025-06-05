# **The Aging Epigenome: Integrative Analyses Reveal Functional Overlap with Alzheimer’s Disease**
Wei Zhang, David Lukacsovich, Juan I. Young, Lissette Gomez, Michael A. Schmidt, Brian W. Kunkle, Xi Chen, Eden R. Martin, Lily Wang

### Description

This github repository includes scripts used for the analyses in the above manuscript. 

**Background**  
Aging is the strongest risk factor for Alzheimer’s disease (AD), yet the role of age-associated DNA methylation (DNAm) changes in blood and their relevance to AD remains poorly understood. We aimed to identify blood-based methylation differences associated with aging in older adults and explore their functional relevance, cross-tissue concordance, and genetic colocalization with AD loci.

**Methods**      
We performed a meta-analysis of DNAm samples from 475 dementia-free subjects aged over 65 years across two independent cohorts, which included 282 subjects from the Framingham Heart Study at Exam 9 and 193 subjects from the Alzheimer’s Disease Neuroimaging Initiative (earliest visit with available DNAm data). All DNAm profiles were measured using the Illumina HumanMethylation EPIC array. Linear regression models, adjusting sex and immune cell-type proportions, were used to assess associations between DNAm and chronological age. Cohort-specific results were corrected for genomic inflation using the bacon method and combined using inverse-variance fixed-effects meta-analysis. Differentially methylated regions were identified using comb-p and coMethDMR approaches. Integrative analyses included overlap with gene expression quantitative trait methylation analysis, pathways enrichment, and colocalization with loci nominated by genome-wide association studies for dementia. We also compared results with recent DNAm studies in AD and performed brain-to-blood correlation analyses using 69 matched samples pairs.

**Results**        
At 5% false discovery rate, we identified 3758 CpGs and 556 differentially methylated regions (DMRs) consistently associated with aging in both cohorts. The majority of the hypermethylated DNAm differences were located in gene promoter regions, while the majority of the hypomethylated DNAm differences were in distal regions. Pathway enrichment analyses highlighted immune response, metabolic regulation, and synaptic plasticity. Approximately 45% of significant CpGs harbored methylation quantitative trait variants, and colocalization analysis revealed 32 regions where shared genetic variants influenced both methylation and dementia risk. Roughly one-third of aging-associated CpGs overlapped with AD-associated methylation changes in independent studies. Furthermore, we prioritized 9 aging-associated CpGs, located in promoter regions of PDE1B, ELOVL2, PODXL2, and other genomic regions, that showed strong positive blood-to-brain methylation concordance, as well as association with AD or AD neuropathology in independent studies.

**Conclusions**   
Our findings provided insights into the functional overlap between the aging process and AD, and nominated promising blood-based biomarkers for future AD research. 


### 1. Study cohorts, Preprocessing of DNA methylation data

The pre-processing procedures were previously described[^1].

| Data                   | Link                                                         |
| ---------------------- | ------------------------------------------------------------ |
| ADNI                   | [Link to the script](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/ADNI/preprocessing) |
| Framingham Heart Study | [Link to the script](https://github.com/TransBioInfoLab/blood-dnam-and-incident-dementia/blob/main/code/Framingham/preprocessing) |

### 2. Single cohor, Meta analysis, and **Pathway analysis**

| File                 | Link |
|----------------------|-------------|
| code/analysis/association_test.Rmd | [Link to the script](https://github.com/TransBioInfoLab/AD-Aging-blood-sample-analysis/blob/main/code/analysis/association_test.Rmd) |

### 3. Combp and cometh dmr results

| File                 | Link |
|----------------------|-------------|
| code/dmr/combp_annotation.R | [Link to the script](https://github.com/TransBioInfoLab/AD-Aging-blood-sample-analysis/blob/main/code/dmr/combp_annotation.R) |
| code/dmr/coMethDMR.Rmd | [Link to the script](https://github.com/TransBioInfoLab/AD-Aging-blood-sample-analysis/blob/main/code/dmr/coMethDMR.Rmd) |

### 4. **Integrative analyses with gene expression**

| File                             | Link                                                         |
| -------------------------------- | ------------------------------------------------------------ |
| code/dnam_vs_rna/dnam_vs_rna.Rmd | [Link to the script](https://github.com/TransBioInfoLab/AD-Aging-blood-sample-analysis/blob/main/code/dnam_vs_rna/dnam_vs_rna.Rmd) |

### 5. Correlate expression in brain and blood samples

| File                 | Link |
|----------------------|-------------|
| code/brain_blood_corr/brain_blood_corr.Rmd | [Link to the script](https://github.com/TransBioInfoLab/AD-Aging-blood-sample-analysis/blob/main/code/brain_blood_corr/brain_blood_corr.Rmd) |

### 6. **Validation using independent datasets**

| File                 | Link |
|----------------------|-------------|
| code/compare_study/compare_study_miami_ad.Rmd | [Link to the script](https://github.com/TransBioInfoLab/AD-Aging-blood-sample-analysis/blob/main/code/compare_study/compare_study_miami_ad.Rmd) |

### 7. Assessment Plots

| File                 | Link |
|----------------------|-------------|
| code/analysis/manhattan_plot.R | [Link to the script](https://github.com/TransBioInfoLab/AD-Aging-blood-sample-analysis/blob/main/code/analysis/manhattan_plot.R) |

## For reproducible research

To perform the analysis, begin by installing the packages found in `session_info.R` ([Link to the script](https://github.com/TransBioInfoLab/AD-Aging-blood-sample-analysis/blob/main/code/session_info.R)). Then, load the auxiliary functions from folder `code/utility` ([Link to the file](https://github.com/TransBioInfoLab/AD-Aging-blood-sample-analysis/tree/main/code/utility)). Follow the sequence provided in the Description to conduct the analysis.

## Acknowledgement

Data used in the preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design and implementation of ADNI and/or provided data but did not participate in the analysis or writing of this report. A complete listing of ADNI investigators can be found at: http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf

## Reference

[^1]: Zhang W, Young JI, Gomez L, et al. Blood DNA methylation signature for incident dementia: Evidence from longitudinal cohorts. *Alzheimer's Dement.* 2025; 21:e14496. https://doi.org/10.1002/alz.14496
