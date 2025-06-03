# **The Aging Epigenome: Integrative Analyses Reveal Functional and Genetic Links to Alzheimer’s Disease**
Wei Zhang, David Lukacsovich, Juan I. Young, Lissette Gomez, Michael A. Schmidt, Brian W. Kunkle, Xi Chen, Eden R. Martin, Lily Wang

## Citing this repository

To be added

### Description

This github repository includes scripts used for the analyses in the above manuscript. 

**Background**  

Aging is the strongest risk factor for Alzheimer’s disease (AD), yet the role of age-associated DNA methylation (DNAm) changes in blood and their relevance to AD remains poorly understood. We aimed to identify blood-based methylation differences associated with aging in older adults and explore their functional relevance, cross-tissue concordance, and genetic colocalization with AD loci. 

**Methods**      

 We conducted a meta-analysis of DNAm data from 475 dementia-free participants aged over 65 years across two independent cohorts: the Framingham Heart Study (n=282) and the Alzheimer’s Disease Neuroimaging Initiative (n=193). DNAm levels were measured using the Illumina HumanMethylation EPIC array. Linear regression models, adjusting sex and immune cell-type proportions, were used to assess associations with chronological age. Cohort-specific results were corrected for genomic inflation using the bacon method and combined using inverse-variance fixed-effects meta-analysis. Differentially methylated regions were identified using comb-p and coMethDMR approaches. Functional annotation included overlap with gene expression quantitative trait methylation analysis, pathways enrichment, and colocalization with loci nominated by genome-wide association studies for dementia. We also compared results with recent DNAm studies in AD and performed brain-to-blood correlation analyses using 69 matched samples pairs.

**Results**        

We identified 3758 CpGs and 556 differentially methylated regions (DMRs) significantly associated with aging. Promoter-centered hypermethylation was the dominant signature. Pathway enrichment analyses highlighted immune response, metabolic regulation, and pathways involved in neuronal signaling and Alzheimer’s disease. Notably, 73 aging-associated CpGs were significantly correlated with expression of nearby genes involved in interferon signaling, lipid metabolism, and energy homeostasis. Approximately 45% of significant CpGs harbored methylation quantitative trait variants, and colocalization analysis revealed 32 regions where shared genetic variants influenced both methylation and dementia risk, including genes such as ABI3 and DLG4. Roughly one-third of aging-associated CpGs overlapped with Alzheimer’s disease-related methylation changes in independent studies. Furthermore, 23 CpGs demonstrated significant positive correlations between blood and brain methylation, suggesting cross-tissue consistency. 

**Conclusions**

Our findings indicate that age-associated DNA methylation changes in blood predominantly occur in gene promoters and converge on pathways and genetic variants implicated in Alzheimer’s disease. These results highlight potential biomarkers of aging and dementia risk and provide insights into the shared biological mechanisms that underlie aging and Alzheimer’s disease.es for AD.

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
| code/dmr/combp_annotation.R | [Link to the script](https://github.com/TransBioInfoLab/AD-aging-blood-samples-analysis/blob/main/code/analysis/dmr/combp_annotation.R) |
| code/dmr/coMethDMR.Rmd | [Link to the script](https://github.com/TransBioInfoLab/AD-aging-blood-samples-analysis/blob/main/code/dmr/coMethDMR.Rmd) |

### 4. **Integrative analyses with gene expression**

| File                             | Link                                                         |
| -------------------------------- | ------------------------------------------------------------ |
| code/dnam_vs_rna/dnam_vs_rna.Rmd | [Link to the script](https://github.com/TransBioInfoLab/AD-aging-blood-samples-analysis/blob/main/code/dnam_vs_rna/dnam_vs_rna.Rmd) |

### 5. Correlate expression in brain and blood samples

| File                 | Link |
|----------------------|-------------|
| code/brain_blood_corr/brain_blood_corr.Rmd | [Link to the script](https://github.com/TransBioInfoLab/AD-aging-blood-samples-analysis/blob/main/code/brain_blood_corr/brain_blood_corr.Rmd) |

### 6. **Validation using independent datasets**

| File                 | Link |
|----------------------|-------------|
| code/compare_study/compare_study_miami_ad.Rmd | [Link to the script](https://github.com/TransBioInfoLab/AD-aging-blood-samples-analysis/blob/main/code/compare_study/compare_study_miami_ad.Rmd) |

### 7. Assessment Plots

| File                 | Link |
|----------------------|-------------|
| code/analysis/manhattan_plot.R | [Link to the script](https://github.com/TransBioInfoLab/AD-aging-blood-samples-analysis/blob/main/code/analysis/manhattan_plot.R) |

## For reproducible research

To perform the analysis, begin by installing the packages found in `session_info.R` ([Link to the script](https://github.com/TransBioInfoLab/AD-aging-blood-samples-analysis/blob/main/code/session_info.R)). Then, load the auxiliary functions from folder `code/utility` ([Link to the file](https://github.com/TransBioInfoLab/AD-aging-blood-samples-analysis/tree/main/code/utility)). Follow the sequence provided in the Description to conduct the analysis.

## Acknowledgement

Data used in the preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design and implementation of ADNI and/or provided data but did not participate in the analysis or writing of this report. A complete listing of ADNI investigators can be found at: http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf

## Reference

- [^1]: Zhang W, Young JI, Gomez L, et al. Blood DNA methylation signature for incident dementia: Evidence from longitudinal cohorts. *Alzheimer's Dement.* 2025; 21:e14496. https://doi.org/10.1002/alz.14496
