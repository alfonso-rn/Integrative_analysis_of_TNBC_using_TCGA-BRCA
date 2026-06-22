# Integrative analysis of TNBC using TCGA-BRCA
This repository contains the R workflow developed for the integrative analysis of triple-negative breast cancer (TNBC) using TCGA-BRCA data. The pipeline combines clinical,transcriptomic, mutational and cis-regulatory information to compare TNBC primary tumour samples against non-TNBC primary breast tumour samples.
The workflow uses TCGAbiolinks to retrieve and process TCGA-BRCA clinical, RNA-seq and somatic mutation data, and RBioGateway to query cis-regulatory modules (CRMs), transcription factors and regulatory annotations. Additional analyses include differential expression analysis, mutational profiling, TFBS disruption prediction and functional enrichment.

## Workflow overview
The analysis must be executed sequentially because each script generates intermediate objects that are required by the following steps.

### 00_Environment_configuration
Creates the project folder structure and installs the required R packages.

### 01_Identification_of_TNBC_patients_within_TCGA_BRCA.R
Downloads TCGA-BRCA clinical data and identifies clinically defined TNBC patients using ER, PR and HER2 status.

### 02_Prepare_TCGA_BRCA_expression_data.R
Downloads TCGA-BRCA RNA-seq STAR-Counts data from primary tumour samples, preprocesses the expression matrix, normalizes it, filters low-expression genes and separates samples into TNBC and non-TNBC groups.

### 03_Differential_expression_analysis.R
Performs differential expression analysis between TNBC and non-TNBC samples using TCGAbiolinks, with FDR < 0.01 and |log2FC| >= 1. It also generates summary tables and plots.

### 04_Somatic_mutations_in_TNBC_associated_DEGs.R
Downloads TCGA-BRCA masked somatic mutation data, compares the mutational profiles of TNBC and non-TNBC samples, and selects non-truncating intronic mutations located in differentially expressed genes.

### 05_Identification_of_affected_CRMs.R
Uses RBioGateway to identify CRMs overlapping intronic mutations in DEGs. CRMs are filtered by breast-tissue activity and prioritised when the mutated DEG, the CRM-located gene and the predicted CRM target gene match.

### 06_Identification_of_affected_TFBS.R
Uses motifbreakR to predict the effect of intronic SNPs on transcription factor binding sites. Strong TFBS effects are overlapped with the prioritised intragenic CRMs.

### 07_Functional_enrichment_analysis.R
Performs functional enrichment analysis of upregulated and downregulated DEGs with gprofiler2, and integrates enriched terms with genes associated with prioritised intragenic CRMs.
