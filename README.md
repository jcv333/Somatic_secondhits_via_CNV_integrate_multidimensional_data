# Code to determine second hits of each germline genetic variant in cohort patients utilizing Copy Number Variants (CNVs) data, followed by the integration of multidimensional data (germline, tumor and clinical information).
## Repository created by Jose Camacho Valenzuela.

### Disclaimer.
This repository provides an example of a mock Rscript which is part of a workflow from unpublished work. The script was anonymized and intended to show part of the workflow.

## Steps of the workflow (written as `Rscript`):

### 1) Determine second hit of each germline variant based on their genomic coordinates using CNVs data.
Evaluate whether the genomic coordinate of a given germline variant is present within a segment of the CNV data to assess whether there is LOH for the variant in question. This analysis is scaled up to every germline variant identified in the germline exome of the patient. A germline variant with a second hit via LOH is considered a form of biallelic inactivation.

### 2) Determine somatic loss of both copies of the HRD-related genes <i> BARD1, BRCA1, BRCA2, RAD51C, RAD51D</i> and <i> PALB2 </i>.
Evaluate whether the genomic coordinates of each of these genes are present within a segment of the CNV data to assess whether there is homozygous deletion of the gene in question. A full loss of both copies of these genes is considered a form of biallelic inactivation.

### 3) Determine second hit of each germline variant based on somatic point mutations.
Evaluate whether there is a somatic inactivating point mutation in the second allele of a germline-mutated gene, which is considered a form of biallelic inactivation.

### 4) Integration of multidimensional data.
Integrate:
- Germline data (information of each variant).
- Tumor data (second hits either as LOH or point mutations, CNV status of the 6 HRD-related genes, methylation in <i>BRCA1/RAD51C</i>, and Signature 3).
- Clinical data (age of diagnosis, sex, ethnicity, clinical subtype, among others).

