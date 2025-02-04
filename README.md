# Pipeline for the 16S-metabarcoding characterization in Pantanos de Villa - Lima - Peru
This computational pipeline is designed for comprehensive microbiome analysis of wetland samples using 16S rRNA gene sequencing. Developed to process complex environmental microbiome data, the pipeline provides a robust, reproducible workflow for taxonomic characterization and ecological insights in wetland microbial communities.

## 03/02/2025 - Finished Cutadapt
- Could not find a reverse complement primer in the evaluated sequences using the `grep -iE` command
- Could find the 5'-3' oriented initial primers in all sequences.
- Cutadapt was used for these primers only.
- Any other treatment will have to be done using quality plot info from DADA2.

## 04/02/2025 - Installing the conda environment
```
conda env create -f dada2_env.yaml
conda activate dada2_env
``` 
