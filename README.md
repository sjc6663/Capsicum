# Capsicum

This repository holds code, R scripts, data, and information for the published manuscript: 

Bierly, S., Van Syoc, E., Westphalen, M., Miles, A., Gaeta, N., Felix, T., Hristov, A., Ganda, E. 2024. "Alterations of rumen and fecal microbiome in growing beef and dairy steers fed rumen-protected *Capsicum* oleoresin." Journal of Animal Science https://doi.org/10.1093/jas/skae014


Analysis was done using the dada2 pipeline and downstream analysis was performed with R in RStudio. Below is a breakdown of what each file is relevant for related to analysis. 

## Data

### 1

------------------
## Scripts

### bash scripts folder
This folder contains all scripts necessary for analysis in the dada2 pipeline. This is a combination of bash and R scripts. 

### generate-phyloseq.R
This R script is utilized to generate a phyloseq object from the output asv and tax tables from dada2. 

### decontam.R
This R script is utilized to remove putative contaminants in samples based on taxa identified in positive and negative controls sequenced with samples. 


