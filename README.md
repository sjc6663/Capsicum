# Capsicum

This repository holds code, R scripts, data, and information for the published manuscript: 

Bierly, S., Van Syoc, E., Westphalen, M., Miles, A., Gaeta, N., Felix, T., Hristov, A., Ganda, E. 2024. "Alterations of rumen and fecal microbiome in growing beef and dairy steers fed rumen-protected *Capsicum* oleoresin." Journal of Animal Science https://doi.org/10.1093/jas/skae014


Analysis was done using the dada2 pipeline and downstream analysis was performed with R in RStudio. Below is a breakdown of what each file is relevant for related to analysis. 

## Data

### rawstats.txt
This file contains information about the raw reads. Obtained by running seqkit stats in miniconda3. 

### metadata.csv
This file contains all the metadata for all samples. 

### asv-table.RData
This file contains the asv count information for all samples. This is an output file from dada2.

### tax-table.RData
This file contains the taxa information for all samples. This is an output file from dada2. 

### phyloseq-fecal-samples.RData
This is a phyloseq object of only fecal samples for both breeds. 

### phyloseq-rumen-samples.RData
This is a phyloseq object of only rumen samples for both breeds. 

## Scripts

### bash scripts folder
This folder contains all scripts necessary for analysis in the dada2 pipeline. This is a combination of bash and R scripts. 

### generate-phyloseq.R
This R script is utilized to generate a phyloseq object from the output asv and tax tables from dada2. 

### decontam.R
This R script is utilized to remove putative contaminants in samples based on taxa identified in positive and negative controls sequenced with samples.

### cross-over-effects.R
This R script tests for cross over effects of both alpha and beta diversity based on treatment order. 

### figure2.R
This R script generates the relative abundance figure for fecal samples (Figure 2). 

### figure3.R
This R script generates the relative abundance figure for rumen samples (Figure 3). 

### rel-abund-percents.R
This R script generates percentages of taxa that are changed based on treatment. 

### beta-diversity-figure4.R
This R script generates p-values for beta diversity of fecal samples by treatment and generates the Principal Coordinate Analysis (PCA) plot for fecal samples (Figure 4).

### beta-diversity-figure5.R
This R script generates p-values for beta diversity of rumen samples by treatment and time and generates the PCA plot for rumen samples for treatment and time (Figure 5).

### alpha-diversity-figure6.R
This R script generates p-values for alpha diversity of fecal samples by treatment and generates the boxplot for fecal samples (Figure 6). 

### alpha-diversity-figure7.R
This R script generates p-values for alpha diversity of rumen samples by treatment and generates the boxplot for fecal samples (Figure 7). 

### supplemental-figure1.R
This R script generates the lineplot for alpha diversity of rumen samples by time (Supplemental Figure 1). 

### DRA-ALDEx-Time.R
This R script uses ALDEx2 to compare taxa based on time (H0 vs. H2 and H0 vs. H18) for rumen samples of both breeds. 

### DRA-ALDEx-treatment.R
This R script uses ALDEx2 to compare taxa between a Control Group and a *Capsicum* Group based on treatment for both fecal and rumen samples. 

