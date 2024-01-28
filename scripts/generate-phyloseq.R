### Create a phyloseq object from dada2 data

# SAC 3/15/2022

## ---- setup ----
# load packages that we need
library(phyloseq)

## if you do not have tidyverse install it
# install.packages("tidyverse")

# load packages
library(tidyverse)

# load ASV table
load("dada2-data/asv-table.RData")

# load tax table
load("dada2-data/tax-table.RData")

## ---- create phyloseq object ----

# matrices in R
dim(seqtab.nochim)
dim(tax)

# create phyloseq object
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = FALSE),
               tax_table(tax))

# look at the phyloseq object
ps
# look at number of taxa
ntaxa(ps)
# get taxa ranks
rank_names(ps)

# access the data "slots" with @
head(ps@tax_table)
head(ps@otu_table)

# fix ASV names
### from dada2 tutorial: fix ASV names
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

# rerun code to look at tax table
head(ps@tax_table)
head(ps@otu_table)

# get taxa names
head(taxa_names(ps))
head(sample_names(ps))


## ---- load in metadata to phyloseq object ----

dat <- readxl::read_xlsx("metadata/Extraction_IDs.xlsx",
                         sheet = "Metadata")

# look at metadata
head(dat)
str(dat)

# fix sample names to get only sample ID
names <- str_extract(sample_names(ps), "C(\\d){1,3}")

# check the top of the object
head(names)

# check the bottom of the object
tail(names)

# change the sample names to NAMES
sample_names(ps) <- names

# did it work?
head(sample_names(ps))

# format our data to add to phyloseq
sampdf <- dat %>% 
  column_to_rownames(var = "Novogene-Sample")

head(dat)
cat(dat)

dat

# add to phyloseq
sample_data(ps) <- sampdf

# did it work?
ps

## ---- save our output ----

# re-name our phyloseq
psraw <- ps

# save as RImage
save(psraw, file = "ps-obj/phyloseq-raw.RData")

load("ps-obj/phyloseq-raw.RData")

psraw

## ---- take a look at what's in our phyloseq ----
get_taxa_unique(ps, "Kingdom")
get_taxa_unique(ps, "Phylum")

# visualize
library(microViz)

ps <- tax_fix(ps)

ps
