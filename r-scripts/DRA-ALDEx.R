# load packages
require(tidyverse)
require(phyloseq)
library(BiocManager)
BiocManager::install("ALDEx2")
library(ALDEx2)
library(microViz)

# load phyloseq of counts
load("ps-obj/phyloseq-fecal-samples-angus-counts.RData")
load("ps-obj/phyloseq-fecal-samples-holstein-counts.RData")

# aggregate at genus level for counts
A_fecal_counts <- A_fecal_counts %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

H_fecal_counts <- H_fecal_counts %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()


# transpose OTU to make ALDEx2 happy
H_pscount3<-t(otu_table(H_fecal_counts))
A_pscount3 <- t(otu_table(A_fecal_counts))

# run AlDEx2 function
H_aldex2_da <- ALDEx2::aldex(data.frame(H_pscount3), phyloseq::sample_data(H_fecal_counts)$Treatment, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")
A_aldex2_da <- ALDEx2::aldex(data.frame(A_pscount3), phyloseq::sample_data(A_fecal_counts)$Treatment, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")



# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)
# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(pscount2))
taxa_info <- taxa_info %>% rownames_to_column(var = “OTU”)
# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type=“MW”, test=“wilcox”, called.cex = 1, cutoff = 0.05)
# make a table of significant corrected p-values
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = “OTU”) %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)
# add in previously formed taxa information to complete the table
sig_aldex2 <- left_join(sig_aldex2, taxa_info)