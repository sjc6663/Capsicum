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
load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")

# aggregate at genus level for counts
A_fecal_counts <- A_fecal_counts %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

H_fecal_counts <- H_fecal_counts %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

A_rumen_counts <- A_rumen_counts %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

H_rumen_counts <- H_rumen_counts %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

A_rumen_counts <- A_rumen_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

H_rumen_counts <- H_rumen_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

# transpose OTU to make ALDEx2 happy
H_pscount3<-t(otu_table(H_fecal_counts))
A_pscount3 <- t(otu_table(A_fecal_counts))
AR_pscount3 <- t(otu_table(A_rumen_counts))
HR_pscount3 <- t(otu_table(H_rumen_counts))

HQ_pscount3 <- t(otu_table(H_quick))
HL_pscount3 <- t(otu_table(H_long))

# run AlDEx2 function
H_aldex2_da <- ALDEx2::aldex(data.frame(H_pscount3), phyloseq::sample_data(H_fecal_counts)$Treatment, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")
A_aldex2_da <- ALDEx2::aldex(data.frame(A_pscount3), phyloseq::sample_data(A_fecal_counts)$Treatment, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

AR_aldex2_da <- ALDEx2::aldex(data.frame(AR_pscount3), phyloseq::sample_data(A_rumen_counts)$SampleBinary, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")
HR_aldex2_da <- ALDEx2::aldex(data.frame(HR_pscount3), phyloseq::sample_data(H_rumen_counts)$SampleBinary, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

HQ_aldex2_da <- ALDEx2::aldex(data.frame(HQ_pscount3), phyloseq::sample_data(H_quick)$Hour, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")
HL_aldex2_da <- ALDEx2::aldex(data.frame(HL_pscount3), phyloseq::sample_data(H_long)$Hour, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")



# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
AR_sig_aldex2 <- AR_aldex2_da %>%
  filter(wi.eBH < 0.05)

HR_sig_aldex2 <- HR_aldex2_da %>%
  filter(wi.eBH < 0.05)

HQ_sig_aldex2 <- HQ_aldex2_da %>%
  filter(wi.eBH < 0.05)

HL_sig_aldex2 <- HL_aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
AR_taxa_info <- data.frame(tax_table(A_rumen_counts))
HQ_taxa_info <- data.frame(tax_table(H_quick))


AR_taxa_info <- AR_taxa_info %>% rownames_to_column(var = "OTU")
HQ_taxa_info <- HQ_taxa_info %>% rownames_to_column(var = "OTU")


# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(AR_aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

ALDEx2::aldex.plot(HR_aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

ALDEx2::aldex.plot(HQ_aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)


# make a table of significant corrected p-values
AR_sig_aldex2 <- AR_aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

HQ_sig_aldex2 <- HQ_aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
AR_sig_aldex2 <- left_join(AR_sig_aldex2, AR_taxa_info)

HQ_sig_aldex2 <- left_join(HQ_sig_aldex2, HQ_taxa_info)

