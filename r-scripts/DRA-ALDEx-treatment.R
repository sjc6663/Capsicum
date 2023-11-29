# load packages
require(tidyverse)
require(phyloseq)
library(BiocManager)
#BiocManager::install("ALDEx2")
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

A_fecal_counts <- A_fecal_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

H_fecal_counts <- H_fecal_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

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

# run AlDEx2 function
H_aldex2_da <- ALDEx2::aldex(data.frame(H_pscount3), phyloseq::sample_data(H_fecal_counts)$SampleBinary, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")
A_aldex2_da <- ALDEx2::aldex(data.frame(A_pscount3), phyloseq::sample_data(A_fecal_counts)$SampleBinary, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")



AR_aldex2_da <- ALDEx2::aldex(data.frame(AR_pscount3), phyloseq::sample_data(A_rumen_counts)$SampleBinary, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")
HR_aldex2_da <- ALDEx2::aldex(data.frame(HR_pscount3), phyloseq::sample_data(H_rumen_counts)$SampleBinary, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")



# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
AR_sig_aldex2 <- AR_aldex2_da %>%
  filter(wi.eBH < 0.05) # no significantly different taxa

HR_sig_aldex2 <- HR_aldex2_da %>%
  filter(wi.eBH < 0.05) # no significantly different taxa

H_sig_aldex2 <- H_aldex2_da %>%
  filter(wi.eBH < 0.05) # no significantly different taxa

A_sig_aldex2 <- A_aldex2_da %>%
  filter(wi.eBH < 0.05) # no significantly different taxa

# setup tax table to be able to merge
AR_taxa_info <- data.frame(tax_table(A_rumen_counts))


AR_taxa_info <- AR_taxa_info %>% rownames_to_column(var = "OTU")


# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(AR_aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)



# make a table of significant corrected p-values
AR_sig_aldex2 <- AR_aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)


# add in previously formed taxa information to complete the table
AR_sig_aldex2 <- left_join(AR_sig_aldex2, AR_taxa_info)


## DRA with regards to time ----

# subset out H0~H2 and H0~H18

quick <- c("H0", "H2")
long <- c("H0", "H18")

hquick <- subset_samples(
  H_rumen_counts,
  Hour %in% quick)

hlong <- subset_samples(
  H_rumen_counts,
  Hour %in% long
)

aquick <- subset_samples(
  A_rumen_counts,
  Hour %in% quick
)

along <- subset_samples(
  A_rumen_counts,
  Hour %in% long
)

# transpose OTU to make ALDEx2 happy
hq3 <-t(otu_table(hquick))
hl3 <- t(otu_table(hlong))
aq3 <- t(otu_table(aquick))
al3 <- t(otu_table(along))

# run AlDEx2 function
hq_aldex2_da <- ALDEx2::aldex(data.frame(hq3), phyloseq::sample_data(hquick)$Hour, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")
hl_aldex2_da <- ALDEx2::aldex(data.frame(hl3), phyloseq::sample_data(hlong)$Hour, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")
aq_aldex2_da <- ALDEx2::aldex(data.frame(aq3), phyloseq::sample_data(aquick)$Hour, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")
al_aldex2_da <- ALDEx2::aldex(data.frame(al3), phyloseq::sample_data(along)$Hour, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")


# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
hq_sig_aldex2 <- hq_aldex2_da %>%
  filter(wi.eBH < 0.05) # one (1) significantly different taxa

hl_sig_aldex2 <- hl_aldex2_da %>%
  filter(wi.eBH < 0.05) # one (1) significantly different taxa

aq_sig_aldex2 <- aq_aldex2_da %>%
  filter(wi.eBH < 0.05) # no significantly different taxa

al_sig_aldex2 <- al_aldex2_da %>%
  filter(wi.eBH < 0.05) # no significantly different taxa

# setup tax table to be able to merge
hq_taxa_info <- data.frame(tax_table(hquick))


hq_taxa_info <- hq_taxa_info %>% rownames_to_column(var = "OTU")


# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(hq_aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)



# make a table of significant corrected p-values
hq_sig_aldex2 <- hq_aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)


# add in previously formed taxa information to complete the table
hq_sig_aldex2 <- left_join(hq_sig_aldex2, hq_taxa_info) # Lactobacillus

# ---------------------------------------------------------------
  
# setup tax table to be able to merge
hl_taxa_info <- data.frame(tax_table(hlong))


hl_taxa_info <- hl_taxa_info %>% rownames_to_column(var = "OTU")


# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(hl_aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)



# make a table of significant corrected p-values
hl_sig_aldex2 <- hl_aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)


# add in previously formed taxa information to complete the table
hl_sig_aldex2 <- left_join(hl_sig_aldex2, hl_taxa_info) # Corynebacterium
