## DRA with regards to time ----
library(microViz)
library(phyloseq)
library(ALDEx2)
library(dplyr)
library(tibble)

set.seed(81299)

# load phyloseq of counts
load("ps-obj/phyloseq-fecal-samples-angus-counts.RData")
load("ps-obj/phyloseq-fecal-samples-holstein-counts.RData")
load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")

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
  filter(wi.eBH < 0.05) # no significantly different taxa

aq_sig_aldex2 <- aq_aldex2_da %>%
  filter(wi.eBH < 0.05) # no significantly different taxa

al_sig_aldex2 <- al_aldex2_da %>%
  filter(wi.eBH < 0.05) # no significantly different taxa

# setup tax table to be able to merge
hq_taxa_info <- data.frame(tax_table(hquick))

hq_taxa_info <- hq_taxa_info %>% rownames_to_column(var = "OTU")
hq_sig_aldex2 <- hq_sig_aldex2 %>% rownames_to_column(var = "OTU")

# make a table of significant corrected p-values
hq_sig_aldex2 <- hq_aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)


# add in previously formed taxa information to complete the table
hq_sig_aldex2 <- left_join(hq_sig_aldex2, hq_taxa_info) # Lactobacillus


# look at plot to see if positive/negative effect size corresponds to which treatment level
#ALDEx2::aldex.plot(hq_aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)


# -------------------------------------------------------------

ggplot(data=hq_sig_aldex2, aes(x=Phylum, y=(effect), col=Genus, label=Genus)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  #scale_color_manual() +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black") +
  theme(legend.position = "none")
