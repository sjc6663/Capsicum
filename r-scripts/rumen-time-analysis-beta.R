## Rumen Time Beta Analysis
## SAB - Last Edited 03-02-2023

# load packages
library(BiocManager)
library(phyloseq)
library(vegan)
library(patchwork)

# load phyloseq objects
load("ps-obj/phyloseq-rumen-samples-angus-relabund.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-relabund.RData")

# check if we have NAs
anyNA(tax_table(H_rumen_rel)[,"Phylum"])
anyNA(tax_table(A_rumen_rel)[,"Phylum"])

# tax fix our phyloseq object
A_rumen_rel <- tax_fix(A_rumen_rel)
H_rumen_rel <- tax_fix(H_rumen_rel)


A_rumen_rel <- A_rumen_rel %>% tax_fix(unknowns = c("Incertae Sedis"))
H_rumen_rel <- H_rumen_rel %>% tax_fix(unknowns = c("Incertae Sedis"))

# put the time series in order
sample_data(A_rumen_rel)$"Treatment" <- factor(sample_data(A_rumen_rel)$"Hour", 
                                                  levels = c("H0", "H2", "H6", "H12", "H18"))
sample_data(H_rumen_rel)$"Treatment" <- factor(sample_data(H_rumen_rel)$"Hour", 
                                               levels = c("H0", "H2", "H6", "H12", "H18"))

# clr transform phyloseq objects at Genus level
H_trans <- H_rumen_rel %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

A_trans <- A_rumen_rel %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
H_dm <- phyloseq::distance(H_trans, method = "euclidean")
A_dm <- phyloseq::distance(A_trans, method = "euclidean")

#ADONIS test
vegan::adonis2(H_dm ~ phyloseq::sample_data(H_trans)$Hour) # p = 0.001***
vegan::adonis2(A_dm ~ phyloseq::sample_data(A_trans)$Hour) # p = 0.002**

##  PCA plot - fecal - holstein
A <- H_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_colour_brewer(palette = "BrBG") +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("A") +
  theme(legend.position = "none") +
  labs(caption = "")
A
ggsave(filename = "plots/PCA-microViz-holstein-rumen.pdf", dpi = 600)

B <- A_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_colour_brewer(palette = "BrBG") +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("B") +
 # theme(legend.position = "none") +
  labs(caption = "")
B
ggsave(filename = "plots/PCA-microViz-holstein-rumen.pdf", dpi = 600)

A|B
ggsave(filename = "plots/PCA-rumen-all.pdf", dpi = 600)
