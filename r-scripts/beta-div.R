## BETA DIVERSITY
# SAB 11-03-2022

## ---- setup ----
# load packages
library(BiocManager)
library(phyloseq)
library(vegan)
library(patchwork)
library(pairwiseAdonis)

# load phyloseq objects
load("ps-obj/phyloseq-fecal-samples-angus-relabund.RData")
load("ps-obj/phyloseq-fecal-samples-holstein-relabund.RData")
load("ps-obj/phyloseq-rumen-samples-angus-relabund.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-relabund.RData")

# check if we have NAs
anyNA(tax_table(H_fecal_rel)[,"Phylum"])
anyNA(tax_table(H_rumen_rel)[,"Phylum"])
anyNA(tax_table(A_fecal_rel)[,"Phylum"])
anyNA(tax_table(A_rumen_rel)[,"Phylum"])

# tax fix our phyloseq object
A_fecal_rel <- tax_fix(A_fecal_rel)
A_rumen_rel <- tax_fix(A_rumen_rel)
H_fecal_rel <- tax_fix(H_fecal_rel)
H_rumen_rel <- tax_fix(H_rumen_rel)

A_fecal_rel <- A_fecal_rel %>% tax_fix(unknowns = c("Incertae Sedis"))
A_rumen_rel <- A_rumen_rel %>% tax_fix(unknowns = c("Incertae Sedis"))
H_fecal_rel <- H_fecal_rel %>% tax_fix(unknowns = c("Incertae Sedis"))
H_rumen_rel <- H_rumen_rel %>% tax_fix(unknowns = c("Incertae Sedis"))

# transform
H_transform_fecal <- microbiome::transform(H_fecal_rel, 'clr')
H_transform_rumen <- microbiome::transform(H_rumen_rel, 'clr')
A_transform_fecal <- microbiome::transform(A_fecal_rel, 'clr')
A_transform_rumen <- microbiome::transform(A_rumen_rel, 'clr')

## ---- PERMANOVA @ Genus Level ----
# clr transform phyloseq objects at Genus level
H_trans_fecal <- H_fecal_rel %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

A_trans_fecal <- A_fecal_rel %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

A_trans_rumen <- A_rumen_rel %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

H_trans_rumen <- H_rumen_rel %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
H_fecal_dist_matrix_trans <- phyloseq::distance(H_trans_fecal, method = "euclidean")
A_fecal_dist_matrix_trans <- phyloseq::distance(A_trans_fecal, method = "euclidean")
A_rumen_dist_matrix_trans <- phyloseq::distance(A_trans_rumen, method = "euclidean")
H_rumen_dist_matrix_trans <- phyloseq::distance(H_trans_rumen, method = "euclidean")

#ADONIS test
vegan::adonis2(H_fecal_dist_matrix_trans ~ phyloseq::sample_data(H_trans_fecal)$Treatment) # p = 0.001
vegan::adonis2(A_fecal_dist_matrix_trans ~ phyloseq::sample_data(A_trans_fecal)$Treatment) # p = 0.384
vegan::adonis2(H_rumen_dist_matrix_trans ~ phyloseq::sample_data(H_trans_rumen)$Treatment)
vegan::adonis2(A_rumen_dist_matrix_trans ~ phyloseq::sample_data(A_trans_rumen)$Treatment)

# Pairwise Adonis for Holstein Samples
H_df <- pairwise.adonis(H_fecal_dist_matrix_trans, sample_data(H_trans_fecal)$Treatment)
H_df

#write.csv(H_df, file = "tables/H-pairwise-adonis-hour.csv")


## ---- PCA Plots Holstein and Angus ----
##  PCA plot - fecal - holstein
A <- H_fecal_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Treatment", shape = "Treatment", plot_taxa = 1:5, size = 2) +
  scale_colour_brewer(palette = "BrBG") +
  stat_ellipse(aes(group = Treatment, color = Treatment)) +
  theme_classic() +
  ggtitle("A") +
  theme(legend.position = "none") +
  labs(caption = "")

ggsave(filename = "plots/PCA-microViz-holstein-fecal.pdf", dpi = 600)

##  PCA plot - fecal - angus
B <- A_fecal_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Treatment", shape = "Treatment", plot_taxa = 1:5, size = 2) +
  scale_colour_brewer(palette = "BrBG") +
  stat_ellipse(aes(group = Treatment, color = Treatment)) + 
  theme_classic() +
  ggtitle("B") + 
  labs(caption = "")

ggsave(filename = "plots/PCA-microViz-angus-fecal.pdf", dpi = 600)

##  PCA plot - fecal - treatment
C <- H_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Treatment", shape = "Treatment", plot_taxa = 1:5, size = 2) +
  scale_colour_brewer(palette = "BrBG") +
  stat_ellipse(aes(group = Treatment, color = Treatment)) + 
  theme_classic() +
  ggtitle("C") +
  theme(legend.position = "none") +
  labs(caption = "")

ggsave(filename = "plots/PCA-microViz-rumen-holstein.pdf", dpi = 600)

##  PCA plot - fecal - treatment
D <- A_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Treatment", shape = "Treatment", plot_taxa = 1:5, size = 2) +
  scale_colour_brewer(palette = "BrBG") +
  stat_ellipse(aes(group = Treatment, color = Treatment)) + 
  theme_classic() +
  ggtitle("D") +
  labs(caption = "")

D
ggsave(filename = "plots/PCA-microViz-angus-rumen.pdf", dpi = 600)

# (A|B)/(C|D)

all <- ggarrange(A, B, C, D,
                 ncol = 2, nrow = 2, 
                 common.legend = TRUE, legend = "bottom")

all

ggsave(filename = "plots/PCA-all-samples.pdf", dpi = 600, width = 12, height = 16)

#A|B
#ggsave(filename = "plots/PCA-all-fecal-samples.pdf", dpi = 600)
