## Rumen Time Beta Analysis
## SAB - Last Edited 03-02-2023

# ---- setup ----
# load packages
library(BiocManager)
library(phyloseq)
library(vegan)
library(patchwork)
library(microViz)
library(ggplot2)
library(dplyr)
library(stringr)

set.seed(081299)

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
sample_data(A_rumen_rel)$"Hour" <- factor(sample_data(A_rumen_rel)$"Hour", 
                                          levels = c("H0", "H2", "H6", "H12", "H18"))
sample_data(H_rumen_rel)$"Hour" <- factor(sample_data(H_rumen_rel)$"Hour", 
                                          levels = c("H0", "H2", "H6", "H12", "H18"))

H_rumen_rel <- H_rumen_rel %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

H_con <- subset_samples(H_rumen_rel,
                        SampleBinary == "Control"
)

H_cap <- subset_samples(H_rumen_rel,
                        SampleBinary == "Capsicum")

# give samples a column to descripe capsicum or no capsicum
A_rumen_rel <- A_rumen_rel %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

A_con <- subset_samples(A_rumen_rel,
                        SampleBinary == "Control")

A_cap <- subset_samples(A_rumen_rel,
                        SampleBinary == "Capsicum")

### ---- BETA DIVERSITY ----

# clr transform phyloseq objects at Genus level
H_trans <- H_rumen_rel %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

A_trans <- A_rumen_rel %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

H_trans_con <- H_con %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

H_trans_cap <- H_cap %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

A_trans_con <- A_con %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

A_trans_cap <- A_cap %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
H_dm <- phyloseq::distance(H_trans, method = "euclidean")
A_dm <- phyloseq::distance(A_trans, method = "euclidean")

H_con_dm <- phyloseq::distance(H_trans_con, method = "euclidean")
H_cap_dm <- phyloseq::distance(H_trans_cap, method = "euclidean")
A_con_dm <- phyloseq::distance(A_trans_con, method = "euclidean")
A_cap_dm <- phyloseq::distance(A_trans_cap, method = "euclidean")

#ADONIS test
vegan::adonis2(H_dm ~ phyloseq::sample_data(H_trans)$Treatment * phyloseq::sample_data(H_trans)$Hour) # R2 = 0.109, P = 1.00, ns
vegan::adonis2(A_dm ~ phyloseq::sample_data(A_trans)$Treatment * phyloseq::sample_data(A_trans)$Hour) # R2 = 0.113, P = 0.99, ns

vegan::adonis2(H_dm ~ phyloseq::sample_data(H_trans)$Hour) # R2 = 0.12, F = 2.63, P = 0.001
vegan::adonis2(A_dm ~ phyloseq::sample_data(A_trans)$Hour) # R2 = 0.08, F = 1.72, P = 0.005


## Pairwise Adonis By Hour ----

# pairwise adonis test
H_df <- pairwise.adonis(H_dm, sample_data(H_trans)$Hour)
H_df

write.csv(H_df, file = "tables/H-pairwise-adonis-hour.csv")

A_df <- pairwise.adonis(A_dm, sample_data(A_trans)$Hour)
A_df
write.csv(A_df, file = "tables/A-pairwise-adonis-hour.csv")

## plots ----
A <- H_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Treatment", shape = "Hour", plot_taxa = 1:5, size = 5, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Treatment, color = Treatment)) +
  theme_classic() +
  ggtitle("A") +
  theme(legend.position = "none") +
  labs(caption = "R2 = 0.109, P = 1.00") +
  theme(text = element_text(size = 20))

B <- A_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Treatment", shape = "Hour", plot_taxa = 1:5, size = 5, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Treatment, color = Treatment)) +
  theme_classic() +
  ggtitle("B") +
  #theme(legend.position = "none") +
  labs(caption = "R2 = 0.113, P = 0.99") +
  theme(text = element_text(size = 20))

c <- H_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", shape = "Treatment", plot_taxa = 1:5, size = 5, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("C") +
  theme(legend.position = "none") +
  labs(caption = "R2 = 0.12, P = 0.001***") +
  theme(text = element_text(size = 20))

d <- A_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", shape = "Treatment", plot_taxa = 1:5, size = 5, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("D") +
  #theme(legend.position = "none") +
  labs(caption = "R2 = 0.08, P = 0.005***") +
  theme(text = element_text(size = 20))

(A|B)/(c|d)

ggsave(filename = "plots/paper/figure-5.tiff", dpi = 300, width = 10, height = 10)