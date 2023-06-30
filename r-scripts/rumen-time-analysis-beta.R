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
  ord_plot(color = "Treatment", shape = "Hour", plot_taxa = 1:5, size = 4, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Treatment, color = Treatment)) +
  theme_classic() +
  ggtitle("A") +
  #theme(legend.position = "none") +
  labs(caption = "")

B <- A_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Treatment", shape = "Hour", plot_taxa = 1:5, size = 4, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Treatment, color = Treatment)) +
  theme_classic() +
  ggtitle("B") +
  #theme(legend.position = "none") +
  labs(caption = "")

c <- H_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", shape = "Treatment", plot_taxa = 1:5, size = 4, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("C") +
  #theme(legend.position = "none") +
  labs(caption = "")

d <- A_rumen_rel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", shape = "Treatment", plot_taxa = 1:5, size = 4, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("D") +
  #theme(legend.position = "none") +
  labs(caption = "")

(A|B)/(c|d)

ggsave(filename = "plots/PCA-time-stuff.pdf", dpi = 600, width = 12, height = 12)

# control capsicum time plots ----
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


## ---- Setup Treatments to be Control or Capsicum ----
A_rumen_rel <- A_rumen_rel %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

H_rumen_rel <- H_rumen_rel %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

H_con <- subset_samples(
  H_rumen_rel,
  SampleBinary == "Control"
)

H_cap <- subset_samples(
  H_rumen_rel,
  SampleBinary == "Capsicum"
)

# PCA plot for all hours - Holstein Control # these plots show the change over time for control steers and capsicum steers
conh <- H_con %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 4, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("A") +
  theme(legend.position = "none") +
  labs(caption = "")
conh

# PCA plot for all hours - Holstein Capsicum
caph <- H_cap %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 4, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("B") +
  theme(legend.position = "none") +
  labs(caption = "")
caph

conh|caph
hplots <- ggarrange(conh, caph, common.legend = TRUE, legend = "right")
hplots​

pdf("hplots.pdf")
hplots
dev.off()

ggsave(plot = hplots, filename = "plots/H-cap-con.pdf", dpi = 600, width = 12, height = 6)

A_con <- subset_samples(
  A_rumen_rel,
  SampleBinary == "Control"
)

A_cap <- subset_samples(
  A_rumen_rel,
  SampleBinary == "Capsicum"
)

# PCA plot for all hours - Holstein Control
cona <- A_con %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 4, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("C") +
  theme(legend.position = "none") +
  labs(caption = "")
cona

# PCA plot for all hours - Holstein Capsicum
capa <- A_cap %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 4, tax_lab_style = tax_lab_style(type = "text", size = 3, fontface = "bold.italic", check_overlap = TRUE)) +
  scale_color_manual(values =c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1", "#496e00")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("D") +
  # theme(legend.position = "none") +
  labs(caption = "")
capa

cona|capa

(conh|caph)/(cona|capa)

aplots <- ggarrange(cona, capa, common.legend = TRUE, legend = "right")

ggsave(plot = aplots, filename = "plots/PCA-angus-rumen-all.pdf", dpi = 600, width = 12, height = 6)

all <- ggarrange(conh, caph, cona, capa,
                              ncol = 2, nrow = 2, 
                              common.legend = TRUE, legend = "bottom")

all

ggsave(plot = all, filename = "plots/PCA-hours-all.pdf", dpi = 600, width = 12, height = 15)

## ---- H0 vs H2 (Holstein) ----

# define the quick hour segment
quickHour <- c("H0", "H2")

# define the long hour segment
longHour <- c("H0", "H18")

# subset the samples
H_quick <- subset_samples(
  H_rumen_rel, 
  Hour %in% quickHour)

H_trans_q <- H_quick %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
H_dmq <- phyloseq::distance(H_trans_q, method = "euclidean")

#ADONIS test
vegan::adonis2(H_dmq ~ phyloseq::sample_data(H_trans_q)$Hour) # R2 = 0.11, F = 3.80, P = 0.004**, significant difference between H0 and H2

# PCA Plot
HQ <- H_quick %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("HQ") +
  # theme(legend.position = "none") +
  labs(caption = "")
HQ

## ---- H0 vs H18 (Holstein) ----

H_long <- subset_samples(
  H_rumen_rel, 
  Hour %in% longHour
  )

H_trans_l <- H_long %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
H_dml <- phyloseq::distance(H_trans_l, method = "euclidean")

#ADONIS test
vegan::adonis2(H_dml ~ phyloseq::sample_data(H_trans_l)$Hour) # R2 = 0.102, F = 3.43, P = 0.005**, significant difference between H0 and H18

# PCA Plot
HL <- H_long %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("HL") +
  # theme(legend.position = "none") +
  labs(caption = "")
HL
## ---- H0 vs H2 (Angus) ----

A_quick <- subset_samples(
  A_rumen_rel, 
  Hour %in% quickHour
  )

A_trans_q <- A_quick %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
A_dmq <- phyloseq::distance(A_trans_q, method = "euclidean")

#ADONIS test
vegan::adonis2(A_dmq ~ phyloseq::sample_data(A_trans_q)$Hour) # R2 = 0.083, F = 2.71, P = 0.008**, significant difference between H0 and H2

# PCA Plot
AQ <- A_quick %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("AQ") +
  # theme(legend.position = "none") +
  labs(caption = "")
AQ
## ---- H0 vs H18 (Angus) ----

A_long <- subset_samples(
  A_rumen_rel, 
  Hour %in% longHour
  )

A_trans_l <- A_long %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
A_dml <- phyloseq::distance(A_trans_l, method = "euclidean")

#ADONIS test
vegan::adonis2(A_dml ~ phyloseq::sample_data(A_trans_l)$Hour) # R2 = 0.076, F = 2.47, P = 0.015*, significant difference between H0 and H18

# PCA Plot
AL <- A_long %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("AL") +
  # theme(legend.position = "none") +
  labs(caption = "")
AL
## ---- H0 vs H12 (Angus) ----

A_med <- subset_samples(
  A_rumen_rel,
  Hour == c("H0", "H12")
)

A_trans_m <- A_med %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
A_dmm <- phyloseq::distance(A_trans_m, method = "euclidean")

#ADONIS test
vegan::adonis2(A_dmm ~ phyloseq::sample_data(A_trans_m)$Hour) # R2 = 0.15, F = 1.95, P = 0.031*, significant difference between H0 and H12

# PCA Plot
AM <- A_med %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("AM") +
  # theme(legend.position = "none") +
  labs(caption = "")
AM
## ---- H0 vs H6 (Angus) ----

A_2 <- subset_samples(
  A_rumen_rel,
  Hour == c("H0", "H6")
)

A_trans_2 <- A_2 %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
A_dm2 <- phyloseq::distance(A_trans_2, method = "euclidean")

#ADONIS test
vegan::adonis2(A_dm2 ~ phyloseq::sample_data(A_trans_2)$Hour) # R2 = 0.165, F = 2.57, P = 0.006**, significant difference between H0 and H6

# PCA Plot
A6 <- A_2 %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("A6") +
  # theme(legend.position = "none") +
  labs(caption = "")
A6

## ---- H0 vs H12 (Holstein) ----

H_med <- subset_samples(
  H_rumen_rel,
  Hour == c("H0", "H12")
)

H_trans_m <- H_med %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
H_dmm <- phyloseq::distance(H_trans_m, method = "euclidean")

#ADONIS test
vegan::adonis2(H_dmm ~ phyloseq::sample_data(H_trans_m)$Hour) # R2 = 0.09, F = 1.64, P = 0.093, no significant difference between H0 and H12

# PCA Plot
HM <- H_med %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("HM") +
  # theme(legend.position = "none") +
  labs(caption = "")
HM

## ---- H0 vs H6 (Holstein) ----

H_2 <- subset_samples(
  H_rumen_rel,
  Hour == c("H0", "H6")
)

H_trans_2 <- H_2 %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

# generate distance matrix
H_dm2 <- phyloseq::distance(H_trans_2, method = "euclidean")

#ADONIS test
vegan::adonis2(H_dm2 ~ phyloseq::sample_data(H_trans_2)$Hour) # R2 = 0.20, F = 4.00, P = 0.002**, significant difference between H0 and H6

# PCA Plot
H6 <- H_2 %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("H6") +
  # theme(legend.position = "none") +
  labs(caption = "")
H6