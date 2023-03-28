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



### ---- BETA DIVERSITY ----

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

# PCA plot for all hours - Holstein Control
conh <- H_con %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr", rank = "Genus") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c( "#028571", "#80cdc1", "#dfc27d", "#a76119", "#543005")) +
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
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c( "#028571", "#80cdc1", "#dfc27d", "#a76119", "#543005")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("B") +
  # theme(legend.position = "none") +
  labs(caption = "")
caph

conh|caph

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
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c( "#028571", "#80cdc1", "#dfc27d", "#a76119", "#543005")) +
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
  ord_plot(color = "Hour", plot_taxa = 1:5, size = 2) +
  scale_color_manual(values = c( "#028571", "#80cdc1", "#dfc27d", "#a76119", "#543005")) +
  stat_ellipse(aes(group = Hour, color = Hour)) +
  theme_classic() +
  ggtitle("D") +
  # theme(legend.position = "none") +
  labs(caption = "")
capa

cona|capa

(conh|caph)/(cona|capa)

ggsave(filename = "plots/PCA-hours-all.pdf", dpi = 600)

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




### ---- ALPHA DIVERSITY ----

load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")

# tax fix our phyloseq object
A_rumen_counts <- tax_fix(A_rumen_counts)
H_rumen_counts <- tax_fix(H_rumen_counts)

# individual tax fix
A_rumen_counts <- A_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
H_rumen_counts <- H_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))

## ---- create data frames for H0 and H2 (Holstein) ----
adivH <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_quick, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_quick)$Treatment,
  "Hour" = phyloseq::sample_data(H_quick)$Hour
)
head(adivH)
# write.csv(adivH,"tables/holstein-fecal.csv", row.names = FALSE)

H_quick <- subset_samples(
  H_rumen_counts,
  Hour == c("H0", "H2")
)

testH <- aov(Shannon ~ Hour, data = adivH)
summary(testH) # p = 0.11


sample_data(H_quick)$"Hour" <- factor(sample_data(H_quick)$"Hour", 
                                               levels = c("H0", "H2"))

A <- plot_richness(H_quick, x="Hour", measures=c("Shannon"), title = "A", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  theme_classic()
 # theme(legend.position = "none")
A

## ---- create data frames for H0 and H18 (Holstein) ----
adivHL <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_long, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_long)$Treatment,
  "Hour" = phyloseq::sample_data(H_long)$Hour
)
head(adivHL)
# write.csv(adivH,"tables/holstein-fecal.csv", row.names = FALSE)

H_long <- subset_samples(
  H_rumen_counts,
  Hour == c("H0", "H18")
)

testHL <- aov(Shannon ~ Hour, data = adivHL)
summary(testHL) # p = 0.502


sample_data(H_long)$"Hour" <- factor(sample_data(H_long)$"Hour", 
                                      levels = c("H0", "H18"))

B <- plot_richness(H_long, x="Hour", measures=c("Shannon"), title = "B", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  theme_classic()
# theme(legend.position = "none")
B

## ---- create data frames for H0 and H2 (Angus) ----
adivA <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_quick, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(A_quick)$Treatment,
  "Hour" = phyloseq::sample_data(A_quick)$Hour
)
head(adivA)
# write.csv(adivH,"tables/holstein-fecal.csv", row.names = FALSE)

A_quick <- subset_samples(
  A_rumen_counts,
  Hour == c("H0", "H2")
)

testA <- aov(Shannon ~ Hour, data = adivA)
summary(testA) # p = 0.509


sample_data(A_quick)$"Hour" <- factor(sample_data(A_quick)$"Hour", 
                                      levels = c("H0", "H2"))

C <- plot_richness(A_quick, x="Hour", measures=c("Shannon"), title = "C", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) +
  theme_classic()
# theme(legend.position = "none")
C

## ---- create data frames for H0 and H18 (Angus) ----
adivAL <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_long, measures = "Shannon"),
  "SampleBinary" = phyloseq::sample_data(A_long)$SampleBinary,
  "Hour" = phyloseq::sample_data(A_long)$Hour
)
head(adivAL)
# write.csv(adivH,"tables/holstein-fecal.csv", row.names = FALSE)

A_long <- subset_samples(
  A_rumen_counts,
  Hour == c("H0", "H18")
)

testAL <- aov(Shannon ~ Hour, data = adivAL)
summary(testAL) # p = 0.819


sample_data(A_long)$"Hour" <- factor(sample_data(A_long)$"Hour", 
                                     levels = c("H0", "H18"))

D <- plot_richness(A_long, x="Hour", measures=c("Shannon"), title = "D", color = "SampleBinary") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571")) +
  theme_classic()
# theme(legend.position = "none")
D

## ---- Setup Treatments to be Control or Capsicum ----
A_rumen_counts <- A_rumen_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

H_rumen_counts <- H_rumen_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 
