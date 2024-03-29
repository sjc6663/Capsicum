# Repeated Measures - Alpha Diversity Rumen Time Analysis
# April 11, 2023 - SAB

# setup ----
# load phyloseq objects
load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")
load("ps-obj/phyloseq-rumen-samples-only-counts.RData")

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(tidyr)
library(dplyr)
library(rstatix)

# create data frames ----
adivA <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(A_rumen_counts)$Treatment,
  "Hour" = phyloseq::sample_data(A_rumen_counts)$Hour,
  "SteerID" = phyloseq::sample_data(A_rumen_counts)$Steer.ID
)

adivH <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_rumen_counts)$Treatment,
  "Hour" = phyloseq::sample_data(H_rumen_counts)$Hour,
  "SteerID" = phyloseq::sample_data(H_rumen_counts)$Steer.ID
)

testHr <- aov(Shannon ~ Treatment, data = adivH)
summary(testHr)

testAr <- aov(Shannon ~ Treatment, data = adivA)
summary(testAr)

# repeated measures anova ----
# https://m-clark.github.io/docs/mixedModels/anovamixed.html 

anovaA <- aov(Shannon ~ Treatment*Hour + Error(SteerID), data = adivA)
summary(anovaA) #  F = 1.090, P = 0.386, Df = 12

anovaH <- aov(Shannon ~ Treatment*Hour + Error(SteerID), data = adivH)
summary(anovaH) #  F = 0.969, P = 0.488, Df = 12

# hour came up significant so I checked just for curiositys sake
testH <- aov(Shannon ~ Hour, data = adivH)
summary(testH) # p = 0.0044**

testA <- aov(Shannon ~ Hour, data = adivA)
summary(testA) # p = 0.492

TukeyHSD(testH)

## box plots -----
sample_data(H_rumen_counts)$"Hour" <- factor(sample_data(H_rumen_counts)$"Hour", 
                                             levels = c("H0", "H2", "H6", "H12", "H18")) 
sample_data(H_rumen_counts)$"Treatment" <- factor(sample_data(H_rumen_counts)$"Treatment", 
                                                  levels = c("Control", "RPC5", "RPC10", "RPC15"))

A <- plot_richness(H_rumen_counts, x="Hour", measures=c("Shannon"), title = "A", color = "Treatment") + 
  geom_boxplot() + 
  geom_jitter() +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1")) + 
  theme_classic() + 
  theme(legend.position = "none")


sample_data(A_rumen_counts)$"Hour" <- factor(sample_data(A_rumen_counts)$"Hour", 
                                             levels = c("H0", "H2", "H6", "H12", "H18"))
sample_data(A_rumen_counts)$"Treatment" <- factor(sample_data(A_rumen_counts)$"Treatment", 
                                                  levels = c("Control", "RPC5", "RPC10", "RPC15"))

B <- plot_richness(A_rumen_counts, x="Hour", measures=c("Shannon"), title = "B", color = "Treatment") + 
  geom_boxplot() + 
  geom_jitter() +
  scale_color_manual(values = c("#c5d280", "#ffc1b0", "#fdde9c", "#80cdc1")) + 
  theme_classic() 
 # theme(legend.position = "none")

A|B

ggsave(filename = "plots/paper/figure-7.tiff", dpi = 300, width = 8, height = 7)


