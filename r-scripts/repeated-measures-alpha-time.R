# Repeated Measures - Alpha Diversity Rumen Time Analysis
# April 11, 2023 - SAB

# load phyloseq objects
load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(tidyr)
library(dplyr)

# get Shannon diversity
adivA <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(A_rumen_counts)$Treatment,
  "Hour" = phyloseq::sample_data(A_rumen_counts)$Hour,
  "SteerID" = phyloseq::sample_data(A_rumen_counts)$Steer.ID
)

sample_data(A_rumen_counts)$"Hour" <- factor(sample_data(A_rumen_counts)$"Hour", 
                                                 levels = c("H0", "H2", "H6", "H12", "H18"))
sample_data(A_rumen_counts)$"Treatment" <- factor(sample_data(A_rumen_counts)$"Treatment", 
                                                  levels = c("Control", "RPC5", "RPC10", "RPC15"))

B <- plot_richness(A_rumen_counts, x="Hour", measures=c("Shannon"), title = "B", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) + 
  theme_classic() 
  theme(legend.position = "none")

adivA <- rownames_to_column(adivA, var = "SampleID")

# repeated measures anova

# mod <- lmer(Shannon ~ Hour*Treatment + (1|SteerID), data = adivA)
# Anova(mod, test = "F", type = "III")

anovaA <- aov(Shannon ~ Treatment*Hour + Error(SteerID), data = adivA)
summary(anovaA) #  F = 1.090, P = 0.386, Df = 12

adivH <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_rumen_counts)$Treatment,
  "Hour" = phyloseq::sample_data(H_rumen_counts)$Hour,
  "SteerID" = phyloseq::sample_data(H_rumen_counts)$Steer.ID
)

anovaH <- aov(Shannon ~ Treatment*Hour + Error(SteerID), data = adivH)
summary(anovaH) #  F = 0.969, P = 0.488, Df = 12

sample_data(H_rumen_counts)$"Hour" <- factor(sample_data(H_rumen_counts)$"Hour", 
                                             levels = c("H0", "H2", "H6", "H12", "H18")) 
sample_data(H_rumen_counts)$"Treatment" <- factor(sample_data(H_rumen_counts)$"Treatment", 
                                                  levels = c("Control", "RPC5", "RPC10", "RPC15"))

A <- plot_richness(H_rumen_counts, x="Hour", measures=c("Shannon"), title = "A", color = "Treatment") + 
  geom_boxplot() + 
  scale_color_manual(values = c("#a76119", "#028571", "#dfc27d", "#80cdc1")) + 
  theme_classic() + 
  theme(legend.position = "none")

testH <- aov(Shannon ~ Hour, data = adivH)
summary(testH) # p = 0.0865

A|B

ggsave(filename = "plots/alpha-div-rumen-all.pdf", dpi = 600)

# break up individual to test ====
RPC5 <- subset_samples(A_rumen_counts,
                       Treatment == "RPC5")

adiv5 <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(RPC5, measures = "Shannon"),
  "Hour" = phyloseq::sample_data(RPC5)$Hour,
  "SteerID" = phyloseq::sample_data(RPC5)$Steer.ID
)

test5 <- aov(Shannon ~ Hour, data = adiv5)
summary(test5) # p = 0.25


RPC15 <- subset_samples(A_rumen_counts,
                       Treatment == "RPC15")

adiv15 <- data.frame(
  # "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(RPC15, measures = "Shannon"),
  "Hour" = phyloseq::sample_data(RPC15)$Hour,
  "SteerID" = phyloseq::sample_data(RPC15)$Steer.ID
)

test15 <- aov(Shannon ~ Hour, data = adiv15)
summary(test15) # p = 0.374

