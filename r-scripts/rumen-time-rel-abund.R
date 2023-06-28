# Relative Abundance Time Plots - Rumen
# SAB 3/16/2023

## ---- setup ----
# load necessary packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(plyr)
library(patchwork)


load("ps-obj/phyloseq-rumen-samples-only-relabund.RData")

# check if we have NAs
anyNA(tax_table(rumen_rel)[,"Phylum"])

# tax fix our phyloseq object
rumen_rel <- tax_fix(rumen_rel)


rumen_rel <- rumen_rel %>% tax_fix(unknowns = c("Incertae Sedis"))

sample_data(rumen_rel)$"Hour" <- factor(sample_data(rumen_rel)$"Hour", 
                                        levels = c("H0", "H2", "H6", "H12", "H18"))

## Separate out by treatment ----

RC <- subset_samples(
  rumen_rel,
  Treatment == "Control"
)

R5 <- subset_samples(
  rumen_rel,
  Treatment == "RPC5"
)

R10 <- subset_samples(
  rumen_rel,
  Treatment == "RPC10"
)

R15 <- subset_samples(
  rumen_rel,
  Treatment == "RPC15"
)

## plot ----
# The palette with grey:
cbPalette <- c("#003C30", "#8C510A",  "#35978F", "#BF812D", "#DFC27D", "#80CDC1", "#999999", "#01665E",  "#F6E8C3", "#543005", "#F5F5F5")

pC <- ggplot(data = psmelt(RC), mapping = aes_string(x = "Hour", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  scale_fill_manual(values = cbPalette) +
  labs(x = NULL, y = "Relative Abundance",
       title = "Control") +
  theme_bw() +
  facet_grid(~Breed) +
  theme(legend.position = "none")


ggsave(filename = "plots/rel-abund-control-rumen.pdf", dpi = 600)

p5 <- ggplot(data = psmelt(R5), mapping = aes_string(x = "Hour", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  scale_fill_manual(values = cbPalette) +
  labs(x = NULL, y = "Relative Abundance",
       title = "RPC5") +
  theme_bw() +
  facet_grid(~Breed) 


ggsave(filename = "plots/rel-abund-rpc5-rumen.pdf", dpi = 600)

p10 <- ggplot(data = psmelt(R10), mapping = aes_string(x = "Hour", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  scale_fill_manual(values = cbPalette) +
  labs(x = NULL, y = "Relative Abundance",
       title = "RPC10") +
  theme_bw() +
  facet_grid(~Breed) +
  theme(legend.position = "none")


ggsave(filename = "plots/rel-abund-rpc10-rumen.pdf", dpi = 600)

p15 <- ggplot(data = psmelt(R15), mapping = aes_string(x = "Hour", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  scale_fill_manual(values = cbPalette) +
  labs(x = NULL, y = "Relative Abundance",
       title = "RPC15") +
  theme_bw() +
  facet_grid(~Breed) + 
  theme(legend.position = "none")


ggsave(filename = "plots/rel-abund-rpc15-rumen.pdf", dpi = 600)

(pC|p5)/(p10|p15)

ggsave(filename = "plots/rel-abund-rumen-all.jpeg", dpi = 600, width = 12, height = 12)
