## Taxa Relative Abundance
# SAB 11/9/2022

## ---- setup ----
# load necessary packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(plyr)
library(patchwork)

# load phyloseq objects
psrel <- readRDS("ps-decontam-relabund.rds")
pscount <- readRDS("ps-decontam-filtered-counts.rds")

psrel

# merge everything at the Genus level so we don't have duplicates
psgrel <- tax_glom(psrel, "Genus")
psgcount <- tax_glom(pscount, "Genus")

## ---- separate out Rumen and Fecal Samples from relative abundance set ----
fecal_rel <- subset_samples(
  psgrel, 
  Sample.Type == "Fecal"
)

fecal_rel
sample_data(fecal_rel)

rumen_rel <- subset_samples(
  psgrel,
  Sample.Type == "Rumen"
)

rumen_rel
sample_data(rumen_rel)

# save these separate phyloseq objects
save(fecal_rel, file = "ps-obj/phyloseq-fecal-samples-only-relabund.RData")
save(rumen_rel, file = "ps-obj/phyloseq-rumen-samples-only-relabund.RData")

## ---- Tax Fix ----
### From here on out you have to run doubles, so one for fecal and one for rumen ###

load("ps-obj/phyloseq-fecal-samples-only-relabund.RData")
# load("ps-obj/phyloseq-rumen-samples-only-relabund.RData")

# check if we have NAs
anyNA(tax_table(fecal_rel)[,"Phylum"])
# anyNA(tax_table(rumen_rel)[,"Phylum"])

# tax fix our phyloseq object
fecal_rel <- tax_fix(fecal_rel)
# rumen_rel <- tax_fix(rumen_rel)

fecal_rel <- fecal_rel %>% tax_fix(unknowns = c("Incertae Sedis"))
# rumen_rel <- rumen_rel %>% tax_fix(unknowns = c("Incertae Sedis"))

## Separate out by treatment ----

FC <- subset_samples(
  fecal_rel,
  Treatment == "Control"
)

F5 <- subset_samples(
  fecal_rel,
  Treatment == "RPC5"
)

F10 <- subset_samples(
 fecal_rel,
  Treatment == "RPC10"
)

F15 <- subset_samples(
  fecal_rel,
  Treatment == "RPC15"
)

## ---- Plot in ggplot2 ----
# The palette with grey:
cbPalette <- c("#003C30", "#8C510A",  "#35978F", "#BF812D", "#DFC27D", "#80CDC1", "#999999", "#01665E",  "#F6E8C3", "#543005", "#F5F5F5")

pC <- ggplot(data = psmelt(FC), mapping = aes_string(x = "Period", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  scale_fill_manual(values = cbPalette) +
  labs(x = NULL, y = "Relative Abundance",
       title = "Control") +
  theme_bw() +
  facet_grid(~Breed) 
  #theme(legend.position = "none")

pC

ggsave(plot = pC, filename = "plots/rel-abund-control-fecal.pdf", dpi = 600)

p5 <- ggplot(data = psmelt(F5), mapping = aes_string(x = "Period", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  scale_fill_manual(values = cbPalette) +
  labs(x = NULL, y = "Relative Abundance",
       title = "RPC5") +
  theme_bw() +
  facet_grid(~Breed) 

p5

ggsave(plot = p5, filename = "plots/rel-abund-rpc5-fecal.pdf", dpi = 600)

p10 <- ggplot(data = psmelt(F10), mapping = aes_string(x = "Period", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  scale_fill_manual(values = cbPalette) +
  labs(x = NULL, y = "Relative Abundance",
       title = "RPC10") +
  theme_bw() +
  facet_grid(~Breed) 
  #theme(legend.position = "none")

p10

ggsave(plot = p10, filename = "plots/rel-abund-rpc10-fecal.pdf", dpi = 600)

p15 <- ggplot(data = psmelt(F15), mapping = aes_string(x = "Period", y = "Abundance")) +
  geom_bar(aes(fill = Phylum), stat = "identity", position = "fill") +
  scale_fill_manual(values = cbPalette) +
  labs(x = NULL, y = "Relative Abundance",
       title = "RPC15") +
  theme_bw() +
  facet_grid(~Breed) 
  # theme(legend.position = "none")

p15

ggsave(plot = p15, filename = "plots/rel-abund-rpc15-fecal.pdf", dpi = 600)

(pC|p5)/(p10|p15)

ggsave(filename = "plots/rel-abund-fecal-all.pdf", dpi = 600, width = 12, height = 12)
