## Taxa Relative Abundance
# SAB 11/9/2022

## ---- setup ----
# load necessary packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(plyr)
library(patchwork)
library(tibble)

# load phyloseq objects
psrel <- readRDS("ps-obj/ps-decontam-relabund.rds")
#pscount <- readRDS("ps-decontam-filtered-counts.rds")

psrel

# merge everything at the Genus level so we don't have duplicates
psgrel <- tax_glom(psrel, "Genus")
#psgcount <- tax_glom(pscount, "Genus")

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

# Control Color Palette ----
hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Genus"

# Sort phyloseq at lower, and then higher ranks
FC2 <- FC %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" genuses will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(FC2)[, hueRank]),
  shade = as.vector(tt_get(FC2)[, shadeRank]),
  counts = taxa_sums(otu_get(FC2))
)

hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )

matrix <- c("#578767", "#ffc2b0", "#00534B", "#4D4D4D", "#78aa88", "#ffcca0", "#2B7A70", "#969696", "#99c086", "#ffd699", "#56A397", "#C3C3C3", "#c6d281", "#fdde9c", "#80CDC1", "#E6E6E6")


hierarchicalPalMatrix <- matrix(matrix, nrow = 4, ncol = 4, byrow = TRUE)

hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))

hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))

tax_palette_plot(hierarchicalPal) +
  theme(axis.text.y.left = element_text(family = "mono"))

# RPC5 Color Palette ----
hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Genus"

# Sort phyloseq at lower, and then higher ranks
F52 <- F5 %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" genuses will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(F52)[, hueRank]),
  shade = as.vector(tt_get(F52)[, shadeRank]),
  counts = taxa_sums(otu_get(F52))
)

hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )

matrix <- c("#578767", "#ffc2b0", "#00534B", "#4D4D4D", "#78aa88", "#ffcca0", "#2B7A70", "#969696", "#99c086", "#ffd699", "#56A397", "#C3C3C3", "#c6d281", "#fdde9c", "#80CDC1", "#E6E6E6")


hierarchicalPalMatrix <- matrix(matrix, nrow = 4, ncol = 4, byrow = TRUE)

hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))

hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))

tax_palette_plot(hierarchicalPal) +
  theme(axis.text.y.left = element_text(family = "mono"))

# RPC10 Color Palette ----
hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Genus"

# Sort phyloseq at lower, and then higher ranks
F102 <- F10 %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" genuses will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(F102)[, hueRank]),
  shade = as.vector(tt_get(F102)[, shadeRank]),
  counts = taxa_sums(otu_get(F102))
)

hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )

matrix <- c("#578767", "#ffc2b0", "#00534B", "#4D4D4D", "#78aa88", "#ffcca0", "#2B7A70", "#969696", "#99c086", "#ffd699", "#56A397", "#C3C3C3", "#c6d281", "#fdde9c", "#80CDC1", "#E6E6E6")


hierarchicalPalMatrix <- matrix(matrix, nrow = 4, ncol = 4, byrow = TRUE)

hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))

hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))

tax_palette_plot(hierarchicalPal) +
  theme(axis.text.y.left = element_text(family = "mono"))

# RPC15 Color Palette ----
hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Genus"

# Sort phyloseq at lower, and then higher ranks
F152 <- F15 %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" genuses will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(F152)[, hueRank]),
  shade = as.vector(tt_get(F152)[, shadeRank]),
  counts = taxa_sums(otu_get(F152))
)

hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )

matrix <- c("#578767", "#ffc2b0", "#00534B", "#4D4D4D", "#78aa88", "#ffcca0", "#2B7A70", "#969696", "#99c086", "#ffd699", "#56A397", "#C3C3C3", "#c6d281", "#fdde9c", "#80CDC1", "#E6E6E6")


hierarchicalPalMatrix <- matrix(matrix, nrow = 4, ncol = 4, byrow = TRUE)

hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))

hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))

tax_palette_plot(hierarchicalPal) +
  theme(axis.text.y.left = element_text(family = "mono"))

# Control Plot ----
A <- FC2 %>%
  ps_get() %>%
  tax_mutate("Phylum: Genus" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Genus", sample_order = "asis",
    tax_order = "asis",
    merge_other = TRUE, x = "Period", # x argument available since microViz 0.9.7,
    n_taxa = length(hierarchicalPal),
    palette = hierarchicalPal,
  ) +
  facet_wrap(
    facets = vars(Breed), labeller = as_labeller(~ paste("", .))
  ) +
  theme_bw() + # slightly clearer axes for facets
  labs(title = "Control", x = "Period", y = "Relative Abundance", fill = "Genus") +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent()
  ) +
  theme(legend.position = "none")

ggsave(plot = pC, filename = "plots/rel-abund-control-fecal.pdf", dpi = 600)


# RPC5 Plot ----
B <- F52 %>%
  ps_get() %>%
  tax_mutate("Phylum: Genus" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Genus", sample_order = "asis",
    tax_order = "asis",
    merge_other = TRUE, x = "Period", # x argument available since microViz 0.9.7,
    n_taxa = length(hierarchicalPal),
    palette = hierarchicalPal
  ) +
  facet_wrap(
    facets = vars(Breed), labeller = as_labeller(~ paste("", .))
  ) +
  theme_bw() + # slightly clearer axes for facets
  labs(title = "RPC5", x = "Period", y = "Relative Abundance", fill = "Phylum: Genus") +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent()
  ) 

ggsave(plot = p5, filename = "plots/rel-abund-rpc5-fecal.pdf", dpi = 600)

# RPC10 Plot ----

C <- F102 %>%
  ps_get() %>%
  tax_mutate("Phylum: Genus" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Genus", sample_order = "asis",
    tax_order = "asis",
    merge_other = TRUE, x = "Period", # x argument available since microViz 0.9.7,
    n_taxa = length(hierarchicalPal),
    palette = hierarchicalPal
  ) +
  facet_wrap(
    facets = vars(Breed), labeller = as_labeller(~ paste("", .))
  ) +
  theme_bw() + # slightly clearer axes for facets
  labs(title = "RPC10", x = "Period", y = "Relative Abundance", fill = "Genus") +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent()
  ) +
  theme(legend.position = "none")


ggsave(plot = p10, filename = "plots/rel-abund-rpc10-fecal.pdf", dpi = 600)

# RPC15 Plot ----

D <- F152 %>%
  ps_get() %>%
  tax_mutate("Phylum: Genus" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Genus", sample_order = "asis",
    tax_order = "asis",
    merge_other = TRUE, x = "Period", # x argument available since microViz 0.9.7,
    n_taxa = length(hierarchicalPal),
    palette = hierarchicalPal
  ) +
  facet_wrap(
    facets = vars(Breed), labeller = as_labeller(~ paste("", .))
  ) +
  theme_bw() + # slightly clearer axes for facets
  labs(title = "RPC15", x = "Period", y = "Relative Abundance", fill = "Genus") +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent()
  ) +
  theme(legend.position = "none")

ggsave(plot = p15, filename = "plots/rel-abund-rpc15-fecal.pdf", dpi = 600)

plots <- (A|B)/(C|D)
plots

ggsave(filename = "plots/rel-abund-fecal-all.jpeg", dpi = 600, width = 18, height = 12)
