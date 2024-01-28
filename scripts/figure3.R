# Relative Abundance Time Plots - Rumen
# SAB 3/16/2023

## ---- setup ----
# load necessary packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(plyr)
library(patchwork)
library(microshades)
library(viridis)
library(microbiome)
library(stringr)
library(tibble)


load("ps-obj/phyloseq-rumen-samples-only-relabund.RData")

# check if we have NAs
anyNA(tax_table(rumen_rel)[,"Phylum"])

# tax fix our phyloseq object
rumen_rel <- tax_fix(rumen_rel)


rumen_rel <- rumen_rel %>% tax_fix(unknowns = c("Incertae Sedis"))

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
# Control Color Palette ----
hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Genus"

# Sort phyloseq at lower, and then higher ranks
RC2 <- RC %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" genuses will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(RC2)[, hueRank]),
  shade = as.vector(tt_get(RC2)[, shadeRank]),
  counts = taxa_sums(otu_get(RC2))
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

#tax_palette_plot(hierarchicalPal) +
 # theme(axis.text.y.left = element_text(family = "mono"))

# RPC5 Color Palette ----
hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Genus"

# Sort phyloseq at lower, and then higher ranks
R52 <- R5 %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" genuses will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(R52)[, hueRank]),
  shade = as.vector(tt_get(R52)[, shadeRank]),
  counts = taxa_sums(otu_get(R52))
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

#tax_palette_plot(hierarchicalPal) +
 # theme(axis.text.y.left = element_text(family = "mono"))

# RPC10 Color Palette ----
hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Genus"

# Sort phyloseq at lower, and then higher ranks
R102 <- R10 %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" genuses will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(R102)[, hueRank]),
  shade = as.vector(tt_get(R102)[, shadeRank]),
  counts = taxa_sums(otu_get(R102))
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

#tax_palette_plot(hierarchicalPal) +
 # theme(axis.text.y.left = element_text(family = "mono"))

# RPC15 Color Palette ----
hueRank <- "Phylum"
hueRankPlural <- "Phyla"
shadeRank <- "Genus"

# Sort phyloseq at lower, and then higher ranks
R152 <- R15 %>%
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" genuses will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(R152)[, hueRank]),
  shade = as.vector(tt_get(R152)[, shadeRank]),
  counts = taxa_sums(otu_get(R152))
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

#tax_palette_plot(hierarchicalPal) +
 # theme(axis.text.y.left = element_text(family = "mono"))

# Control Plot ----

sample_data(RC2)$"Hour" <- factor(sample_data(RC2)$"Hour", 
                                        levels = c("H0", "H2", "H6", "H12", "H18"))

A <- RC2 %>%
  ps_get() %>%
  tax_mutate("Phylum: Genus" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Genus", sample_order = "asis",
    tax_order = "asis",
    merge_other = TRUE, x = "Hour", # x argument available since microViz 0.9.7,
    n_taxa = length(hierarchicalPal),
    palette = hierarchicalPal,
  ) +
  facet_wrap(
    facets = vars(Breed), labeller = as_labeller(~ paste("", .))
  ) +
  theme_bw() + # slightly clearer axes for facets
  labs(title = "Control", x = "Hour", y = "Relative Abundance", fill = "Genus") +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent()
  ) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 12))
  
#ggsave(filename = "plots/rel-abund-control-rumen.pdf", dpi = 600, width = 12, height = 10)


# RPC5 Plot ----

sample_data(R52)$"Hour" <- factor(sample_data(R52)$"Hour", 
                                  levels = c("H0", "H2", "H6", "H12", "H18"))

B <- R52 %>%
  ps_get() %>%
  tax_mutate("Phylum: Genus" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Genus", sample_order = "asis",
    tax_order = "asis",
    merge_other = TRUE, x = "Hour", # x argument available since microViz 0.9.7,
    n_taxa = length(hierarchicalPal),
    palette = hierarchicalPal
  ) +
  facet_wrap(
    facets = vars(Breed), labeller = as_labeller(~ paste("", .))
  ) +
  theme_bw() + # slightly clearer axes for facets
  labs(title = "RPC5", x = "Hour", y = "Relative Abundance", fill = "Phylum: Genus") +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent()
  ) +
  theme(text = element_text(size = 12))

#ggsave(filename = "plots/rel-abund-rpc5-rumen.pdf", dpi = 600)

# RPC10 Plot ----
sample_data(R102)$"Hour" <- factor(sample_data(R102)$"Hour", 
                                  levels = c("H0", "H2", "H6", "H12", "H18"))

C <- R102 %>%
  ps_get() %>%
  tax_mutate("Phylum: Genus" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Genus", sample_order = "asis",
    tax_order = "asis",
    merge_other = TRUE, x = "Hour", # x argument available since microViz 0.9.7,
    n_taxa = length(hierarchicalPal),
    palette = hierarchicalPal
  ) +
  facet_wrap(
    facets = vars(Breed), labeller = as_labeller(~ paste("", .))
  ) +
  theme_bw() + # slightly clearer axes for facets
  labs(title = "RPC10", x = "Hour", y = "Relative Abundance", fill = "Genus") +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent()
  ) + 
  theme(legend.position = "none") +
  theme(text = element_text(size = 12))

#ggsave(filename = "plots/rel-abund-rpc10-rumen.pdf", dpi = 600)

# RPC15 Plot ----
sample_data(R152)$"Hour" <- factor(sample_data(R152)$"Hour", 
                                   levels = c("H0", "H2", "H6", "H12", "H18"))

D <- R152 %>%
  ps_get() %>%
  tax_mutate("Phylum: Genus" = hierarchicalPalInfo$Taxa, .keep = "none") %>%
  comp_barplot(
    tax_level = "Phylum: Genus", sample_order = "asis",
    tax_order = "asis",
    merge_other = TRUE, x = "Hour", # x argument available since microViz 0.9.7,
    n_taxa = length(hierarchicalPal),
    palette = hierarchicalPal
  ) +
  facet_wrap(
    facets = vars(Breed), labeller = as_labeller(~ paste("", .))
  ) +
  theme_bw() + # slightly clearer axes for facets
  labs(title = "RPC15", x = "Hour", y = "Relative Abundance", fill = "Genus") +
  scale_y_continuous(
    expand = expansion(add = c(0, 0.1)), # axis starts exactly at 0
    labels = scales::label_percent()
  ) + 
  theme(legend.position = "none") +
  theme(text = element_text(size = 12))

#ggsave(filename = "plots/rel-abund-rpc15-rumen.pdf", dpi = 600)


plots <- (A|B)/(C|D)
plots

ggsave(filename = "plots/paper/figure-3.tiff", dpi = 300, width = 11, height = 10)
