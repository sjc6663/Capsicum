## Differential Relative Abundance with MicroViz
# 3-28-2023 SAB

library(microViz)
library(corncob)
library(dplyr)
library(ggplot2)
library(skimr)

load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")

A_rumen_counts <- tax_fix(A_rumen_counts)


A_rumen_counts <- A_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))

A_rumen_counts <- A_rumen_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

A_rumen_counts %>% 
  samdat_tbl() %>% 
  dplyr::mutate(across(where(is.character), as.factor)) %>% 
  skimr::skim()

A_rumen_counts <- A_rumen_counts %>% 
  ps_mutate(
    Cap = ifelse(SampleBinary == "Capsicum", yes = 1, no = 0)
  )
sample_data(A_rumen_counts)

lm_models <- A_rumen_counts %>% 
  tax_fix() %>% 
  tax_prepend_ranks() %>% 
  # it makes sense to perform the compositional transformation BEFORE filtering
  tax_transform("compositional", rank = "Genus", keep_counts = TRUE) %>% 
  tax_filter(min_prevalence = 0.1, undetected = 0, use_counts = TRUE) %>% 
  tax_transform(
    trans = "log2", chain = TRUE, zero_replace = "halfmin"
  ) %>% 
  taxatree_models(
    type = lm, 
    ranks = NULL, # uses every rank available except the first
    variables = c("Cap")
  )

lm_models
lm_stats <- taxatree_models2stats(lm_models)
lm_stats
lm_stats %>% taxatree_stats_get()

lm_stats <- taxatree_stats_p_adjust(
  data = lm_stats, method = "BH", grouping = "rank"
)

lm_stats %>% taxatree_stats_get()


lm_stats %>% taxatree_plots(
  node_size_range = c(1,3), var_renamer = toupper
) %>% 
  patchwork::wrap_plots(
    ncol = 2, guides = "collect"
  )

set.seed(123) # label position 
key <- taxatree_plotkey(
  data = lm_stats, 
  taxon_renamer = function(x) stringr::str_remove(x, "[PFG]: "),
  # 2 lines of conditions below, for filtering taxa to be labelled
  rank == "Phylum" | rank == "Genus" & prevalence > 0.25, 
  !grepl("Kingdom", taxon)
) + 
  # add a bit more space for the longer labels by expanding the x axis
  scale_x_continuous(expand = expansion(mult = 0.2))
# all phyla are labelled, and all genera with a prevalence of over 0.2
# except for any taxa whose names (partly) match "Kingdom" 
# (i.e. an unclassified taxon)

key


lm_stats %>% 
  taxatree_label(
    rank == "Genus", p.value < 0.05 | prevalence > 0.5, estimate > 0
  ) %>% 
  taxatree_plots() %>% 
  .[[1]] %>% # show only the first plot
  taxatree_plot_labels(
    taxon_renamer = function(x) stringr::str_remove(x, "G: "),
    fun = ggrepel::geom_label_repel, x_nudge = 0.7, hjust = 0.5, size = 2
  ) 




