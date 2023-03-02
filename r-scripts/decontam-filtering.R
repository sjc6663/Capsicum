## Decontamination and filtering
# Emily Van Syoc 11/21 - adapted by Sophia 6/22 - adapted by Stephanie 6/23/22

## ---- setup ----
# load packages
require(tidyverse)
require(phyloseq)
require(dada2)
# install latest version of microViz
#devtools::install_github("david-barnett/microViz@0.7.1")
require(microViz)
#require(BiocManager)
#BiocManager::install("decontam")
require(decontam)
require(ggpubr)

# get data - pull in raw phyloseq
load("data-objects-redo/phyloseq-raw.RData")

## ---- basic info ----

ps <- psraw

ntaxa(ps) # 35028 taxa

get_taxa_unique(ps, "Phylum") # 48 different phyla including NA
get_taxa_unique(ps, "Kingdom") # 4 kingdoms incling NA- need to filter this first

# how many controls at each step
ps %>% samdat_tbl() %>% group_by(Steer.ID) %>% summarize(n())

# ---- filter for Bacteria & remove NA phyla ----

# get only bacteria
psf <- subset_taxa(ps, Kingdom == "Bacteria")

# remove NA phylum
psf1 <- subset_taxa(psf, !is.na(Phylum) & !Phylum %in% c("", "NA"))

# validate
psf2 <- tax_fix(psf1)
psf3 <- phyloseq_validate(psf2, remove_undetected = TRUE)

## ---- filter for relative abundance ----

# transform to relative abundance
pst <- transform_sample_counts(psf3, function(x) x / sum(x) )

# remove taxa with total relative abundance less than 10e-5
psr <- filter_taxa(pst, function(x) mean(x) > 1e-5, TRUE) #2381 taxa

## ---- decontam ----

## workflow from: https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# add sampling variables
psd <- psr %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Steer.ID,"Control"), true = "Control", false = "Sample")
  ) 

#check

sample_data(psd)

## remove positive controls for this 
## update - keep and remove contaminants like a sample
#psnopos <- psd %>% ps_filter(!Treatment %in% c("PCP2", "PCP3", "PCP1", "PCISO1"))

## inspect library sizes 
sampdf <- data.frame(sample_data(psd))
sampdf$LibrarySize <- sample_sums(psd)
sampdf <- sampdf[order(sampdf$LibrarySize), ]
sampdf$Index <- seq(nrow(sampdf))
ggplot(data = sampdf, aes(x = Index, y = LibrarySize, color = SampleBinary)) + geom_point()

## controls are scattered throughout

## use negative controls to identify contaminants
# add distinguisher for positive controls and treat as a sample
psp <- psd %>% ps_mutate(SampleBinary = case_when(
  SampleBinary == "Control" & str_detect(Steer.ID, "Pos") ~ "Sample",
  SampleBinary == "Control" & !str_detect(Steer.ID, "Pos") ~ "Control",
  SampleBinary == "Sample" ~ "Sample"
))
pds <- psp %>% ps_mutate(is.neg = if_else(SampleBinary == "Control", true = TRUE, false = FALSE))
contamdf.prev <- isContaminant(pds, method="prevalence", neg="is.neg", threshold = 0.5) # more aggressive threshold is 0.5
table(contamdf.prev$contaminant) # this identifies 690 contaminants and 1691 'true' taxa

## identify the contaminant taxa and remove them from the phyloseq object
psdecon <- prune_taxa(!contamdf.prev$contaminant, pds)

## remove the positive controls and save the decontaminated phyloseq object
pssave <- ps_filter(psdecon, SampleBinary == "Sample") %>% 
  # remove positive controls
  ps_filter(!str_detect(Steer.ID, "Pos")) %>% 
  ps_select(-c(is.neg, SampleBinary))

factor(sample_data(pssave)$Steer.ID)

#this leaves 1691 taxa

# write
saveRDS(pssave, "data-objects-redo/ps-decontam-filtered.rds")

## ---- negative controls ----

# what do raw controls look like
neg <- ps %>% # ps_filter(Steer.ID != "ID" | Steer.ID != "1-CB") %>% 
  # get negatives
  ps_filter(str_detect(Steer.ID, "Neg")) %>% 
  tax_fix()

factor(sample_data(neg)$Steer.ID)

ntaxa(neg) # 2110 taxa

# barplot
neg %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 10,  
    merge_other = FALSE 
    # set merge_other = TRUE (the default) to remove outlines from inside "other" category
  ) +
  #facet_wrap("Treatment", scales = "free") +
  labs(x = NULL, y = NULL)

## filtered neg controls
negf <- subset_taxa(ps, taxa_names(ps) %in% taxa_names(psr)) %>% 
  #ps_filter(Treatment != "Control" | Treatment != "Probiotic") %>% 
  # get negatives
  ps_filter(str_detect(Steer.ID, "Neg")) %>% 
  tax_fix()

ntaxa(negf) # 1268
get_taxa_unique(negf, "Genus")

# barplot
sample_names(negf) <- c("NC1", "NC2", "NC3", "NC4")

p <- negf %>%
  ps_mutate(Name = c("NC1", "NC2", "NC3", "NC4")) %>% 
  comp_barplot(
    tax_level = "Phylum", #n_taxa = 5,
    #tax_order = "name",
    label = "Name",
    sample_order = c("NC1", "NC2", "NC3", "NC4"),
    merge_other = TRUE,
    bar_width = 0.9,
    # set merge_other = TRUE (the default) to remove outlines from inside "other" category
  ) +
  labs(y = "Relative Abundance") +
  theme(
    text = element_text(size = 16))


# hacky fix to make taxa names alphabetical and remove "other"
p$data <- filter(p$data, !unique == "other")
p$data$top %>% head()
p$data$top <- as.character(p$data$top)
p$data$top[p$data$top == "other"] <- NA
p <- p + guides(fill = guide_legend(title = "Phylum", reverse = FALSE))
p 

p$data$unique <- as.character(p$data$unique)
p$data$unique[p$data$unique == "other"] <- NA
p


# save
ggsave(filename = "neg-control-barplot.png", plot = p, dpi = 600)

# save
saveRDS(negf, file = "data-objects-redo/ps-negcontrols-filtered.rds")


## ----- positive control -----

# "raw" positive controls
pos <- ps %>%
  ps_filter(str_detect(Steer.ID, "Pos")) %>% 
  tax_fix()

pos %>% 
  comp_barplot(
    tax_level = "Phylum", n_taxa = 8
  ) +
  facet_wrap(~Treatment, scales = "free")

## get positive controls after removing contaminants and filtering for relative abundance
pspos <- psdecon %>% 
  ps_filter(str_detect(Steer.ID, "Pos")) 
# remove PCR control - added incorrectly
# ps_filter(!Treatment == "PCR_PC")

# there are 1117 unique taxa
ntaxa(pspos)

pspos <- tax_fix(pspos)
# visualize in relative abundance barplot
pspos %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 3,  
    merge_other = TRUE,
    #label = 'Treatment',
    bar_width = 0.9
    # set merge_other = TRUE (the default) to remove outlines from inside "other" category
  )

# barplot
sample_names(pspos) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")
p <- pspos %>%
  ps_mutate(Name = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")) %>% 
  comp_barplot(
    tax_level = "Phylum", n_taxa = 3,
    label = "Name",
    sample_order = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7"),
    merge_other = TRUE,
    bar_width = 0.9,
    # set merge_other = TRUE (the default) to remove outlines from inside "other" category
  ) 

# save
ggsave(filename = "pos-control-barplot.png", plot = last_plot(), dpi = 600)

## ---- build positive Zymo control ----

## control used in DNA extraction: ZymoBiomics Community Standard Lot ZRC190633, Cat D6300

# strain info for 16S here: https://files.zymoresearch.com/pdf/d6300-_zymobiomics_microbial_community_standard_v1-1-3.pdf

# build taxonomy - only need Phylum & Genus
tax <- data.frame(
  #row.names = c("Kingdom", "Genus"),
  p.aeru = c("Proteobacteria", "Pseudomonas"),
  e.coli = c("Proteobacteria", "Escherichia"),
  s.ent = c("Proteobacteria", "Salmonella"),
  l.ferm = c("Firmicutes", "Lactobacillus"),
  e.fae = c("Firmicutes", "Enterococcus"),
  s.aur = c("Firmicutes", "Staphylococcus"),
  l.mon = c("Firmicutes", "Listeria"),
  b.sub = c("Firmicutes", "Bacillus")
  
)  %>% t()

colnames(tax) <- c("Phylum", "Genus")


zymo <- phyloseq(otu_table(
  data.frame(
    row.names = "ZymoComm",
    p.aeru = 0.042,
    e.coli = 0.101,
    s.ent = 0.104,
    l.ferm = 0.184,
    e.fae = 0.099,
    s.aur = 0.155,
    l.mon = 0.141,
    b.sub = 0.174
  ), taxa_are_rows = FALSE
),
tax_table(
  tax
  
),
# have to do this for microViz...
sample_data(data.frame(
  row.names = "ZymoComm",
  fakedata = 1
))
)

# plot
z<- comp_barplot(zymo, "Phylum",
                 bar_width = 0.5)

# plot with both
ggarrange(p, z, common.legend = TRUE, legend = "right")



ggsave(filename = "pos-with-zymo.png", dpi = 600)

## ---- find taxa in pos not in zymo ----
get_taxa_unique(pspos, "Phylum")
get_taxa_unique(pspos, "Genus")
get_taxa_unique(zymo, "Phylum")
get_taxa_unique(zymo, "Genus")

pstruepos <- subset_taxa(pspos, Phylum == "Proteobacteria" | Phylum == "Firmicutes") 
pstruepos <- subset_taxa(pspos, Genus == "Pseudomonas"| Genus == "Escherichia-Shigella" | Genus == "Salmonella" | Genus == "Lactobacillus" | Genus == "Enterococcus" | Genus == "Staphylococcus" | Genus == "Listeria" | Genus == "Bacillus")

#create ps of those dropped from pos control
psposrm <- subset_taxa(pspos, !taxa_names(pspos) %in% taxa_names(pstruepos)) #these are contams and need filtered out

get_taxa_unique(pstruepos, "Phylum")
get_taxa_unique(pstruepos, "Genus")

get_taxa_unique(psposrm, "Phylum")
get_taxa_unique(psposrm, "Genus")

psposrm #1112 taxa are pos control contaminants 

#plot again

p2 <- pstruepos %>%
  ps_mutate(Name = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7")) %>% 
  comp_barplot(
    tax_level = "Phylum", n_taxa = 3,
    label = "Name",
    sample_order = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7"),
    merge_other = TRUE,
    bar_width = 0.9,
    # set merge_other = TRUE (the default) to remove outlines from inside "other" category
  ) 

ggarrange(p2, z, common.legend = TRUE, legend = "right")

#save 
ggsave(filename = "posfilt-with-zymo.png", dpi = 600)


## ---- get phyloseq of counts ----

## we need the count data with the taxa present in the decontaminated and filtered PS
pscount <- subset_taxa(ps, taxa_names(ps) %in% taxa_names(pssave)) %>% 
  # remove controls
  ps_filter(!str_detect(Steer.ID, "Control")) %>% 
  #ps_filter(Sample_Type != "air") 
  ps_filter(!str_detect(Steer.ID, "1-CB"))

## ---- filter out positive control contaminants from both counts and rel abund phyloseqs ----
pscount2 <- subset_taxa(pscount, !taxa_names(pscount) %in% taxa_names(psposrm)) # down to 579 taxa
psrelabund  <- subset_taxa(pssave, !taxa_names(pssave) %in% taxa_names(psposrm)) # down to 579 taxa

psrelabund2 <- ps_filter(psrelabund, !str_detect(Steer.ID, "1-CB"))

factor(sample_data(pscount2)$Steer.ID)
factor(sample_data(psrelabund2)$Steer.ID)
sample_data(pscount2)

# save
saveRDS(pscount2, file = "ps-decontam-filtered-counts.rds") # use this file for alpha diversity
saveRDS(psrelabund2, file = "ps-decontam-relabund.rds") # use this file for beta diversity and clr transform

## ---- Relative abundance barplot ----

## all samples
pscount2 %>% 
  tax_fix() %>%
  comp_barplot(
    tax_level = "Phylum", n_taxa = 10,
    merge_other = TRUE,
  )+
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  facet_wrap(~ Sample.Type)

save.image(file = "data-objects-redo/working_decontam.RData")


pssave


