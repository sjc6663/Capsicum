# Checking for Cross Over Effects
# 8/19/2022 SJC

# one number will be shannon diversity, beta diversity, diff rel abund, etc
## ---- setup ----
# load libraries and data
library(phyloseq)
library(microViz)
library(tidyverse)
library(vegan)
library(dplyr)

load("data-objects-redo/phyloseq-fecal-counts.RData")
load("data-objects-redo/phyloseq-rumen-counts.RData")

# check if we have NAs
anyNA(tax_table(fecal_counts)[,"Phylum"])
anyNA(tax_table(rumen_counts)[,"Phylum"])

# tax fix our phyloseq object
fecal_counts <- tax_fix(fecal_counts)
rumen_counts <- tax_fix(rumen_counts)

fecal_counts <- fecal_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
rumen_counts <- rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))

# load in effect coding table
dat <- readxl::read_xlsx("metadata/Extraction_IDs.xlsx",
                         sheet = "effect-coding")

dat

head(H_fecal_counts@sam_data)

# renamed Steer ID to match, rename Treatment so it DOESN'T match, change Period to match
dat1 <- dat %>% rename(Steer.ID = Steer, Treatment.1 = Treatment) %>% 
  mutate(Period = paste0("P", Period))

# join with old sample data - FECAL
newsampledat <- samdat_tbl(fecal_counts) %>% 
  dplyr::select(Steer.ID, .sample_name, Breed, Period) %>% 
  left_join(dat1) %>% 
  column_to_rownames(var = ".sample_name")

# join with old sample data - RUMEN
newsampledatR <- samdat_tbl(rumen_counts) %>% 
  dplyr::select(Steer.ID, .sample_name, Breed, Period) %>% 
  left_join(dat1) %>% 
  column_to_rownames(var = ".sample_name")



# create a new sample data
sample_data(fecal_counts) = newsampledat
sample_data(rumen_counts) = newsampledatR

## ---- cross over effects = ALPHA DIVERSITY ----
# load alpha diversity adivf and adivr tables - in 04-alpha-diversity.R script
#Generate a data.frame with adiv measures (richness types and what groups we want to compare)
adivf <- data.frame(
  "Observed" = phyloseq::estimate_richness(fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(fecal_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(fecal_counts)$Treatment, 
  "Breed" = phyloseq::sample_data(fecal_counts)$Breed)
head(adivf)
adivr <- data.frame(
  "Observed" = phyloseq::estimate_richness(rumen_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(rumen_counts)$Treatment,
  "Breed" = phyloseq::sample_data(rumen_counts)$Breed)
head(adivr)
# this gives us shannon and observed alpha diversity for each sample

df <- rownames_to_column(adivf, var = ".sample_name")
both <- df %>% merge(newsampledat)
both

dfr <- rownames_to_column(adivr, var = ".sample_name")
bothr <- dfr %>% merge(newsampledatR)

hist(both$Shannon)
mod <- lm(Shannon ~ Treatment + x1 + x2, data = both)
summary(mod) # Px1 = 1, Px2 = 1

hist(bothr$Shannon)
modR <- lm(Shannon ~ Treatment + x1 + x2, data = bothr)
summary(modR) #Px1 = 1, Px2 = 1

# check for normality
hist(resid(mod))
qqnorm(resid(mod))
qqline(resid(mod))

hist(resid(modR))
qqnorm(resid(modR))
qqline(resid(modR))


# checked for crossover effect but there was none (in neither fecal nor rumen) so we took it out

## ---- cross over effects = BETA DIVERSITY ----

transform_fecal <- microbiome::transform(fecal_counts, 'clr')
sample_data(transform_fecal)

# generate distance matrix
fecal_dist_matrix <- phyloseq::distance(transform_fecal, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dm <- vegan::adonis2(fecal_dist_matrix ~ phyloseq::sample_data(transform_fecal)$Treatment.1 + phyloseq::sample_data(transform_fecal)$x1 + phyloseq::sample_data(transform_fecal)$x2)
mod_dm # Px1 = 0.81, Px2 = 0.93

# there is no significance in fecal samples for treatment order

# RUMEN
sample_data(rumen_counts)
transform_rumen <- microbiome::transform(rumen_counts, 'clr')
sample_data(transform_rumen)

# generate distance matrix
rumen_dist_matrix <- phyloseq::distance(transform_rumen, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dmr <- vegan::adonis2(rumen_dist_matrix ~ phyloseq::sample_data(transform_rumen)$Treatment.1 + phyloseq::sample_data(transform_rumen)$x1 + phyloseq::sample_data(transform_rumen)$x2)
mod_dmr # Px1 = 0.69, Px2 = 0.47

# there is no signficance in rumen samples for treatment order


## ---- cross over effects = ALPHA DIVERSITY + HOLSTEIN ----
H_fecal_counts <- subset_samples(fecal_counts,
                                 Breed == "Holstein")

H_rumen_counts <- subset_samples(rumen_counts,
                                 Breed == "Holstein")

# create data frame 
adivHf <- data.frame(
  "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_fecal_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_fecal_counts)$Treatment,
  "x1" = phyloseq::sample_data(H_fecal_counts)$x1,
  "x2" = phyloseq::sample_data(H_fecal_counts)$x2)
head(adivHf)

# create data frame 
adivHr <- data.frame(
  "Observed" = phyloseq::estimate_richness(H_rumen_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_rumen_counts)$Treatment,
  "x1" = phyloseq::sample_data(H_rumen_counts)$x1,
  "x2" = phyloseq::sample_data(H_rumen_counts)$x2)
head(adivHr)

# check for crossover effects 
modHf <- lm(Shannon ~ Treatment + x1 + x2, data = adivHf)
summary(modHf) # no effects, Px1 = 0.72, Px2 = 0.56

modHr <- lm(Shannon ~ Treatment + x1 + x2, data = adivHr)
summary(modHr) # no effects, Px1 = 0.56, Px2 = 0.013 (remember our threshold for significance is 0.0125 after adjusting for multiple comparisons) 

## ---- setup for Angus Counts ----
A_fecal_counts <- subset_samples(fecal_counts,
                                 Breed == "Angus")

A_rumen_counts <- subset_samples(rumen_counts,
                                 Breed == "Angus")

## ---- cross over effects = ALPHA DIVERSITY + ANGUS ----
# create data frame 
adivAf <- data.frame(
  "Observed" = phyloseq::estimate_richness(A_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_fecal_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(A_fecal_counts)$Treatment,
  "x1" = phyloseq::sample_data(A_fecal_counts)$x1,
  "x2" = phyloseq::sample_data(A_fecal_counts)$x2)
head(adivAf)

# create data frame 
adivAr <- data.frame(
  "Observed" = phyloseq::estimate_richness(A_rumen_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(A_rumen_counts)$Treatment,
  "x1" = phyloseq::sample_data(A_rumen_counts)$x1,
  "x2" = phyloseq::sample_data(A_rumen_counts)$x2)
head(adivAr)

# check for crossover effects 
modAf <- lm(Shannon ~ Treatment + x1 + x2, data = adivAf)
summary(modAf) # no effects, Px1 = 0.23, Px2 = 0.19

modAr <- lm(Shannon ~ Treatment + x1 + x2, data = adivAr)
summary(modAr) # no effects, Px1 = 0.13, Px2 = 0.12




## ---- setup for Holstein Rel Abund ----
# use phyloseq objects from Alpha Diversity because they already have the x1 and x2 in them
# clr transform data
H_trans_fecal <- H_fecal_counts %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

H_trans_rumen <- H_rumen_counts %>% 
  tax_transform(trans = "clr", rank = "Genus") %>% 
  ps_get()

## ---- cross over effects = BETA DIVERSITY + HOLSTEIN ----

# FECAL

# generate distance matrix
H_fecal_dist_matrix <- phyloseq::distance(H_trans_fecal, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dfH <- vegan::adonis2(H_fecal_dist_matrix ~ phyloseq::sample_data(H_trans_fecal)$Treatment.1 + phyloseq::sample_data(H_trans_fecal)$x1 + phyloseq::sample_data(H_trans_fecal)$x2)
mod_dfH # no effects, Px1 = 0.43, Px2 = 0.44


# RUMEN

# generate distance matrix
H_rumen_dist_matrix <- phyloseq::distance(H_trans_rumen, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dfH <- vegan::adonis2(H_rumen_dist_matrix ~ phyloseq::sample_data(H_trans_rumen)$Treatment.1 + phyloseq::sample_data(H_trans_rumen)$x1 + phyloseq::sample_data(H_trans_rumen)$x2)
mod_dfH # no effects, Px1 = 0.67, Px2 = 0.98


## ---- setup for Angus Rel Abund ----
# use phyloseq objects from Alpha Diversity because they already have the x1 and x2 in them
# clr transform data
A_transform_fecal <- microbiome::transform(A_fecal_counts, 'clr')

A_transform_rumen <- microbiome::transform(A_rumen_counts, 'clr')


## ---- cross over effects = BETA DIVERSITY + ANGUS ----

# FECAL

# generate distance matrix
A_fecal_dist_matrix <- phyloseq::distance(A_transform_fecal, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dfH <- vegan::adonis2(A_fecal_dist_matrix ~ phyloseq::sample_data(A_transform_fecal)$Treatment.1 + phyloseq::sample_data(A_transform_fecal)$x1 + phyloseq::sample_data(A_transform_fecal)$x2)
mod_dfH # no effects, Px1 = 0.54, Px2 = 0.09

# RUMEN

# generate distance matrix
A_rumen_dist_matrix <- phyloseq::distance(A_transform_rumen, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dfH <- vegan::adonis2(A_rumen_dist_matrix ~ phyloseq::sample_data(A_transform_rumen)$Treatment.1 + phyloseq::sample_data(A_transform_rumen)$x1 + phyloseq::sample_data(A_transform_rumen)$x2)
mod_dfH # no effects, Px1 = 0.43, Px2 = 0.52
