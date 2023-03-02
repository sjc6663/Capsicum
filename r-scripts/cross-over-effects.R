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

load("ps-obj/phyloseq-fecal-samples-only-counts.RData")
load("ps-obj/phyloseq-rumen-samples-only-counts.RData")

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
summary(mod)

hist(bothr$Shannon)
modR <- lm(Shannon ~ Treatment + x1 + x2, data = bothr)
summary(modR)

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

OTUf <- as(sample_data(transform_fecal), "matrix")
OTUdf <- as.data.frame(OTUf)
head(OTUdf)
OTUdf1 <- rownames_to_column(OTUdf, var = ".sample_name")
head(OTUdf1)
head(fecal_dist_matrix)

# generate distance matrix
fecal_dist_matrix <- phyloseq::distance(transform_fecal, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dm <- vegan::adonis2(fecal_dist_matrix ~ Treatment.1 + x1 + x2, data = OTUdf1)
mod_dm

# RUMEN
sample_data(rumen_counts)
transform_rumen <- microbiome::transform(rumen_counts, 'clr')
sample_data(transform_rumen)

OTUr <- as(sample_data(transform_rumen), "matrix")
OTUdr <- as.data.frame(OTUr)
head(OTUdr)
OTUdr1 <- rownames_to_column(OTUdr, var = ".sample_name")
head(OTUdr1)

# generate distance matrix
rumen_dist_matrix <- phyloseq::distance(transform_rumen, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dmr <- vegan::adonis2(rumen_dist_matrix ~ Treatment.1 + x1 + x2, data = OTUdr1)
mod_dmr

# there is a signficance in rumen samples for treatment order



## ---- setup for Holstein Counts ----
# load in Holstein Fecal and Rumen Counts
load("ps-obj/phyloseq-fecal-samples-holstein-counts.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")

# load in effect coding table
dat <- readxl::read_xlsx("metadata/Extraction_IDs.xlsx",
                         sheet = "effect-coding")

head(dat)

# renamed Steer ID to match, rename Treatment so it DOESN'T match, change Period to match
dat1 <- dat %>% rename(Steer.ID = Steer, Treatment.1 = Treatment) %>% 
  mutate(Period = paste0("P", Period))

# join with old sample data - FECAL
Hf_newsampledat <- samdat_tbl(H_fecal_counts) %>% 
  dplyr::select(Steer.ID, .sample_name, Breed, Period) %>% 
  left_join(dat1) %>% 
  column_to_rownames(var = ".sample_name")

# create a new sample data
sample_data(H_fecal_counts) = Hf_newsampledat

# add new sample data to phyloseq object
head(sample_data(H_fecal_counts))

# join with old sample data - RUMEN
Hr_newsampledat <- samdat_tbl(H_rumen_counts) %>% 
  dplyr::select(Steer.ID, .sample_name, Breed, Period) %>% 
  left_join(dat1) %>% 
  column_to_rownames(var = ".sample_name")

# create a new sample data
sample_data(H_rumen_counts) = Hr_newsampledat

head(sample_data(H_rumen_counts))


## ---- cross over effects = ALPHA DIVERSITY + HOLSTEIN ----
# create data frame 
adivHf <- data.frame(
  "Observed" = phyloseq::estimate_richness(H_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_fecal_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_fecal_counts)$Treatment)
head(adivHf)

# create data frame 
adivHr <- data.frame(
  "Observed" = phyloseq::estimate_richness(H_rumen_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(H_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(H_rumen_counts)$Treatment)
head(adivHr)

# change rownames to column 
dHf <- rownames_to_column(adivHf, var = ".sample_name")
dHr <- rownames_to_column(adivHr, var = ".sample_name")

# merge data frames together to get our variables with our Shannon Diveristy values
bothHf <- dHf %>% merge(Hf_newsampledat)
bothHr <- dHr %>% merge(Hr_newsampledat)

# check for crossover effects 
modHf <- lm(Shannon ~ Treatment + x1 + x2, data = bothHf)
summary(modHf) # no effects

modHr <- lm(Shannon ~ Treatment + x1 + x2, data = bothHr)
summary(modHr) # no effects

## ---- setup for Angus Counts ----
# load in Angus Fecal and Rumen Counts
load("ps-obj/phyloseq-fecal-samples-angus-counts.RData")
load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")

# load in effect coding table
dat <- readxl::read_xlsx("metadata/Extraction_IDs.xlsx",
                         sheet = "effect-coding")

head(dat)

# renamed Steer ID to match, rename Treatment so it DOESN'T match, change Period to match
dat1 <- dat %>% rename(Steer.ID = Steer, Treatment.1 = Treatment) %>% 
  mutate(Period = paste0("P", Period))

# join with old sample data - FECAL
Af_newsampledat <- samdat_tbl(A_fecal_counts) %>% 
  dplyr::select(Steer.ID, .sample_name, Breed, Period) %>% 
  left_join(dat1) %>% 
  column_to_rownames(var = ".sample_name")

# create a new sample data
sample_data(A_fecal_counts) = Af_newsampledat

# add new sample data to phyloseq object
head(sample_data(A_fecal_counts))

# join with old sample data - RUMEN
Ar_newsampledat <- samdat_tbl(A_rumen_counts) %>% 
  dplyr::select(Steer.ID, .sample_name, Breed, Period) %>% 
  left_join(dat1) %>% 
  column_to_rownames(var = ".sample_name")

# create a new sample data
sample_data(A_rumen_counts) = Ar_newsampledat

head(sample_data(A_rumen_counts))
## ---- cross over effects = ALPHA DIVERSITY + ANGUS ----
# create data frame 
adivAf <- data.frame(
  "Observed" = phyloseq::estimate_richness(A_fecal_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_fecal_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(A_fecal_counts)$Treatment)
head(adivAf)

# create data frame 
adivAr <- data.frame(
  "Observed" = phyloseq::estimate_richness(A_rumen_counts, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(A_rumen_counts, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(A_rumen_counts)$Treatment)
head(adivAr)

# change rownames to column 
dAf <- rownames_to_column(adivAf, var = ".sample_name")
dAr <- rownames_to_column(adivAr, var = ".sample_name")

# merge data frames together to get our variables with our Shannon Diveristy values
bothAf <- dAf %>% merge(Af_newsampledat)
bothAr <- dAr %>% merge(Ar_newsampledat)

# check for crossover effects 
modAf <- lm(Shannon ~ Treatment + x1 + x2, data = bothAf)
summary(modAf) # no effects

modAr <- lm(Shannon ~ Treatment + x1 + x2, data = bothAr)
summary(modAr) # no effects




## ---- setup for Holstein Rel Abund ----
# use phyloseq objects from Alpha Diversity because they already have the x1 and x2 in them
# clr transform data
H_transform_fecal <- microbiome::transform(H_fecal_counts, 'clr')
sample_data(H_transform_fecal)

H_transform_rumen <- microbiome::transform(H_rumen_counts, 'clr')
sample_data(H_transform_rumen)

## ---- cross over effects = BETA DIVERSITY + HOLSTEIN ----

# FECAL
# pull out OTU table as data frame
OTUHf <- as(sample_data(H_transform_fecal), "matrix")
OTUdHf <- as.data.frame(OTUHf)
head(OTUdHf)
OTUdHf1 <- rownames_to_column(OTUdHf, var = ".sample_name")
head(OTUdHf1)

# generate distance matrix
H_fecal_dist_matrix <- phyloseq::distance(H_transform_fecal, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dmHf <- vegan::adonis2(H_fecal_dist_matrix ~ Treatment.1 + x1 + x2, data = OTUdHf1)
mod_dmHf # no effect

# RUMEN
# pull out OTU table as data frame
OTUHr <- as(sample_data(H_transform_rumen), "matrix")
OTUdHr <- as.data.frame(OTUHr)
head(OTUdHr)
OTUdHr1 <- rownames_to_column(OTUdHr, var = ".sample_name")
head(OTUdHr1)

# generate distance matrix
H_rumen_dist_matrix <- phyloseq::distance(H_transform_rumen, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dmHr <- vegan::adonis2(H_rumen_dist_matrix ~ Treatment.1 + x1 + x2, data = OTUdHr1)
mod_dmHr # x1 is significant so we are seeing a significant effect on beta diversity in rumen for the order

## ---- setup for Angus Rel Abund ----
# use phyloseq objects from Alpha Diversity because they already have the x1 and x2 in them
# clr transform data
A_transform_fecal <- microbiome::transform(A_fecal_counts, 'clr')
sample_data(A_transform_fecal)

A_transform_rumen <- microbiome::transform(A_rumen_counts, 'clr')
sample_data(A_transform_rumen)

## ---- cross over effects = BETA DIVERSITY + ANGUS ----

# FECAL
# pull out OTU table as data frame
OTUAf <- as(sample_data(A_transform_fecal), "matrix")
OTUdAf <- as.data.frame(OTUAf)
head(OTUdAf)
OTUdAf1 <- rownames_to_column(OTUdAf, var = ".sample_name")
head(OTUdAf1)

# generate distance matrix
A_fecal_dist_matrix <- phyloseq::distance(A_transform_fecal, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dmAf <- vegan::adonis2(A_fecal_dist_matrix ~ Treatment.1 + x1 + x2, data = OTUdAf1)
mod_dmAf # no effect

# RUMEN
# pull out OTU table as data frame
OTUAr <- as(sample_data(A_transform_rumen), "matrix")
OTUdAr <- as.data.frame(OTUAr)
head(OTUdAr)
OTUdAr1 <- rownames_to_column(OTUdAr, var = ".sample_name")
head(OTUdAr1)

# generate distance matrix
A_rumen_dist_matrix <- phyloseq::distance(A_transform_rumen, method = "euclidean")

# check the significance with treatment and crossover identifiers
mod_dmAr <- vegan::adonis2(A_rumen_dist_matrix ~ Treatment.1 + x1 + x2, data = OTUdAr1)
mod_dmAr # x1 and x2 are significant so we are seeing a significant effect on beta diversity in rumen for the order

