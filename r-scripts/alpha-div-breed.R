# Shannon's Diversity Index for Breeds
# 11-14-22 SAB

## ---- load phyloseq objects ----
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

# glom to Genus level
psgf <- fecal_counts %>% tax_glom("Genus")
psgr <- rumen_counts %>% tax_glom("Genus")

#Generate a data.frame with adiv measures (richness types and what groups we want to compare)
adivf <- data.frame(
  "Observed" = phyloseq::estimate_richness(psgf, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(psgf, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(psgf)$Treatment, 
  "Breed" = phyloseq::sample_data(psgf)$Breed)
head(adivf)
adivr <- data.frame(
  "Observed" = phyloseq::estimate_richness(psgr, measures = "Observed"),
  "Shannon" = phyloseq::estimate_richness(psgr, measures = "Shannon"),
  "Treatment" = phyloseq::sample_data(psgr)$Treatment,
  "Breed" = phyloseq::sample_data(psgr)$Breed)
head(adivr)

# fit a linear model to test for assumptions of parametric tests
modf <- lm(Observed ~ Treatment, data = adivf)
modr <- lm(Observed ~ Treatment, data = adivr)

# plot the redisuals
hist(resid(modf))
hist(resid(modr))

# plot a normality plot 
qqnorm(resid(modf))

# add a normal line
qqline(resid(modf))

# repeat for rumen
qqnorm(resid(modr))
qqline(resid(modr))

## ---- Shannon ANOVA with Holstein and Angus ----
testFB <- aov(Shannon ~ Breed, data = adivf)
summary(testFB)

testRB <- aov(Shannon ~ Breed, data = adivr)
summary(testRB)
