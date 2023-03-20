# Differential Relative Abundance - Capsicum vs. No Capsicum (CNTRL)
# SAB 03-20-2023

# setup ----
# load libraries and data
library(corncob)
library(phyloseq)
library(magrittr)
library(microViz)
library(BiocManager)
# BiocManager::install("edgeR")
library(edgeR)
# BiocManager::install("DESeq2")
library(DESeq2)
library(knitr)
library(dplyr)
library(stringr)

# load phyloseq objects
load("ps-obj/phyloseq-fecal-samples-angus-counts.RData")
load("ps-obj/phyloseq-fecal-samples-holstein-counts.RData")
load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")

# check if we have NAs
anyNA(tax_table(A_fecal_counts)[,"Phylum"])
anyNA(tax_table(H_fecal_counts)[,"Phylum"])
anyNA(tax_table(A_rumen_counts)[,"Phylum"])
anyNA(tax_table(H_rumen_counts)[,"Phylum"])

# tax fix our phyloseq object
A_fecal_counts <- tax_fix(A_fecal_counts)
H_fecal_counts <- tax_fix(H_fecal_counts)
H_rumen_counts <- tax_fix(H_rumen_counts)
A_rumen_counts <- tax_fix(A_rumen_counts)

A_fecal_counts <- A_fecal_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
H_fecal_counts <- H_fecal_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
A_rumen_counts <- A_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
H_rumen_counts <- H_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))

# identify samples by Cap or No Cap ----
# give samples a column to descripe capsicum or no capsicum
A_rumen_counts <- A_rumen_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

A_fecal_counts <- A_fecal_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

H_rumen_counts <- H_rumen_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 

H_fecal_counts <- H_fecal_counts %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(Treatment,"Control"), true = "Control", false = "Capsicum")
  ) 


## DESeq - Holstein Rumen ----
# add pseudo count of 1 
otu_table(H_rumen_counts) <- otu_table(H_rumen_counts) + 1

# phyloseq-DEseq
ddsrH <- phyloseq_to_deseq2(H_rumen_counts, ~ SampleBinary)

# perform DESeq
dsrH <- DESeq(ddsrH, test = "Wald", fitType = "parametric")

# get results
resrH = results(dsrH, cooksCutoff = FALSE)
alpha = 0.05
sigtabrH = resrH[which(resrH$padj < alpha), ]
sigtabrH = cbind(as(sigtabrH, "data.frame"), as(tax_table(H_rumen_counts)[rownames(sigtabrH), ], "matrix"))
write.csv(sigtabrH, file = "tables/rumen-H.csv")

## plot
# Phylum order
x = tapply(sigtabrH$log2FoldChange, sigtabrH$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabrH$Phylum = factor(as.character(sigtabrH$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabrH$log2FoldChange, sigtabrH$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabrH$Genus = factor(as.character(sigtabrH$Genus), levels=names(x))
sigtabrH

RH <- ggdotchart(sigtabrH, y = "log2FoldChange", x = "Genus",
                  color = "Phylum", shape = "Phylum",
                  sorting = "descending",
                  add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                  dot.size = 5,
                  # add log2 fold change in the dot
                  label = round(sigtabrH$log2FoldChange, 2),
                  font.label = list(color = "black", size = 11, hjust = 1.25),
                  rotate = TRUE,
                  ylab = "Log2 Fold Change",
                  xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("A")

ggsave("plots/rumen-holstein-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 8)

## DESeq - Holstein Fecal ----
# add pseudo count of 1 
otu_table(H_fecal_counts) <- otu_table(H_fecal_counts) + 1

# phyloseq-DEseq
ddsHf <- phyloseq_to_deseq2(H_fecal_counts, ~ SampleBinary)

# perform DESeq
dsHf <- DESeq(ddsHf, test = "Wald", fitType = "parametric")

# get results
resHf = results(dsHf, cooksCutoff = FALSE)
alpha = 0.05
sigtabHf = resHf[which(resHf$padj < alpha), ]
sigtabHf = cbind(as(sigtabHf, "data.frame"), as(tax_table(H_fecal_counts)[rownames(sigtabHf), ], "matrix"))
sigtabHf
write.csv(sigtabHf, file = "tables/fecal-H.csv")

## plot
# Phylum order
x = tapply(sigtabHf$log2FoldChange, sigtabHf$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabHf$Phylum = factor(as.character(sigtabHf$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabHf$log2FoldChange, sigtabHf$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabHf$Genus = factor(as.character(sigtabHf$Genus), levels=names(x))
sigtabHf

FH <- ggdotchart(sigtabHf, y = "log2FoldChange", x = "Genus",
           color = "Phylum", shape = "Phylum",
           sorting = "descending",
           add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
           dot.size = 5,
           # add log2 fold change in the dot
           label = round(sigtabHf$log2FoldChange, 2),
           font.label = list(color = "black", size = 11, hjust = 1.25),
           rotate = TRUE,
           ylab = "Log2 Fold Change",
           xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("C")

ggsave("plots/fecal-holstein-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 8)

# DESeq - RUMEN ANGUS ----
# add pseudo count of 1 
otu_table(A_rumen_counts) <- otu_table(A_rumen_counts) + 1

# phyloseq-DEseq
ddsrA <- phyloseq_to_deseq2(A_rumen_counts, ~ SampleBinary)

# perform DESeq
dsrA <- DESeq(ddsrA, test = "Wald", fitType = "parametric")

# get results
resrA = results(dsrA, cooksCutoff = FALSE)
alpha = 0.05
sigtabrA = resrA[which(resrA$padj < alpha), ]
sigtabrA = cbind(as(sigtabrA, "data.frame"), as(tax_table(A_rumen_counts)[rownames(sigtabrA), ], "matrix"))
sigtabrA
write.csv(sigtabrA, file = "tables/rumen-A.csv")

## plot
# Phylum order
x = tapply(sigtabrA$log2FoldChange, sigtabrA$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabrA$Phylum = factor(as.character(sigtabrA$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabrA$log2FoldChange, sigtabrA$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabrA$Genus = factor(as.character(sigtabrA$Genus), levels=names(x))
sigtabrA

RA <- ggdotchart(sigtabrA, y = "log2FoldChange", x = "Genus",
           color = "Phylum", shape = "Phylum",
           sorting = "descending",
           add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
           dot.size = 5,
           # add log2 fold change in the dot
           label = round(sigtabrA$log2FoldChange, 2),
           font.label = list(color = "black", size = 11, hjust = 1.25),
           rotate = TRUE,
           ylab = "Log2 Fold Change",
           xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("B")

ggsave("plots/rumen-angus-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 8)

# DEseq - FECAL ANGUS ----
# add pseudo count of 1 
otu_table(A_fecal_counts) <- otu_table(A_fecal_counts) + 1

# phyloseq-DEseq
ddsaf <- phyloseq_to_deseq2(A_fecal_counts, ~ SampleBinary)

# perform DESeq
dsaf <- DESeq(ddsaf, test = "Wald", fitType = "parametric")

# get results
resaf = results(dsaf, cooksCutoff = FALSE)
alpha = 0.05
sigtabaf = resrH[which(resaf$padj < alpha), ]
sigtabaf = cbind(as(sigtabaf, "data.frame"), as(tax_table(A_fecal_counts)[rownames(sigtabaf), ], "matrix"))
write.csv(sigtabaf, file = "tables/fecal-A.csv")

## plot
# Phylum order
x = tapply(sigtabaf$log2FoldChange, sigtabaf$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabaf$Phylum = factor(as.character(sigtabaf$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabaf$log2FoldChange, sigtabaf$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabaf$Genus = factor(as.character(sigtabaf$Genus), levels=names(x))
sigtabaf

AF <- ggdotchart(sigtabaf, y = "log2FoldChange", x = "Genus",
           color = "Phylum", shape = "Phylum",
           sorting = "descending",
           add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
           dot.size = 5,
           # add log2 fold change in the dot
           label = round(sigtabaf$log2FoldChange, 2),
           font.label = list(color = "black", size = 11, hjust = 1.25),
           rotate = TRUE,
           ylab = "Log2 Fold Change",
           xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("D")

ggsave("plots/rumen-holstein-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 8)

library(patchwork)

(RH|RA)/(FH|AF)

ggsave("plots/diff-rel-abund-all.pdf", dpi = 600, width = 24, height = 20)
