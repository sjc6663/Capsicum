## Differential Relative Abundance - Time (0v2, 0v18)
## SAB 03/20/2023

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
library(patchwork)

# load phyloseq objects
load("ps-obj/phyloseq-rumen-samples-angus-counts.RData")
load("ps-obj/phyloseq-rumen-samples-holstein-counts.RData")

# check if we have NAs
anyNA(tax_table(A_rumen_counts)[,"Phylum"])
anyNA(tax_table(H_rumen_counts)[,"Phylum"])

# tax fix our phyloseq object
H_rumen_counts <- tax_fix(H_rumen_counts)
A_rumen_counts <- tax_fix(A_rumen_counts)

A_rumen_counts <- A_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))
H_rumen_counts <- H_rumen_counts %>% tax_fix(unknowns = c("Incertae Sedis"))

## subset groups by time ----
# define the quick hour segment
quickHour <- c("H0", "H2")

# define the long hour segment
longHour <- c("H0", "H18")

A2 <- subset_samples(
  A_rumen_counts,
  Hour %in% quickHour
)

A18 <- subset_samples(
  A_rumen_counts,
  Hour %in% longHour
)

H2 <- subset_samples(
  H_rumen_counts,
  Hour %in% quickHour
)

H18 <- subset_samples(
  H_rumen_counts,
  Hour %in% longHour
)

## DESeq Holstein 0v2 ----
# add pseudo count of 1 
otu_table(H2) <- otu_table(H2) + 1

# phyloseq-DEseq
ddsh2 <- phyloseq_to_deseq2(H2, ~ Hour)

# perform DESeq
dsh2 <- DESeq(ddsh2, test = "Wald", fitType = "parametric")

# get results
resh2 = results(dsh2, cooksCutoff = FALSE)
# we are going to set this to be more stringent because we got A LOT of taxa at 0.05.
alpha = 0.01
sigtabh2 = resh2[which(resh2$padj < alpha), ]
sigtabh2 = cbind(as(sigtabh2, "data.frame"), as(tax_table(H2)[rownames(sigtabh2), ], "matrix"))
sigtabh2
write.csv(sigtabh2, file = "tables/rumen-H-0v2.csv")

## plot
# Phylum order
x = tapply(sigtabh2$log2FoldChange, sigtabh2$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabh2$Phylum = factor(as.character(sigtabh2$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabh2$log2FoldChange, sigtabh2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabh2$Genus = factor(as.character(sigtabh2$Genus), levels=names(x))
sigtabh2

# ph2 <- 
  ggdotchart(sigtabh2, y = "log2FoldChange", x = "Genus",
                 color = "Phylum",
                 sorting = "descending",
                 add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                 dot.size = 5,
                 # add log2 fold change in the dot
                 label = round(sigtabh2$log2FoldChange, 2),
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

ph2

ggsave("plots/H-0v2-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 8)

## DESeq Holstein 0v18 ----
# add pseudo count of 1 
otu_table(H18) <- otu_table(H18) + 1

# phyloseq-DEseq
ddsh18 <- phyloseq_to_deseq2(H18, ~ Hour)

# perform DESeq
dsh18 <- DESeq(ddsh18, test = "Wald", fitType = "parametric")

# get results
resh18 = results(dsh18, cooksCutoff = FALSE)
# we are going to set this to be more stringent because we got A LOT of taxa at 0.05.
alpha = 0.01
sigtabh18 = resh18[which(resh18$padj < alpha), ]
sigtabh18 = cbind(as(sigtabh18, "data.frame"), as(tax_table(H18)[rownames(sigtabh18), ], "matrix"))
write.csv(sigtabh18, file = "tables/rumen-H-0v18.csv")

## plot
# Phylum order
x = tapply(sigtabh18$log2FoldChange, sigtabh18$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabh18$Phylum = factor(as.character(sigtabh18$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabh18$log2FoldChange, sigtabh18$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabh18$Genus = factor(as.character(sigtabh18$Genus), levels=names(x))
sigtabh18

# pH18 <- 
  ggdotchart(sigtabh18, y = "log2FoldChange", x = "Genus",
                 color = "Phylum",
                 sorting = "descending",
                 add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                 dot.size = 5,
                 # add log2 fold change in the dot
                 label = round(sigtabh18$log2FoldChange, 2),
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

pH18

ggsave("plots/rumen-H-0v18-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 8)

# DESeq - Angus 0v2 ----
# add pseudo count of 1 
otu_table(A2) <- otu_table(A2) + 1

# phyloseq-DEseq
ddsa2 <- phyloseq_to_deseq2(A2, ~ Hour)

# perform DESeq
dsa2 <- DESeq(ddsa2, test = "Wald", fitType = "parametric")

# get results
resa2 = results(dsa2, cooksCutoff = FALSE)
# we are going to set this to be more stringent because we got A LOT of taxa at 0.05.
alpha = 0.01
sigtaba2 = resa2[which(resa2$padj < alpha), ]
sigtaba2 = cbind(as(sigtaba2, "data.frame"), as(tax_table(A2)[rownames(sigtaba2), ], "matrix"))
write.csv(sigtaba2, file = "tables/rumen-A-0v2.csv")

## plot
# Phylum order
x = tapply(sigtaba2$log2FoldChange, sigtaba2$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtaba2$Phylum = factor(as.character(sigtaba2$Phylum), levels=names(x))
# Genus order
x = tapply(sigtaba2$log2FoldChange, sigtaba2$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtaba2$Genus = factor(as.character(sigtaba2$Genus), levels=names(x))
sigtaba2

# pA2 <- 
  ggdotchart(sigtaba2, y = "log2FoldChange", x = "Genus",
                   color = "Phylum",
                   sorting = "descending",
                   add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                   dot.size = 5,
                   # add log2 fold change in the dot
                   label = round(sigtaba2$log2FoldChange, 2),
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

pA2

ggsave("plots/rumen-A-0v2-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 8)

# DESeq - Angus 0v18 ----
# add pseudo count of 1 
otu_table(A18) <- otu_table(A18) + 1

# phyloseq-DEseq
ddsa18 <- phyloseq_to_deseq2(A18, ~ Hour)

# perform DESeq
dsa18 <- DESeq(ddsa18, test = "Wald", fitType = "parametric")

# get results
resa18 = results(dsa18, cooksCutoff = FALSE)
# we are going to set this to be more stringent because we got A LOT of taxa at 0.05.
alpha = 0.01
sigtaba18 = resa18[which(resa18$padj < alpha), ]
sigtaba18 = cbind(as(sigtaba18, "data.frame"), as(tax_table(A18)[rownames(sigtaba18), ], "matrix"))
sigtaba18
write.csv(sigtaba18, file = "tables/rumen-A-0v18.csv")

## plot
# Phylum order
x = tapply(sigtaba18$log2FoldChange, sigtaba18$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtaba18$Phylum = factor(as.character(sigtaba18$Phylum), levels=names(x))
# Genus order
x = tapply(sigtaba18$log2FoldChange, sigtaba18$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtaba18$Genus = factor(as.character(sigtaba18$Genus), levels=names(x))
sigtaba18

# pA18 <- 
  ggdotchart(sigtaba18, y = "log2FoldChange", x = "Genus",
                   color = "Phylum",
                   sorting = "descending",
                   add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                   dot.size = 5,
                   # add log2 fold change in the dot
                   label = round(sigtaba18$log2FoldChange, 2),
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

pA18

ggsave("plots/rumen-A-0v18-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 8)

(ph2|pH18)/(pA2|pH18)

ggsave("plots/diff-rel-abund-time-all.pdf", dpi = 600, width = 12, height = 20)

ph2|pH18
ggsave("plots/diff-rel-abund-holstein-all.pdf", dpi = 600, width = 24, height = 8)

pA2|pH18
ggsave("plots/diff-rel-abund-angus-all.pdf", dpi = 600, width = 24, height = 8)

