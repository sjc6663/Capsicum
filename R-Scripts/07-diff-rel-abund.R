# Differential Relative Abundace using DESeq
# 11/7/2022 SAB
# Tutorial: https://www.yanh.org/2021/01/01/microbiome-r/ 

## ---- setup ----
# load libraries and data
library(corncob)
library(phyloseq)
library(magrittr)
library(microViz)
library(edgeR)
library(DESeq2)
library(knitr)

load("alpha-diversity/phyloseq-fecal-samples-angus-counts.RData")
load("alpha-diversity/phyloseq-fecal-samples-holstein-counts.RData")
load("alpha-diversity/phyloseq-rumen-samples-angus-counts.RData")
load("alpha-diversity/phyloseq-rumen-samples-holstein-counts.RData")

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

## ---- subset samples by each treatment ----
# RUMEN - HOLSTEIN
H_rumen_counts5 <- subset_samples(
  H_rumen_counts, 
  Treatment == c("RPC5", "Control")
)

H_rumen_counts10 <-subset_samples(
  H_rumen_counts, 
  Treatment == c("RPC10", "Control")
)

H_rumen_counts15 <- subset_samples(
  H_rumen_counts, 
  Treatment == c("RPC15", "Control")
)

# RUMEN - ANGUS

A_rumen_counts5 <- subset_samples(
  A_rumen_counts, 
  Treatment == c("RPC5", "Control")
)

A_rumen_counts10 <- subset_samples(
  A_rumen_counts, 
  Treatment == c("RPC10", "Control")
)

A_rumen_counts15 <- subset_samples(
  A_rumen_counts, 
  Treatment == c("RPC15", "Control")
)

# FECAL - HOLSTEIN
H_fecal_counts5 <- subset_samples(
  H_fecal_counts, 
  Treatment == c("RPC5", "Control")
)

H_fecal_counts10 <- subset_samples(
  H_fecal_counts, 
  Treatment == c("RPC10", "Control")
)

H_fecal_counts15 <- subset_samples(
  H_fecal_counts, 
  Treatment == c("RPC15", "Control")
)

# FECAL - ANGUS
A_fecal_counts5 <- subset_samples(
  A_fecal_counts, 
  Treatment == c("RPC5", "Control")
)

A_fecal_counts10 <- subset_samples(
  A_fecal_counts, 
  Treatment == c("RPC10", "Control")
)

A_fecal_counts15 <- subset_samples(
  A_fecal_counts, 
  Treatment == c("RPC15", "Control")
)

## ---- RUMEN, Holstein RPC5 ----
# add pseudo count of 1 
otu_table(H_rumen_counts5) <- otu_table(H_rumen_counts5) + 1

# phyloseq-DEseq
ddsrH5 <- phyloseq_to_deseq2(H_rumen_counts5, ~ Treatment)

# perform DESeq
dsrH5 <- DESeq(ddsrH5, test = "Wald", fitType = "parametric")

# get results
resrH5 = results(dsrH5, cooksCutoff = FALSE)
alpha = 0.05
sigtabrH5 = resrH5[which(resrH5$padj < alpha), ]
sigtabrH5 = cbind(as(sigtabrH5, "data.frame"), as(tax_table(H_rumen_counts5)[rownames(sigtabrH5), ], "matrix"))
write.csv(sigtabrH5, file = "differential-rel-abund/Rumen/H5.csv")

## plot
# Phylum order
x = tapply(sigtabrH5$log2FoldChange, sigtabrH5$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabrH5$Phylum = factor(as.character(sigtabrH5$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabrH5$log2FoldChange, sigtabrH5$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabrH5$Genus = factor(as.character(sigtabrH5$Genus), levels=names(x))
sigtabrH5
ggplot(sigtabrH5, aes(x=log2FoldChange, y=Genus, color=Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + 
  scale_color_brewer(palette = "BrBG") +
  theme_classic()

Hr5_A <- ggplot(data = sigtabrH5, 
                aes(x = log2FoldChange, y = Genus, color = Phylum, shape = Phylum)) + geom_point(size = 6) + geom_vline(xintercept = 0) + 
  #  geom_bar(stat = "identity", position=position_dodge(), width = 0.1, fill = "light grey", colour = "light grey", linetype = 1) +
  theme_classic() +
  ggtitle("A")

ggsave("differential-rel-abund/Rumen/holstein-rpc5-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

rH5 <- ggdotchart(sigtabrH5, y = "log2FoldChange", x = "Genus",
                  color = "Phylum", shape = "Phylum",
                  sorting = "descending",
                  add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                  dot.size = 5,
                  # add log2 fold change in the dot
                  label = round(sigtabrH5$log2FoldChange, 2),
                  font.label = list(color = "black", size = 11, vjust = -1.5),
                  rotate = TRUE,
                  ylab = "Log2 Fold Change",
                  xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("")

rH5

ggsave("differential-rel-abund/Rumen/holstein-rpc5-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)


# use patchwork to plot multiple graphs onto one graph (if that makes sense)
## ---- RUMEN, Holstein RPC10 ----
# add pseudo count of 1 
otu_table(H_rumen_counts10) <- otu_table(H_rumen_counts10) + 1

# phyloseq-DEseq
ddsrH10 <- phyloseq_to_deseq2(H_rumen_counts10, ~ Treatment)

# perform DESeq
dsrH10 <- DESeq(ddsrH10, test = "Wald", fitType = "parametric")

# get results
resrH10 = results(dsrH10, cooksCutoff = FALSE)
alpha = 0.05
sigtabrH10 = resrH10[which(resrH10$padj < alpha), ]
sigtabrH10 = cbind(as(sigtabrH10, "data.frame"), as(tax_table(H_rumen_counts10)[rownames(sigtabrH10), ], "matrix"))
# THERE ARE NONE THAT ARE SIGNFICANT SO WE CAN'T PLOT IT
sigtabrH10


## ---- RUMEN, Holstein RPC15 ----
# add pseudo count of 1 
otu_table(H_rumen_counts15) <- otu_table(H_rumen_counts15) + 1

# phyloseq-DEseq
ddsrH15 <- phyloseq_to_deseq2(H_rumen_counts15, ~ Treatment)

# perform DESeq
dsrH15 <- DESeq(ddsrH15, test = "Wald", fitType = "parametric")

# get results
resrH15 = results(dsrH15, cooksCutoff = FALSE)
alpha = 0.05
sigtabrH15 = resrH15[which(resrH15$padj < alpha), ]
sigtabrH15
sigtabrH15 = cbind(as(sigtabrH15, "data.frame"), as(tax_table(H_rumen_counts15)[rownames(sigtabrH15), ], "matrix"))
# THERE ARE NONE THAT ARE SIGNFICANT SO WE CANNOT PLOT IT 


## ---- RUMEN, Angus RPC5 ----
# add pseudo count of 1 
otu_table(A_rumen_counts5) <- otu_table(A_rumen_counts5) + 1

# phyloseq-DEseq
ddsrA5 <- phyloseq_to_deseq2(A_rumen_counts5, ~ Treatment)

# perform DESeq
dsrA5 <- DESeq(ddsrA5, test = "Wald", fitType = "parametric")

# get results
resrA5 = results(dsrA5, cooksCutoff = FALSE)
alpha = 0.05
sigtabrA5 = resrA5[which(resrA5$padj < alpha), ]
sigtabrA5
sigtabrA5 = cbind(as(sigtabrA5, "data.frame"), as(tax_table(A_rumen_counts5)[rownames(sigtabrA5), ], "matrix"))
write.csv(sigtabrA5, file = "differential-rel-abund/Rumen/A5.csv")

## plot
# Phylum order
x = tapply(sigtabrA5$log2FoldChange, sigtabrA5$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabrA5$Phylum = factor(as.character(sigtabrA5$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabrA5$log2FoldChange, sigtabrA5$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabrA5$Genus = factor(as.character(sigtabrA5$Genus), levels=names(x))
rA5 <- ggplot(sigtabrA5, aes(x=log2FoldChange, y=Genus, color=Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_brewer(palette = "BrBG") +
  ggtitle("A") +
  theme_classic()

rA5

ggsave("differential-rel-abund/Rumen/angus-rpc5-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

rA5 <- ggdotchart(sigtabrA5, y = "log2FoldChange", x = "Genus",
                  color = "Phylum", shape = "Phylum",
                  sorting = "descending",
                  add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                  dot.size = 5,
                  # add log2 fold change in the dot
                  label = round(sigtabrA5$log2FoldChange, 2),
                  font.label = list(color = "black", size = 11, vjust = -1.5),
                  rotate = TRUE,
                  ylab = "Log2 Fold Change",
                  xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("A")

## ---- RUMEN, Angus RPC10 ----
# add pseudo count of 1 
otu_table(A_rumen_counts10) <- otu_table(A_rumen_counts10) + 1

# phyloseq-DEseq
ddsrA10 <- phyloseq_to_deseq2(A_rumen_counts10, ~ Treatment)

# perform DESeq
dsrA10 <- DESeq(ddsrA10, test = "Wald", fitType = "parametric")

# get results
resrA10 = results(dsrA10, cooksCutoff = FALSE)
alpha = 0.05
sigtabrA10 = resrA10[which(resrA10$padj < alpha), ]
sigtabrA10
sigtabrA10 = cbind(as(sigtabrA10, "data.frame"), as(tax_table(A_rumen_counts10)[rownames(sigtabrA10), ], "matrix"))
write.csv(sigtabrA10, file = "differential-rel-abund/Rumen/A10.csv")

## plot
# Phylum order
x = tapply(sigtabrA10$log2FoldChange, sigtabrA10$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabrA10$Phylum = factor(as.character(sigtabrA10$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabrA10$log2FoldChange, sigtabrA10$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabrA10$Genus = factor(as.character(sigtabrA10$Genus), levels=names(x))
rA10 <- ggplot(sigtabrA10, aes(x=log2FoldChange, y=Genus, color=Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_brewer(palette = "BrBG") +
  ggtitle("B") +
  theme_classic()

ggsave("differential-rel-abund/Rumen/angus-rpc10-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

rA10 <- ggdotchart(sigtabrA10, y = "log2FoldChange", x = "Genus",
                   color = "Phylum", shape = "Phylum",
                   sorting = "descending",
                   add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                   dot.size = 5,
                   # add log2 fold change in the dot
                   label = round(sigtabrA10$log2FoldChange, 2),
                   font.label = list(color = "black", size = 11, vjust = -1.5),
                   rotate = TRUE,
                   ylab = "Log2 Fold Change",
                   xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("B")

## ---- RUMEN, Angus RPC15 ----
# add pseudo count of 1 
otu_table(A_rumen_counts15) <- otu_table(A_rumen_counts15) + 1

# phyloseq-DEseq
ddsrA15 <- phyloseq_to_deseq2(A_rumen_counts15, ~ Treatment)

# perform DESeq
dsrA15 <- DESeq(ddsrA15, test = "Wald", fitType = "parametric")

# get results
resrA15 = results(dsrA15, cooksCutoff = FALSE)
alpha = 0.05
sigtabrA15 = resrA15[which(resrA15$padj < alpha), ]
sigtabrA15
sigtabrA15 = cbind(as(sigtabrA15, "data.frame"), as(tax_table(A_rumen_counts15)[rownames(sigtabrA15), ], "matrix"))
write.csv(sigtabrA15, file = "differential-rel-abund/Rumen/A15.csv")

## plot
# Phylum order
x = tapply(sigtabrA15$log2FoldChange, sigtabrA15$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabrA15$Phylum = factor(as.character(sigtabrA15$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabrA15$log2FoldChange, sigtabrA15$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabrA15$Genus = factor(as.character(sigtabrA15$Genus), levels=names(x))
rA15

ggdotchart(sigtabrA15, aes(x=log2FoldChange, y=Genus, color = Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_brewer(palette = "BrBG") +
  ggtitle("C") +
  theme_classic()

rA15

ggsave("differential-rel-abund/Rumen/angus-rpc15-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

rA15 <- ggdotchart(sigtabrA15, y = "log2FoldChange", x = "Genus",
                   color = "Phylum", shape = "Phylum",
                   sorting = "descending",
                   add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                   dot.size = 5,
                   # add log2 fold change in the dot
                   label = round(sigtabrA15$log2FoldChange, 2),
                   font.label = list(color = "black", size = 11, vjust = -1.5),
                   rotate = TRUE,
                   ylab = "Log2 Fold Change",
                   xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("C")

# combine all the graphs together into one
library(patchwork)
(rA5 / rA10) | rA15

ggsave("differential-rel-abund/Rumen/angus-all.pdf", dpi = 600, width = 20, height = 16)

## ---- FECAL, Holstein RPC5 ----
# add pseudo count of 1 
otu_table(H_fecal_counts5) <- otu_table(H_fecal_counts5) + 1

# phyloseq-DEseq
ddsfH5 <- phyloseq_to_deseq2(H_fecal_counts5, ~ Treatment)

# perform DESeq
dsfH5 <- DESeq(ddsfH5, test = "Wald", fitType = "parametric")

# get results
resfH5 = results(dsfH5, cooksCutoff = FALSE)
alpha = 0.05
sigtabfH5 = resfH5[which(resfH5$padj < alpha), ]
sigtabfH5
sigtabfH5 = cbind(as(sigtabfH5, "data.frame"), as(tax_table(H_fecal_counts5)[rownames(sigtabfH5), ], "matrix"))
write.csv(sigtabfH5, file = "differential-rel-abund/Fecal/H5.csv")

## plot
# Phylum order
x = tapply(sigtabfH5$log2FoldChange, sigtabfH5$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabfH5$Phylum = factor(as.character(sigtabfH5$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabfH5$log2FoldChange, sigtabfH5$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabfH5$Genus = factor(as.character(sigtabfH5$Genus), levels=names(x))
fH5 <- ggplot(sigtabfH5, aes(x=log2FoldChange, y=Genus, color=Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_brewer(palette = "BrBG") +
  ggtitle("A") + 
  theme_classic()

ggsave("differential-rel-abund/Fecal/holstein-rpc5-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

fH5 <- ggdotchart(sigtabfH5, y = "log2FoldChange", x = "Genus",
                  color = "Phylum", shape = "Phylum",
                  sorting = "descending",
                  add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                  dot.size = 5,
                  # add log2 fold change in the dot
                  label = round(sigtabfH5$log2FoldChange, 2),
                  font.label = list(color = "black", size = 11, vjust = -1.5),
                  rotate = TRUE,
                  ylab = "Log2 Fold Change",
                  xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("A")

## ---- FECAL, Holstein RPC10 ----
# add pseudo count of 1 
otu_table(H_fecal_counts10) <- otu_table(H_fecal_counts10) + 1

# phyloseq-DEseq
ddsfH10 <- phyloseq_to_deseq2(H_fecal_counts10, ~ Treatment)

# perform DESeq
dsfH10 <- DESeq(ddsfH10, test = "Wald", fitType = "parametric")

# get results
resfH10 = results(dsfH10, cooksCutoff = FALSE)
alpha = 0.05
sigtabfH10 = resfH10[which(resfH10$padj < alpha), ]
sigtabfH10
sigtabfH10 = cbind(as(sigtabfH10, "data.frame"), as(tax_table(H_fecal_counts10)[rownames(sigtabfH10), ], "matrix"))
write.csv(sigtabfH10, file = "differential-rel-abund/Fecal/H10.csv")

## plot
# Phylum order
x = tapply(sigtabfH10$log2FoldChange, sigtabfH10$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabfH10$Phylum = factor(as.character(sigtabfH10$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabfH10$log2FoldChange, sigtabfH10$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabfH10$Genus = factor(as.character(sigtabfH10$Genus), levels=names(x))
fH10 <- ggplot(sigtabfH10, aes(x=log2FoldChange, y=Genus, color=Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_brewer(palette = "BrBG") +
  ggtitle("B") +
  theme_classic()

ggsave("differential-rel-abund/Fecal/holstein-rpc10-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

fH10 <- ggdotchart(sigtabfH10, y = "log2FoldChange", x = "Genus",
                   color = "Phylum", shape = "Phylum",
                   sorting = "descending",
                   add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                   dot.size = 5,
                   # add log2 fold change in the dot
                   label = round(sigtabfH10$log2FoldChange, 2),
                   font.label = list(color = "black", size = 11, vjust = -1.5),
                   rotate = TRUE,
                   ylab = "Log2 Fold Change",
                   xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("B")

## ---- FECAL, Holstein RPC15 ----
# add pseudo count of 1 
otu_table(H_fecal_counts15) <- otu_table(H_fecal_counts15) + 1

# phyloseq-DEseq
ddsfH15 <- phyloseq_to_deseq2(H_fecal_counts15, ~ Treatment)

# perform DESeq
dsfH15 <- DESeq(ddsfH15, test = "Wald", fitType = "parametric")

# get results
resfH15 = results(dsfH15, cooksCutoff = FALSE)
alpha = 0.05
sigtabfH15 = resfH15[which(resfH15$padj < alpha), ]
sigtabfH15
sigtabfH15 = cbind(as(sigtabfH15, "data.frame"), as(tax_table(H_fecal_counts15)[rownames(sigtabfH15), ], "matrix"))
write.csv(sigtabfH15, file = "differential-rel-abund/Fecal/H15.csv")

## plot
# Phylum order
x = tapply(sigtabfH15$log2FoldChange, sigtabfH15$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabfH15$Phylum = factor(as.character(sigtabfH15$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabfH15$log2FoldChange, sigtabfH15$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabfH15$Genus = factor(as.character(sigtabfH15$Genus), levels=names(x))
fH15 <- ggplot(sigtabfH15, aes(x=log2FoldChange, y=Genus, color=Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_brewer(palette = "BrBG") +
  ggtitle("C") +
  theme_classic()
fH15
ggsave("differential-rel-abund/Fecal/holstein-rpc15-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

fH15 <- ggdotchart(sigtabfH15, y = "log2FoldChange", x = "Genus",
                   color = "Phylum", shape = "Phylum",
                   sorting = "descending",
                   add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                   dot.size = 5,
                   # add log2 fold change in the dot
                   label = round(sigtabfH15$log2FoldChange, 2),
                   font.label = list(color = "black", size = 11, vjust = -1.5),
                   rotate = TRUE,
                   ylab = "Log2 Fold Change",
                   xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("C")

(fH5 / fH15) | fH10
ggsave("differential-rel-abund/Fecal/holstein-all.pdf", dpi = 600, width = 16, height = 12)

## ---- FECAL, Angus RPC5 ----
# add pseudo count of 1 
otu_table(A_fecal_counts5) <- otu_table(A_fecal_counts5) + 1

# phyloseq-DEseq
ddsfA5 <- phyloseq_to_deseq2(A_fecal_counts5, ~ Treatment)

# perform DESeq
dsfA5 <- DESeq(ddsfA5, test = "Wald", fitType = "parametric")

# get results
resfA5 = results(dsfA5, cooksCutoff = FALSE)
alpha = 0.05
sigtabfA5 = resfA5[which(resfA5$padj < alpha), ]
sigtabfA5
sigtabfA5 = cbind(as(sigtabfA5, "data.frame"), as(tax_table(A_fecal_counts5)[rownames(sigtabfA5), ], "matrix"))
write.csv(sigtabfA5, file = "differential-rel-abund/Fecal/A5.csv")

## plot
# Phylum order
x = tapply(sigtabfA5$log2FoldChange, sigtabfA5$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabfA5$Phylum = factor(as.character(sigtabfA5$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabfA5$log2FoldChange, sigtabfA5$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabfA5$Genus = factor(as.character(sigtabfA5$Genus), levels=names(x))
fA5 <- ggplot(sigtabfA5, aes(x=log2FoldChange, y=Genus, color=Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_brewer(palette = "BrBG") +
  ggtitle("A") +
  theme_classic()

ggsave("differential-rel-abund/Fecal/angus-rpc5-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

fA5 <- ggdotchart(sigtabfA5, y = "log2FoldChange", x = "Genus",
                  color = "Phylum", shape = "Phylum",
                  sorting = "descending",
                  add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                  dot.size = 5,
                  # add log2 fold change in the dot
                  label = round(sigtabfA5$log2FoldChange, 2),
                  font.label = list(color = "black", size = 11, vjust = -1.5),
                  rotate = TRUE,
                  ylab = "Log2 Fold Change",
                  xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("A")

## ---- FECAL, Angus RPC10 ----
# add pseudo count of 1 
otu_table(A_fecal_counts10) <- otu_table(A_fecal_counts10) + 1

# phyloseq-DEseq
ddsfA10 <- phyloseq_to_deseq2(A_fecal_counts10, ~ Treatment)

# perform DESeq
dsfA10 <- DESeq(ddsfA10, test = "Wald", fitType = "parametric")

# get results
resfA10 = results(dsfA10, cooksCutoff = FALSE)
alpha = 0.05
sigtabfA10 = resfA10[which(resfA10$padj < alpha), ]
sigtabfA10
sigtabfA10 = cbind(as(sigtabfA10, "data.frame"), as(tax_table(A_fecal_counts10)[rownames(sigtabfA10), ], "matrix"))
write.csv(sigtabfA10, file = "differential-rel-abund/Fecal/A10.csv")

## plot
# Phylum order
x = tapply(sigtabfA10$log2FoldChange, sigtabfA10$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabfA10$Phylum = factor(as.character(sigtabfA10$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabfA10$log2FoldChange, sigtabfA10$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabfA10$Genus = factor(as.character(sigtabfA10$Genus), levels=names(x))
fA10 <- ggplot(sigtabfA10, aes(x=log2FoldChange, y=Genus, color=Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_brewer(palette = "BrBG") +
  ggtitle("B") +
  theme_classic()

ggsave("differential-rel-abund/Fecal/angus-rpc10-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

fA10 <- ggdotchart(sigtabfA10, y = "log2FoldChange", x = "Genus",
                   color = "Phylum", shape = "Phylum",
                   sorting = "descending",
                   add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                   dot.size = 5,
                   # add log2 fold change in the dot
                   label = round(sigtabfA10$log2FoldChange, 2),
                   font.label = list(color = "black", size = 11, vjust = -1.5),
                   rotate = TRUE,
                   ylab = "Log2 Fold Change",
                   xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("B")

## ---- FECAL, Angus RPC15 ----
# add pseudo count of 1 
otu_table(A_fecal_counts15) <- otu_table(A_fecal_counts15) + 1

# phyloseq-DEseq
ddsfA15 <- phyloseq_to_deseq2(A_fecal_counts15, ~ Treatment)

# perform DESeq
dsfA15 <- DESeq(ddsfA15, test = "Wald", fitType = "parametric")

# get results
resfA15 = results(dsfA15, cooksCutoff = FALSE)
alpha = 0.05
sigtabfA15 = resfA15[which(resfA15$padj < alpha), ]
sigtabfA15
sigtabfA15 = cbind(as(sigtabfA15, "data.frame"), as(tax_table(A_fecal_counts15)[rownames(sigtabfA15), ], "matrix"))
write.csv(sigtabfA15, file = "differential-rel-abund/Fecal/A15.csv")

## plot
# Phylum order
x = tapply(sigtabfA15$log2FoldChange, sigtabfA15$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabfA15$Phylum = factor(as.character(sigtabfA15$Phylum), levels=names(x))
# Genus order
x = tapply(sigtabfA15$log2FoldChange, sigtabfA15$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabfA15$Genus = factor(as.character(sigtabfA15$Genus), levels=names(x))
fA15 <- ggplot(sigtabfA15, aes(x=log2FoldChange, y=Genus, color=Phylum, shape=Phylum)) + geom_point(size=6) + geom_vline(xintercept = 0) +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) +
  scale_color_brewer(palette = "BrBG") +
  ggtitle("C") + 
  theme_classic()

ggsave("differential-rel-abund/Fecal/angus-rpc15-diff-rel-abund-deseq.pdf", dpi = 600, width = 12, height = 7)

fA15 <- ggdotchart(sigtabfA15, y = "log2FoldChange", x = "Genus",
                   color = "Phylum", shape = "Phylum",
                   sorting = "descending",
                   add = "segments", add.params = list(size = 2, alpha = 0.5, linetype = "dashed"),
                   dot.size = 5,
                   # add log2 fold change in the dot
                   label = round(sigtabfA15$log2FoldChange, 2),
                   font.label = list(color = "black", size = 11, vjust = -1.5),
                   rotate = TRUE,
                   ylab = "Log2 Fold Change",
                   xlab = "Genus") +
  geom_hline(yintercept = 0, linetype = 1) +
  scale_y_continuous(expand = expansion(add = 0.5)) +
  # remove legend label
  labs(color = "Phylum", shape = "Phylum") + 
  scale_color_brewer(palette = "BrBG") + 
  ggtitle("C")

fA5 | (fA10 / fA15)
ggsave("differential-rel-abund/Fecal/angus-all.pdf", dpi = 600, width = 17, height = 14)
