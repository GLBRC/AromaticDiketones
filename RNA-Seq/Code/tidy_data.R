# Read in data tables
# Tidy
# Save as .Rdata files in the project directory

# Source: Kevin's processed files in his lab notebook
# Raw data files: /mnt/bigdata/processed_data/Alex/
# See Kevin's lab notebook for processing notes


library(tidyverse)
library(readxl)
#library(edgeR)

# Kevin's code for testing differential expression in edgeR:
# all <- read.table(file = "HTseq_values.txt", header = TRUE, sep = "\t", row.names = 1)
# targets <- read.table(file = "targets.txt", header = TRUE, sep = "\t")
# aromatic <- factor(targets$Condition)
# rep <- factor(targets$Rep)
# design <- model.matrix(~0 + aromatic + rep)
# dge <- DGEList(counts = all)
# dge <- calcNormFactors(dge)
# x<-estimateGLMCommonDisp(dge, design)
# y<-estimateGLMTrendedDisp(x, design)
# z<-estimateGLMTagwiseDisp(y, design)
# fit<-glmQLFit(z, design)
# my.contrasts <- makeContrasts(GvsCtrl = aromaticG - aromaticGlu, HvsCtrl = aromaticH - aromaticGlu, SvsCtrl = aromaticS - aromaticGlu, levels = design)
# GvsCtrl <- glmQLFTest(fit, contrast = my.contrasts[,"GvsCtrl"])
# FDR <- p.adjust(GvsCtrl$table$PValue, method = "BH")
# table <- cbind(GvsCtrl$table, FDR)
# write.table(table, file = "G_vs_Control_Results_DE_EdgeR.txt", quote = FALSE, sep = "\t")

# I'll start with his processed tables and tidy them up

RPKM <- read_excel("/Users/Alex/Downloads/RPKM.results_linear.xlsx", col_names = TRUE)
Gdiffs <- read.table("/Users/Alex/Downloads/G_vs_Control_Results_DE_EdgeR.txt", header = T)
Sdiffs <- read.table("/Users/Alex/Downloads/S_vs_Control_Results_DE_EdgeR.txt", header = T)
Hdiffs <- read.table("/Users/Alex/Downloads/H_vs_Control_Results_DE_EdgeR.txt", header = T)

RPKM2 <- gather(RPKM, Sample, Count, S_2:G_3)
condition <- c()
condition[which(RPKM2$Sample == "S_1" | RPKM2$Sample == "S_2" | RPKM2$Sample == "S_3")] <- "syringic_acid"
condition[which(RPKM2$Sample == "G_1" | RPKM2$Sample == "G_2" | RPKM2$Sample == "G_3")] <- "vanillic_acid"
condition[which(RPKM2$Sample == "H_1" | RPKM2$Sample == "H_2" | RPKM2$Sample == "H_3")] <- "p-hydroxybenzoic_acid"
condition[which(RPKM2$Sample == "Glu_1" | RPKM2$Sample == "Glu_2" | RPKM2$Sample == "Glu_3")] <- "glucose"
RPKM2$Condition <- condition

# The differential expression tests are already long, but need to be combined into a single dataset

Gdiffs$Condition <- "VAvsGlu"
Sdiffs$Condition <- "SAvsGlu"
Hdiffs$Condition <- "HAvsGlu"

diffex <- rbind(Gdiffs, Sdiffs, Hdiffs)

# Save as Rdata files. Note that I'll grab the gene info from my other RNASeq experiment.

saveRDS(RPKM2, "GSH_RPKM_normalize_read_counts.rds")
saveRDS(diffex, "GSH_Differential_expression_testing.rds")
