# Giketone paper figures list

# Set path to AromaticDiketones folder here
setwd("/Users/Alex/Desktop/AromaticDiketones/")

# Load packages and data
library(tidyverse)
library(cowplot)
library(reshape2)

RPKM <- readRDS("RNA-Seq/Data/RPKM_normalized_read_counts.rds")
diff_ex <- readRDS("RNA-Seq/Data/Differential_expression_testing_results.rds")
gene_info <- readRDS("RNA-Seq/Data/Gene_information_Apr2019.rds")
gene_names <- read_csv("RNA-Seq/Data/gene_names.csv")
growth_data <- read_csv("Growth_curves/ligLNDO-deletion_growth_curves_cleaned.csv")
enzyme_data <- read_csv("Enzyme_kinetics/ligLND_enzymes_2020-07-09.csv")
GDK_tests <- read_csv("Growth_curves/Diketone_growth_and_HPLC.csv")
fix_locus_tags <- readRDS("RNA-Seq/Data/gene_names2locus_tags.rds")

fix_locus_tags$to_change <- as.character(fix_locus_tags$to_change)
fix_locus_tags$to_replace <- as.character(fix_locus_tags$to_replace)

for(i in 1:dim(fix_locus_tags)[1]){
  z <- which(RPKM$Gene == fix_locus_tags$to_change[i])
  RPKM$Gene[z] <- fix_locus_tags$to_replace[i]
  y <- which(diff_ex$Gene == fix_locus_tags$to_change[i])
  diff_ex$Gene[y] <- fix_locus_tags$to_replace[i]
}

diff_ex_control <- diff_ex[grep("Ctrl", diff_ex$Comparison), ]
substrate <- c()
substrate[which(diff_ex_control$Comparison == "PCAvsCtrl")] <- "PCA"
substrate[which(diff_ex_control$Comparison == "VAvsCtrl")] <- "vanillic acid"
substrate[which(diff_ex_control$Comparison == "VvsCtrl")] <- "vanillin"
substrate[which(diff_ex_control$Comparison == "FAvsCtrl")] <- "ferulic acid"
substrate[which(diff_ex_control$Comparison == "UNvsCtrl")] <- "GDK intermediate"
substrate[which(diff_ex_control$Comparison == "GDKvsCtrl")] <- "G-diketone"
diff_ex_control$Substrate <- substrate

condition <- c()
condition[grep("Glu_", RPKM$Sample)] <- "glucose"
condition[grep("PCA_", RPKM$Sample)] <- "PCA"
condition[grep("VA_", RPKM$Sample)] <- "vanillic acid"
condition[grep("V_", RPKM$Sample)] <- "vanillin"
condition[grep("FA_", RPKM$Sample)] <- "ferulic acid"
condition[grep("UN_", RPKM$Sample)] <- "GDK intermediate"
condition[grep("GDK_", RPKM$Sample)] <- "G-diketone"
RPKM$Condition <- condition
RPKM$Condition <- factor(RPKM$Condition, levels = c("glucose", "PCA", "vanillic acid", "vanillin", "ferulic acid", "GDK intermediate", "G-diketone"))

#WT on diketone only
GDK_tests$substrate <- NULL
long_GDK <- melt(GDK_tests, id.vars = c("culture", "hours"))

averaged_GDK <- aggregate(value ~ hours + variable,  data = long_GDK, mean)
sd_GDK <- aggregate(value ~ hours + variable,  data = long_GDK, sd)
averaged_GDK$error <- sd_GDK$value

density <- ggplot(averaged_GDK[which(averaged_GDK$variable == "klett"), ], aes(x = hours, y = value)) + geom_pointrange(aes(ymin = value - error, ymax = value + error), size = 0.5, alpha = 0.5) + theme_bw() + labs(x = "Hours", y = "Klett units") + scale_color_brewer(palette = "Paired") + theme(legend.position = "none")  + geom_line()

concentrations <- ggplot(data = averaged_GDK[which(averaged_GDK$variable == "G-diketone"), ], aes(x = hours, y = value)) + geom_pointrange(aes(ymin = value - error, ymax = value + error), size = 0.5, alpha = 0.5) + geom_line() + theme_bw() + labs(x = "Hours", y = "G-diketone \n mmol/L")

intermediates <- ggplot(data = averaged_GDK[which(averaged_GDK$variable == "GP"), ], aes(x = hours, y = value)) + geom_pointrange(aes(ymin = value - error, ymax = value + error), size = 0.5, alpha = 0.5) + geom_line() + theme_bw() + labs(x = "Hours", y = "Intermediate \n relative peak area")

GDK_plots <- plot_grid(density, concentrations, intermediates, nrow = 3, labels = c("A.", "B.", "C."))

save_plot("Plots/GDK_plots.pdf", GDK_plots, base_height = 6.1, base_aspect_ratio = 1)


# How many genes expressed in each condition? 

diff_ex_glucose <- diff_ex_control %>% filter(FDR < 0.01)
# Reorder factors to match the pathway
diff_ex_glucose$Comparison[which(diff_ex_glucose$Comparison == "FAvsCtrl" | diff_ex_glucose$Comparison == "PCAvsCtrl" | diff_ex_glucose$Comparison == "VAvsCtrl" | diff_ex_glucose$Comparison == "VvsCtrl")] <- "BaseGvsCtrl"
diff_ex_glucose <- diff_ex_glucose %>% mutate(Comparison = fct_relevel(Comparison, c("BaseGvsCtrl", "UNvsCtrl", "GDKvsCtrl")))

gene_info <- gene_info %>% filter(locus_tag %in% diff_ex_glucose$Gene)
pathway_order <- c()
for(i in 1:dim(gene_info)[1]){
  gene <- gene_info$locus_tag[i]
  diffs <- diff_ex_glucose[which(diff_ex_glucose$Gene == gene), ]
  if("BaseGvsCtrl" %in% diffs$Comparison){
    pathway_order[i] <- "BaseG"
  }else if ("UNvsCtrl" %in% diffs$Comparison){
    pathway_order[i] <- "GP-1"
  }else if ("GDKvsCtrl" %in% diffs$Comparison){
    pathway_order[i] <- "GDK"
  }
}
gene_info$Condition <- pathway_order
diff_ex_glucose <- inner_join(diff_ex_glucose, gene_info, by = c("Gene" = "locus_tag"))
diff_ex_glucose <- diff_ex_glucose %>% select(Gene, logFC, FDR, Comparison, Condition, Product)

make_legend <- ggplot(diff_ex_glucose, aes(x = Comparison, fill = Condition)) + geom_bar() + theme_bw() + labs(x = NULL, y = "Genes significantly differentially expressed from glucose") + scale_x_discrete(labels=c("Base G", "G-diketone", "GP-1")) + scale_fill_brewer(labels=c("Base G", "GP-1", "G-diketone"), name = "Substrate", palette = "Paired")

add_legend <- get_legend(make_legend)

diff_ex_glucose$Condition <- factor(diff_ex_glucose$Condition, levels = c("BaseG", "GP-1", "GDK"))
add_increase <- ggplot(diff_ex_glucose[which(diff_ex_glucose$logFC > 0), ], aes(x = Comparison, fill = Condition)) + geom_bar() + theme_bw() + labs(x = NULL, y = "# of genes with significantly different transcript abundance", title = "A. Increased") + scale_x_discrete(labels=c("Base G", "GP-1", "G-diketone")) + scale_fill_brewer(labels=c("Base G", "GP-1", "G-diketone"), name = "Substrate", palette = "Paired") + theme(legend.position = "none")  + scale_y_continuous(limits = c(0,325), expand = c(0, 0))

add_decrease <- ggplot(diff_ex_glucose[which(diff_ex_glucose$logFC < 0), ], aes(x = Comparison, fill = Condition)) + geom_bar() + theme_bw() + labs(x = NULL, y = NULL, title = "B. Decreased")  + scale_x_discrete(labels=c("Base G", "GP-1", "G-diketone")) + scale_fill_brewer(labels=c("Base G", "GP-1", "G-diketone"), name = "Substrate", palette = "Paired") + theme(legend.position = "none") + scale_y_continuous(limits = c(0,325), expand = c(0, 0))

additive_fig <- plot_grid(add_increase, add_decrease, add_legend, nrow = 1, rel_widths = c(2.1, 2, 0.5))

save_plot("Plots/Additive_gene_counts.pdf", additive_fig, base_height = 4, base_aspect_ratio = 2.5)

# While we're here, make the diketone gene expression datasets

# table1 <- diff_ex_glucose[which(diff_ex_glucose$Condition == "GDK"), ]
# table1 <- table1[order(table1$logFC, decreasing = T), ]
# table1 <- table1[, c(1, 6, 2, 3)]
# colnames(table1) <- c("Gene", "Product", "log2 Fold Change from Glucose", "FDR-corrected p-value")
# write.csv(table1, "RNA-Seq/Data/GDK_genes.csv")
# 
# # Table 2
# table2<- diff_ex_glucose[which(diff_ex_glucose$logFC < 0 & diff_ex_glucose$Condition == "GDK"), ]
# table2 <- table2[order(table2$logFC), ]
# table2 <- table2[, c(1, 5, 2)]
# write.csv(table2, "RNA-Seq/Data/GDK_downreg_genes.csv")

# Figure 3 - Expression of known aromatic genes by substrate 

gene_RPKM <- RPKM[which(RPKM$Gene %in% gene_names$Locus_tag | RPKM$Gene %in% gene_names$Codename),]
gene_RPKM$Name <- gene_names$Codename[match(gene_RPKM$Gene, gene_names$Locus_tag)]
gene_RPKM$Name[which(is.na(gene_RPKM$Name) == T)] <- gene_RPKM$Gene[which(is.na(gene_RPKM$Name) == T)]
gene_RPKM <- gene_RPKM[grep("xyl", gene_RPKM$Name, invert = T), ]

colorkey <- c()
colorkey[which(gene_RPKM$Name == "ligA"| gene_RPKM$Name == "ligB" | gene_RPKM$Name == "ligC"| gene_RPKM$Name == "ligI"| gene_RPKM$Name == "ligJ"| gene_RPKM$Name == "ligK"| gene_RPKM$Name == "ligU")] <- "group1"
colorkey[which(gene_RPKM$Name == "ligM"| gene_RPKM$Name == "metF"| gene_RPKM$Name == "folD"| gene_RPKM$Name == "purU" | gene_RPKM$Name == "ligV" | gene_RPKM$Name == "ferA"| gene_RPKM$Name == "ferB")] <- "group2"
colorkey[which(gene_RPKM$Name == "ligE"| gene_RPKM$Name == "ligF" | gene_RPKM$Name == "naGST" | gene_RPKM$Name == "baeB" | gene_RPKM$Name == "baeA" | gene_RPKM$Name == "desC"| gene_RPKM$Name == "desD"| gene_RPKM$Name == "desA" |gene_RPKM$Name == "ligL"| gene_RPKM$Name == "ligN"| gene_RPKM$Name == "ligD"| gene_RPKM$Name == "ligO")] <- "group3"
gene_RPKM$colorkey <- colorkey

gene_RPKM$Name <- factor(gene_RPKM$Name, levels = c("ligA", "ligB", "ligC", "ligI", "ligJ", "ligU", "ligK", "ligM", "folD", "metF", "purU", "ligV", "ferA", "ferB", "naGST", "baeA", "baeB", "desA", "desC", "desD", "ligE", "ligF", "ligL", "ligN", "ligD", "ligO"))

diff_ex_control$Name <- gene_names$Codename[match(diff_ex_control$Gene, gene_names$Locus_tag)]
diff_ex_control$Name[which(is.na(diff_ex_control$Name) == T)] <- diff_ex_control$Gene[which(is.na(diff_ex_control$Name) == T)]
#diff_ex_control$Name <- factor(diff_ex_control$Name, levels = levels(gene_RPKM$Name))

#sig_diffs <- diff_ex_control %>% filter(Name %in% readcounts$Gene & FDR < 0.01 & substrate == carbon)
sigs <- c()
for(i in 1:dim(gene_RPKM)[1]){
  if(gene_RPKM$Condition[i] == "glucose"){
    sigs[i] <- NA
  }else{
    test <- diff_ex_control$FDR[which(diff_ex_control$Name == gene_RPKM$Name[i] & diff_ex_control$Substrate == gene_RPKM$Condition[i])]
    if(length(test) > 1){
      sigs[i] <- mean(test)
    }else{
      sigs[i] <- test
    }
  }
}

gene_RPKM$sig <- sigs

plot_genes <- function(carbon){
  readcounts <- gene_RPKM[which(gene_RPKM$Condition == carbon), ]
  
  ggplot(readcounts, aes(x = Name, y = log(RPKM), fill = colorkey)) + geom_boxplot(alpha = 0.5) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") + labs(x = "", y = "log2 RPKM", title = paste(carbon)) + scale_fill_manual(values = c("red", "goldenrod", "dodgerblue")) + geom_point(data = readcounts[which(readcounts$sig < 0.05), ], aes(x = Name, y = 9.4), size = 5, color = "black", shape = "*") + geom_point(data = readcounts[which(readcounts$sig < 0.01), ], aes(x = Name, y = 9.7), size = 5, color = "black", shape = "*") + geom_point(data = readcounts[which(readcounts$sig < 0.01), ], aes(x = Name, y = 10), size = 5, color = "black", shape = "*") + background_grid(major = "xy", minor = "none") + panel_border(colour = "black", size = 1, linetype = 1,remove = FALSE) + ylim(c(0, 10))
}

pca_plot <- plot_genes("PCA")
va_plot <- plot_genes("vanillic acid")
v_plot <- plot_genes("vanillin")
fa_plot <- plot_genes("ferulic acid")
gp_plot <- plot_genes("GDK intermediate")
gdk_plot <- plot_genes("G-diketone")

save_plot("Plots/PCA_boxplots.pdf", pca_plot, base_height = 3, base_aspect_ratio = 2)
save_plot("Plots/VA_boxplots.pdf", va_plot, base_height = 3, base_aspect_ratio = 2)
save_plot("Plots/V_boxplots.pdf", v_plot, base_height = 3, base_aspect_ratio = 2)
save_plot("Plots/FA_boxplots.pdf", fa_plot, base_height = 3, base_aspect_ratio = 2)
save_plot("Plots/GP_boxplots.pdf", gp_plot, base_height = 3, base_aspect_ratio = 2)
save_plot("Plots/GDK_boxplots.pdf", gdk_plot, base_height = 3, base_aspect_ratio = 2)


# Enzyme data
averaged_data2 <- aggregate(Diketone ~ Sample + Time,  data = enzyme_data, mean)
sd_data2 <- aggregate(Diketone ~ Sample + Time,  data = enzyme_data, sd)
averaged_data2$error <- sd_data2$Diketone


enzyme_key <- c()
enzyme_key[which(enzyme_data$Sample == "Control+NADH" | enzyme_data$Sample == "Control-NADH")] <- "No enzyme"
enzyme_key[which(enzyme_data$Sample == "ligD-NADH" | enzyme_data$Sample == "ligD+NADH")] <- "ligD"
enzyme_key[which(enzyme_data$Sample == "ligN-NADH" | enzyme_data$Sample == "ligN+NADH")] <- "ligN"
enzyme_key[which(enzyme_data$Sample == "ligL-NADH" | enzyme_data$Sample == "ligL+NADH")] <- "ligL"
enzyme_data$Enzyme <- enzyme_key

enzyme_key <- c()
enzyme_key[which(averaged_data2$Sample == "Control+NADH" | averaged_data2$Sample == "Control-NADH")] <- "No enzyme"
enzyme_key[which(averaged_data2$Sample == "ligD-NADH" | averaged_data2$Sample == "ligD+NADH")] <- "ligD"
enzyme_key[which(averaged_data2$Sample == "ligN-NADH" | averaged_data2$Sample == "ligN+NADH")] <- "ligN"
enzyme_key[which(averaged_data2$Sample == "ligL-NADH" | averaged_data2$Sample == "ligL+NADH")] <- "ligL"
averaged_data2$Enzyme <- enzyme_key

NADH_key <- c()
NADH_key[which(averaged_data2$Sample == "Control+NADH" | averaged_data2$Sample == "ligD-NADH" | averaged_data2$Sample == "ligN-NADH" | averaged_data2$Sample == "ligL-NADH")] <- "-"
NADH_key[which(averaged_data2$Sample == "Control-NADH" | averaged_data2$Sample == "ligD+NADH" | averaged_data2$Sample == "ligN+NADH" | averaged_data2$Sample == "ligL+NADH")] <- "+"
averaged_data2$NADH <- NADH_key

enzyme_curve <- ggplot(averaged_data2[which(averaged_data2$Enzyme != "No enzyme"), ], aes(x = Time, y = Diketone, color = Enzyme, shape = NADH)) + geom_pointrange(aes(ymin = Diketone - error, ymax = Diketone + error), size = 0.5, alpha = 0.5) + theme_bw() + labs(title = "", x = "Hours", y = "mmol/L diketone") + scale_color_manual(labels=c("LigD", "LigN", "LigL"), name = "Enzyme", values = c("dodgerblue", "orange", "darkgrey"))  + geom_line()

save_plot("Plots/enzyme_curve.pdf", enzyme_curve, base_height = 3, base_aspect_ratio = 1.5)

# Figure S3

# Growth and enzyme data on ligN with the diketone
growth_data$`Culture number` <- as.character(growth_data$`Culture number`)
long_data <- melt(growth_data[, c(1, 5:dim(growth_data)[2])])

long_data$variable <- as.POSIXct(long_data$variable, format = c("%m/%d/%y %H:%M"))
long_data$Hours <- (long_data$variable - as.POSIXct("2019-10-07 15:25:00"))/(60*60)
long_data$Strain <- growth_data$Strain[match(long_data$`Culture number`, growth_data$`Culture number`)]
long_data$Substrate <- growth_data$Carbon[match(long_data$`Culture number`, growth_data$`Culture number`)]

averaged_data <- aggregate(value ~ Strain + Hours + Substrate,  data = long_data, mean)
sd_data <- aggregate(value ~ Strain + Hours + Substrate,  data = long_data, sd)
averaged_data$error <- sd_data$value

# Diketone data

growth_curve <- ggplot(averaged_data[which(averaged_data$Substrate == "G-diketone + glucose"), ], aes(x = Hours, y = value, color = Strain)) + geom_pointrange(aes(ymin = value - error, ymax = value + error), size = 0.5, alpha = 0.5) + theme_bw() + labs(x = "Hours", y = "Klett units") + scale_color_brewer(palette = "Paired") + labs(title = "1 mmol/L G-diketone + 1mmol/L glucose") + geom_line(aes(linetype = Strain))

# Glucose data


growth_curve2 <- ggplot(averaged_data[which(averaged_data$Substrate == "Glucose"), ], aes(x = Hours, y = value, color = Strain)) + geom_pointrange(aes(ymin = value - error, ymax = value + error), size = 0.5, alpha = 0.5) + theme_bw() + labs(x = "Hours", y = "Klett units") + scale_color_brewer(palette = "Paired") + labs(title = "2 mmol/L glucose") + theme(legend.position = "none")  + geom_line(aes(linetype = Strain))

combined_growth_curve <- plot_grid(growth_curve2, growth_curve, nrow = 1, labels = c("A.", "B."), rel_widths = c(1, 1.2))

save_plot("Plots/growth_curve.pdf", combined_growth_curve, base_height = 3, base_aspect_ratio = 3)

