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
enzyme_data2 <- read_csv("Enzyme_kinetics/GP-HPLC-Samplekey_2021-07-18.csv")
GDK_tests <- read_csv("Growth_curves/Diketone_growth_and_HPLC.csv")
fix_locus_tags <- readRDS("RNA-Seq/Data/gene_names2locus_tags.rds")

# Some of the locus tags don't fit the overall pattern - for example, many are 4-letter codes. Change here.
fix_locus_tags$to_change <- as.character(fix_locus_tags$to_change)
fix_locus_tags$to_replace <- as.character(fix_locus_tags$to_replace)

for(i in 1:dim(fix_locus_tags)[1]){
  z <- which(RPKM$Gene == fix_locus_tags$to_change[i])
  RPKM$Gene[z] <- fix_locus_tags$to_replace[i]
  y <- which(diff_ex$Gene == fix_locus_tags$to_change[i])
  diff_ex$Gene[y] <- fix_locus_tags$to_replace[i]
}

# Clean up the differential expression dataset - I'm only interested in comparisons to glucose control
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

#########

# Figure 1 - growth experiment

GDK_tests$substrate <- NULL
long_GDK <- melt(GDK_tests, id.vars = c("culture", "hours"))

averaged_GDK <- aggregate(value ~ hours + variable,  data = long_GDK, mean)
sd_GDK <- aggregate(value ~ hours + variable,  data = long_GDK, sd)
averaged_GDK$error <- sd_GDK$value

density <- ggplot(averaged_GDK[which(averaged_GDK$variable == "klett"), ], aes(x = hours, y = value)) + geom_pointrange(aes(ymin = value - error, ymax = value + error), size = 0.5, alpha = 0.5) + theme_bw() + labs(x = "Hours", y = "Klett units") + scale_color_brewer(palette = "Paired") + theme(legend.position = "none")  + geom_line()

concentrations <- ggplot(data = averaged_GDK[which(averaged_GDK$variable == "G-diketone" | averaged_GDK$variable == "GP-1" | averaged_GDK$variable == "threo-GD"), ], aes(x = hours, y = value, color = variable)) + geom_pointrange(aes(ymin = value - error, ymax = value + error), size = 0.5, alpha = 0.5) + geom_line() + theme_bw() + labs(x = "Hours", y = "mmol/L") + theme(legend.title = element_blank(), legend.position = "bottom") 

GDK_plots <- plot_grid(density, concentrations, nrow = 2, labels = c("A", "B"), rel_heights = c(1, 1.2))

save_plot("Plots/GDK_plots.pdf", GDK_plots, base_height = 5, base_aspect_ratio = 1)
save_plot("Plots/GDK_plots.png", GDK_plots, base_height = 5, base_aspect_ratio = 1)

############
# Figure S3
# How many genes expressed in each condition? 

diff_ex_glucose <- diff_ex_control %>% filter(FDR < 0.01)

# Rename comparisons to substrates
substrates <- c()
substrates[which(diff_ex_glucose$Comparison == "FAvsCtrl")] <- "ferulic acid"
substrates[which(diff_ex_glucose$Comparison == "PCAvsCtrl")] <- "PCA"
substrates[which(diff_ex_glucose$Comparison == "VAvsCtrl")] <- "vanillic acid"
substrates[which(diff_ex_glucose$Comparison == "VvsCtrl")] <- "vanillin"
substrates[which(diff_ex_glucose$Comparison == "UNvsCtrl")] <- "GP-1"
substrates[which(diff_ex_glucose$Comparison == "GDKvsCtrl")] <- "G-diketone"
diff_ex_glucose$Substrate <- substrates

diff_ex_glucose$Substrate <- factor(diff_ex_glucose$Substrate, levels = c("PCA", "vanillic acid", "vanillin", "ferulic acid", "GP-1", "G-diketone"))

add_increase <- ggplot(diff_ex_glucose[which(diff_ex_glucose$logFC > 0), ], aes(x = Substrate)) + geom_bar(stat = "count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) + labs(x = NULL, y = "# of genes with significantly\ndifferent transcript abundance", title = "A. Increased")  + scale_y_continuous(limits = c(0,325), expand = c(0, 0)) + scale_x_discrete(labels=c("PCA", "vanillic acid", "vanillin", "ferulic acid", "GP-1", "G-diketone"))

add_decrease <- ggplot(diff_ex_glucose[which(diff_ex_glucose$logFC < 0), ], aes(x = Substrate)) + geom_bar(stat = "count") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) + labs(x = NULL, y = " ", title = "B. Decreased")  + scale_y_continuous(limits = c(0,325), expand = c(0, 0)) + scale_x_discrete(labels=c("PCA", "vanillic acid", "vanillin", "ferulic acid", "GP-1", "G-diketone")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

additive_fig <- plot_grid(add_increase, add_decrease, nrow = 1)

save_plot("Plots/Additive_gene_counts_2021-05-24.pdf", additive_fig, base_height = 4, base_aspect_ratio = 2.5)

##########
# Figure 2 - Expression of known aromatic genes by substrate 

# Fetch known aromatic gnes
gene_RPKM <- RPKM[which(RPKM$Gene %in% gene_names$Locus_tag | RPKM$Gene %in% gene_names$Codename),]
gene_RPKM$Name <- gene_names$Codename[match(gene_RPKM$Gene, gene_names$Locus_tag)]
gene_RPKM$Name[which(is.na(gene_RPKM$Name) == T)] <- gene_RPKM$Gene[which(is.na(gene_RPKM$Name) == T)]
gene_RPKM <- gene_RPKM[grep("xyl", gene_RPKM$Name, invert = T), ]

# Colorcode by steps in metabolism
colorkey <- c()
colorkey[which(gene_RPKM$Name == "ligA"| gene_RPKM$Name == "ligB" | gene_RPKM$Name == "ligC"| gene_RPKM$Name == "ligI"| gene_RPKM$Name == "ligJ"| gene_RPKM$Name == "ligK"| gene_RPKM$Name == "ligU")] <- "group1"
colorkey[which(gene_RPKM$Name == "ligM"| gene_RPKM$Name == "metF"| gene_RPKM$Name == "folD"| gene_RPKM$Name == "purU" | gene_RPKM$Name == "ligV" | gene_RPKM$Name == "ferA"| gene_RPKM$Name == "ferB")] <- "group2"
colorkey[which(gene_RPKM$Name == "ligE"| gene_RPKM$Name == "ligF" | gene_RPKM$Name == "naGST" | gene_RPKM$Name == "baeB" | gene_RPKM$Name == "baeA" | gene_RPKM$Name == "desC"| gene_RPKM$Name == "desD"| gene_RPKM$Name == "desA" |gene_RPKM$Name == "ligL"| gene_RPKM$Name == "ligN"| gene_RPKM$Name == "ligD"| gene_RPKM$Name == "ligO")] <- "group3"
gene_RPKM$colorkey <- colorkey

# Set order for plotting
gene_RPKM$Name <- factor(gene_RPKM$Name, levels = c("ligA", "ligB", "ligC", "ligI", "ligJ", "ligU", "ligK", "ligM", "folD", "metF", "purU", "ligV", "ferA", "ferB", "naGST", "baeA", "baeB", "desA", "desC", "desD", "ligE", "ligF", "ligL", "ligN", "ligD", "ligO"))

diff_ex_control$Name <- gene_names$Codename[match(diff_ex_control$Gene, gene_names$Locus_tag)]
diff_ex_control$Name[which(is.na(diff_ex_control$Name) == T)] <- diff_ex_control$Gene[which(is.na(diff_ex_control$Name) == T)]

# Get significance levels
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

# Plot data for each substrate

plot_genes <- function(carbon){
  readcounts <- gene_RPKM[which(gene_RPKM$Condition == carbon), ]
  
  ggplot(readcounts, aes(x = Name, y = log(RPKM), fill = colorkey)) + geom_boxplot(alpha = 0.5)  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position = "none") + labs(x = "", y = "log2 RPKM", title = paste(carbon)) + scale_fill_manual(values = c("red", "goldenrod", "dodgerblue")) + geom_point(data = readcounts[which(readcounts$sig < 0.05), ], aes(x = Name, y = 9.4), size = 5, color = "black", shape = "*") + geom_point(data = readcounts[which(readcounts$sig < 0.01), ], aes(x = Name, y = 9.7), size = 5, color = "black", shape = "*") + geom_point(data = readcounts[which(readcounts$sig < 0.01), ], aes(x = Name, y = 10), size = 5, color = "black", shape = "*") + panel_border(colour = "black", size = 1, linetype = 1,remove = FALSE) + ylim(c(0, 10))
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

#######
# Figure 3 - 24 hour enzyme experiment

# Get averages and standard error for each timepoint
averaged_data2 <- aggregate(Diketone ~ Sample + Time,  data = enzyme_data, mean)
sd_data2 <- aggregate(Diketone ~ Sample + Time,  data = enzyme_data, sd)
averaged_data2$error <- sd_data2$Diketone

# Rename conditions
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

# Make the plot
enzyme_curve <- ggplot(averaged_data2[which(averaged_data2$Enzyme != "No enzyme"), ], aes(x = Time, y = Diketone, color = Enzyme, shape = NADH)) + geom_pointrange(aes(ymin = Diketone - error, ymax = Diketone + error), size = 0.5, alpha = 0.5) + theme_bw() + labs(title = "", x = "Hours", y = "mmol/L G-diketone") + scale_color_manual(labels=c("LigD", "LigN", "LigL"), name = "Enzyme", values = c("dodgerblue", "orange", "darkgrey"))  + geom_line() + theme(legend.position = "none")

# Do the 2nd GP panel
averaged_data3 <- aggregate(`GP (mM)` ~ `Sample Name` + Hours,  data = enzyme_data2, mean)
sd_data3 <- aggregate(`GP (mM)` ~ `Sample Name` + Hours,  data = enzyme_data2, sd)
averaged_data3$error <- sd_data3$`GP (mM)`

# Rename conditions
enzyme_key <- c()
enzyme_key[which(enzyme_data2$`Sample Name` == "No enzyme+NADH" | enzyme_data2$`Sample Name` == "Control-NADH")] <- "No enzyme"
enzyme_key[which(enzyme_data2$`Sample Name` == "LigD-NADH" | enzyme_data2$`Sample Name` == "LigD+NADH")] <- "LigD"
enzyme_key[which(enzyme_data2$`Sample Name` == "LigN-NADH" | enzyme_data2$`Sample Name` == "LigN+NADH")] <- "LigN"
enzyme_key[which(enzyme_data2$`Sample Name` == "LigL-NADH" | enzyme_data2$`Sample Name` == "LigL+NADH")] <- "LigL"
enzyme_data2$Enzyme <- enzyme_key

enzyme_key <- c()
enzyme_key[which(averaged_data3$`Sample Name` == "No enzyme+NADH" | averaged_data3$`Sample Name` == "Control-NADH")] <- "No enzyme"
enzyme_key[which(averaged_data3$`Sample Name` == "LigD-NADH" | averaged_data3$`Sample Name` == "LigD+NADH")] <- "LigD"
enzyme_key[which(averaged_data3$`Sample Name` == "LigN-NADH" | averaged_data3$`Sample Name` == "LigN+NADH")] <- "LigN"
enzyme_key[which(averaged_data3$`Sample Name` == "LigL-NADH" | averaged_data3$`Sample Name` == "LigL+NADH")] <- "LigL"
averaged_data3$Enzyme <- enzyme_key

NADH_key <- c()
NADH_key[which(averaged_data3$`Sample Name` == "No enzyme+NADH" | averaged_data3$`Sample Name` == "LigD-NADH" | averaged_data3$`Sample Name` == "LigN-NADH" | averaged_data3$`Sample Name` == "LigL-NADH")] <- "-"
NADH_key[which(averaged_data3$`Sample Name` == "No enzyme-NADH" | averaged_data3$`Sample Name` == "LigD+NADH" | averaged_data3$`Sample Name` == "LigN+NADH" | averaged_data3$`Sample Name` == "LigL+NADH")] <- "+"
averaged_data3$NADH <- NADH_key

averaged_data3$error[which(is.na(averaged_data3$error) == T)] <- 0

# Make the plot
enzyme_curve2 <- ggplot(averaged_data3[which(averaged_data3$Enzyme != "No enzyme"), ], aes(x = Hours, y = `GP (mM)`, color = Enzyme, shape = NADH)) + geom_pointrange(aes(ymin = `GP (mM)` - error, ymax = `GP (mM)` + error), size = 0.5, alpha = 0.5) + theme_bw() + labs(title = "", x = "Hours", y = "mmol/L GP-1") + scale_color_manual( name = "Enzyme", values = c("dodgerblue", "darkgrey", "orange"))  + geom_line()

legend <- get_legend(enzyme_curve2)

enzyme_curve2 <- ggplot(averaged_data3[which(averaged_data3$Enzyme != "No enzyme"), ], aes(x = Hours, y = `GP (mM)`, color = Enzyme, shape = NADH)) + geom_pointrange(aes(ymin = `GP (mM)` - error, ymax = `GP (mM)` + error), size = 0.5, alpha = 0.5) + theme_bw() + labs(title = "", x = "Hours", y = "mmol/L GP-1") + scale_color_manual( name = "Enzyme", values = c("dodgerblue", "darkgrey", "orange"))  + geom_line() + theme(legend.position = "none")


all_enzyme_curves <- plot_grid(enzyme_curve, enzyme_curve2, legend, rel_widths = c(1, 1, 0.25), nrow = 1, labels = c("A", "B", ""))

save_plot("Plots/enzyme_curve2.pdf", all_enzyme_curves, base_height = 3, base_aspect_ratio = 2.5)
save_plot("Plots/enzyme_curve2.png", all_enzyme_curves, base_height = 3, base_aspect_ratio = 2.5)

###########
# Figure S5 - deletion mutant growth data

# Convert growth curve data to long format and fix dates
growth_data$`Culture number` <- as.character(growth_data$`Culture number`)
long_data <- melt(growth_data[, c(1, 5:dim(growth_data)[2])])

long_data$variable <- as.POSIXct(long_data$variable, format = c("%m/%d/%y %H:%M"))
long_data$Hours <- (long_data$variable - as.POSIXct("2019-10-07 15:25:00"))/(60*60)
long_data$Strain <- growth_data$Strain[match(long_data$`Culture number`, growth_data$`Culture number`)]
long_data$Substrate <- growth_data$Carbon[match(long_data$`Culture number`, growth_data$`Culture number`)]

# Get means and SD as in the Fig 3 code
averaged_data <- aggregate(value ~ Strain + Hours + Substrate,  data = long_data, mean)
sd_data <- aggregate(value ~ Strain + Hours + Substrate,  data = long_data, sd)
averaged_data$error <- sd_data$value

# Plot growth on diketone + glucose

growth_curve <- ggplot(averaged_data[which(averaged_data$Substrate == "G-diketone + glucose"), ], aes(x = Hours, y = value, color = Strain)) + geom_pointrange(aes(ymin = value - error, ymax = value + error), size = 0.5, alpha = 0.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + labs(x = "Hours", y = "Klett units") + scale_color_brewer(palette = "Paired") + labs(title = "1 mmol/L G-diketone + 1mmol/L glucose") + geom_line(aes(linetype = Strain)) + ylim(0, 225) + xlim(0, 150)

# Plot growth on glucose

growth_curve2 <- ggplot(averaged_data[which(averaged_data$Substrate == "Glucose"), ], aes(x = Hours, y = value, color = Strain)) + geom_pointrange(aes(ymin = value - error, ymax = value + error), size = 0.5, alpha = 0.5) + labs(x = "Hours", y = "Klett units") + scale_color_brewer(palette = "Paired") + labs(title = "2 mmol/L glucose") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "none")  + geom_line(aes(linetype = Strain)) + xlim(0, 150) + ylim(0, 225)

combined_growth_curve <- plot_grid(growth_curve2, growth_curve, nrow = 1, labels = c("A.", "B."), rel_widths = c(1, 1.2))

save_plot("Plots/growth_curve.pdf", combined_growth_curve, base_height = 3, base_aspect_ratio = 3)

