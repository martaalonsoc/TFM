##################### MAKE CUSTOM PCA PLOTS FOR OLINK DATA #####################
## This script is used to make custom PCA plots for the exploratory part of the
## Olink data analysis. Samples are colored according to their disease status.
## The second part is used to obtain the plots that represent the top 20 
## proteins mostly contributing to PC1 and PC2 or PC2 and PC3.

## INPUT: NPX data frame in long format

## Load libraries
library(OlinkAnalyze)
library(dplyr)
library(ggplot2)
library(stringr)
library(xlsx)
library(gplots)
library(gridExtra)
library(enrichR)
library(RColorBrewer,quietly=Tq, warn.conflicts=F, verbose=F)
library(plotrix)
library(openxlsx)
library(FactoMineR)
library(factoextra)
library(ggbiplot)
library(tidyr)
library(ggfortify)
library(pals)


## Define analysis
exp <- "CardiovascularPlasma"
# exp <- "InflammationPlasma"

## Define input data
if (exp == "CardiovascularPlasma") {
  dataFile        <- "Path/To/Input/Data/CardiovascularPlasma.xlsx"
  outPutFolder    <- "Path/To/Output/Folder/CardiovascularPlasma"
  assaysPerSample <- 92 
} else if (exp == "InflammationPlasma") {
  dataFile        <- "Path/To/Input/Data/InflammationPlasma.xlsx"
  outPutFolder    <- "Path/To/Output/Folder/InflammationPlasma"
  assaysPerSample <- 368
}

## Load data dataframe
data <- openxlsx::read.xlsx(xlsxFile = dataFile, sheet = 1)
data$Disease <- as.factor(data$Disease)

## Compute basic PCA using Olink's default function
pca_plot <- data %>% olink_pca_plot(df = ., label_samples = FALSE, byPanel = FALSE, color_g = "Disease", x_val = 1, y_val = 2, n_loadings = 20) # x_val and y_val are the PCs that are represented in the X and Y axes, respectively; n_loadings = 20 to plot top 20 proteins mostly contributing to PCs 1 and 2

## Extract the plot data
pca_plot_data <- pca_plot[[1]]$data

## Add sample ID data to the plot data
pca_plot_data["SampleID"] <- rownames(pca_plot_data)

## Merge plot data with general Olink and phenotype data
merged_data <- merge(pca_plot_data, data, by.x = "SampleID", by.y = "SampleID", all.x = FALSE)
merged_data$Disease <- as.factor(merged_data$Disease)

## Generate custom PCA plot
if (exp == "CardiovascularPlasma") {
  # Define colors for the custom PCA plot
  color_mapping <- c("COV mild" = "#FF9999", "COV moderate" = "#FF3333", "CT no fever" = "darkgrey", "KD acute" = "lightblue", "KD recovery" = "#3399FF", "MIS-C acute" = "#99FF99", "MIS-C recovery" = "#006633")
  legend_labels <- c("COV mild", "COV moderate", "CT no fever", "KD acute", "KD recovery", "MIS-C acute", "MIS-C recovery")
  
  # Generate custom PCA plot
  custom_pca_plot <- ggplot(merged_data, aes(x = PCX, y = PCY, color = Disease, shape = Disease)) +
    geom_point(size = 3) +
    labs(x = pca_plot[[1]]$labels$x, y = pca_plot[[1]]$labels$y, color = "Disease") +
    scale_color_manual(name = "Disease", values = color_mapping, labels = legend_labels) + 
    scale_shape_manual(values = c("COV mild" = 16, "COV moderate" = 16, "CT no fever" = 17, "KD acute" = 15, "KD recovery" = 15, "MIS-C acute" = 16, "MIS-C recovery" = 16)) +  
    theme_classic()

} else if (exp == "InflammationPlasma") {
  # Define colors for the custom PCA plot
  color_mapping <- c("COV mild" = "#FF9999", "COV moderate" = "#FF3333", "CT no fever" = "darkgrey", "CT with fever" = "black", "KD acute" = "lightblue", "KD recovery" = "#3399FF", "MIS-C acute" = "#99FF99", "MIS-C recovery" = "#006633")
  legend_labels <- c("COV mild", "COV moderate", "CT no fever", "CT with fever", "KD acute", "KD recovery", "MIS-C acute", "MIS-C recovery")
  
  # Generate custom PCA plot
  custom_pca_plot <- ggplot(merged_data, aes(x = PCX, y = PCY, color = Disease, shape = Disease)) +
    geom_point(size = 3) +
    labs(x = pca_plot[[1]]$labels$x, y = pca_plot[[1]]$labels$y, color = "Disease") +
    scale_color_manual(name = "Disease", values = color_mapping, labels = legend_labels) + 
    scale_shape_manual(values = c("COV mild" = 16, "COV moderate" = 16, "CT no fever" = 17, "CT with fever" = 17, "KD acute" = 15, "KD recovery" = 15, "MIS-C acute" = 16, "MIS-C recovery" = 16)) +  
    theme_classic()
  
}

## Save the plot as a PNG file
ggsave(paste0(outPutFolder, "Olink_pca_plot.png"), custom_pca_plot, width = 8, height = 7)


############ TOP 20 PROTEINS MOSTLY CONTRIBUTING TO PC1/2 AND PC2/3 #############

## Transform data re-formatted as columns (proteins) and rows (samples)
data2 <- data %>%
  select(SampleID, Assay, NPX, Disease) %>%
  pivot_wider(names_from = Assay, values_from = NPX)

## Compute PCA
pca_res <- prcomp(data2[, -c(1,2)], scale. = TRUE, center = TRUE)

## PC contribution plot (PCs 1 & 2) - USING PROTEINS THAT RESULT FROM APPLYING THE n_loadings = 20 MODIFIER TO olink_pca_plot() FUNCTION
interested_olinkIDs <- pca_plot[[1]][["layers"]][[2]][["data"]][["variables"]]

# Define the top 20 proteins that contribute most to the PCs 1 and 2 according to the results of olink_pca_plot(n_loadings = 20)
if (exp == "CardiovascularPlasma") {
  interested_proteins <- c("GP6", "SELP", "PECAM-1", "SORT1", "CASP-3", "PDGF subunit A", "JAM-A", "CXCL1", "TNF-R2", "TNFRSF14", "NEMO", "TNF-R1", "CD84", "ST2", "IL6", "PAI", "IL-18BP", "U-PAR", "IL-27", "IL2-RA")
} else if (exp == "InflammationPlasma") {
  interested_proteins <- c("IL12RB1", "DBNL", "CASP2", "TBC1D5", "CD4", "STX8", "PDLIM7", "EIF4G1", "NUB1", "MAP2K6", "NUDC", "DNAJA2", "CRKL", "SAMD9L", "DFFA", "PGF", "SKAP2", "LTBR", "CRIM1", "IKBKG")
}

pca_var <- fviz_pca_var(pca_res, axes = c(1, 2), select.var = list(name = interested_proteins), repel=TRUE, col.circle = NA, labelsize = 3, arrowsize = 0.25, axes.linesize = 0.25) +
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), panel.grid = element_blank(), element_text(size = 3), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8)) +
  labs(x = pca_plot[[1]]$labels$x, y = pca_plot[[1]]$labels$y, title = " ") 

ggsave(paste0(outPutFolder, "Olink_pca_top20_PC1_PC2.png"), pca_var, bg = "white", width = 7, height = 7)


## PC contribution plot (PCs 2 & 3) - USING PROTEINS THAT RESULT FROM APPLYING THE n_loadings = 20 MODIFIER TO olink_pca_plot() FUNCTION

interested_olinkIDs <- pca_plot[[1]][["layers"]][[2]][["data"]][["variables"]]

# Define the top 20 proteins that contribute most to the PCs 2 and 3 according to the results of olink_pca_plot(n_loadings = 20)
if (exp == "InflammationPlasma") {
  interested_proteins <- c("B4GALT1", "CCL23", "IL6", "TNFSF11", "LAMA4", "TNFSF13", "CDON", "LILRB4", "BSG", "TNF", "OSCAR", "CSF1", "FSTL3", "CXCL10", "HGF", "CXCL9", "SMOC2", "SIRPB1", "IL1RN", "IL15")
}

pca_var <- fviz_pca_var(pca_res, axes = c(2, 3), select.var = list(name = interested_proteins), repel=TRUE, col.circle = NA, labelsize = 3, arrowsize = 0.25, axes.linesize = 0.25) +
  theme(plot.title = element_text(hjust = 0.5), panel.background = element_blank(), panel.grid = element_blank(), element_text(size = 3), 
        axis.title = element_text(size = 10), axis.text = element_text(size = 8)) +
  labs(x = pca_plot[[1]]$labels$x, y = pca_plot[[1]]$labels$y, title = " ") 

ggsave(paste0(outPutFolder, "Olink_pca_top20_PC2_PC3.png"), pca_var, bg = "white", width = 7, height = 7)

