######################## OLINK DATA ANALYSIS PIPELINE ##########################
## This script contains the analysis part for the Olink data. It is used to 
## obtain custom volcano plots and heatmaps (ComplexHeatmap) for the Differential 
## Expression analysis part. It also performs enrichment analysis against different
## databases and it does the representation of the results.
## The second part of the script is used to obtain the Venn diagrams (and their
## information) that show the overlap and differences between the obtained DEP 
## between the contrasts of interest (Ctrl no fever vs KD Acute, Ctrl no fever vs 
## MIS-C Acute and KD Acute vs MIS-C Acute).
## The third part of the script is used to obtain boxplots of the NPX values of
## the proteins of interest grouping the samples by Disease status.

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
library(ComplexHeatmap)
library(ggvenn)

source("./doEnrichRanalyses.R")


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

## Significance thresholds for contrasts
usePval <- TRUE # False to use adjusted p-value
maxPval <- 0.05

## EnrichR parameters
enrichRdbs <- c("GO_Biological_Process_2018","WikiPathways_2016","KEGG_2016")
maxEnrichRpVal <- 0.05
allDataGO      <- data.frame()
runPairedTtest <- TRUE
EnrichR <- "RUN"

## Load data dataframe
data <- openxlsx::read.xlsx(xlsxFile = dataFile, sheet = 1)
data$Disease <- as.factor(data$Disease)

## Run paired t-test (represent results as volcano plots and ComplexHeatmaps) and Enrichment analysis (represent results as barplots)
if (runPairedTtest) {
  varToColourList <- c("Sex","AgeInterval","Disease")
  Cond <- sort(unique(data$custom2))
  nConds <- length(Cond)

  # Lists to store the results of each contrast so that they are later stored in different sheets of an xlsx file
  PairedTtest.dfListALL <- list()
  PairedTtest.dfListDEP <- list()
  
  depPairedTtest.df <- data.frame()
  for (i in 1:(nConds-1)){
    for (j in (i+1):nConds){
      sel <- data.frame()
      sel <- which(data$custom2 == Cond[i])
      sel <- c(sel,which(data$custom2 == Cond[j]))
      selCont <- data[sel,]
      selCont$custom2 <- factor(selCont$custom2)
      selCont$Subject.ID <- factor(selCont$Subject.ID)
      
      # Skip contrasts in which the number of cases is too low
      if (sum(table(selCont$custom2)/assaysPerSample > 1) < 2){
        next
      }
      
      # Skip contrasts in which the number of samples per group is 1 or smaller
      if (nrow(data %>% filter(custom2 == Cond[i]) %>% distinct(Subject.ID)) <= 1) {
        next
      }
      
      if (nrow(data %>% filter(custom2 == Cond[j]) %>% distinct(Subject.ID)) <= 1) {
        next
      }
      
      # Define contrasts in which a paired t-test must be performed (when comparing acute and recovery samples of the same patient)
      if (paste0(Cond[i], "_vs_", Cond[j]) == "MISC_AC_vs_MISC_REC") {
        ttest_results <- olink_ttest(df = selCont, variable = "custom2", pair_id = "Subject.ID")
      } else if (paste0(Cond[i], "_vs_", Cond[j]) == "Kaw_AC_vs_Kaw_REC") {
        ttest_results <- olink_ttest(df = selCont, variable = "custom2", pair_id = "Subject.ID")
      } else {
        ttest_results <- olink_ttest(df = selCont, variable = "custom2")
      }
      
      contName <- paste0(Cond[i], "_vs_", Cond[j])
      PairedTtest.dfListALL[[contName]] <- ttest_results
      openxlsx::write.xlsx(PairedTtest.dfListALL, paste0(outPutFolder, "/TTEST_PAIRED_ALL_CONTRASTS.xlsx"))
      
      
      if (usePval == FALSE) {
        contrastResults <- ttest_results[ttest_results$Adjusted_pval <= maxPval,]
      } else {
        contrastResults <- ttest_results[ttest_results$p.value <= maxPval,]
      }
      
      # Select names of proteins to show
      label_sel <- as.vector(as.data.frame(ttest_results[ttest_results$Adjusted_pval <= maxPval,"OlinkID"])$OlinkID)
      
      if (nrow(contrastResults) > 0) {
        depPairedTtest.df[paste0(Cond[i], "_vs_", Cond[j]),"DEP"] <- nrow(contrastResults)
        print(paste0(Cond[i], "_vs_", Cond[j], " ", nrow(contrastResults)))
        PairedTtest.dfListDEP[[contName]] <- contrastResults
        openxlsx::write.xlsx(PairedTtest.dfListDEP, paste0(outPutFolder, "/DEP_TTEST_PAIRED_ALL_CONTRASTS.xlsx"))
      }
      
      ### Volcano plot
      png(paste0(outPutFolder, "/PLOTS/", Cond[i], "_vs_", Cond[j], "_volcano.png"), 800, 800)
      
      p <- olink_volcano_plot(p.val_tbl = ttest_results,
                              x_lab = paste0(Cond[i], " - ", Cond[j]),
                              olinkid_list = label_sel)
      
      volcano_plot <- p + scale_color_manual(values = c("turquoise", "red"), labels = c("Adjusted p >= 0.05", "Adjusted p < 0.05")) +
        theme(axis.title = element_text(size = 14),  # Size of axes' titles
              axis.text = element_text(size = 12),   # Size of axes' values
              legend.title = element_blank(),        # Remove legend title
              legend.text = element_text(size = 14)) # Size of legend values
      
      print(volcano_plot)
      dev.off()
      
      
      ### Heatmap
      if (nrow(contrastResults) > 1) {
        
        label_sel_heatmap <- as.vector(contrastResults$OlinkID)
        
        if (((Cond[i] == "Ctrl_FEV" && Cond[j] == "Ctrl_noF") || (Cond[j] == "Ctrl_FEV" && Cond[i] == "Ctrl_noF"))){
          next
        } else {
          dataHeatmap <- selCont[selCont$OlinkID %in% label_sel_heatmap,]
          
          dataHeatmap <- dataHeatmap[,c("SampleID", "Assay", "NPX", varToColourList, "AgeK", "AgeM","custom2", "Disease", "Age")]
          
          # Order samples by Disease and within each Disease group order them by Age (ascending)
          dataHeatmap <- dataHeatmap[order(dataHeatmap$Disease, dataHeatmap$Age), ]
          
          dataHeatmapWide <- reshape(dataHeatmap, v.names = "NPX", idvar = "SampleID", timevar = "Assay", direction = "wide")
          row.names(dataHeatmapWide)<-dataHeatmapWide$SampleID
          
          NPXcolumns <- grep("NPX",names(dataHeatmapWide))
          NPXWide <- dataHeatmapWide[,NPXcolumns]
          names(NPXWide) <- gsub(paste0("NPX","."),"",names(NPXWide))
          
          outPlot <- paste0(outPutFolder, "/PLOTS/", Cond[i], "_vs_", Cond[j], "_by_Disease_heatmap.png")
          png(outPlot, width=25,height=25,units="cm", res=600)
          
          
          # Define custom colors for each contrast
          custom_colors <- list(
            "CovMild_vs_CovMod" = c("CovMild" = "#FF9999", "CovMod" = "#FF3333"),
            "CovMild_vs_Ctrl_noF" = c("CovMild" = "#FF9999", "Ctrl_noF" = "darkgrey"),
            "CovMild_vs_Ctrl_FEV" = c("CovMild" = "#FF9999", "Ctrl_FEV" = "black"),
            "CovMild_vs_Kaw_AC" = c("CovMild" = "#FF9999", "Kaw_AC" = "lightblue"),
            "CovMild_vs_Kaw_REC" = c("CovMild" = "#FF9999", "Kaw_REC" = "#3399FF"),
            "CovMild_vs_MISC_AC" = c("CovMild" = "#FF9999", "MISC_AC" = "#99FF99"),
            "CovMild_vs_MISC_REC" = c("CovMild" = "#FF9999", "MISC_REC" = "#006633"),
            "CovMod_vs_Ctrl_noF" = c("CovMod" = "#FF3333", "Ctrl_noF" = "darkgrey"),
            "CovMod_vs_Ctrl_FEV" = c("CovMod" = "#FF3333", "Ctrl_FEV" = "black"),
            "CovMod_vs_Kaw_AC" = c("CovMod" = "#FF3333", "Kaw_AC" = "lightblue"),
            "CovMod_vs_Kaw_REC" = c("CovMod" = "#FF3333", "Kaw_REC" = "#3399FF"),
            "CovMod_vs_MISC_AC" = c("CovMod" = "#FF3333", "MISC_AC" = "#99FF99"),
            "CovMod_vs_MISC_REC" = c("CovMod" = "#FF3333", "MISC_REC" = "#006633"),
            "Ctrl_noF_vs_Ctrl_FEV" = c("Ctrl_noF" = "darkgrey", "Ctrl_FEV" = "black"),
            "Ctrl_noF_vs_Kaw_AC" = c("Ctrl_noF" = "darkgrey", "Kaw_AC" = "lightblue"),
            "Ctrl_noF_vs_Kaw_REC" = c("Ctrl_noF" = "darkgrey", "Kaw_REC" = "#3399FF"),
            "Ctrl_noF_vs_MISC_AC" = c("Ctrl_noF" = "darkgrey", "MISC_AC" = "#99FF99"),
            "Ctrl_noF_vs_MISC_REC" = c("Ctrl_noF" = "darkgrey", "MISC_REC" = "#006633"),
            "Ctrl_FEV_vs_Kaw_AC" = c("Ctrl_FEV" = "black", "Kaw_AC" = "lightblue"),
            "Ctrl_FEV_vs_Kaw_REC" = c("Ctrl_FEV" = "black", "Kaw_REC" = "#3399FF"),
            "Ctrl_FEV_vs_MISC_AC" = c("Ctrl_FEV" = "black", "MISC_AC" = "#99FF99"),
            "Ctrl_FEV_vs_MISC_REC" = c("Ctrl_FEV" = "black", "MISC_REC" = "#006633"),
            "Kaw_AC_vs_Kaw_REC" = c("Kaw_AC" = "lightblue", "Kaw_REC" = "#3399FF"),
            "Kaw_AC_vs_MISC_AC" = c("Kaw_AC" = "lightblue", "MISC_AC" = "#99FF99"),
            "Kaw_AC_vs_MISC_REC" = c("Kaw_AC" = "lightblue", "MISC_REC" = "#006633"),
            "Kaw_REC_vs_MISC_AC" = c("Kaw_REC" = "#3399FF", "MISC_AC" = "#99FF99"),
            "Kaw_REC_vs_MISC_REC" = c("Kaw_REC" = "#3399FF", "MISC_REC" = "#006633"),
            "MISC_AC_vs_MISC_REC" = c("MISC_AC" = "#99FF99", "MISC_REC" = "#006633")
          )
          
          # Get the levels being compared
          comparison <- paste0(Cond[i], "_vs_", Cond[j])
          
          # Get the colors for the specific comparison
          if (comparison %in% names(custom_colors)) {
            Disease_color_vector <- custom_colors[[comparison]]
          } else {
            # Default to rainbow palette if specific comparison not found
            rainbowPalette <- rainbow(length(levels(dataHeatmapWide$Disease)))
            Disease_color_vector <- setNames(rainbowPalette, levels(dataHeatmapWide$Disease))
          }
          
          # Age color
          dataHeatmapWide$Age <- as.numeric(as.character(dataHeatmapWide$Age))
          dataHeatmapWide$Age <- round(dataHeatmapWide$Age, 0)
          dataHeatmapWide$Age <- as.factor(dataHeatmapWide$Age)
          num_levels <- length(unique(dataHeatmapWide$Age))
          RdBu_palette <- colorRampPalette(c("red", "lightyellow", "blue"))(length(unique(dataHeatmapWide$Age)))
          Age_color_vector <- setNames(RdBu_palette, sort(unique(dataHeatmapWide$Age), decreasing = TRUE))
          
          # Sex color
          Sex_color_vector <- brewer.pal(8, "Dark2")[seq_along(levels(dataHeatmapWide$Sex))]
          names(Sex_color_vector) <- levels(dataHeatmapWide$Sex)
          
          # Annotation
          ha <- HeatmapAnnotation(df = data.frame(Disease = dataHeatmapWide$custom2, Age = dataHeatmapWide$Age, Sex = dataHeatmapWide$Sex),
                                  col = list(Disease = Disease_color_vector, Age = Age_color_vector, Sex = Sex_color_vector),
                                  annotation_height = 1, height = unit(6, "cm"),
                                  annotation_legend_param = list(Disease = list(direction = "horizontal", labels = unique(dataHeatmapWide$Disease)),
                                                                 Age = list(direction = "horizontal", labels = rev(names(Age_color_vector))),
                                                                 Sex = list(direction = "horizontal", labels = levels(dataHeatmapWide$Sex))))
          
          
          # Plot heatmap
          hm <- ComplexHeatmap::Heatmap(t(as.matrix(scale(NPXWide))), name = " ", cluster_columns = FALSE, cluster_rows = TRUE, show_row_names = TRUE,
                                        show_column_names = FALSE, show_heatmap_legend = TRUE, top_annotation = ha, row_dend_reorder = TRUE,
                                        heatmap_legend_param = list(direction = "vertical"), row_names_gp = gpar(fontsize = 5),
                                        heatmap_width = unit(10, "cm"), heatmap_height = unit(15, "cm"))
          
          draw(hm, annotation_legend_side = "right")
          
          dev.off()
        }
      }
      
      
      ### Enrichment analysis
      if (nrow(contrastResults) >= 3 & EnrichR == "RUN") {
        dataGO <- doEnrichRanalyses(as.character(contrastResults$Assay)
                                     , paste0(outPutFolder, "/ENRICHR/", Cond[i], "_vs_", Cond[j])
                                     , enrichRdbs
                                     , maxEnrichRpVal)
        if (!is.null(dataGO)) {
          if (nrow(allDataGO) == 0) {
            allDataGO <- dataGO[,c("Term","Adjusted.P.value")]
            names(allDataGO) <- c("Term",paste0(Cond[i], "_vs_", Cond[j]))
          } else {
            newNames <- c(names(allDataGO),paste0(Cond[i], "_vs_", Cond[j]))
            allDataGO <- merge(allDataGO,dataGO[,c("Term","Adjusted.P.value")],by="Term", all=T)
            names(allDataGO) <-  newNames         
          }
        }
      }
    }
  }
}



################################# VENN DIAGRAM #################################

## Input data: list of differentially expressed protein obtained running the previous section
if (exp == "CardiovascularPlasma") {
  dataFile        <- "Path/To/File/CardiovascularPlasma/DEP_TTEST_PAIRED_ALL_CONTRASTS.xlsx"
  outPutFolder    <- "Path/To/Output/Folder/CardiovascularPlasma"
  assaysPerSample <- 92
  
} else if (exp == "InflammationPlasma") {
  dataFile        <- "Path/To/File/InflammationPlasma/DEP_TTEST_PAIRED_ALL_CONTRASTS.xlsx"
  outPutFolder    <- "Path/To/Output/Folder/InflammationPlasma"
  assaysPerSample <- 368
}


## Load data of the DEP for the contrasts of interest
CtrlKD_DEP_data <- openxlsx::read.xlsx(xlsxFile = dataFile, sheet = "Ctrl_noF_vs_Kaw_AC")
CtrlMISC_DEP_data <- openxlsx::read.xlsx(xlsxFile = dataFile, sheet = "Ctrl_noF_vs_MISC_AC")
KDMISC_DEP_data <- openxlsx::read.xlsx(xlsxFile = dataFile, sheet = "Kaw_AC_vs_MISC_AC")

## Merge the data of the DEP for the contrasts of interest
DEP_data <- list(
  Ctrl_KD = CtrlKD_DEP_data$Assay,
  Ctrl_MISC = CtrlMISC_DEP_data$Assay,
  KD_MISC = KDMISC_DEP_data$Assay
)

## Plot and save Venn Diagram
names(DEP_data) <- c(paste0("CT no fever vs KD Acute (", length(DEP_data$Ctrl_KD), " DEP)"), paste0("CT no fever vs MIS-C Acute (", length(DEP_data$Ctrl_MISC), " DEP)"), paste0("KD Acute vs MIS-C Acute (", length(DEP_data$KD_MISC), " DEP)"))
venn <- ggvenn(DEP_data, fill_color = c("#3399FF", "hotpink", "springgreen3"), stroke_size = 0.5, set_name_size = 4.5, text_size = 5,
               set_name_color = c("#3399FF", "hotpink", "springgreen3"), show_percentage = FALSE)
ggsave(paste0(outPutFolder, "/PLOTS/Ctrl_KD_MISC_venn.png"), venn, width = 8, height = 8, bg = "white")

## Obtain information of the DE proteins in each case and save the resulting file
venn_info <-  ggvenn(DEP_data, fill_color = c("#3399FF", "hotpink", "springgreen3"), stroke_size = 0.5, set_name_size = 4.5, text_size = 5,
                     set_name_color = c("#3399FF", "hotpink", "springgreen3"), show_elements = TRUE, label_sep = ",")

proteins <- venn_info[["layers"]][[4]][["data"]][["text"]]
group <- venn_info[["layers"]][[4]][["data"]][["name"]]
venn_info_list <- as.data.frame(proteins, group)
rownames(venn_info_list) <- c("Ctrl_KD", "Ctrl_MISC", "KD_MISC", "Ctrl_KD_Ctrl_MISC", "Ctrl_KD_KD_MISC", "Ctrl_MISC_KD_MISC", "All")
openxlsx::write.xlsx(venn_info_list, paste0(outPutFolder, "/Ctrl_KD_MISC_venn.xlsx"), rowNames = TRUE)




############### BOXPLOTS for specific proteins grouped by Disease ###############

plotProteins <- TRUE

## Define proteins to plot
if (exp == "CardiovascularPlasma") {
  proteins_boxplots <- c("PD-L2", "CCL24") # List of proteins to plot
} else if (exp == "InflammationPlasma") {
  proteins_boxplots <- c("GBP2", "IL17C", "CCL17") # List of proteins to plot
}


if (plotProteins) {
  for (protein in proteins_boxplots){
    
    # Filter data for the specific protein
    protein_data <- data[data$Assay == protein, ]
    protein_data$Disease <- as.character(protein_data$Disease)
    protein_data$Disease <- as.factor(protein_data$Disease)
    
    # Define custom colors
    disease_colors <- list(
      "CT no fever" = "darkgrey",
      "KD acute" = "lightblue",
      "KD recovery" = "#3399FF",
      "MIS-C acute" = "#99FF99",
      "MIS-C recovery" = "#006633"
    )
    
    # Create the boxplot
    box <- ggplot(protein_data, aes(x = Disease, y = NPX)) +
      geom_boxplot(outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
      geom_jitter(width = 0.2, size = 1.5, aes(colour = Disease)) +
      scale_color_manual(values = disease_colors) +
      labs(title = paste0("NPX Values for ", protein),
           x = "Disease Group",
           y = "NPX") +
      theme_minimal() +
      theme(legend.position = "none",
            plot.title = element_text(size = 15, hjust = 0.5),
            axis.title.x = element_text(size = 15, margin = margin(t = 10)),
            axis.title.y = element_text(size = 15, margin = margin(r = 10)),
            axis.text.x = element_text(size = 13, angle = 45, hjust = 1),
            axis.text.y = element_text(size = 12))
    
    # Save the boxplot
    ggsave(paste0(outPutFolder, "/BOXPLOTS/", protein, "_boxplot.png"), box, width = 8, height = 8, bg = "white")
  }
}

