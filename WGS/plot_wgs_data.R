########################### WGS DATA REPRESENTATION ############################
## This script is used to represent the results from the WGS data analysis. In 
## The first part allows the obtention of the segregation histograms for all 
## detected variants in the cohort and for the pre-filtered variants.
## The second part of the script is used to represent the PCA plots used to 
## assess for possible population stratification (with and without the control
## IBS population from the 1000 Genomes Project).


library(ggplot2)
library(ggrepel)
library(dplyr)


##### Segregation histograms

## Input: CSV file with one variant per row and information about the number of
## families in which segregation of that variant occurs (segregation meaning
## that case and control siblings have different genotypes at that locus)

# Define input data
data <- read.csv("Path/To/File/out_all_variants_information.csv") # For all variants
# data <- read.csv("Path/To/File/out_selected_variants_information.csv") # For pre-filtered variants

# Create a new column num_segr_total
data$num_segr_total <- data$num_segr_ek + data$num_segr_misc

# Count the number of variants for each family count
family_counts <- table(data$num_segr_total)

# Convert the table to a data frame and rename columns
plot_data <- as.data.frame(family_counts)
colnames(plot_data) <- c("count", "frequency")

# Create a data frame with all possible values from 0 to 21
all_counts <- data.frame(count = 0:20)

# Merge with the original plot_data2 to retain counts where available
plot_data <- merge(all_counts, plot_data, by = "count", all.x = TRUE)

# Fill missing frequencies with 0
plot_data$frequency[is.na(plot_data$frequency)] <- 0

# Create a bar plot
p <- ggplot(plot_data, aes(x = factor(count), y = frequency)) +
  geom_bar(stat = "identity", fill = "#00AFBB") +
  labs(x = "Number of Families in which segregation occurs",
       y = "Number of Variants") +
  scale_x_discrete(labels = 0:20) +
  theme_minimal(base_size = 12)

# Define output file
ggsave("Path/To/File/ALLvar_distrib_histogram.png", p, width = 7, height = 7, bg = "white") # For all variants
# ggsave("Path/To/File/Filteredvar_distrib_histogram.png", p, width = 7, height = 7, bg = "white") # For pre-filtered variants



##### Population stratification PCA plots

## Input: eigenvectors and eigenvalues files obtained after running PLINK v1.9 PCA

### Not including IBS population

# Read eigenvectors file
eigenvec <- read.table("Path/To/File/CoKid.norm_pca.eigenvec", header=FALSE) 

# Color by individual status (KD, MIS-C and Control samples in different colors)
# Load phenotype data
pheno_data <- read.table("Path/To/File/pheno_data.csv", header=TRUE)

# Merge eigenvec and pheno_data based on ID
merged_data <- merge(eigenvec, pheno_data, by.x = "V1", by.y = "ID", all.x = TRUE)
merged_data$GRUPO <- as.factor(merged_data$GRUPO)

## Proportion of variance explained by the first 2 PCs
# Read eigenvalues file
eigenval <- read.table("Path/To/File/CoKid.norm_pca.eigenval", header=FALSE)

# Compute % of variance explained per PC as:
# % Variance explained = [(Eigenvalue of PC)/(Sum of all Eigenvalues)]*100
all_var_explained <- c()

for(i in 1:nrow(eigenval)) {
  row <- eigenval[i,]
  var_explained <- (row/sum(eigenval))*100
  all_var_explained <- c(all_var_explained, var_explained)
}

# Add it to the eigenvalues table
eigenval["var_explained"] <- all_var_explained

x_percent <- round(eigenval[1,"var_explained"],2)
y_percent <- round(eigenval[2,"var_explained"],2)


color_mapping <- c("C" = "black", "KD" = "red", "MISC" = "blue")
legend_labels <- c("Control", "KD", "MIS-C")

# Plot the first two principal components coloring by individual status
pca_plot_colored <- ggplot(merged_data, aes(x = V3, y = V4, label = V1, color = GRUPO)) +
  geom_point(aes(shape = STATUS), size = 2) +
  labs(x = paste("PC1 (", x_percent, "%)", sep = ""), y = paste("PC2 (", y_percent, "%)", sep = ""), color = "Individual Status") +
  scale_shape_manual(name = "Individual Status", values = c("Control" = 17, "KD" = 15, "MIS-C" = 16)) +
  scale_color_manual(name = "Individual Status", values = color_mapping, labels = legend_labels) +
  theme_classic()

# Save the plot as a PNG file
ggsave("Path/To/File/CoKid_populationstrat_pca.png", pca_plot_colored)


### Including IBS population

# Read eigenvectors file
eigenvec <- read.table("Path/To/File/CoKid.norm_IBS.merged_pca.eigenvec", header=FALSE) 

# Color by individual status (KD, MIS-C, Control and 1KG IBS samples in different colors)
# Load phenotype data
CoKid_pheno_data <- read.table("Path/To/File/pheno_data.csv", header=TRUE, sep = ",")

# Merge eigenvec and pheno_data based on CoKid ID
merged_data <- merge(eigenvec, CoKid_pheno_data, by.x = "V1", by.y = "ID", all.x = TRUE)

# Add information about IBS population to the merged dataframe
merged_data <- merged_data %>%
  mutate(STATUS = ifelse(substr(V1, 1, 2) == "HG", "IBS", STATUS))

merged_data <- merged_data[!is.na(merged_data$STATUS), ]

merged_data$STATUS <- as.factor(merged_data$STATUS)

## Proportion of variance explained by the first 2 PCs
# Read eigenvalues file
eigenval <- read.table("Path/To/File/CoKid.norm_IBS.merged_pca.eigenval", header=FALSE)

# Compute % of variance explained per PC as:
# % Variance explained = [(Eigenvalue of PC)/(Sum of all Eigenvalues)]*100
all_var_explained <- c()

for(i in 1:nrow(eigenval)) {
  row <- eigenval[i,]
  var_explained <- (row/sum(eigenval))*100
  all_var_explained <- c(all_var_explained, var_explained)
}

# Add it to the eigenvalues table
eigenval["var_explained"] <- all_var_explained

x_percent <- round(eigenval[1,"var_explained"],2)
y_percent <- round(eigenval[2,"var_explained"],2)

color_mapping <- c("Control" = "black", "IBS" = "tomato", "KD" = "skyblue1", "MIS-C" = "#99FF99")
legend_labels <- c("Control", "IBS", "KD", "MIS-C")

# Plot the first two principal components coloring by individual status
pca_plot_colored <- ggplot(merged_data, aes(x = V3, y = V4, label = V1, color = STATUS)) +
  geom_point(aes(shape = STATUS), size = 2) +
  labs(x = paste("PC1 (", x_percent, "%)", sep = ""), y = paste("PC2 (", y_percent, "%)", sep = ""), color = "Individual Status") +
  scale_shape_manual(name = "Individual Status", values = c("Control" = 17, "IBS" = 17, "KD" = 15, "MIS-C" = 16)) +
  scale_color_manual(name = "Individual Status", values = color_mapping, labels = legend_labels) +
  theme_classic()

# Save the plot as a PNG file
ggsave("Path/To/File/CoKid_IBS_populationstrat_pca.png", pca_plot_colored, width = 7, height = 7)

