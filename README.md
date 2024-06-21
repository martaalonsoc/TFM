# README

This repository contains the code used for the analysis of the Master's Thesis entitled "Multi-omics characterization of Kawasaki disease versus Multisystem inflammatory syndrome in children (MIS-C)" by Marta Alonso-Caubilla et al.

It is structured as follows:
- Proteomics: this folder contains the code used to analyze proteomic data.
    - PlasmaProteinAnalysis.R: script used to perform differential expression analysis and enrichment analysis for Olink data. It represents the results as volcano plots, heatmaps, barplots, venn diagrams and boxplots.
    - doEnrichRanalyses.R: generic function to make EnrichR analyses.
    - PlasmaProteinPCA.R: script used to make PCA plots for Olink data and to obtain plots that represent the top 20 proteins mostly contributing to PC1 and PC2 or PC2 and PC3.
- WGS: this folder contains the code used to analyze WGS data.
    - Germline_SNV: this folder contains the code used for the germline variant calling process.
    - variant_filtering_prioritization.py: script used to obtain a table with the IDs and additional information of the filtered variants.
    - vcfETL.py: class to parse a VCF file.
    - vcfETLExceptions.py: exceptions when parsing a VCF file.
    - build_selected_variants_vcf.sh: script used to obtain a filtered VCF file.
    - plink_pca.sh: script used to compute the PCA of the WGS data using PLINK v1.9.
    - plot_wgs_data.R: script used to represent the results from the WGS data analysis (histograms of variant segregation in families and population stratification PCA plots).
    - plink_association_fisher.sh: script used to perform standard case/control association analyses with the filtered variants using Fisher's exact test to generate significance using PLINK v1.9.
    - MAGMA.sh: script that summarizes the commands used to perform gene-based analysis and gene-set analysis (GSA) using MAGMA v1.10.
