## This script summarizes the commands used to perform gene-based analysis and gene-set analysis (GSA) using MAGMA (v1.10)

## Step 1.- Annotation of pre-filtered variants to genes.
## This is done based on NCBI [38] gene definitions, mapping variants to a protein-coding gene if they are located in the 
## transcription region of that gene, or within 5 kilobase upstream or 3 kilobase downstream of the transcription region.

## Input: BIM file with the pre-filtered variants and reference gene location file
## Output: .genes.annot file with each row corresponding to a gene, containing the gene ID, a specification of the gene's location, and a list of SNP IDs of SNPs mapped to that gene

magma --annotate window=5,3 --snp-loc Path/To/File/FilteredVariants.bim --gene-loc ./NCBI38/NCBI38.gene.loc --out FilteredVariants


## Step 2.- Multi-model gene analysis.
## In this step the gene p-values and other gene-level metrics are computed.

## Input: pre-filtered variants data in binary PLINK format; phenotype data for each sample and .genes.annot file generated in the previous step 
## Output: .genes.out file with the gene analysis results; .genes.raw file with the same results plus gene correlations

magma --gene-model multi=all --bfile Path/To/File/FilteredVariants --pheno file=Path/To/File/pheno_data.txt --gene-annot ./FilteredVariants.genes.annot --out FilteredVariants


## Step 3.- Gene-set analysis.
## In this step basic gene-set analysis is performed.

## Input: .genes.raw file generated in the previous step and gene-set definition file (each row corresponds to a gene-set definition; the first value on each row is the name of the gene-set, followed by the IDs of genes in that gene-set)
## Output: .gsa.out file with the gene-set analysis results

magma --gene-results ./FilteredVariants.genes.raw --set-annot Path/To/File/GO_BP_EntrezID.txt --out FilteredVariants