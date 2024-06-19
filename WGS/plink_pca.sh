## This script is used to compute the PCA of the WGS data using PLINK v1.9

input_prefix="Path/To/File/CoKid.norm_IBS.merged" # Prefix of input PLINK binary files
output_prefix="Path/To/File/CoKid.norm_IBS.merged_pca"

plink --bfile ${input_prefix} --pca --out ${output_prefix}