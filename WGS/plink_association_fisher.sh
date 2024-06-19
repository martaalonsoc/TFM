## This script is used to perform standard case/control association analyses with the 19,971 filtered single variants using Fisher's exact test to generate significance using PLINK v1.9

input_prefix="Path/To/File/FilteredVariants_KD_families" # Prefix of input PLINK binary files; KD families
# input_prefix="Path/To/File/FilteredVariants_MISC_families" # Prefix of input PLINK binary files; MIS-C families

plink --file ${input_prefix} --fisher