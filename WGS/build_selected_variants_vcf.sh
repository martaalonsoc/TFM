## This script is used to obtain a filtered VCF file in which only the variants of interest that are contained in a txt file are included in the final VCF file

txt_file="Path/To/File/var_id.txt"
vcf_file="Path/To/File/input.norm.vcf"
output_file="Path/To/File/Filtered.vcf"

# Iterate through each line in the VCF file
while IFS= read -r line || [[ -n "$line" ]]; do
    # Check if the line is a header line (starts with "#") or a variant line
    if [[ $line == \#* ]]; then
        # If it's a header line, just copy it to the output file
        echo "$line" >> "$output_file"

    else
        # Extract relevant fields from the VCF file (columns 1, 2, 4, and 5)
        vcf_fields=$(echo "$line" | awk '{print $1"_"$2"_"$4"_"$5}')

        # Check if the extracted fields exist in the txt file
        if grep -qw "$vcf_fields" "$txt_file"; then
            # If the variant is present, append the line to the output file
            echo "$line" >> "$output_file"
        fi
    fi
done < "$vcf_file"
