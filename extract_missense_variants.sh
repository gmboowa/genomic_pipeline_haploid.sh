#!/bin/bash

# Ensure correct usage
if [ "$#" -ne 4 ] || [ "$1" != "-i" ] || [ "$3" != "-o" ]; then
    echo "Usage: $0 -i <input_vcf> -o <output_vcf>"
    exit 1
fi

# Input and output VCF files
VCF_FILE="$2"
OUTPUT_FILE="$4"

# Check if input file exists
if [ ! -f "$VCF_FILE" ]; then
    echo "Error: File '$VCF_FILE' not found!"
    exit 1
fi

# Extract header and missense variants
grep "^#" "$VCF_FILE" > "$OUTPUT_FILE"
grep -i "missense_variant" "$VCF_FILE" >> "$OUTPUT_FILE"

echo "Missense variants extracted to: $OUTPUT_FILE"
