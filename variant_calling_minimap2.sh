#!/bin/bash

# Ensure correct argument usage
if [[ "$#" -ne 4 ]] || [[ "$1" != "--ref" ]] || [[ "$3" != "-s" ]]; then
    echo "Usage: $0 --ref <reference.fasta> -s <sample_list.txt>"
    exit 1
fi

# Command-line arguments
REF_GENOME="$2"      # Reference genome file
SAMPLE_LIST="$4"     # List of draft genome assemblies

# Check if required tools are installed
for TOOL in minimap2 samtools bcftools; do
    if ! command -v $TOOL &> /dev/null; then
        echo "Error: $TOOL is not installed. Install it before proceeding."
        exit 1
    fi
done

# Index the reference genome
echo "=== Step 1: Indexing Reference Genome ==="
samtools faidx "$REF_GENOME"

# Create output directory
mkdir -p variant_results

# Loop through each sample in the list
while read -r SAMPLE; do
    if [[ ! -f "$SAMPLE" ]]; then
        echo "❌ ERROR: File $SAMPLE not found! Skipping..."
        continue
    fi

    SAMPLE_NAME=$(basename "$SAMPLE" .fasta)  # Extract filename without extension
    echo "Processing: $SAMPLE_NAME"

    # Step 2: Align draft genome to reference
    echo "=== Step 2: Aligning $SAMPLE_NAME to Reference ==="
    minimap2 -ax asm5 "$REF_GENOME" "$SAMPLE" | samtools view -bS -o "variant_results/${SAMPLE_NAME}.bam"

    # Step 3: Sort BAM file
    echo "=== Step 3: Sorting BAM File for $SAMPLE_NAME ==="
    samtools sort -o "variant_results/${SAMPLE_NAME}.sorted.bam" "variant_results/${SAMPLE_NAME}.bam"

    # Step 4: Index the sorted BAM file
    samtools index "variant_results/${SAMPLE_NAME}.sorted.bam"

    # Step 5: Call variants using BCFtools (Haploid mode)
    echo "=== Step 5: Calling Variants for $SAMPLE_NAME (Haploid Mode) ==="
    bcftools mpileup -Ou -f "$REF_GENOME" "variant_results/${SAMPLE_NAME}.sorted.bam" | \
    bcftools call --ploidy 1 -mv -Ov -o "variant_results/${SAMPLE_NAME}.vcf"

    echo "✅ Variants saved for $SAMPLE_NAME: variant_results/${SAMPLE_NAME}.vcf"

    # Cleanup intermediate files
    rm -f "variant_results/${SAMPLE_NAME}.bam"

done < "$SAMPLE_LIST"

echo "🎉 Variant Calling Completed! Results are in 'variant_results/'"