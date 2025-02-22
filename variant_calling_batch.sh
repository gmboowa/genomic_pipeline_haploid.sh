#!/bin/bash

# Ensure correct argument usage
if [ "$#" -ne 4 ] || [ "$1" != "--ref" ] || [ "$3" != "-s" ]; then
    echo "Usage: $0 --ref <Reference_Genome.fasta> -s <Sample_List.txt>"
    exit 1
fi

# Command-line arguments
REF_GENOME="$2"    # Path to reference genome
SAMPLE_LIST="$4"   # File containing sample IDs
THREADS=8          # Number of CPU threads

# Ensure required tools are installed
for TOOL in bwa samtools bcftools; do
    if ! command -v $TOOL &> /dev/null; then
        echo "Error: $TOOL is not installed. Install it before proceeding."
        exit 1
    fi
done

# Read sample list and process each sample
while read -r SAMPLE_ID; do
    echo "=== Processing Sample: ${SAMPLE_ID} ==="
    
    READ1="${SAMPLE_ID}_1.clean_1.fastq.gz"  # Path to forward read
    READ2="${SAMPLE_ID}_2.clean_2.fastq.gz"  # Path to reverse read

    # Ensure the reads exist
    if [[ ! -f "$READ1" || ! -f "$READ2" ]]; then
        echo "Error: Missing read files for ${SAMPLE_ID}. Skipping..."
        continue
    fi

    # Create an output directory for the sample
    mkdir -p "$SAMPLE_ID"

    echo "=== Step 1: Indexing the reference genome ==="
    bwa index "$REF_GENOME"
    samtools faidx "$REF_GENOME"

    echo "=== Step 2: Aligning reads to reference genome ==="
    bwa mem -t "$THREADS" "$REF_GENOME" "$READ1" "$READ2" | samtools view -bS -o "${SAMPLE_ID}.unsorted.bam"

    echo "=== Step 3: Sorting BAM by query name (Required for Fixmate) ==="
    samtools sort -n -@ "$THREADS" -o "${SAMPLE_ID}.queryname.bam" "${SAMPLE_ID}.unsorted.bam"

    echo "=== Step 4: Fixing mate information ==="
    samtools fixmate -m "${SAMPLE_ID}.queryname.bam" "${SAMPLE_ID}.fixmate.bam"

    echo "=== Step 5: Sorting BAM by coordinate (Required for Markdup) ==="
    samtools sort -@ "$THREADS" -o "${SAMPLE_ID}.positionsort.bam" "${SAMPLE_ID}.fixmate.bam"

    echo "=== Step 6: Marking duplicates ==="
    samtools markdup -@ "$THREADS" "${SAMPLE_ID}.positionsort.bam" "${SAMPLE_ID}.dedup.bam"

    echo "=== Step 7: Indexing the BAM file ==="
    samtools index "${SAMPLE_ID}.dedup.bam"

    echo "=== Step 8: Calling variants using BCFtools ==="
    bcftools mpileup -Ou -f "$REF_GENOME" "${SAMPLE_ID}.dedup.bam" | \
    bcftools call --ploidy 1 -mv -Oz -o "${SAMPLE_ID}.vcf.gz"

    echo "=== Step 9: Filtering Variants ==="
    bcftools filter -s LowQual -e "QUAL<20" "${SAMPLE_ID}.vcf.gz" -Oz -o "${SAMPLE_ID}.filtered.vcf.gz"
    bcftools index "${SAMPLE_ID}.filtered.vcf.gz"

    echo "=== Step 10: Converting Filtered VCF to Uncompressed Format ==="
    bcftools view "${SAMPLE_ID}.filtered.vcf.gz" -o "${SAMPLE_ID}/${SAMPLE_ID}.filtered.vcf"

    if [ ! -f "${SAMPLE_ID}/${SAMPLE_ID}.filtered.vcf" ]; then
        echo "Error: Unzipped filtered VCF file not created!"
        exit 1
    fi

    echo "=== Cleaning Up Intermediate Files ==="
    rm -f "${SAMPLE_ID}.unsorted.bam" \
          "${SAMPLE_ID}.queryname.bam" \
          "${SAMPLE_ID}.fixmate.bam" \
          "${SAMPLE_ID}.positionsort.bam" \
          "${SAMPLE_ID}.dedup.bam" \
          "${SAMPLE_ID}.dedup.bam.bai" \
          "${SAMPLE_ID}.vcf.gz" \
          "${SAMPLE_ID}.filtered.vcf.gz" \
          "${SAMPLE_ID}.filtered.vcf.gz.csi"

    echo "=== Variant Calling Completed for Sample: ${SAMPLE_ID} ==="
    echo "Final output stored in: ${SAMPLE_ID}/${SAMPLE_ID}.filtered.vcf"

done < "$SAMPLE_LIST"

echo "=== Batch Variant Calling Pipeline Completed Successfully ==="
