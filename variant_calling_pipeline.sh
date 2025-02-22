#!/bin/bash

# Ensure correct argument usage
if [ "$#" -ne 8 ] || [ "$1" != "--ref" ] || [ "$3" != "-1" ] || [ "$5" != "-2" ] || [ "$7" != "-o" ]; then
    echo "Usage: $0 --ref <Reference_Genome.fasta> -1 <Read1.fastq.gz> -2 <Read2.fastq.gz> -o <SampleName>"
    exit 1
fi

# Command-line arguments
REF_GENOME="$2"    # Path to reference genome
READ1="$4"         # Path to forward read
READ2="$6"         # Path to reverse read
SAMPLE_NAME="$8"   # Sample name for output
THREADS=8          # Number of CPU threads

# Create a directory for final output
mkdir -p "$SAMPLE_NAME"

# Ensure required tools are installed
for TOOL in bwa samtools bcftools; do
    if ! command -v $TOOL &> /dev/null; then
        echo "Error: $TOOL is not installed. Install it before proceeding."
        exit 1
    fi
done

echo "=== Step 1: Indexing the reference genome ==="
bwa index "$REF_GENOME"
samtools faidx "$REF_GENOME"

echo "=== Step 2: Aligning reads to reference genome ==="
bwa mem -t "$THREADS" "$REF_GENOME" "$READ1" "$READ2" | samtools view -bS -o "${SAMPLE_NAME}.unsorted.bam"

echo "=== Step 3: Sorting BAM by query name (Required for Fixmate) ==="
samtools sort -n -@ "$THREADS" -o "${SAMPLE_NAME}.queryname.bam" "${SAMPLE_NAME}.unsorted.bam"

echo "=== Step 4: Fixing mate information ==="
samtools fixmate -m "${SAMPLE_NAME}.queryname.bam" "${SAMPLE_NAME}.fixmate.bam"

echo "=== Step 5: Sorting BAM by coordinate (Required for Markdup) ==="
samtools sort -@ "$THREADS" -o "${SAMPLE_NAME}.positionsort.bam" "${SAMPLE_NAME}.fixmate.bam"

echo "=== Step 6: Marking duplicates ==="
samtools markdup -@ "$THREADS" "${SAMPLE_NAME}.positionsort.bam" "${SAMPLE_NAME}.dedup.bam"

echo "=== Step 7: Indexing the BAM file ==="
samtools index "${SAMPLE_NAME}.dedup.bam"

echo "=== Step 8: Calling variants using BCFtools ==="
bcftools mpileup -Ou -f "$REF_GENOME" "${SAMPLE_NAME}.dedup.bam" | \
bcftools call --ploidy 1 -mv -Oz -o "${SAMPLE_NAME}.vcf.gz"

echo "=== Step 9: Filtering Variants ==="
bcftools filter -s LowQual -e "QUAL<20" "${SAMPLE_NAME}.vcf.gz" -Oz -o "${SAMPLE_NAME}.filtered.vcf.gz"
bcftools index "${SAMPLE_NAME}.filtered.vcf.gz"

echo "=== Step 10: Converting Filtered VCF to Uncompressed Format ==="
bcftools view "${SAMPLE_NAME}.filtered.vcf.gz" -o "${SAMPLE_NAME}/${SAMPLE_NAME}.filtered.vcf"

if [ ! -f "${SAMPLE_NAME}/${SAMPLE_NAME}.filtered.vcf" ]; then
    echo "Error: Unzipped filtered VCF file not created!"
    exit 1
fi

echo "=== Cleaning Up Intermediate Files ==="
rm -f "${SAMPLE_NAME}.unsorted.bam" \
      "${SAMPLE_NAME}.queryname.bam" \
      "${SAMPLE_NAME}.fixmate.bam" \
      "${SAMPLE_NAME}.positionsort.bam" \
      "${SAMPLE_NAME}.dedup.bam" \
      "${SAMPLE_NAME}.dedup.bam.bai" \
      "${SAMPLE_NAME}.vcf.gz" \
      "${SAMPLE_NAME}.filtered.vcf.gz" \
      "${SAMPLE_NAME}.filtered.vcf.gz.csi"

echo "=== Variant Calling Pipeline Completed Successfully ==="
echo "Final output stored in: ${SAMPLE_NAME}/${SAMPLE_NAME}.filtered.vcf"