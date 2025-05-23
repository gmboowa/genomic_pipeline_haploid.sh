# Genomic pipeline for analysis of pathogen sequence data

### Reference data retrieval, variant calling, custom database setup, annotation & missense variant filtering

This repository provides a **comprehensive genomic analysis pipeline** for **haploid organisms**, integrating **reference data retrieval, variant calling, annotation, & missense variant filtering**.  
It also supports **custom SnpEff database setup** for genome annotation.

---

## Features
| Step # | Task Description                                        |
|--------|---------------------------------------------------------|
| ✔️ 1   | Download sequencing data from SRA                       |
| ✔️ 2   | Human DNA sequence removal using hostile                |
| ✔️ 3   | Retrieve reference genome (FASTA format)                |
| ✔️ 4   | Perform haploid variant calling using paired-end data   |
| ✔️ 5   | Download and process GenBank files                      |
| ✔️ 6   | Convert GenBank to GFF3 format                          |
| ✔️ 7   | Set up a custom SnpEff annotation database              |
| ✔️ 8   | Annotate variants and extract missense variants         |

---

## Prerequisites

Ensure you have the following installed:
- [`fastq-dump`](https://github.com/ncbi/sra-tools) (**SRA Toolkit**)
- [`Biopython`](https://biopython.org/) (`pip install biopython`)
- [`BWA`](http://bio-bwa.sourceforge.net/), [`Samtools`](http://www.htslib.org/), [`BCFtools`](http://www.htslib.org/)  
- [`SnpEff`](https://pcingola.github.io/SnpEff/) (for variant annotation)  
- `genbank_to_gff3.py` (for **GenBank to GFF3 conversion**)  
- `extract_missense_variants` (**for missense variant filtering**)  

---

## Step-by-step workflow

### Download sequence data from SRA

```bash

fastq-dump --gzip --split-3 SRR1735032

```

### Human DNA sequence removal

```bash

*Illumina*

python3 ~~/hostile_clean_human.py --fasta GCF_009914755.4 -i ~/fastq_samples_list.txt -o Clean

```
**Illumina_fastq_samples_list.txt

| Forward Read (R1)                             | Reverse Read (R2)                             |
|-----------------------------------------------|-----------------------------------------------|
| `~/39907_1_83_R1.fastq.gz`                    | `~/39907_1_83_R2.fastq.gz`                    |
| `~/39907_1_84_R1.fastq.gz`                    | `~/39907_1_84_R2.fastq.gz`                    |
| `~/39907_1_85_R1.fastq.gz`                    | `~/39907_1_85_R2.fastq.gz`                    |

```
*Nanopore*

python3 ~~/hostile_clean_ont_human_minimap2.py --fasta GCF_009914755.4 -i ~/ONT_fastq_samples_list.txt -o Clean

```
**ONT_fastq_samples_list.txt

| Read (R)                                   | 
|--------------------------------------------|
| `~/39907_1_83.fastq.gz`                    | 
| `~/39907_1_84.fastq.gz`                    | 
| `~/39907_1_85.fastq.gz`                    | 


### Download reference genome (FASTA)

```bash

python download_fasta.py -i "KR781608.1"

```

### Perform haploid variant calling



**Single sample**

```bash

bash ~/snpEff/variant_calling_pipeline --ref ~/snpEff/EVD/Makona-SLE.fa -1 ~/SRR1735032_1.fastq.gz -2 ~/SRR1735032_2.fastq.gz -o SRR1735032

```


 ** For a batch run, create a text file (samples.txt) containing one sample ID per line, e.g.:**

   SRR1735032
   SRR1735033
   SRR1735034

 **a for a batch run, execute the variant_calling_batch.sh script with your reference genome and sample list:**

```bash

bash variant_calling_batch.sh --ref reference.fasta -s samples.txt

```

### Download GenBank file

```bash
python download_genbank.py -i "Ebola virus zaire"

python download_genbank.py -i "KR781608.1"

```

### Convert GenBank to GFF3

```bash

python ~/genbank_to_gff3.py -i ~/KR781608.1.gb --output ~/KR781608.1.gff


 
```


### Configure SnpEff for custom databases

Edit the `snpEff.config` file and add:

```

Ebola_virus_zaire.genome : Ebola virus zaire

```

### Create custom genome directories**

Create the necessary directories inside the `data/` directory:

```

...Ebola_virus_zaire

```
Each directory should contain:

```
genes.gbk
sequence.gff
sequence.fa
```

### Build custom SnpEff database

```bash

java -Xmx8g -jar ~/snpEff/snpEff.jar build -genbank -v Ebola_virus_zaire

```

### Annotate variants

```bash

java -Xmx8g -jar ~/snpEff/snpEff.jar Ebola_virus_zaire ~/SRR1735032.filtered.vcf > ~/SRR1735032.output.ann.vcf

```

### Extract missense variants

```bash
extract_missense_variants -i ~/SRR1735032.output.ann.vcf -o ~/SRR1735032.missense.vcf

```

---

## Expected output

After running the pipeline, the key output files include:
- **Variant calling output:** `SRR1735032.filtered.vcf`
- **Annotated variants:** `SRR1735032.output.ann.vcf`
- **Filtered missense variants:** `SRR1735032.missense.vcf`



## Variant calling using Minimap2 (for large draft assemblies)

If you are working with long-read assemblies (PacBio/Nanopore), use Minimap2 for variant detection.

## Install Minimap2

```bash

conda install -c bioconda minimap2

```

## Align draft genome to reference

```bash

minimap2 -ax asm5 reference.fasta draft_genome.fasta > alignment.sam

Convert to BAM & Sort

samtools view -bS alignment.sam | samtools sort -o alignment_sorted.bam

samtools index alignment_sorted.bam

```

## Call variants

```bash

bcftools mpileup -Ou -f reference.fasta alignment_sorted.bam | bcftools call -mv -Oz -o variants.vcf.gz

Output: variants.vcf.gz (VCF file containing SNPs and structural variants).

--ref <reference.fasta> → Reference genome.

-s sample_list.txt → A text file containing paths to multiple draft genome assemblies.

```
## Create a sample_list.txt with paths to draft assemblies

Example (sample_list.txt):

```bash

/path/to/draft_genome1.fasta
/path/to/draft_genome2.fasta
/path/to/draft_genome3.fasta

```

## Run the script

```bash

bash variant_calling_minimap2.sh --ref reference.fasta -s sample_list.txt

```

Output:

```bash

Sorted BAM files: variant_results/*.sorted.bam

Final Decompressed VCFs: variant_results/*.vcf (for each draft genome)

```

## Summary of methods

Approach	best for	tool	output

Whole-genome alignment (WGA) comparing multiple draft assemblies	nucmer (MUMmer4)	variants.snps

Read-Mapping Approach	If you have sequencing reads	BWA + BCFtools	variants.vcf.gz

Long-Read Draft Assemblies	Nanopore/PacBio assemblies	Minimap2 + BCFtools	variants.vcf.gz

## Now you can perform variant calling on your draft genome assemblies!

You could also perform annotation (SnpEff) or visualization, this can follow earlier steps described

---

##  Troubleshooting

- Ensure all **paths are correct** before running the scripts.
  
- If **SnpEff fails to build a custom database**, check that `genes.gbk`, `sequence.gff`, and `sequence.fa` are correctly formatted.
  
- If **variant calling fails**, verify that **BWA, Samtools, and BCFtools** are installed and correctly linked.

---

## Citation

If you use this workflow in your research, please cite:
- **SnpEff:** Cingolani et al., 2012.
- **BCFtools:** Li et al., 2009.
- **Biopython:** Cock et al., 2009.

---


