Comprehensive Genomic Workflow
Reference Data Retrieval, Variant Calling, Custom Database Setup, Annotation & Missense Variant Filtering
This repository provides an automated genomic analysis pipeline that integrates reference data retrieval, variant calling, annotation, and filtering. The workflow is tailored for haploid genomes and supports custom database building for SnpEff-based annotation.

🚀 Features
Download sequencing data from SRA.
Retrieve reference genome in FASTA format.
Perform variant calling using paired-end sequencing data.
Download and process GenBank files.
Convert GenBank to GFF3 format.
Set up a custom SnpEff annotation database.
Annotate variants and extract missense variants.
📥 Prerequisites
Ensure you have the following installed:

fastq-dump (from SRA Toolkit)
Python 3 with Biopython
bwa, samtools, bcftools (for variant calling)
snpEff (for variant annotation)
genbank_to_gff3.py (for GenBank to GFF3 conversion)
extract_missense_variants (for filtering missense variants)
📌 Step-by-Step Workflow
1️⃣ Download Sequence Data from SRA
bash
Copy
Edit
fastq-dump --gzip --split-3 SRR1735032
2️⃣ Download Reference Genome (FASTA)
bash
Copy
Edit
python download_fasta.py -i "KR781608.1"
3️⃣ Perform Variant Calling
bash
Copy
Edit
bash /Users/gerald/snpEff/variant_calling_pipeline --ref ~/snpEff/EVD/Makona-SLE.fa -1 ~/SRR1735032_1.fastq.gz -2 ~/SRR1735032_2.fastq.gz -o SRR1735032
4️⃣ Download GenBank File
bash
Copy
Edit
python download_genbank.py -i "Ebola virus zaire"
python download_genbank.py -i "KR781608.1"
5️⃣ Convert GenBank to GFF3
bash
Copy
Edit
python ~/genbank_to_gff3.py -i ~/KR781608.1.gb --output ~/KR781608.1.gff
6️⃣ Configure SnpEff for Custom Databases
Edit the snpEff.config file and add:

yaml
Copy
Edit
Pseudomonas_aeruginosa_PPF-1.genome : Pseudomonas aeruginosa PPF-1
Ebola_virus_zaire.genome : Ebola virus zaire
7️⃣ Create Custom Genome Directories
Create the necessary directories inside the data/ directory:

Copy
Edit
...Pseudomonas_aeruginosa_PPF-1 
...Ebola_virus_zaire
Each directory should contain:

pgsql
Copy
Edit
genes.gbk
sequence.gff
sequence.fa
8️⃣ Build Custom SnpEff Database
bash
Copy
Edit
java -Xmx8g -jar ~/snpEff/snpEff.jar build -genbank -v Ebola_virus_zaire
9️⃣ Annotate Variants
bash
Copy
Edit
java -Xmx8g -jar ~/snpEff/snpEff.jar Ebola_virus_zaire ~/SRR1735032.filtered.vcf > ~/SRR1735032.output.ann.vcf
🔟 Extract Missense Variants
bash
Copy
Edit
extract_missense_variants -i ~/SRR1735032.output.ann.vcf -o ~/SRR1735032.missense.vcf
📌 Expected Output
After running the pipeline, the key output files include:

Variant calling output: SRR1735032.filtered.vcf
Annotated variants: SRR1735032.output.ann.vcf
Filtered missense variants: SRR1735032.missense.vcf
🛠️ Troubleshooting
Ensure all paths are correct before running the scripts.
If SnpEff fails to build a custom database, check that genes.gbk, sequence.gff, and sequence.fa are correctly formatted.
If variant calling fails, verify that BWA, Samtools, and BCFtools are installed and correctly linked.
📝 Citation
If you use this workflow in your research, please cite:

SnpEff: Cingolani et al., 2012.
BCFtools: Li et al., 2009.
Biopython: Cock et al., 2009.
