# 🧬 Genomic Pipeline for Haploid Pathogens

### *Reference data retrieval, variant calling, custom database setup, annotation & missense variant filtering*

This repository provides a **comprehensive genomic analysis pipeline** for **haploid organisms**, integrating **reference data retrieval, variant calling, annotation, and missense variant filtering**.  
It also supports **custom SnpEff database setup** for genome annotation.

---

## 🚀 Features

✔️ **Download sequencing data** from SRA  
✔️ **Retrieve reference genome (FASTA format)**  
✔️ **Perform haploid variant calling** using paired-end sequencing data  
✔️ **Download and process GenBank files**  
✔️ **Convert GenBank to GFF3 format**  
✔️ **Set up a custom SnpEff annotation database**  
✔️ **Annotate variants and extract missense variants**  

---

## 📥 Prerequisites

Ensure you have the following installed:
- [`fastq-dump`](https://github.com/ncbi/sra-tools) (**SRA Toolkit**)
- [`Biopython`](https://biopython.org/) (`pip install biopython`)
- [`BWA`](http://bio-bwa.sourceforge.net/), [`Samtools`](http://www.htslib.org/), [`BCFtools`](http://www.htslib.org/)  
- [`SnpEff`](https://pcingola.github.io/SnpEff/) (for variant annotation)  
- `genbank_to_gff3.py` (for **GenBank to GFF3 conversion**)  
- `extract_missense_variants` (for **missense variant filtering**)  

---

## 📌 Step-by-Step Workflow

### **1️⃣ Download Sequence Data from SRA**

```

fastq-dump --gzip --split-3 SRR1735032

```

### **2️⃣ Download Reference Genome (FASTA)**

```
python download_fasta.py -i "KR781608.1"
```

### **3️⃣ Perform Haploid Variant Calling**

```

**Single sample**
bash /Users/gerald/snpEff/variant_calling_pipeline --ref ~/snpEff/EVD/Makona-SLE.fa -1 ~/SRR1735032_1.fastq.gz -2 ~/SRR1735032_2.fastq.gz -o SRR1735032

```


 **For a batch run, create a text file (samples.txt) containing one sample ID per line, e.g.:**

   SRR1735032
   SRR1735033
   SRR1735034

 **a For a batch run, execute the variant_calling_batch.sh script with your reference genome and sample list:**

```
bash variant_calling_batch.sh --ref reference.fasta -s samples.txt

```

### **4️⃣ Download GenBank File**

```
python download_genbank.py -i "Ebola virus zaire"
python download_genbank.py -i "KR781608.1"

```

### **5️⃣ Convert GenBank to GFF3**

```
python ~/genbank_to_gff3.py -i ~/KR781608.1.gb --output ~/KR781608.1.gff

```
Rename genbank.gbk to genes.gbk in the .../data directory
 
```


### **6️⃣ Configure SnpEff for Custom Databases**
Edit the `snpEff.config` file and add:

```

Ebola_virus_zaire.genome : Ebola virus zaire

```

### **7️⃣ Create Custom Genome Directories**
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

### **8️⃣ Build Custom SnpEff Database**

```
java -Xmx8g -jar ~/snpEff/snpEff.jar build -genbank -v Ebola_virus_zaire

```

### **9️⃣ Annotate Variants**

```
java -Xmx8g -jar ~/snpEff/snpEff.jar Ebola_virus_zaire ~/SRR1735032.filtered.vcf > ~/SRR1735032.output.ann.vcf

```

### **🔟 Extract Missense Variants**

```
extract_missense_variants -i ~/SRR1735032.output.ann.vcf -o ~/SRR1735032.missense.vcf

```

---

## 📌 Expected Output

After running the pipeline, the key output files include:
- **Variant calling output:** `SRR1735032.filtered.vcf`
- **Annotated variants:** `SRR1735032.output.ann.vcf`
- **Filtered missense variants:** `SRR1735032.missense.vcf`


```
## Variant Calling Using Minimap2 (For Large Draft Assemblies)

If you are working with long-read assemblies (PacBio/Nanopore), use Minimap2 for variant detection.

Install Minimap2

conda install -c bioconda minimap2

## Align Draft Genome to Reference

minimap2 -ax asm5 reference.fasta draft_genome.fasta > alignment.sam
Convert to BAM & Sort

samtools view -bS alignment.sam | samtools sort -o alignment_sorted.bam

samtools index alignment_sorted.bam

## Call Variants

bcftools mpileup -Ou -f reference.fasta alignment_sorted.bam | bcftools call -mv -Oz -o variants.vcf.gz

Output: variants.vcf.gz (VCF file containing SNPs and structural variants).

--ref <reference.fasta> → Reference genome.

-s sample_list.txt → A text file containing paths to multiple draft genome assemblies.


## Create a sample_list.txt with Paths to Draft Assemblies

Example (sample_list.txt):

/path/to/draft_genome1.fasta
/path/to/draft_genome2.fasta
/path/to/draft_genome3.fasta

## Run the Script

bash variant_calling_minimap2.sh --ref reference.fasta -s sample_list.txt

Output:

Sorted BAM files: variant_results/*.sorted.bam

Final Decompressed VCFs: variant_results/*.vcf (for each draft genome)

## 📌 Summary of Methods

Approach	Best For	Tool	Output

Whole-Genome Alignment (WGA)	Comparing multiple draft assemblies	nucmer (MUMmer4)	variants.snps

Read-Mapping Approach	If you have sequencing reads	BWA + BCFtools	variants.vcf.gz

Long-Read Draft Assemblies	Nanopore/PacBio assemblies	Minimap2 + BCFtools	variants.vcf.gz

## 🚀 Now you can perform variant calling on your draft genome assemblies!

  You could also perform annotation (SnpEff) or visualization, this can follow earlier steps described

---

## 🛠️ Troubleshooting

- Ensure all **paths are correct** before running the scripts.
  
- If **SnpEff fails to build a custom database**, check that `genes.gbk`, `sequence.gff`, and `sequence.fa` are correctly formatted.
  
- If **variant calling fails**, verify that **BWA, Samtools, and BCFtools** are installed and correctly linked.

---

## 📝 Citation

If you use this workflow in your research, please cite:
- **SnpEff:** Cingolani et al., 2012.
- **BCFtools:** Li et al., 2009.
- **Biopython:** Cock et al., 2009.

---


