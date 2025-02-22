# 🧬 Genomic Pipeline for Haploid Organisms
### *Reference Data Retrieval, Variant Calling, Custom Database Setup, Annotation & Missense Variant Filtering*

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
```bash
fastq-dump --gzip --split-3 SRR1735032
```

### **2️⃣ Download Reference Genome (FASTA)**
```bash
python download_fasta.py -i "KR781608.1"
```

### **3️⃣ Perform Haploid Variant Calling**
```bash
bash /Users/gerald/snpEff/variant_calling_pipeline --ref ~/snpEff/EVD/Makona-SLE.fa -1 ~/SRR1735032_1.fastq.gz -2 ~/SRR1735032_2.fastq.gz -o SRR1735032
```

### **4️⃣ Download GenBank File**
```bash
python download_genbank.py -i "Ebola virus zaire"
python download_genbank.py -i "KR781608.1"
```

### **5️⃣ Convert GenBank to GFF3**
```bash
python ~/genbank_to_gff3.py -i ~/KR781608.1.gb --output ~/KR781608.1.gff
```

### **6️⃣ Configure SnpEff for Custom Databases**
Edit the `snpEff.config` file and add:
```plaintext
Pseudomonas_aeruginosa_PPF-1.genome : Pseudomonas aeruginosa PPF-1
Ebola_virus_zaire.genome : Ebola virus zaire
```

### **7️⃣ Create Custom Genome Directories**
Create the necessary directories inside the `data/` directory:
```plaintext
...Pseudomonas_aeruginosa_PPF-1 
...Ebola_virus_zaire
```
Each directory should contain:
```plaintext
genes.gbk
sequence.gff
sequence.fa
```

### **8️⃣ Build Custom SnpEff Database**
```bash
java -Xmx8g -jar ~/snpEff/snpEff.jar build -genbank -v Ebola_virus_zaire
```

### **9️⃣ Annotate Variants**
```bash
java -Xmx8g -jar ~/snpEff/snpEff.jar Ebola_virus_zaire ~/SRR1735032.filtered.vcf > ~/SRR1735032.output.ann.vcf
```

### **🔟 Extract Missense Variants**
```bash
extract_missense_variants -i ~/SRR1735032.output.ann.vcf -o ~/SRR1735032.missense.vcf
```

---

## 📌 Expected Output
After running the pipeline, the key output files include:
- **Variant calling output:** `SRR1735032.filtered.vcf`
- **Annotated variants:** `SRR1735032.output.ann.vcf`
- **Filtered missense variants:** `SRR1735032.missense.vcf`

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


