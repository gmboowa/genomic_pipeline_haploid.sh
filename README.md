# Comprehensive genomic workflow - reference data retrieval, variant calling, custom database setup, annotation & missense variant filtering

Download from SRA

$ fastq-dump --gzip --split-3 SRR1735032

Download reference fasta

$ python download_fasta.py -i "KR781608.1"


Variant calling

$ bash /Users/gerald/snpEff/variant_calling_pipeline --ref ~/snpEff/EVD/Makona-SLE.fa -1 ~/SRR1735032_1.fastq.gz -2 ~/SRR1735032_2.fastq.gz -o SRR1735032

Download Genbank

$ python download_genbank.py -i "Ebola virus zaire"

$ python download_genbank.py -i "KR781608.1"

Download genbank_to_gff3.py

$ python ~/genbank_to_gff3.py -i ~/KR781608.1.gb --output ~/KR781608.1.gff

Edit the snpEff.config to add the following lines as needed.

Pseudomonas_aeruginosa_PPF-1.genome : Pseudomonas aeruginosa PPF-1

Ebola_virus_zaire.genome : Ebola virus zaire

Create the following directory with ..data/....

...Pseudomonas_aeruginosa_PPF-1 
...Ebola_virus_zaire

In each directory add 
....genes.gbk
....sequence.gff
....sequence.fa

Build custom database

$ java -Xmx8g -jar ~/snpEff/snpEff.jar build -genbank -v Ebola_virus_zaire

Annotate

$ java -Xmx8g -jar ~/snpEff/snpEff.jar Ebola_virus_zaire  ~/SRR1735032.filtered.vcf > ~/SRR1735032.output.ann.vcf

Extract only missense_variants

$ extract_missense_variants -i ~/SRR1735032.output.ann.vcf -o ~/SRR1735032.missense.vcf




