#!/usr/bin/env python3
# encoding: utf-8

"""
genbank_to_gff3.py

Converts a GenBank (.gbk) file to GFF3 format.

- Requires: Biopython
- Usage: python3 genbank_to_gff3.py -i input.gbk -o output.gff3
"""

import argparse
from Bio import SeqIO

def genbank_to_gff3(input_gbk, output_gff3):
    """
    Converts a GenBank file to GFF3 format.
    """
    with open(output_gff3, "w") as gff:
        gff.write("##gff-version 3\n")  # GFF3 header

        for record in SeqIO.parse(input_gbk, "genbank"):
            seqid = record.id  # Contig name (or accession number)

            for feature in record.features:
                if feature.qualifiers:
                    start = int(feature.location.start) + 1  # Convert to 1-based
                    end = int(feature.location.end)
                    strand = "+" if feature.strand == 1 else "-"

                    # Extract feature type (gene, CDS, etc.)
                    feature_type = feature.type

                    # Extract gene name, locus tag, or product name
                    gene_name = feature.qualifiers.get("gene", 
                                feature.qualifiers.get("locus_tag", 
                                feature.qualifiers.get("product", ["unknown"])))[0]

                    # Debug print to check extracted names
                    print(f"Processing feature: {feature_type}, Location: {start}-{end}, Gene: {gene_name}")

                    # GFF3 fields: seqid, source, type, start, end, score, strand, phase, attributes
                    gff_line = f"{seqid}\tGenBank\t{feature_type}\t{start}\t{end}\t.\t{strand}\t.\tID={gene_name};Name={gene_name}\n"
                    gff.write(gff_line)

    print(f"✅ Conversion complete: {output_gff3}")

def main():
    parser = argparse.ArgumentParser(description="Convert a GenBank file to GFF3 format.")
    parser.add_argument("-i", "--input", required=True, help="Path to input GenBank (.gbk) file.")
    parser.add_argument("-o", "--output", required=True, help="Path to output GFF3 (.gff3) file.")
    args = parser.parse_args()

    genbank_to_gff3(args.input, args.output)

if __name__ == "__main__":
    main()
