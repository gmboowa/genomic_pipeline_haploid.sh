#!/usr/bin/env python3
"""
hostile_clean_paired_human.py - Paired-end cleaning with T2T-HLA removal
"""

import argparse
import subprocess
import re
import zipfile
import shutil
import sys
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Paired-end read cleaning with human/swine host removal',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   
    parser.add_argument('--fasta', required=True,
                        help='NCBI accession (GCF_/GCA_) or path to FASTA')
    parser.add_argument('-i', '--input', required=True,
                        help='File list containing paired FASTQ paths')
    parser.add_argument('-o', '--output', required=True,
                        help='Main output directory')
    parser.add_argument('--index-dir', default='bowtie2_index',
                        help='Directory for Bowtie2 indexes')
    parser.add_argument('--reference-dir', default='references',
                        help='Directory for downloaded references')
    parser.add_argument('--remove-human', action='store_true',
                        help='Include T2T-CHM13v2.0 (human-t2t-hla) reference')
   
    return parser.parse_args()

def validate_dependencies():
    """Check for required executables"""
    for cmd in ['datasets', 'bowtie2', 'hostile']:
        if not shutil.which(cmd):
            raise SystemExit(f"Missing dependency: {cmd}. Install with 'conda install -c bioconda ncbi-datasets-cli bowtie2 hostile'")

def extract_accession(input_str):
    """Extract valid NCBI accession from input string"""
    match = re.search(r'(GC[AF]_\d+)', input_str)
    if not match:
        raise SystemExit(f"Invalid accession format: {input_str}")
    return match.group(1).split('.')[0]

def find_existing_reference(accession, reference_dir):
    """Check for existing reference files"""
    ref_dir = Path(reference_dir)
    patterns = [
        f"{accession}*_genomic.fna",
        f"{accession}*.fna",
        f"{accession}*.fa",
        f"{accession}*.fasta"
    ]
   
    for pattern in patterns:
        matches = list(ref_dir.rglob(pattern))
        if matches:
            return matches[0]
    return None

def get_reference(fasta_input, reference_dir):
    """Get reference path with smart checking"""
    local_path = Path(fasta_input)
    if local_path.exists():
        print(f"Using local reference: {local_path.resolve()}")
        return local_path

    base_accession = extract_accession(fasta_input)
    existing_ref = find_existing_reference(base_accession, reference_dir)
    if existing_ref:
        print(f"Using existing reference: {existing_ref.resolve()}")
        return existing_ref

    print(f"Downloading reference: {base_accession}")
    ref_dir = Path(reference_dir)
    ref_dir.mkdir(parents=True, exist_ok=True)
    zip_path = ref_dir / f"{base_accession}.zip"
   
    try:
        subprocess.run(
            ['datasets', 'download', 'genome', 'accession', base_accession,
             '--include', 'genome', '--filename', str(zip_path)],
            check=True
        )
       
        with zipfile.ZipFile(zip_path, 'r') as zf:
            zf.extractall(ref_dir)
       
        new_ref = find_existing_reference(base_accession, ref_dir)
        if not new_ref:
            raise SystemExit(f"Download failed for {base_accession}")
        print(f"Successfully downloaded reference: {new_ref.resolve()}")
        return new_ref
   
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"Failed to download {base_accession}: {e}")
    finally:
        if zip_path.exists():
            zip_path.unlink()

def build_bowtie_index(fasta_path, index_dir):
    """Build Bowtie2 index if missing"""
    index_dir = Path(index_dir)
    index_dir.mkdir(parents=True, exist_ok=True)
    index_base = index_dir / "reference"
   
    if not any(index_dir.glob("*.bt2")):
        print(f"Building Bowtie2 index at: {index_base.resolve()}")
        subprocess.run(
            ['bowtie2-build', str(fasta_path), str(index_base)],
            check=True
        )
    else:
        print(f"Using existing Bowtie2 index in: {index_dir.resolve()}")
    return index_base

def validate_paired_fastq(fq1, fq2):
    """Validate paired FASTQ files"""
    for path in [fq1, fq2]:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"FASTQ file not found: {p.resolve()}")
        if not p.is_file():
            raise ValueError(f"Path is not a file: {p.resolve()}")
        if p.stat().st_size == 0:
            raise ValueError(f"Empty FASTQ file: {p.resolve()}")
    return Path(fq1).resolve(), Path(fq2).resolve()

def process_samples(input_list, index_base, output_root):
    """Clean paired samples using hostile"""
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    processed = 0
    failed_samples = []
   
    with open(input_list) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            try:
                fq1, fq2 = line.split()
                fastq1, fastq2 = validate_paired_fastq(fq1, fq2)
                sample_id = fastq1.stem.split('_')[0]
                sample_out = output_root / f"{sample_id}_cleaned"
                sample_out.mkdir(parents=True, exist_ok=True)
               
                print(f"\nProcessing sample: {sample_id}")
                print(f"Input FASTQs: {fastq1} | {fastq2}")
                print(f"Output directory: {sample_out.resolve()}")
               
                cmd = [
                    'hostile', 'clean',
                    '--fastq1', str(fastq1),
                    '--fastq2', str(fastq2),
                    '--index', str(index_base),
                    '--aligner', 'bowtie2',
                    '-o', str(sample_out)
                ]
               
                subprocess.run(cmd, check=True)
                processed += 1
                print(f"Successfully processed: {sample_id}")
               
            except Exception as e:
                print(f"❌ Failed to process {line}: {str(e)}")
                failed_samples.append((line, str(e)))
                continue
   
    if failed_samples:
        print("\nFailed samples:")
        for sample, error in failed_samples:
            print(f"- {sample}: {error}")
   
    return processed

def main():
    args = parse_arguments()
    validate_dependencies()
   
    try:
        # Handle primary reference
        ref_path = get_reference(args.fasta, args.reference_dir)

        # Add human reference if requested
        if args.remove_human:
            print("\nIncluding T2T-CHM13v2.0 human reference")
            human_accession = 'GCA_009914755.4'  # T2T-CHM13v2.0
            human_ref = get_reference(human_accession, args.reference_dir)
            combined_path = Path(args.reference_dir) / "combined_host_human.fasta"
            
            # Combine references
            print(f"Combining references: {ref_path.name} + {human_ref.name}")
            with open(combined_path, 'wb') as outfile:
                for source in [ref_path, human_ref]:
                    with open(source, 'rb') as infile:
                        outfile.write(infile.read())
            ref_path = combined_path

        # Build Bowtie2 index
        index_base = build_bowtie_index(ref_path, args.index_dir)
       
        # Process samples
        print(f"\n{' PAIRED-END PROCESSING STARTED ':=^60}")
        processed = process_samples(args.input, index_base, args.output)
        print(f"\nCompleted {processed} paired samples successfully")
       
    except Exception as e:
        print(f"\n❌ Pipeline failed: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()