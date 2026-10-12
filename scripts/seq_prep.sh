#!/usr/bin/env bash

set -euo pipefail

# set up project path
PROJECT_DIR="$HOME/cfDNA_fragmentomics"

FASTQ_DIR="$PROJECT_DIR/data/fastq"
FASTQ_URL="https://drive.google.com/drive/folders/16BYCUUCDDjyOKV6b5gw6kMDCTQiaHXYT"

REFSEQ_DIR="$PROJECT_DIR/data/refseq"

# download sequencing data (from Google Drive)
mkdir -p "$FASTQ_DIR"

gdown "$FASTQ_URL" -O "$FASTQ_DIR"
for f in "$FASTQ_DIR"/*.fastq.gz; do
    echo "Checking $f"
    gzip -t "$f" || exit 1
done

# download reference genome and bwa index files
mkdir -p "$REFSEQ_DIR"

wget -P "$REFSEQ_DIR" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai" \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.26_GRCh38/GRCh38_major_release_seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz"

gunzip "$REFSEQ_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
tar -xzf "$REFSEQ_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz" -C "$REFSEQ_DIR"
rm "$REFSEQ_DIR/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bwa_index.tar.gz"

for file in "$REFSEQ_DIR/"GCA_000001405.15_GRCh38_no_alt_analysis_set*; do
    [[ -f "$file" ]] || continue

    base_name=$(basename "$file")
    suffix="${base_name#GCA_000001405.15_GRCh38_no_alt_analysis_set}"
    mv -- "$file" "$REFSEQ_DIR/GRCh38_no_alt$suffix"
done
