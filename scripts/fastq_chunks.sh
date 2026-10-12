#!/usr/bin/env bash

set -euo pipefail

# set up project path
PROJECT_DIR="$HOME/cfDNA_fragmentomics"

FASTQ_DIR="$PROJECT_DIR/data/trimmed_fastq"
CHUNK_DIR="$PROJECT_DIR/data/chunk_trimmed_fastq"

# divide sequencing data into chunks
mkdir -p "$CHUNK_DIR"

for r1 in "$FASTQ_DIR/"*_R1_trimmed.fastq.gz; do
    sample=$(basename "$r1" _R1_trimmed.fastq.gz)
    r2="$FASTQ_DIR/${sample}_R2_trimmed.fastq.gz"

    seqkit split2 \
        -1 "$r1" \
        -2 "$r2" \
        -s 500000 \
        -O "$CHUNK_DIR"
done
