#!/usr/bin/env bash

set -euo pipefail

# set up script path
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# set up project path
PROJECT_DIR="$HOME/cfDNA_fragmentomics"

REFSEQ_DIR="$PROJECT_DIR/data/refseq"
BAM_DIR="$PROJECT_DIR/result/bam"
TMP_DIR="$PROJECT_DIR/result/bam/tmp"

FRAGMENT_DIR="$PROJECT_DIR/result/fragments"

# reconstruct fragments
mkdir -p "$FRAGMENT_DIR" "$TMP_DIR"

for bam in "$BAM_DIR/markup/"*.marked.bam; do
    sample=$(basename "$bam" .marked.bam)
    output="$FRAGMENT_DIR/${sample}.fragments.tsv"

    samtools collate -O -u \
        -T "$TMP_DIR/${sample}.collate" \
        "$bam" \
        | python "$SCRIPT_DIR/fragment_reconstruction.py" \
            - "$REFSEQ_DIR/GRCh38_no_alt.fna" "$output"
done
