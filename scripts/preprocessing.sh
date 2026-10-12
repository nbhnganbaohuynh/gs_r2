#!/usr/bin/env bash

set -euo pipefail

# set up project path
PROJECT_DIR="$HOME/cfDNA_fragmentomics"

FASTQ_DIR="$PROJECT_DIR/data/fastq"
TRIMMED_DIR="$PROJECT_DIR/data/trimmed_fastq"

PREQC_DIR="$PROJECT_DIR/result/preprocessing/pre_fastq"
FASTP_DIR="$PROJECT_DIR/result/preprocessing/fastp"
POSTQC_DIR="$PROJECT_DIR/result/preprocessing/post_fastq"

# fastqc (pre-trimming and filtering)
mkdir -p "$PREQC_DIR"

fastqc "$FASTQ_DIR/"*.fastq.gz --outdir "$PREQC_DIR" --threads 8
multiqc "$PREQC_DIR" --outdir "$PREQC_DIR"

# fastp
mkdir -p "$TRIMMED_DIR" "$FASTP_DIR"

for r1 in "$FASTQ_DIR/"*_R1.fastq.gz; do
    sample=$(basename "$r1" _R1.fastq.gz)
    r2="$FASTQ_DIR/${sample}_R2.fastq.gz"

    fastp \
        -i "$r1" \
        -I "$r2" \
        -o "$TRIMMED_DIR/${sample}_R1_trimmed.fastq.gz" \
        -O "$TRIMMED_DIR/${sample}_R2_trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --html "$FASTP_DIR/${sample}_fastp.html" \
        --json "$FASTP_DIR/${sample}_fastp.json" \
        --thread "$(nproc)"
done

# fastqc (post-trimming and filtering)
mkdir -p "$POSTQC_DIR"

fastqc "$TRIMMED_DIR/"*.fastq.gz --outdir "$POSTQC_DIR" --threads 8
multiqc "$POSTQC_DIR" --outdir "$POSTQC_DIR"
