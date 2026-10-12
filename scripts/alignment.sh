#!/usr/bin/env bash

set -euo pipefail

# set up project path
PROJECT_DIR="$HOME/cfDNA_fragmentomics"

REFSEQ_DIR="$PROJECT_DIR/data/refseq"
FASTQ_DIR="$PROJECT_DIR/data/trimmed_fastq"

MINIBWA_DIR="$PROJECT_DIR/result/minibwa_logs"
BAM_DIR="$PROJECT_DIR/result/bam"

# index reference genome (once only)
mkdir -p "$MINIBWA_DIR" "$BAM_DIR/sort" "$BAM_DIR/qc" "$BAM_DIR/markup"

minibwa index -l "$REFSEQ_DIR/GRCh38_no_alt.fna"

for r1 in "$FASTQ_DIR"/*_R1_trimmed.fastq.gz; do
    sample=$(basename "$r1" _R1_trimmed.fastq.gz)
    r2="$FASTQ_DIR/${sample}_R2_trimmed.fastq.gz"

    # alignment
    minibwa map \
        -t 8 \
        "$REFSEQ_DIR/GRCh38_no_alt.fna" \
        "$FASTQ_DIR/${sample}_R1_trimmed.fastq.gz" \
        "$FASTQ_DIR/${sample}_R2_trimmed.fastq.gz" \
        2> "$MINIBWA_DIR/${sample}.minibwa.log" \
    | samtools sort \
        -@ 8 \
        -o "$BAM_DIR/sort/${sample}.sorted.bam" \
        -

    samtools index "$BAM_DIR/sort/${sample}.sorted.bam"

    # alignment qc
    samtools flagstat \
        "$BAM_DIR/sort/${sample}.sorted.bam" \
        > "$BAM_DIR/qc/${sample}.flagstat.txt"

    samtools stats \
        "$BAM_DIR/sort/${sample}.sorted.bam" \
        > "$BAM_DIR/qc/${sample}.stats.txt"

    samtools idxstats \
        "$BAM_DIR/sort/${sample}.sorted.bam" \
        > "$BAM_DIR/qc/${sample}.idxstats.txt"

    samtools coverage \
        "$BAM_DIR/sort/${sample}.sorted.bam" \
        > "$BAM_DIR/qc/${sample}.coverage.txt"

    # mark duplicates
    picard MarkDuplicates \
        I="$BAM_DIR/sort/${sample}.sorted.bam" \
        O="$BAM_DIR/markup/${sample}.marked.bam" \
        M="$BAM_DIR/qc/${sample}.duplication_metrics.txt"

    samtools index "$BAM_DIR/markup/${sample}.marked.bam"
done

multiqc "$BAM_DIR/qc" --outdir "$BAM_DIR/qc"
