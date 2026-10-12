#!usr/bin/env bash

set -euo pipefail

# set up project path
PROJECT_DIR="$HOME/cfDNA_fragmentomics"

REFSEQ_DIR="$PROJECT_DIR/data/refseq"
FASTQ_DIR="$PROJECT_DIR/data/trimmed_fastq"
CHUNK_DIR="$PROJECT_DIR/data/chunk_fastq"

BWA_DIR="$PROJECT_DIR/result/bwa_logs"
BAM_DIR="$PROJECT_DIR/result/bam"

mkdir -p "$CHUNK_DIR" "$BWA_DIR" "$BAM_DIR/sort" "$BAM_DIR/qc" "$BAM_DIR/markup"

for r1 in "$FASTQ_DIR/"*_R1_trimmed.fastq.gz; do
    sample=$(basename "$r1" _R1_trimmed.fastq.gz)
    r2="$FASTQ_DIR/${sample}_R2_trimmed.fastq.gz"

    # randomly get a chunk of 500,000 reads
    #seqkit sample2 -n 500000 -2 -s 8 "$r1" -o "$CHUNK_DIR/${sample}_R1_chunk.fastq"
    #seqkit sample2 -n 500000 -2 -s 8 "$r2" -o "$CHUNK_DIR/${sample}_R2_chunk.fastq"

    # alignment
    bwa mem \
        -t 8 \
        "$REFSEQ_DIR/GRCh38_no_alt.fna" \
        "$CHUNK_DIR/${sample}_R1_chunk.fastq" \
        "$CHUNK_DIR/${sample}_R2_chunk.fastq" \
        2> "$BWA_DIR/${sample}.bwamem.log" \
    | samtools sort \
        -@ 8 \
        -o "$BAM_DIR/sort/${sample}_chunk.sorted.bam"

    samtools index "$BAM_DIR/sort/${sample}_chunk.sorted.bam"

    # alignment qc
    samtools flagstat \
        "$BAM_DIR/sort/${sample}_chunk.sorted.bam" \
        > "$BAM_DIR/qc/${sample}_chunk.flagstat.txt"

    samtools stats \
        "$BAM_DIR/sort/${sample}_chunk.sorted.bam" \
        > "$BAM_DIR/qc/${sample}_chunk.stats.txt"

    samtools idxstats \
        "$BAM_DIR/sort/${sample}_chunk.sorted.bam" \
        > "$BAM_DIR/qc/${sample}_chunk.idxstats.txt"

    # mark duplicates
    picard MarkDuplicates \
        I="$BAM_DIR/sort/${sample}_chunk.sorted.bam" \
        O="$BAM_DIR/markup/${sample}_chunk.marked.bam" \
        M="$BAM_DIR/markup/${sample}_chunk.duplication_metrics.txt"

    samtools index "$BAM_DIR/markup/${sample}_chunk.marked.bam"
done
