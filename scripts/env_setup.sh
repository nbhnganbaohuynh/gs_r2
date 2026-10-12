#!/usr/bin/env bash

set -euo pipefail

# set up bioconda
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# install gdown
if ! command -v gdown &> /dev/null; then
    pip3 install -U gdown
fi

# install fastqc
if ! command -v fastqc &> /dev/null; then
    conda install fastqc -y
fi

# install multiqc
if ! command -v multiqc &> /dev/null; then
    conda install multiqc -y
fi

# install fastp
if ! command -v fastp &> /dev/null; then
    conda install fastp -y
fi

# install seqkit
if ! command -v seqkit &> /dev/null; then
    conda install seqkit -y
fi

# install minibwa/bwa-mem
if ! command -v minibwa &> /dev/null; then
    conda install minibwa -y
fi

if ! command -v bwa &> /dev/null; then
    conda install bwa -y
fi

# install samtools
if ! command -v samtools &> /dev/null; then
    conda install samtools ncurses -y
fi

# install picard
if ! command -v picard &> /dev/null; then
    conda install picard -y
fi

# install pysam
if ! command -v pysam &> /dev/null; then
    conda install pysam -y
fi

# install numpy
if ! command -v nump &> /dev/null; then
    pip3 install -U numpy
fi

# install pandas
if ! command -v pandas &> /dev/null; then
    pip3 install -U pandas
fi

# install matplotlib
if ! command -v matplotlib &> /dev/null; then
    pip3 install -U matplotlib
fi

# install scipy
if ! command -v scipy &> /dev/null; then
    pip3 install -U scipy
fi

# install statsmodels
if ! command -v statsmodels &> /dev/null; then
    pip3 install -U statsmodels
fi

# install statsmodels
if ! command -v scikit_posthocs &> /dev/null; then
    pip3 install -U scikit_posthocs
fi
