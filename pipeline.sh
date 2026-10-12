#!/usr/bin/env bash

set -euo pipefail

# set up project directory
PROJECT_DIR="$HOME/cfDNA_fragmentomics"
mkdir -p "$PROJECT_DIR"

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/scripts"

# set up conda env
bash "$SCRIPT_DIR/env_setup.sh"

# prepare sequences
bash "$SCRIPT_DIR/seq_prep.sh"

# preprocess sequencing data
bash "$SCRIPT_DIR/preprocessing.sh"

# divide sequencing data into chunks
#bash "$SCRIPT_DIR/fastq_chunks.sh"

# alignment (minibwa)
bash "$SCRIPT_DIR/alignment.sh"

# reconstruct fragments (fragment length and end motif)
bash "$SCRIPT_DIR/fragment_reconstruction.sh"

# perform feature analysis
bash "$SCRIPT_DIR/feature_analysis.sh"
