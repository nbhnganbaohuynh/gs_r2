#!/usr/bin/env bash

set -euo pipefail

# set up script path
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

# set up project path
PROJECT_DIR="$HOME/cfDNA_fragmentomics"

FRAGMENT_DIR="$PROJECT_DIR/result/fragments"
FEATURE_DIR="$PROJECT_DIR/result/features"

# analyze fragment length
mkdir -p "$FEATURE_DIR/fragment_length"

python "$SCRIPT_DIR/fragment_length.py"

# analyze end motifs
mkdir -p "$FEATURE_DIR/end_motif"

python "$SCRIPT_DIR/end_motif.py"
