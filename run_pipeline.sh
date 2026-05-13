#!/bin/bash

set -e 

echo "Starting RNA-seq analysis..."
Rscript scripts/rna_seq_analysis.R

echo "Starting ChIP-seq analysis..."
bash scripts/chip_seq_analysis.sh

echo "All analyses complete."