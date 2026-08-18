#!/bin/bash
# Copy this file to config/<your_species>.sh and edit the four lines below.
# shellcheck disable=SC1091
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.example.sh"

SPECIES="your_species"
GENOME="/path/to/${SPECIES}_genome.fa"
GENE_PREFIX="NITR"
OUT="$UP/results/${SPECIES}"

# Optional RNA-seq. Leave blank to skip step 08a.
RNA_R1=""
RNA_R2=""
STAR_READ_CMD="zcat"
