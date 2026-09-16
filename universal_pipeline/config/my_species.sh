#!/bin/bash
# Copy this file to config/<your_species>.sh and edit the four lines below.
# shellcheck disable=SC1091
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/config.example.sh"

SPECIES="your_species"
GENOME="/path/to/${SPECIES}_genome.fa"
GENE_PREFIX="NITR"
OUT="$UP/results/${SPECIES}"

# Optional RNA. Leave blank to skip 08a; 8b then reports V/I only.
RNA_R1=""
RNA_R2=""
STAR_READ_CMD="zcat"
# PacBio Iso-Seq: path to one FASTA or a directory of FASTAs
RNA_ISOSEQ=""
