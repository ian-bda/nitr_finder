#!/bin/bash
# Shared defaults. Species configs in config/ source this file, then override
# SPECIES, GENOME, GENE_PREFIX, and RNA paths.

UP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- species (required — override in config/<species>.sh) ---
SPECIES="your_species"
GENOME="/path/to/${SPECIES}_genome.fa"
GENE_PREFIX="NITR"                              # GeneIDs become NITR1, NITR2, ...
OUT="$UP/results/${SPECIES}"

# --- queries ---
I_QUERY="$UP/queries/zfish_nitr_pub_i_domain.fa"
SCRIPTS="$UP/scripts"
# Download SMART.hmm from https://smart.embl.de/ (the HMM file, not the web search).
SMART_HMM="${SMART_HMM:-$UP/tools/SMART.hmm}"

# --- optional RNA-seq (leave empty to skip STAR/StringTie) ---
RNA_R1=""
RNA_R2=""
STAR_READ_CMD="zcat"                            # zcat for .gz, cat for uncompressed .fastq

# --- window sizes ---
CLUSTER_FLANK_BP=100000                         # outermost confirmed I ± this (step 06)
FLANK_BP=1000000                                # NITR span ± this for flanking genes (step 10)
N_FLANK=6
MIN_IDENTITY=0.30
TBLASTN_EVALUE="1e-5"
STAR_SA_INDEX_NBASES=13

# --- related proteomes for flanking-gene homology (downloaded on first step 10) ---
PEPDIR="$UP/related_proteomes"
ZF_FAA="$PEPDIR/GCF_000002035.6_GRCz11_protein.faa"
GAR_FAA="$PEPDIR/GCF_000242695.1_LepOcu1_protein.faa"
ZF_FAA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_protein.faa.gz"
GAR_FAA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/242/695/GCF_000242695.1_LepOcu1/GCF_000242695.1_LepOcu1_protein.faa.gz"

# --- tools (leave empty to use PATH after mamba install) ---
# Directories that do not exist are skipped.
BLAST_BIN=""
BEDTOOLS_BIN=""
CONDA_BIN=""
STAR="$(command -v STAR 2>/dev/null || true)"
STRINGTIE="$(command -v stringtie 2>/dev/null || true)"
INTERPROSCAN="$(command -v interproscan.sh 2>/dev/null || true)"
JAVA_HOME="${JAVA_HOME:-}"
MINIPROT="${MINIPROT:-miniprot}"
SAMTOOLS="${SAMTOOLS:-samtools}"

# --- this cluster (Yoder lab) — uncomment if those paths exist ---
# BLAST_BIN="/opt/blast+/2.14.0/bin"
# BEDTOOLS_BIN="/opt/bedtools/2.30.0/bin"
# CONDA_BIN="/home5/ibirchl/miniforge3/bin"
# STAR="/home5/ibirchl/Bioinformatics_tools/STAR-2.7.11b/bin/Linux_x86_64_static/STAR"
# STRINGTIE="/home5/ibirchl/Yoder_Lab/Polypterus_NITR_Project/tools/stringtie"
# INTERPROSCAN="/home5/ibirchl/Bioinformatics_tools/interproscan-5.65-97.0/interproscan.sh"
# JAVA_HOME="/usr/lib/jvm/java-11-openjdk-amd64"
# SMART_HMM="/home5/ibirchl/HMMER/SMART.hmm"

NITR_THREADS="${NITR_THREADS:-8}"
