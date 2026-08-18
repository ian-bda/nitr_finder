#!/bin/bash
# Shared paths and helpers for NITR step scripts.
# Each step: source this file, then call nitr_init.
# Set NITR_CONFIG to an absolute path before running a step.

nitr_init() {
  if [[ -z "${NITR_CONFIG:-}" ]]; then
    echo "ERROR: NITR_CONFIG is not set." >&2
    echo "  export NITR_CONFIG=/absolute/path/to/config/my_species.sh" >&2
    echo "  then:  sbatch --export=ALL,NITR_CONFIG=\$NITR_CONFIG steps/01_makeblastdb.sh" >&2
    echo "  or:    bash steps/03_inspect_round1.sh" >&2
    exit 1
  fi
  if [[ ! -f "$NITR_CONFIG" ]]; then
    echo "ERROR: config not found: $NITR_CONFIG" >&2
    exit 1
  fi
  # shellcheck disable=SC1090
  source "$NITR_CONFIG"

  : "${SPECIES:?Set SPECIES in the species config}"
  : "${GENOME:?Set GENOME in the species config}"

  # This toolkit folder (parent of lib.sh), even if the folder is moved.
  UP="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  SCRIPTS="${SCRIPTS:-$UP/scripts}"
  I_QUERY="${I_QUERY:-$UP/queries/zfish_nitr_pub_i_domain.fa}"
  PEPDIR="${PEPDIR:-$UP/related_proteomes}"
  OUT="${OUT:-$UP/results/$SPECIES}"

  BLASTDB="${BLASTDB:-$OUT/blast_db/ref_db}"
  BLAST="${BLAST:-$OUT/blast}"
  MODELS="${MODELS:-$OUT/gene_models}"
  IGDIR="${IGDIR:-$OUT/ig_domains}"
  RNASEQ="${RNASEQ:-$OUT/rnaseq}"
  ANN="${ANN:-$OUT/annotation}"
  FLANK_DIR="${FLANK_DIR:-$ANN/flanks}"
  LOGDIR="${LOGDIR:-$OUT/logs}"
  SMART_HMM="${SMART_HMM:-$UP/tools/SMART.hmm}"
  SAMTOOLS="${SAMTOOLS:-samtools}"
  MINIPROT="${MINIPROT:-miniprot}"
  CONFIRMED_R1="${CONFIRMED_R1:-$BLAST/confirmed_round1_I_domains.fasta}"
  CONFIRMED_I="${CONFIRMED_I:-$BLAST/confirmed_I_domains.fasta}"
  GENE_PREFIX="${GENE_PREFIX:-NITR}"
  CLUSTER_FLANK_BP="${CLUSTER_FLANK_BP:-100000}"
  FLANK_BP="${FLANK_BP:-1000000}"
  N_FLANK="${N_FLANK:-6}"
  MIN_IDENTITY="${MIN_IDENTITY:-0.30}"
  TBLASTN_EVALUE="${TBLASTN_EVALUE:-1e-5}"
  ZF_FAA="${ZF_FAA:-$PEPDIR/GCF_000002035.6_GRCz11_protein.faa}"
  GAR_FAA="${GAR_FAA:-$PEPDIR/GCF_000242695.1_LepOcu1_protein.faa}"

  mkdir -p "$LOGDIR" "$BLAST" "$OUT/blast_db" "$MODELS" "$IGDIR" "$ANN" "$FLANK_DIR" "$PEPDIR"

  for _p in "${BLAST_BIN:-}" "${BEDTOOLS_BIN:-}" "${CONDA_BIN:-}"; do
    if [[ -n "$_p" && -d "$_p" ]]; then
      export PATH="$_p:$PATH"
    fi
  done
  if [[ -n "${JAVA_HOME:-}" ]]; then
    export JAVA_HOME
    export PATH="$JAVA_HOME/bin:$PATH"
  fi

  THREADS="${SLURM_CPUS_PER_TASK:-${NITR_THREADS:-8}}"

  echo "=============================================="
  echo "Host:    $(hostname)"
  echo "Start:   $(date)"
  echo "Species: $SPECIES"
  echo "Genome:  $GENOME"
  echo "Out:     $OUT"
  echo "Config:  $NITR_CONFIG"
  echo "Threads: $THREADS"
  echo "=============================================="
}

nitr_need() {
  local f
  for f in "$@"; do
    if [[ ! -e "$f" ]]; then
      echo "ERROR: missing $f" >&2
      exit 1
    fi
  done
}

nitr_need_nonzero() {
  local f
  for f in "$@"; do
    if [[ ! -s "$f" ]]; then
      echo "ERROR: missing or empty $f" >&2
      exit 1
    fi
  done
}

nitr_has_rna() {
  [[ -n "${RNA_R1:-}" && -s "${RNA_R1}" && -n "${RNA_R2:-}" && -s "${RNA_R2}" ]]
}
