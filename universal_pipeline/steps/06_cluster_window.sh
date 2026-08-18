#!/bin/bash
#SBATCH --job-name=nitr_06
#SBATCH --partition=standard
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=nitr_%x_%j.out
#SBATCH --error=nitr_%x_%j.err

set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init
nitr_need_nonzero "$CONFIRMED_I" "$GENOME"

if [[ ! -s "${GENOME}.fai" ]]; then
  "$SAMTOOLS" faidx "$GENOME"
fi

echo "Confirmed I: $CONFIRMED_I ($(grep -c '^>' "$CONFIRMED_I") seqs)"
echo "Mode: cluster (outermost I ± ${CLUSTER_FLANK_BP} bp)"

python3 "$SCRIPTS/build_locus_windows.py" \
  "$CONFIRMED_I" \
  -g "${GENOME}.fai" \
  -o "$MODELS/cluster_windows.bed" \
  --flank "$CLUSTER_FLANK_BP" \
  --mode cluster

bedtools getfasta \
  -fi "$GENOME" \
  -bed "$MODELS/cluster_windows.bed" \
  -name+ \
  -fo "$MODELS/cluster_windows.fasta"

echo "Done: $(date)"
grep '^>' "$MODELS/cluster_windows.fasta" || true
cat "$MODELS/cluster_windows.bed"
ls -lh "$MODELS/cluster_windows.fasta"
