#!/bin/bash
#SBATCH --job-name=nitr_08b
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
nitr_need_nonzero "$GENOME" "$MODELS/nitr_loci_summary.tsv"

GTF_ARGS=()
if [[ -s "$RNASEQ/cluster_transcripts.gtf" ]]; then
  GTF_ARGS+=(--gtf "$RNASEQ/cluster_transcripts.gtf")
fi

python3 "$SCRIPTS/stitch_nitr_genomic.py" \
  --genome "$GENOME" \
  --summary "$MODELS/nitr_loci_summary.tsv" \
  "${GTF_ARGS[@]}" \
  -o "$MODELS/nitr_full"

echo "Done: $(date)"
ls -lh "$MODELS"/nitr_full*
cut -f1-12 "$MODELS/nitr_full_summary.tsv" | column -t -s $'\t'
