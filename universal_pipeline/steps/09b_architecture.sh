#!/bin/bash
#SBATCH --job-name=nitr_09b
#SBATCH --partition=standard
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --output=nitr_%x_%j.out
#SBATCH --error=nitr_%x_%j.err

set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init

SIMPLE_FA="$ANN/nitr_proteins_simple.fa"
nitr_need_nonzero \
  "$SIMPLE_FA" \
  "$MODELS/nitr_full_summary.tsv" \
  "$ANN/nitr_ips.tsv" \
  "$ANN/nitr_smart.domtblout" \
  "$ANN/signalp/output.gff3" \
  "$ANN/signalp/output_protein_type.txt" \
  "$ANN/deeptmhmm/TMRs.gff3"

SUM8A_ARGS=()
if [[ -s "$RNASEQ/nitr_summary.tsv" ]]; then
  SUM8A_ARGS+=(--summary8a "$RNASEQ/nitr_summary.tsv")
fi

python3 "$SCRIPTS/annotate_nitr_architecture.py" \
  --proteins "$SIMPLE_FA" \
  --summary8b "$MODELS/nitr_full_summary.tsv" \
  "${SUM8A_ARGS[@]}" \
  --ips "$ANN/nitr_ips.tsv" \
  --domtblout "$ANN/nitr_smart.domtblout" \
  --signalp-gff "$ANN/signalp/output.gff3" \
  --signalp-summary "$ANN/signalp/output_protein_type.txt" \
  --deeptmhmm-gff "$ANN/deeptmhmm/TMRs.gff3" \
  -o "$ANN"

echo "Done: $(date)"
cut -f1-14 "$ANN/nitr_architecture.tsv" | column -t -s $'\t'
