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
if [[ ! -s "$SIMPLE_FA" ]]; then
  nitr_need_nonzero "$MODELS/nitr_full_proteins.fa"
  mkdir -p "$ANN"
  awk '/^>/ {print ">" substr($1,2); next} {print}' "$MODELS/nitr_full_proteins.fa" > "$SIMPLE_FA"
fi
nitr_need_nonzero \
  "$SIMPLE_FA" \
  "$MODELS/nitr_full_summary.tsv"

IPS_ARGS=()
[[ -s "$ANN/nitr_ips.tsv" ]] && IPS_ARGS+=(--ips "$ANN/nitr_ips.tsv")
[[ -s "$ANN/nitr_smart.domtblout" ]] && IPS_ARGS+=(--domtblout "$ANN/nitr_smart.domtblout")

SUM8A_ARGS=()
if [[ -s "$RNASEQ/nitr_summary.tsv" ]]; then
  SUM8A_ARGS+=(--summary8a "$RNASEQ/nitr_summary.tsv")
fi

SIG_ARGS=()
[[ -s "$ANN/signalp/output.gff3" ]] && SIG_ARGS+=(--signalp-gff "$ANN/signalp/output.gff3")
[[ -s "$ANN/signalp/output_protein_type.txt" ]] && SIG_ARGS+=(--signalp-summary "$ANN/signalp/output_protein_type.txt")
[[ -s "$ANN/deeptmhmm/TMRs.gff3" ]] && SIG_ARGS+=(--deeptmhmm-gff "$ANN/deeptmhmm/TMRs.gff3")
[[ -s "$MODELS/nitr_full_models.gtf" ]] && SIG_ARGS+=(--gtf "$MODELS/nitr_full_models.gtf")

python3 "$SCRIPTS/annotate_nitr_architecture.py" \
  --proteins "$SIMPLE_FA" \
  --summary8b "$MODELS/nitr_full_summary.tsv" \
  "${SUM8A_ARGS[@]}" \
  "${IPS_ARGS[@]}" \
  "${SIG_ARGS[@]}" \
  -o "$ANN"

echo "Done: $(date)"
cut -f1-14 "$ANN/nitr_architecture.tsv" | column -t -s $'\t'
