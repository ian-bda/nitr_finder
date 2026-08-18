#!/bin/bash
#SBATCH --job-name=nitr_11
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

nitr_need_nonzero \
  "$GENOME" \
  "$MODELS/nitr_full_models.gtf" \
  "$MODELS/nitr_full_cds.fa" \
  "$ANN/nitr_proteins_simple.fa" \
  "$ANN/nitr_ips.tsv" \
  "$ANN/signalp/output.gff3" \
  "$ANN/signalp/output_protein_type.txt" \
  "$ANN/deeptmhmm/TMRs.gff3" \
  "$ANN/nitr_architecture.tsv"

python3 "$SCRIPTS/build_nitr_exon_csv.py" \
  --genome "$GENOME" \
  --gtf "$MODELS/nitr_full_models.gtf" \
  --cds "$MODELS/nitr_full_cds.fa" \
  --proteins "$ANN/nitr_proteins_simple.fa" \
  --ips "$ANN/nitr_ips.tsv" \
  --signalp-gff "$ANN/signalp/output.gff3" \
  --signalp-summary "$ANN/signalp/output_protein_type.txt" \
  --deeptmhmm-gff "$ANN/deeptmhmm/TMRs.gff3" \
  --architecture "$ANN/nitr_architecture.tsv" \
  -o "$ANN/nitr_exons.csv"

echo "Done: $(date)"
cut -f1-7 "$ANN/nitr_exons.csv" | column -t -s $'\t'
