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
  "$MODELS/nitr_full_models.gtf"

ARCH_ARGS=()
if [[ -s "$ANN/nitr_architecture.tsv" ]]; then
  ARCH_ARGS+=(--architecture "$ANN/nitr_architecture.tsv")
fi

python3 "$SCRIPTS/build_nitr_exon_csv.py" \
  --genome "$GENOME" \
  --gtf "$MODELS/nitr_full_models.gtf" \
  "${ARCH_ARGS[@]}" \
  -o "$ANN/nitr_exons.csv"

echo "Done: $(date)"
cut -f1-7 "$ANN/nitr_exons.csv" | column -t -s $'\t'
