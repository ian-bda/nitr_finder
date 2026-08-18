#!/bin/bash
#SBATCH --job-name=nitr_08
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
nitr_need_nonzero "$GENOME" "$CONFIRMED_I" "$IGDIR/ig_domains.csv"

GTF_ARGS=()
if [[ -s "$RNASEQ/cluster_transcripts.gtf" ]]; then
  GTF_ARGS+=(--gtf "$RNASEQ/cluster_transcripts.gtf")
  echo "Using StringTie GTF: $RNASEQ/cluster_transcripts.gtf"
else
  echo "No StringTie GTF — pairing V/I from Ig catalog only."
fi

python3 "$SCRIPTS/build_nitr_from_transcripts.py" \
  --genome "$GENOME" \
  "${GTF_ARGS[@]}" \
  --i-domains "$CONFIRMED_I" \
  --ig-csv "$IGDIR/ig_domains.csv" \
  --gene-prefix "$GENE_PREFIX" \
  -o "$MODELS/nitr_loci"

# Keep a copy under the historical 8a name so 08b / 09 notes still find it.
mkdir -p "$RNASEQ"
cp -f "$MODELS/nitr_loci_summary.tsv" "$RNASEQ/nitr_summary.tsv"

echo "Done: $(date)"
column -t -s $'\t' "$MODELS/nitr_loci_summary.tsv" | cut -c1-200
