#!/bin/bash
#SBATCH --job-name=nitr_04
#SBATCH --partition=standard
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=nitr_%x_%j.out
#SBATCH --error=nitr_%x_%j.err

set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init
nitr_need_nonzero "$CONFIRMED_R1" "$GENOME"

OUT="$BLAST/round2_hits.tsv"
echo "Query: $CONFIRMED_R1 ($(grep -c '^>' "$CONFIRMED_R1") sequences)"

tblastn \
  -query "$CONFIRMED_R1" \
  -db "$BLASTDB" \
  -evalue "$TBLASTN_EVALUE" \
  -num_threads "$THREADS" \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out "$OUT"

python3 "$SCRIPTS/extract_tblastn_hit_seqs.py" "$OUT" -g "$GENOME"

echo "Done: $(date)"
echo "Hits: $(wc -l < "$OUT")"
ls -lh "$OUT" "${OUT%.tsv}.fasta" "${OUT%.tsv}_aa.fasta"
