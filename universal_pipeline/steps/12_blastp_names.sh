#!/bin/bash
#SBATCH --job-name=nitr_12
#SBATCH --partition=standard
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=nitr_%x_%j.out
#SBATCH --error=nitr_%x_%j.err

set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init
nitr_need_nonzero "$ANN/cluster_genes.tsv" "$ZF_FAA" "$GAR_FAA"

OUTDIR="$FLANK_DIR/blastp"
mkdir -p "$OUTDIR"

python3 "$SCRIPTS/assign_cluster_blast_names.py" \
  --table "$ANN/cluster_genes.tsv" \
  --fasta-out "$OUTDIR/cluster_proteins.fa"

echo "Queries: $(grep -c '^>' "$OUTDIR/cluster_proteins.fa")"

if [[ ! -s "$OUTDIR/zf_prot.pin" && ! -s "$OUTDIR/zf_prot.pdb" ]]; then
  makeblastdb -in "$ZF_FAA" -dbtype prot -title GRCz11 -out "$OUTDIR/zf_prot"
fi
if [[ ! -s "$OUTDIR/gar_prot.pin" && ! -s "$OUTDIR/gar_prot.pdb" ]]; then
  makeblastdb -in "$GAR_FAA" -dbtype prot -title LepOcu1 -out "$OUTDIR/gar_prot"
fi

FMT="6 qseqid sseqid pident length qlen slen qcovs evalue bitscore stitle"

blastp \
  -query "$OUTDIR/cluster_proteins.fa" \
  -db "$OUTDIR/zf_prot" \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -num_threads "$THREADS" \
  -outfmt "$FMT" \
  -out "$OUTDIR/blastp_zebrafish.tsv"

blastp \
  -query "$OUTDIR/cluster_proteins.fa" \
  -db "$OUTDIR/gar_prot" \
  -evalue 1e-5 \
  -max_target_seqs 5 \
  -num_threads "$THREADS" \
  -outfmt "$FMT" \
  -out "$OUTDIR/blastp_gar.tsv"

python3 "$SCRIPTS/assign_cluster_blast_names.py" \
  --table "$ANN/cluster_genes.tsv" \
  --blast-zf "$OUTDIR/blastp_zebrafish.tsv" \
  --blast-gar "$OUTDIR/blastp_gar.tsv" \
  --hits-out "$OUTDIR/cluster_blastp_summary.tsv" \
  -o "$ANN/cluster_genes.tsv"

cp -f "$ANN/nitr_exons.csv" "$OUT/nitr_exons.csv"
cp -f "$ANN/cluster_genes.tsv" "$OUT/cluster_genes.tsv"

echo "Done: $(date)"
cut -f1-6 "$ANN/cluster_genes.tsv" | column -t -s $'\t'
echo
echo "Deliverables:"
echo "  $OUT/nitr_exons.csv"
echo "  $OUT/cluster_genes.tsv"
