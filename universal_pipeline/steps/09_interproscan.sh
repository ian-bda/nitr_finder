#!/bin/bash
#SBATCH --job-name=nitr_09
#SBATCH --partition=standard
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --output=nitr_%x_%j.out
#SBATCH --error=nitr_%x_%j.err

set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init
nitr_need_nonzero "$MODELS/nitr_full_proteins.fa" "$MODELS/nitr_full_summary.tsv"
nitr_need "$INTERPROSCAN" "$SMART_HMM"

mkdir -p "$ANN/ips_tmp"

# InterProScan rejects spaces in FASTA headers
SIMPLE_FA="$ANN/nitr_proteins_simple.fa"
awk '/^>/ {print ">" substr($1,2); next} {print}' "$MODELS/nitr_full_proteins.fa" > "$SIMPLE_FA"
echo "Proteins: $(grep -c '^>' "$SIMPLE_FA")"

echo "=== InterProScan (skip Gene3D — broken HMM index on some installs) ==="
set +e
"$INTERPROSCAN" \
  -i "$SIMPLE_FA" \
  -b "$ANN/nitr_ips" \
  -f TSV,GFF3 \
  -dp \
  --cpu "$THREADS" \
  -T "$ANN/ips_tmp" \
  -appl Pfam,SMART,ProSiteProfiles,ProSitePatterns,HAMAP,NCBIfam,AntiFam
IPS_RC=$?
set -e
if [[ $IPS_RC -ne 0 ]]; then
  echo "WARNING: InterProScan exited $IPS_RC — continuing with SMART hmmscan"
fi

echo "=== SMART hmmscan ==="
hmmscan \
  --cpu "$THREADS" \
  --domtblout "$ANN/nitr_smart.domtblout" \
  "$SMART_HMM" \
  "$SIMPLE_FA" \
  > "$ANN/nitr_smart.out"

echo "Done: $(date)"
echo "Next: run step 09w (SignalP 5.0 + DeepTMHMM web), then 09b."
ls -lh "$ANN"/nitr_ips.tsv "$ANN"/nitr_smart.domtblout "$SIMPLE_FA" 2>/dev/null || true
