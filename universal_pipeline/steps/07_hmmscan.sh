#!/bin/bash
#SBATCH --job-name=nitr_07
#SBATCH --partition=standard
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=nitr_%x_%j.out
#SBATCH --error=nitr_%x_%j.err

set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init
nitr_need_nonzero "$MODELS/cluster_windows.fasta" "$MODELS/cluster_windows.bed" "$SMART_HMM"

echo "=== 1) 6-frame translate ==="
python3 "$SCRIPTS/six_frame_translate.py" \
  "$MODELS/cluster_windows.fasta" \
  --bed "$MODELS/cluster_windows.bed" \
  -o "$IGDIR/cluster_six_frames"

nitr_need_nonzero "$IGDIR/cluster_six_frames.fa"
echo "ORFs: $(grep -c '^>' "$IGDIR/cluster_six_frames.fa")"

echo "=== 2) hmmscan SMART ==="
hmmscan \
  --cpu "$THREADS" \
  --domtblout "$IGDIR/cluster_smart.domtblout" \
  "$SMART_HMM" \
  "$IGDIR/cluster_six_frames.fa" \
  > "$IGDIR/cluster_smart.out"

echo "=== 3) parse Ig domains ==="
python3 "$SCRIPTS/parse_ig_hmmscan.py" \
  --domtblout "$IGDIR/cluster_smart.domtblout" \
  --orf-map "$IGDIR/cluster_six_frames.tsv" \
  --orf-fa "$IGDIR/cluster_six_frames.fa" \
  --window-fasta "$MODELS/cluster_windows.fasta" \
  -o "$IGDIR/ig_domains.csv"

echo "Done: $(date)"
ls -lh "$IGDIR"/ig_domains.csv
cut -f1-9 "$IGDIR/ig_domains.csv" | column -t -s $'\t' | head -40
