#!/bin/bash
# MANUAL checkpoint — SP and TM come from the web tools, not heuristics.
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init

SIMPLE_FA="$ANN/nitr_proteins_simple.fa"
nitr_need_nonzero "$SIMPLE_FA"

mkdir -p "$ANN/signalp" "$ANN/deeptmhmm"

cat <<EOF

================================================================
MANUAL STEP 09w — SignalP 5.0 and DeepTMHMM
================================================================
Upload this FASTA (one protein per NITR):

  $SIMPLE_FA

1) SignalP 5.0  https://services.healthtech.dtu.dk/services/SignalP-5.0/
   Organism: Eukarya
   Download into $ANN/signalp/
     output.gff3
     output_protein_type.txt

   Only SignalP SP calls are trusted. DeepTMHMM often labels an N-terminal
   "signal" on every protein — do NOT treat that as SP.

2) DeepTMHMM    https://dtu.biolib.com/DeepTMHMM
   Download into $ANN/deeptmhmm/
     TMRs.gff3
     predicted_topologies.3line   (optional, kept for records)

Then run step 09b.

================================================================
EOF

echo "Proteins to upload: $(grep -c '^>' "$SIMPLE_FA")"
grep '^>' "$SIMPLE_FA"
echo "Done: $(date)"
