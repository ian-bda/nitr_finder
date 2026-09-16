#!/bin/bash
# Optional checkpoint — SignalP / DeepTMHMM are checks, not exon finders.
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init

SIMPLE_FA="$ANN/nitr_proteins_simple.fa"
if [[ ! -s "$SIMPLE_FA" && -s "$MODELS/nitr_full_proteins.fa" ]]; then
  SIMPLE_FA="$MODELS/nitr_full_proteins.fa"
fi
nitr_need_nonzero "$SIMPLE_FA"

mkdir -p "$ANN/signalp" "$ANN/deeptmhmm"

cat <<EOF

================================================================
OPTIONAL STEP 09w — SignalP 5.0 and DeepTMHMM (checks only)
================================================================
SP / TM / cyto exons are already called in step 8b from RNA or a GT–AG
cassette. Do NOT use these web tools to invent exons.

You may upload this FASTA to see whether SignalP/DeepTMHMM agree:

  $SIMPLE_FA

1) SignalP 5.0  https://services.healthtech.dtu.dk/services/SignalP-5.0/
   Organism: Eukarya
   Download into $ANN/signalp/
     output.gff3
     output_protein_type.txt

2) DeepTMHMM    https://dtu.biolib.com/DeepTMHMM
   Download into $ANN/deeptmhmm/
     TMRs.gff3

DeepTMHMM's N-terminal "signal" is not an SP. Skip this step if you have
no RNA and the cassette already left SP/TM blank — that is intentional.

Then run step 09b (works with or without these files).

================================================================
EOF

echo "Proteins: $(grep -c '^>' "$SIMPLE_FA")"
grep '^>' "$SIMPLE_FA"
echo "Done: $(date)"
