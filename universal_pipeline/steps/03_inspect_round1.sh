#!/bin/bash
# MANUAL checkpoint. Scores cysteine count; you still decide keep/drop.
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init

AA="$BLAST/round1_hits_aa.fasta"
nitr_need_nonzero "$AA"

TSV="$BLAST/round1_i_inspect.tsv"
python3 "$SCRIPTS/flag_i_hits.py" "$AA" -o "$TSV"

KEEP="$BLAST/round1_keep.txt"
cat <<EOF

================================================================
MANUAL STEP 03 — confirm round-1 I domains
================================================================
TBLASTN peptides:  $AA
Cysteine scores:   $TSV
  likely_I  = ≥4 cysteines (typical NITR I)
  inspect   = 2–3 cysteines (allow ~1 missing from the motif)
  unlikely  = 0–1 cysteines (almost never a real I)

1. Open the AA FASTA (or the TSV) and confirm the cysteine motif by eye.
   Canonical NITR I domains usually have several C's; allow ~1 missing.
2. Copy keeper sequence IDs (first token of each header) into:
     $KEEP
   one ID per line.
3. Build the confirmed FASTA:

   python3 $SCRIPTS/flag_i_hits.py \\
     $AA \\
     -o $TSV \\
     --keep $KEEP \\
     --keep-fasta $CONFIRMED_R1

Headers already contain chrom:start-end(strand) — do not strip them.
Then run step 04.
================================================================
EOF

echo "Preview (note / n_cys):"
cut -f1,2,5 "$TSV" | column -t | head -40
echo "Done: $(date)"
