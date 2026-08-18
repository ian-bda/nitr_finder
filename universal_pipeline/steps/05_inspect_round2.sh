#!/bin/bash
# MANUAL checkpoint. Same cysteine check as round 1; this is the final I seed set.
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init

AA="$BLAST/round2_hits_aa.fasta"
nitr_need_nonzero "$AA"

TSV="$BLAST/round2_i_inspect.tsv"
python3 "$SCRIPTS/flag_i_hits.py" "$AA" -o "$TSV"

KEEP="$BLAST/round2_keep.txt"
cat <<EOF

================================================================
MANUAL STEP 05 — confirm round-2 I domains (final seed set)
================================================================
TBLASTN peptides:  $AA
Cysteine scores:   $TSV

Self-blast recovers extra / divergent family members zebrafish missed.
Confirm cysteines again (allow ~1 missing). Then:

   python3 $SCRIPTS/flag_i_hits.py \\
     $AA \\
     -o $TSV \\
     --keep $KEEP \\
     --keep-fasta $CONFIRMED_I

Headers MUST keep chrom:start-end(strand) — step 06 parses them.
Nearby I domains are usually separate genes. Do not merge two I hits
into one sequence.
Then run step 06.
================================================================
EOF

echo "Preview (note / n_cys):"
cut -f1,2,5 "$TSV" | column -t | head -40
echo "Done: $(date)"
