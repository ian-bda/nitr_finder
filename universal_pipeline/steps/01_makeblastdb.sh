#!/bin/bash
#SBATCH --job-name=nitr_01
#SBATCH --partition=standard
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=nitr_%x_%j.out
#SBATCH --error=nitr_%x_%j.err

set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init
nitr_need_nonzero "$GENOME"

echo "=== samtools faidx ==="
if [[ ! -s "${GENOME}.fai" ]]; then
  "$SAMTOOLS" faidx "$GENOME"
fi
ls -lh "${GENOME}.fai"

echo "=== makeblastdb ==="
makeblastdb \
  -in "$GENOME" \
  -dbtype nucl \
  -parse_seqids \
  -title "${SPECIES}_genome" \
  -out "$BLASTDB"

echo "Done: $(date)"
blastdbcmd -info -db "$BLASTDB" | head -20
