#!/bin/bash
#SBATCH --job-name=nitr_08a
#SBATCH --partition=standard
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=nitr_%x_%j.out
#SBATCH --error=nitr_%x_%j.err

set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init
nitr_need_nonzero "$GENOME" "$MODELS/cluster_windows.bed"

if ! nitr_has_rna; then
  echo "No RNA_R1/RNA_R2 in config — skipping STAR/StringTie."
  echo "Step 08 will pair genes from Ig domains only."
  exit 0
fi

nitr_need "$STAR" "$STRINGTIE" "$RNA_R1" "$RNA_R2"
mkdir -p "$RNASEQ"

REGION=$(awk '{print $1":"$2+1"-"$3}' "$MODELS/cluster_windows.bed")
echo "Cluster: $REGION"
echo "Reads: $RNA_R1 + $RNA_R2"

INDEX="$RNASEQ/star_index"
if [[ ! -s "$INDEX/SA" ]]; then
  echo "=== STAR genomeGenerate ==="
  mkdir -p "$INDEX"
  "$STAR" \
    --runMode genomeGenerate \
    --runThreadN "$THREADS" \
    --genomeDir "$INDEX" \
    --genomeFastaFiles "$GENOME" \
    --genomeSAindexNbases "${STAR_SA_INDEX_NBASES:-13}" \
    --limitGenomeGenerateRAM 50000000000
else
  echo "=== STAR index exists, skipping ==="
fi

READ_ARGS=()
if [[ -n "${STAR_READ_CMD:-}" ]]; then
  READ_ARGS+=(--readFilesCommand "$STAR_READ_CMD")
fi

echo "=== STAR align ==="
"$STAR" \
  --runMode alignReads \
  --runThreadN "$THREADS" \
  --genomeDir "$INDEX" \
  --readFilesIn "$RNA_R1" "$RNA_R2" \
  "${READ_ARGS[@]}" \
  --outFileNamePrefix "$RNASEQ/star_" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMstrandField intronMotif \
  --outSAMattributes NH HI AS nM jM jI XS \
  --twopassMode Basic \
  --outFilterMultimapNmax 20 \
  --alignIntronMin 20 \
  --alignIntronMax 20000 \
  --alignMatesGapMax 20000 \
  --limitBAMsortRAM 30000000000

BAM="$RNASEQ/star_Aligned.sortedByCoord.out.bam"
"$SAMTOOLS" index -@ "$THREADS" "$BAM"

echo "=== cluster BAM $REGION ==="
CLUSTER_BAM="$RNASEQ/cluster.bam"
"$SAMTOOLS" view -@ "$THREADS" -b "$BAM" "$REGION" > "$CLUSTER_BAM"
"$SAMTOOLS" index -@ "$THREADS" "$CLUSTER_BAM"
echo "Cluster reads: $("$SAMTOOLS" view -c "$CLUSTER_BAM")"

echo "=== StringTie ==="
"$STRINGTIE" \
  "$CLUSTER_BAM" \
  -o "$RNASEQ/cluster_transcripts.gtf" \
  -p "$THREADS" \
  --conservative \
  -c 1.5 \
  -m 200 \
  -A "$RNASEQ/cluster_gene_abund.tab"

echo "Done: $(date)"
echo "Transcripts: $(awk '$3=="transcript"' "$RNASEQ/cluster_transcripts.gtf" | wc -l)"
