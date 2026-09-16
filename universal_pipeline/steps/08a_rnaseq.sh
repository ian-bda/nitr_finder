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
  echo "No RNA_R1/RNA_R2 or RNA_ISOSEQ in config — skipping mapping/StringTie."
  echo "Step 08 will pair genes from Ig domains only. Step 8b will keep V/I"
  echo "(GT–AG snap only). SP/TM/cyto require RNA exons. No invented Met→stop."
  exit 0
fi

mkdir -p "$RNASEQ"
REGION=$(awk '{print $1":"$2+1"-"$3}' "$MODELS/cluster_windows.bed")
echo "Cluster: $REGION"
CLUSTER_BAM="$RNASEQ/cluster.bam"
MINIMAP2="${MINIMAP2:-minimap2}"

if nitr_has_isoseq; then
  nitr_need "$MINIMAP2" "$STRINGTIE"
  echo "=== PacBio Iso-Seq (minimap2 splice:hq + StringTie -L) ==="
  ISO_FA=()
  if [[ -d "$RNA_ISOSEQ" ]]; then
    while IFS= read -r fa; do
      [[ -n "$fa" ]] && ISO_FA+=("$fa")
    done < <(find "$RNA_ISOSEQ" -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fa.gz' -o -name '*.fasta.gz' \) | LC_ALL=C sort)
  else
    ISO_FA+=("$RNA_ISOSEQ")
  fi
  if [[ ${#ISO_FA[@]} -eq 0 ]]; then
    echo "ERROR: RNA_ISOSEQ set but no FASTA found: $RNA_ISOSEQ" >&2
    exit 1
  fi
  BAMS=()
  i=0
  for fa in "${ISO_FA[@]}"; do
    i=$((i + 1))
    stem="$RNASEQ/isoseq_$(printf '%02d' "$i")"
    echo "Map $fa"
    if [[ "$fa" == *.gz ]]; then
      gzip -dc "$fa" | "$MINIMAP2" -ax splice:hq -uf -t "$THREADS" --secondary=no "$GENOME" - \
        | "$SAMTOOLS" sort -@ "$THREADS" -o "${stem}.bam"
    else
      "$MINIMAP2" -ax splice:hq -uf -t "$THREADS" --secondary=no "$GENOME" "$fa" \
        | "$SAMTOOLS" sort -@ "$THREADS" -o "${stem}.bam"
    fi
    "$SAMTOOLS" index -@ "$THREADS" "${stem}.bam"
    BAMS+=("${stem}.bam")
  done
  MERGED="$RNASEQ/isoseq_all.bam"
  if [[ ${#BAMS[@]} -eq 1 ]]; then
    cp "${BAMS[0]}" "$MERGED"
    cp "${BAMS[0]}.bai" "${MERGED}.bai"
  else
    "$SAMTOOLS" merge -f -@ "$THREADS" "$MERGED" "${BAMS[@]}"
    "$SAMTOOLS" index -@ "$THREADS" "$MERGED"
  fi
  "$SAMTOOLS" view -@ "$THREADS" -b "$MERGED" "$REGION" > "$CLUSTER_BAM"
  "$SAMTOOLS" index -@ "$THREADS" "$CLUSTER_BAM"
  echo "Cluster reads: $("$SAMTOOLS" view -c "$CLUSTER_BAM")"
  echo "=== StringTie -L ==="
  "$STRINGTIE" \
    "$CLUSTER_BAM" \
    -L \
    -o "$RNASEQ/cluster_transcripts.gtf" \
    -p "$THREADS" \
    -c 1 \
    -m 200 \
    -A "$RNASEQ/cluster_gene_abund.tab"
else
  nitr_need "$STAR" "$STRINGTIE" "$RNA_R1" "$RNA_R2"
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
fi

echo "Done: $(date)"
echo "Transcripts: $(awk '$3=="transcript"' "$RNASEQ/cluster_transcripts.gtf" | wc -l)"
