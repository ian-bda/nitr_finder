#!/bin/bash
#SBATCH --job-name=nitr_10
#SBATCH --partition=standard
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=nitr_%x_%j.out
#SBATCH --error=nitr_%x_%j.err

set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "$PIPE/lib.sh"
nitr_init
nitr_need_nonzero \
  "$GENOME" \
  "$MODELS/nitr_full_models.gtf" \
  "$MODELS/nitr_full_proteins.fa"
if [[ ! -s "${GENOME}.fai" ]]; then
  "$SAMTOOLS" faidx "$GENOME"
fi

echo "=== flank window: first NITR to last NITR ± ${FLANK_BP} bp ==="
python3 - "$MODELS/nitr_full_models.gtf" "${GENOME}.fai" "$FLANK_DIR/flank_window.bed" "$FLANK_BP" <<'PY'
import sys
gtf, fai, out, flank = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
chrom, lo, hi = None, 10**18, 0
with open(gtf) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 5 or p[2] != "exon":
            continue
        chrom = p[0]
        lo = min(lo, int(p[3]))
        hi = max(hi, int(p[4]))
chrlen = None
with open(fai) as fh:
    for line in fh:
        name, n = line.split("\t")[0], int(line.split("\t")[1])
        if name == chrom:
            chrlen = n
            break
if chrom is None or chrlen is None:
    raise SystemExit("could not size window from NITR GTF + fai")
start0 = max(0, lo - 1 - flank)
end1 = min(chrlen, hi + flank)
with open(out, "w") as fh:
    fh.write(f"{chrom}\t{start0}\t{end1}\tflank_window\t0\t+\n")
print(f"NITR span {chrom}:{lo}-{hi}")
print(f"Window    {chrom}:{start0+1}-{end1} ({end1-start0} bp)")
PY

REGION=$(awk '{printf "%s:%d-%d", $1, $2+1, $3}' "$FLANK_DIR/flank_window.bed")
echo "samtools region: $REGION"

bedtools getfasta \
  -fi "$GENOME" \
  -bed "$FLANK_DIR/flank_window.bed" \
  -name+ \
  -fo "$FLANK_DIR/flank_window.fasta"

echo "=== related proteomes (zebrafish + spotted gar) ==="
mkdir -p "$PEPDIR"
if [[ ! -s "$ZF_FAA" ]]; then
  wget -O "${ZF_FAA}.gz" "$ZF_FAA_URL"
  gzip -df "${ZF_FAA}.gz"
fi
if [[ ! -s "$GAR_FAA" ]]; then
  wget -O "${GAR_FAA}.gz" "$GAR_FAA_URL"
  gzip -df "${GAR_FAA}.gz"
fi
COMBINED="$FLANK_DIR/zf_gar_proteins.faa"
cat "$ZF_FAA" "$GAR_FAA" > "$COMBINED"
echo "Proteins: $(grep -c '^>' "$COMBINED")"

echo "=== miniprot (do NOT use NCBI annotation of the target species) ==="
"$MINIPROT" \
  -t "$THREADS" \
  --gff \
  -j 2 \
  -G 200k \
  --outc 0.4 \
  --outn 1 \
  -P FLK \
  "$FLANK_DIR/flank_window.fasta" \
  "$COMBINED" \
  > "$FLANK_DIR/miniprot_zf_gar.gff"

echo "miniprot mRNA: $(awk -F'\t' '$3=="mRNA"' "$FLANK_DIR/miniprot_zf_gar.gff" | wc -l)"

ST_ARGS=()
STAR_BAM="$RNASEQ/star_Aligned.sortedByCoord.out.bam"
if [[ -s "$STAR_BAM" && -x "$STRINGTIE" ]]; then
  echo "=== StringTie on expanded window BAM ==="
  "$SAMTOOLS" view -@ "$THREADS" -b "$STAR_BAM" "$REGION" > "$FLANK_DIR/flank.bam"
  "$SAMTOOLS" index -@ "$THREADS" "$FLANK_DIR/flank.bam"
  echo "Window reads: $("$SAMTOOLS" view -c "$FLANK_DIR/flank.bam")"
  "$STRINGTIE" \
    "$FLANK_DIR/flank.bam" \
    -o "$FLANK_DIR/flank_transcripts.gtf" \
    -p "$THREADS" \
    --conservative \
    -c 1.5 \
    -m 200 \
    -A "$FLANK_DIR/flank_gene_abund.tab"
  ST_ARGS+=(--stringtie-gtf "$FLANK_DIR/flank_transcripts.gtf")
else
  echo "No STAR BAM — flanking genes from miniprot only."
fi

python3 "$SCRIPTS/build_flanking_genes.py" \
  --nitr-gtf "$MODELS/nitr_full_models.gtf" \
  --nitr-proteins "$MODELS/nitr_full_proteins.fa" \
  --genome "$GENOME" \
  --miniprot-gff "$FLANK_DIR/miniprot_zf_gar.gff" \
  --window-bed "$FLANK_DIR/flank_window.bed" \
  --proteins "$COMBINED" \
  "${ST_ARGS[@]}" \
  --n-flank "$N_FLANK" \
  --min-identity "$MIN_IDENTITY" \
  -o "$ANN/cluster_genes.tsv"

echo "Done: $(date)"
cut -f1-7 "$ANN/cluster_genes.tsv" | column -t -s $'\t'
echo
echo "Inspect flanks: drop RS-repeat fragments that sit inside a large gene"
echo "(senegalus: fake SRSF next to SON). Then run step 12 to BLASTP-name them."
