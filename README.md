# nitr_finder
# -- Find NITR genes, domain exons, and flanking genes --

**Two tables at the end**

| File | What it is |
|---|---|
| `nitr_exons.csv` | One row per NITR domain exon (SP, V/D1, I/D2, TM, cyto) with genomic start/stop, strand, nucleotide, peptide |
| `cluster_genes.tsv` | NITRs, non-NITRs inside the array, ~6 flanks each side |

```
TBLASTN zebrafish I domains
    → confirm cysteines                         [YOU]
    → TBLASTN self (round 2)
    → confirm cysteines                         [YOU]
    → one cluster window (outermost I ± 100 kb)
    → 6-frame translate + SMART hmmscan         ← Ig exons BLAST missed
    → pair V/I into genes (never two I’s)
    → RNA exons if present (SP→V→I→TM→cyto); else V/I only
    → InterProScan (V vs I)
    → exon table + flanking genes
      (optional SignalP/DeepTMHMM checks only)
```

Canonical architecture: **SP → V → I → TM → cyto**. I-only, V-only, and partial genes are allowed.

---

## Layout

```
universal_pipeline/
  README.md
  config.example.sh
  lib.sh
  requirements.txt
  environment.yml
  queries/zfish_nitr_pub_i_domain.fa
  scripts/                 Python CLIs
  steps/
    00_check_tools.sh      # optional, foreground
    01_makeblastdb.sh
    02_tblastn_round1.sh
    03_inspect_round1.sh   # MANUAL, foreground
    04_tblastn_round2.sh
    05_inspect_round2.sh   # MANUAL, foreground
    06_cluster_window.sh
    07_hmmscan.sh
    08a_rnaseq.sh          # optional
    08_pair_genes.sh
    08b_stitch.sh
    09_interproscan.sh
    09w_web_help.sh        # MANUAL, foreground
    09b_architecture.sh
    10_flanks.sh
    11_nitr_exons.sh
    12_blastp_names.sh
```

---

## Tools

```bash
pip install -r requirements.txt
mamba install -c bioconda blast bedtools samtools hmmer star stringtie miniprot
```

Also: [SMART.hmm](https://smart.embl.de/) (the HMM file, not the website search), [InterProScan](https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html) 5.65 + Java 11.

STAR / StringTie are optional (Illumina). minimap2 is optional (PacBio Iso-Seq).
SignalP / DeepTMHMM are optional **checks**, not how SP/TM exons are found.

Put `SMART.hmm` at `tools/SMART.hmm`, or set `SMART_HMM` in your species config.

---

## Once per species

```bash
cd universal_pipeline          # this folder
cp config/my_species.sh config/your_species.sh
# edit SPECIES, GENOME, GENE_PREFIX, SMART_HMM, RNA_R1/RNA_R2 if you have them
```

```bash
export NITR_CONFIG=$PWD/config/your_species.sh
```

Every `.sh` below reads that config. Leave it exported in your shell for the whole project.

Outputs go to `results/<species>/` unless you set `OUT` in the config.

---

## How to submit (one step at a time)

Compute steps already have `#SBATCH` headers (mem / cpus / time / partition). Edit `--partition=standard` in a step script if your cluster uses a different queue. From this folder:

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/01_makeblastdb.sh
```

Wait until it finishes. Check `nitr_*.out` / `nitr_*.err` (or `results/<species>/logs/` once the job starts writing). Then submit the next script. **Do not** submit 02 before 01 is done.

Manual steps (`03`, `05`, `09w`) are **not** sbatch:

```bash
bash steps/03_inspect_round1.sh
```

To run a compute step on a laptop instead of Slurm:

```bash
bash steps/01_makeblastdb.sh
```

Optional tool check (foreground):

```bash
bash steps/00_check_tools.sh
```

---

## Rules

- Query with **I-domain peptides** and **TBLASTN**, not BLASTP and not full-length NITR proteins.
- Confirm cysteines yourself. The TSV is a triage, not a keep list.
- Cut **one window** from the outermost confirmed I to the other, plus 100 kb.
- **Never** merge two I exons into one gene.
- **Never** translate through an intron. SP/TM/cyto come from spliced RNA exons (Illumina or Iso-Seq). Without RNA, the gene is V and/or I only.
- Do not invent a signal peptide from “nearest Met upstream of V,” and do not search the genome for SP/TM.
- SignalP / DeepTMHMM must not override 8b exons. DeepTMHMM’s N-terminal “signal” is **not** an SP.
- Do not use the target species’ NCBI annotation for flanks. Use zebrafish + spotted gar (miniprot) ± RNA.

---

## Step 1 — BLAST database

`steps/01_makeblastdb.sh` · 16G, 2 cpu, 2 h

```bash
samtools faidx $GENOME

makeblastdb \
  -in $GENOME \
  -dbtype nucl \
  -parse_seqids \
  -title ${SPECIES}_genome \
  -out $OUT/blast_db/ref_db
```

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/01_makeblastdb.sh
```

**Output:** `$GENOME.fai`, `$OUT/blast_db/ref_db.*`

---

## Step 2 — TBLASTN zebrafish I domains vs the genome

`steps/02_tblastn_round1.sh` · 32G, 16 cpu, 12 h

Protein vs DNA. Not BLASTP.

```bash
tblastn \
  -query $I_QUERY \
  -db $OUT/blast_db/ref_db \
  -evalue 1e-5 \
  -num_threads $THREADS \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out $OUT/blast/round1_hits.tsv

python3 scripts/extract_tblastn_hit_seqs.py \
  $OUT/blast/round1_hits.tsv -g $GENOME
```

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/02_tblastn_round1.sh
```

**Output:** `round1_hits.tsv`, `round1_hits.fasta`, `round1_hits_aa.fasta`  
Headers contain `chrom:start-end(strand)` — keep that format.

---

## Step 3 — YOU confirm round-1 I domains  `[MANUAL]`

`steps/03_inspect_round1.sh` · foreground, no sbatch

```bash
python3 scripts/flag_i_hits.py \
  $OUT/blast/round1_hits_aa.fasta \
  -o $OUT/blast/round1_i_inspect.tsv
```

```bash
bash steps/03_inspect_round1.sh
```

| note | meaning |
|---|---|
| `likely_I` | ≥4 cysteines |
| `inspect` | 2–3 (allow ~1 missing) |
| `unlikely` | 0–1 (almost never a real I) |

Open `round1_hits_aa.fasta`. Confirm the motif by eye. Put keeper IDs in `$OUT/blast/round1_keep.txt`, then the script prints:

```bash
python3 scripts/flag_i_hits.py \
  $OUT/blast/round1_hits_aa.fasta \
  -o $OUT/blast/round1_i_inspect.tsv \
  --keep $OUT/blast/round1_keep.txt \
  --keep-fasta $OUT/blast/confirmed_round1_I_domains.fasta
```

Do not strip coordinates from headers. Nearby I’s stay as separate sequences.

---

## Step 4 — TBLASTN round 2 (self)

`steps/04_tblastn_round2.sh` · 32G, 16 cpu, 12 h

Confirmed I peptides vs the same genome.

```bash
tblastn \
  -query $OUT/blast/confirmed_round1_I_domains.fasta \
  -db $OUT/blast_db/ref_db \
  -evalue 1e-5 \
  -num_threads $THREADS \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out $OUT/blast/round2_hits.tsv

python3 scripts/extract_tblastn_hit_seqs.py \
  $OUT/blast/round2_hits.tsv -g $GENOME
```

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/04_tblastn_round2.sh
```

**Output:** `round2_hits.tsv`, `round2_hits.fasta`, `round2_hits_aa.fasta`

---

## Step 5 — YOU confirm round-2 I domains  `[MANUAL]`

`steps/05_inspect_round2.sh` · foreground, no sbatch

Same cysteine check. This is the **final I seed set**.

```bash
bash steps/05_inspect_round2.sh
```

Keepers → `$OUT/blast/round2_keep.txt`, then:

```bash
python3 scripts/flag_i_hits.py \
  $OUT/blast/round2_hits_aa.fasta \
  -o $OUT/blast/round2_i_inspect.tsv \
  --keep $OUT/blast/round2_keep.txt \
  --keep-fasta $OUT/blast/confirmed_I_domains.fasta
```

---

## Step 6 — One cluster window

`steps/06_cluster_window.sh` · 8G, 1 cpu, 1 h

Outermost confirmed I on each scaffold ± **100 kb**.

```bash
python3 scripts/build_locus_windows.py \
  $OUT/blast/confirmed_I_domains.fasta \
  -g ${GENOME}.fai \
  -o $OUT/gene_models/cluster_windows.bed \
  --flank 100000 \
  --mode cluster

bedtools getfasta \
  -fi $GENOME \
  -bed $OUT/gene_models/cluster_windows.bed \
  -name+ \
  -fo $OUT/gene_models/cluster_windows.fasta
```

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/06_cluster_window.sh
```

**Output:** `cluster_windows.bed`, `cluster_windows.fasta`

---

## Step 7 — Six-frame translate + SMART hmmscan

`steps/07_hmmscan.sh` · 16G, 8 cpu, 4 h

Finds **V and I exons TBLASTN missed**. Catalog of domain exons, not genes.

Filters: Ig / V-set / I-set, alignment **>65 aa**, domain i-Evalue **<1e-5**.

```bash
python3 scripts/six_frame_translate.py \
  $OUT/gene_models/cluster_windows.fasta \
  --bed $OUT/gene_models/cluster_windows.bed \
  -o $OUT/ig_domains/cluster_six_frames

hmmscan \
  --cpu $THREADS \
  --domtblout $OUT/ig_domains/cluster_smart.domtblout \
  $SMART_HMM \
  $OUT/ig_domains/cluster_six_frames.fa \
  > $OUT/ig_domains/cluster_smart.out

python3 scripts/parse_ig_hmmscan.py \
  --domtblout $OUT/ig_domains/cluster_smart.domtblout \
  --orf-map $OUT/ig_domains/cluster_six_frames.tsv \
  --orf-fa $OUT/ig_domains/cluster_six_frames.fa \
  --window-fasta $OUT/gene_models/cluster_windows.fasta \
  -o $OUT/ig_domains/ig_domains.csv
```

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/07_hmmscan.sh
```

**Output:** `ig_domains.csv`

---

## Step 8a — RNA-seq (optional)

`steps/08a_rnaseq.sh` · 64G, 16 cpu, 24 h

Skip if you have no FASTQ and no `RNA_ISOSEQ`. The script exits 0. Steps 8 / 8b still run.

- Illumina PE: STAR vs the **full genome**, subset the cluster BAM, StringTie `--conservative`.
- PacBio Iso-Seq FASTA(s): `minimap2 -ax splice:hq`, subset, StringTie `-L`.

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/08a_rnaseq.sh
```

**Output (if RNA):** `rnaseq/cluster_transcripts.gtf`

---

## Step 8 — Pair V and I into genes

`steps/08_pair_genes.sh` · 8G, 1 cpu, 1 h

- One I + nearest same-strand V within 25 kb → one gene
- Unpaired I → I-only
- Unpaired V → V-only
- **Never** two I’s in one gene

If StringTie ran, split a fusion transcript so each exon goes to exactly one gene (nearest same-strand V/I); do not skip it.

```bash
python3 scripts/build_nitr_from_transcripts.py \
  --genome $GENOME \
  --gtf $OUT/rnaseq/cluster_transcripts.gtf \
  --i-domains $OUT/blast/confirmed_I_domains.fasta \
  --ig-csv $OUT/ig_domains/ig_domains.csv \
  --gene-prefix $GENE_PREFIX \
  -o $OUT/gene_models/nitr_loci
```

(Omit `--gtf` if there is no RNA.)

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/08_pair_genes.sh
```

**Output:** `nitr_loci_summary.tsv`  
GeneIDs are `{GENE_PREFIX}1…n` in genomic order.

---

## Step 8b — Transcript exons (or V/I only)

`steps/08b_stitch.sh` · 8G, 1 cpu, 1 h

Same-strand StringTie exons win when they cover that gene’s V/I. Each exon is assigned to one gene: overlap with V/I wins; in the gap between two genes, the nearer Ig edge wins (so an SP between NITR4 and NITR5 belongs to NITR5). Extra 5′ acceptor nt are trimmed so SP–V–I–TM is one ORF. Without RNA, V/I are snapped to GT–AG; SP/TM/cyto are left blank (no Met→stop, no genome SP/TM search).

```bash
python3 scripts/stitch_nitr_genomic.py \
  --genome $GENOME \
  --summary $OUT/gene_models/nitr_loci_summary.tsv \
  --gtf $OUT/rnaseq/cluster_transcripts.gtf \
  -o $OUT/gene_models/nitr_full
```

(Omit `--gtf` if there is no RNA.)

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/08b_stitch.sh
```

**Output:** `nitr_full_proteins.fa`, `nitr_full_cds.fa`, `nitr_full_models.gtf`, `nitr_full_summary.tsv`

---

## Step 9 — InterProScan + SMART

`steps/09_interproscan.sh` · 32G, 8 cpu, 6 h

Labels V vs I. Does **not** label SP/TM.

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/09_interproscan.sh
```

**Output:** `annotation/nitr_proteins_simple.fa`, `nitr_ips.tsv`, `nitr_smart.domtblout`

---

## Step 9w — SignalP / DeepTMHMM `[OPTIONAL WEB CHECK]`

`steps/09w_web_help.sh` · foreground, no sbatch

Skip this if you have no RNA and 8b already left SP/TM blank. Upload `annotation/nitr_proteins_simple.fa` only to see whether the web tools **agree** with 8b.

---

## Step 9b — Architecture

`steps/09b_architecture.sh` · 8G, 1 cpu, 20 min

```bash
python3 scripts/annotate_nitr_architecture.py \
  --proteins $OUT/annotation/nitr_proteins_simple.fa \
  --summary8b $OUT/gene_models/nitr_full_summary.tsv \
  --ips $OUT/annotation/nitr_ips.tsv \
  --domtblout $OUT/annotation/nitr_smart.domtblout \
  --signalp-gff $OUT/annotation/signalp/output.gff3 \
  --signalp-summary $OUT/annotation/signalp/output_protein_type.txt \
  --deeptmhmm-gff $OUT/annotation/deeptmhmm/TMRs.gff3 \
  -o $OUT/annotation
```

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/09b_architecture.sh
```

**Output:** `nitr_architecture.tsv`, `nitr_domain_segments.tsv`

---

## Step 10 — Flanking genes

`steps/10_flanks.sh` · 32G, 8 cpu, 4 h

First NITR to last NITR **± 1 Mb**. miniprot zebrafish + spotted gar. Optional StringTie on that window. Keep ~6 flanks each side. Glance at RS-repeat artifacts next to SON.

Downloads the two proteomes into `related_proteomes/` on the first run.

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/10_flanks.sh
```

**Output:** `annotation/cluster_genes.tsv`

---

## Step 11 — NITR exon CSV

`steps/11_nitr_exons.sh` · 8G, 1 cpu, 20 min

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/11_nitr_exons.sh
```

**Output:** `annotation/nitr_exons.csv`

---

## Step 12 — Name the flanks

`steps/12_blastp_names.sh` · 16G, 8 cpu, 1 h

BLASTP vs zebrafish + gar; rewrite GeneID (RUNX1, DONSON, SON, …).

```bash
sbatch --export=ALL,NITR_CONFIG=$NITR_CONFIG steps/12_blastp_names.sh
```

**Done**

- `$OUT/annotation/nitr_exons.csv`
- `$OUT/annotation/cluster_genes.tsv`


## Python scripts

| Script | Used in |
|---|---|
| `extract_tblastn_hit_seqs.py` | 2, 4 |
| `flag_i_hits.py` | 3, 5 |
| `build_locus_windows.py` | 6 |
| `six_frame_translate.py` | 7 |
| `parse_ig_hmmscan.py` | 7 |
| `build_nitr_from_transcripts.py` | 8 |
| `stitch_nitr_genomic.py` | 8b |
| `annotate_nitr_architecture.py` | 9b |
| `build_flanking_genes.py` | 10 |
| `build_nitr_exon_csv.py` | 11 |
| `assign_cluster_blast_names.py` | 12 |
