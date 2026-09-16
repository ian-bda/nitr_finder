#!/usr/bin/env python3
"""Extract stranded nt + translated protein sequences from TBLASTN outfmt-6 hits.

Writes:
  <prefix>_with_seq.tsv
  <prefix>.fasta          (nucleotide)
  <prefix>_aa.fasta       (protein)
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from Bio.Seq import Seq
from pyfaidx import Fasta

COLS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def blast_sid_to_chrom(sid: str) -> str:
    if sid.startswith("gb|") and sid.endswith("|"):
        return sid[3:-1]
    parts = sid.split("|")
    if len(parts) >= 2 and parts[0] in {"gb", "ref", "emb", "dbj"}:
        return parts[1]
    return sid


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def extract(hits_tsv: Path, genome_fasta: Path, prefix: Path | None = None) -> Path:
    hits_tsv = hits_tsv.resolve()
    if prefix is None:
        prefix = hits_tsv.with_suffix("")
    else:
        prefix = prefix.resolve()

    out_tsv = Path(f"{prefix}_with_seq.tsv")
    out_nt = Path(f"{prefix}.fasta")
    out_aa = Path(f"{prefix}_aa.fasta")

    genome = Fasta(str(genome_fasta), as_raw=True, sequence_always_upper=True)
    rows: list[dict[str, str]] = []

    with hits_tsv.open() as fh, out_nt.open("w") as nt_fh, out_aa.open("w") as aa_fh:
        for i, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) != 12:
                raise SystemExit(f"bad line {i} in {hits_tsv}: expected 12 columns, got {len(parts)}")
            rec = dict(zip(COLS, parts))
            chrom = blast_sid_to_chrom(rec["sseqid"])
            sstart, send = int(rec["sstart"]), int(rec["send"])
            if sstart <= send:
                strand = "+"
                start, end = sstart, send
            else:
                strand = "-"
                start, end = send, sstart

            if chrom not in genome:
                raise SystemExit(f"chrom {chrom} not in genome (from {rec['sseqid']})")

            nt = str(genome[chrom][start - 1 : end])
            if strand == "-":
                nt = str(Seq(nt).reverse_complement())
            aa = str(Seq(nt).translate(to_stop=False))

            # Keep headers filesystem/BLAST-friendly
            qsafe = re.sub(r"[^\w.:()+-]+", "_", rec["qseqid"])[:80]
            hit_id = f"{qsafe}__{chrom}:{start}-{end}({strand})__hit{i}"

            nt_fh.write(f">{hit_id} evalue={rec['evalue']} pident={rec['pident']}\n{wrap(nt)}\n")
            aa_fh.write(f">{hit_id} evalue={rec['evalue']} pident={rec['pident']}\n{wrap(aa)}\n")

            rec.update(
                {
                    "chrom": chrom,
                    "hit_start": str(start),
                    "hit_end": str(end),
                    "strand": strand,
                    "nt_length": str(len(nt)),
                    "nt_sequence": nt,
                    "aa_sequence": aa,
                    "hit_id": hit_id,
                }
            )
            rows.append(rec)

    out_cols = COLS + [
        "chrom",
        "hit_start",
        "hit_end",
        "strand",
        "nt_length",
        "nt_sequence",
        "aa_sequence",
        "hit_id",
    ]
    with out_tsv.open("w") as out:
        out.write("\t".join(out_cols) + "\n")
        for rec in rows:
            out.write("\t".join(rec[c] for c in out_cols) + "\n")

    print(f"Wrote {len(rows)} hits")
    print(f"  {out_tsv}")
    print(f"  {out_nt}")
    print(f"  {out_aa}")
    return out_aa


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("hits_tsv", type=Path, help="TBLASTN outfmt-6 TSV")
    p.add_argument(
        "-g",
        "--genome",
        type=Path,
        required=True,
        help="Reference genome FASTA (indexed)",
    )
    p.add_argument(
        "-o",
        "--prefix",
        type=Path,
        default=None,
        help="Output prefix (default: hits TSV path without .tsv)",
    )
    args = p.parse_args()
    extract(args.hits_tsv, args.genome, args.prefix)


if __name__ == "__main__":
    main()
