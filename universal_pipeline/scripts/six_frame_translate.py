#!/usr/bin/env python3
"""6-frame translate a cluster-window FASTA (same filters as six_frame_universal.pl).

Keeps ORFs >70 aa with <10% X. Writes:
  peptide FASTA (hmmscan input)
  ORF map TSV (window + genomic coords for the parser)
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq

MIN_ORF = 70
MAX_X_FRAC = 0.1
HEADER_COORD_RE = re.compile(r"::(?P<chrom>[^:\s]+):(?P<start>\d+)-(?P<end>\d+)\s*$")


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def load_bed(path: Path) -> list[dict]:
    rows = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            p = line.rstrip("\n").split("\t")
            rows.append(
                {
                    "chrom": p[0],
                    "start": int(p[1]),
                    "end": int(p[2]),
                    "name": p[3] if len(p) > 3 else p[0],
                }
            )
    return rows


def window_meta(record, bed_rows: list[dict] | None, idx: int) -> tuple[str, int, int]:
    """Return (chrom, bed_start_0based, bed_end) for this FASTA record."""
    if bed_rows:
        if len(bed_rows) == 1:
            b = bed_rows[0]
        elif idx < len(bed_rows):
            b = bed_rows[idx]
        else:
            raise SystemExit(f"BED has {len(bed_rows)} intervals but FASTA record {idx}")
        return b["chrom"], b["start"], b["end"]
    m = HEADER_COORD_RE.search(record.description)
    if not m:
        raise SystemExit(
            f"Cannot parse ::chrom:start-end from FASTA header (pass --bed):\n{record.description}"
        )
    return m.group("chrom"), int(m.group("start")), int(m.group("end"))


def emit_orfs(dna: str, frame: int, strand: str) -> list[tuple[int, int, str]]:
    """ORFs as (aa_start_0based_in_frame, aa_end_exclusive, peptide)."""
    offset = (frame - 1) % 3
    frame_dna = dna[offset:]
    frame_dna = frame_dna[: len(frame_dna) - (len(frame_dna) % 3)]
    pep = str(Seq(frame_dna).translate(to_stop=False))
    out = []
    i = 0
    n = len(pep)
    while i < n:
        if pep[i] == "*":
            i += 1
            continue
        j = i
        while j < n and pep[j] != "*":
            j += 1
        seq = pep[i:j]
        if len(seq) > MIN_ORF and (seq.count("X") / len(seq)) < MAX_X_FRAC:
            out.append((i, j, seq))
        i = j
    return out


def translate_window(record, chrom: str, bed_start: int, bed_end: int, orf_i: int):
    dna = str(record.seq).upper()
    win_len = len(dna)
    expected = bed_end - bed_start
    if win_len != expected:
        raise SystemExit(
            f"{record.id}: FASTA length {win_len} != BED span {expected} ({chrom}:{bed_start}-{bed_end})"
        )
    rc = str(Seq(dna).reverse_complement())
    rows = []
    recs = []

    for frame, seq, strand in (
        (1, dna, "+"),
        (2, dna, "+"),
        (3, dna, "+"),
        (4, rc, "-"),
        (5, rc, "-"),
        (6, rc, "-"),
    ):
        offset = (frame - 1) % 3
        for aa0, aa1, pep in emit_orfs(seq, frame, strand):
            orf_i += 1
            dna0 = offset + aa0 * 3
            dna1 = offset + aa1 * 3  # exclusive on this strand's DNA
            if strand == "+":
                win0, win1 = dna0, dna1
            else:
                win0, win1 = win_len - dna1, win_len - dna0
            gstart = bed_start + win0 + 1
            gend = bed_start + win1
            orf_id = f"orf{orf_i:05d}"
            recs.append((orf_id, pep))
            rows.append(
                {
                    "orf_id": orf_id,
                    "chrom": chrom,
                    "strand": strand,
                    "frame": str(frame),
                    "window_start0": str(win0),
                    "window_end": str(win1),
                    "genomic_start": str(gstart),
                    "genomic_end": str(gend),
                    "aa_len": str(len(pep)),
                    "orf_aa": pep,
                }
            )
    return orf_i, recs, rows


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("fasta", type=Path, help="Cluster window FASTA (one seq per window)")
    p.add_argument(
        "-b",
        "--bed",
        type=Path,
        default=None,
        help="Matching BED (0-based) for genomic lift; else parse ::chrom:start-end from header",
    )
    p.add_argument(
        "-o",
        "--prefix",
        type=Path,
        required=True,
        help="Output prefix (writes .fa and .tsv)",
    )
    args = p.parse_args()

    bed_rows = load_bed(args.bed) if args.bed else None
    prefix = args.prefix
    prefix.parent.mkdir(parents=True, exist_ok=True)
    fa_out = Path(f"{prefix}.fa")
    tsv_out = Path(f"{prefix}.tsv")

    orf_i = 0
    all_recs = []
    all_rows = []
    records = list(SeqIO.parse(str(args.fasta), "fasta"))
    if not records:
        raise SystemExit(f"No sequences in {args.fasta}")

    for idx, rec in enumerate(records):
        chrom, bed_start, bed_end = window_meta(rec, bed_rows, idx)
        orf_i, recs, rows = translate_window(rec, chrom, bed_start, bed_end, orf_i)
        all_recs.extend(recs)
        all_rows.extend(rows)

    with fa_out.open("w") as fh:
        for orf_id, pep in all_recs:
            fh.write(f">{orf_id}\n{wrap(pep)}\n")

    cols = [
        "orf_id",
        "chrom",
        "strand",
        "frame",
        "window_start0",
        "window_end",
        "genomic_start",
        "genomic_end",
        "aa_len",
        "orf_aa",
    ]
    with tsv_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(all_rows)

    print(f"Windows: {len(records)}")
    print(f"ORFs >{MIN_ORF} aa: {len(all_recs)}")
    print(f"  {fa_out}")
    print(f"  {tsv_out}")


if __name__ == "__main__":
    main()
