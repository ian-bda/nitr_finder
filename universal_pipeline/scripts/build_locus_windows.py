#!/usr/bin/env python3
"""Build genomic windows BED from confirmed I-domain FASTA headers.

Modes:
  per-locus (old): one ±flank window per I-domain hit
  cluster (default): one window from outermost I domains on each scaffold ±flank

Headers must contain chrom:start-end(strand); the *last* such block is used
(Round-2 style headers).
"""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

COORD_RE = re.compile(
    r"(?<![A-Za-z0-9.])(?P<chrom>[A-Za-z][A-Za-z0-9.]*):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)"
)


def parse_loci(fasta: Path) -> list[dict]:
    loci = []
    seen = set()
    with fasta.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            matches = list(COORD_RE.finditer(header))
            if not matches:
                raise SystemExit(f"No chrom:start-end(strand) in header: {header}")
            m = matches[-1]
            chrom = m.group("chrom")
            start = int(m.group("start"))
            end = int(m.group("end"))
            strand = m.group("strand")
            if start > end:
                start, end = end, start
            key = (chrom, start, end, strand)
            if key in seen:
                continue
            seen.add(key)
            loci.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "strand": strand,
                    "name": header.split()[0],
                }
            )
    return loci


def load_fai(fai: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with fai.open() as fh:
        for line in fh:
            chrom, length, *_ = line.split("\t")
            sizes[chrom] = int(length)
    return sizes


def write_per_locus_bed(loci: list[dict], sizes: dict[str, int], flank: int, out_bed: Path) -> None:
    out_bed.parent.mkdir(parents=True, exist_ok=True)
    with out_bed.open("w") as out:
        for i, loc in enumerate(loci, 1):
            chrom = loc["chrom"]
            if chrom not in sizes:
                raise SystemExit(f"{chrom} not in fai")
            wstart = max(0, loc["start"] - 1 - flank)
            wend = min(sizes[chrom], loc["end"] + flank)
            name = f"locus{i}|{chrom}:{loc['start']}-{loc['end']}({loc['strand']})"
            out.write(f"{chrom}\t{wstart}\t{wend}\t{name}\t0\t{loc['strand']}\n")
    print(f"Wrote {len(loci)} per-locus windows -> {out_bed} (flank={flank} bp)")


def write_cluster_bed(loci: list[dict], sizes: dict[str, int], flank: int, out_bed: Path) -> None:
    """One window per scaffold: min I-start .. max I-end, plus flank."""
    by_chrom: dict[str, list[dict]] = defaultdict(list)
    for loc in loci:
        by_chrom[loc["chrom"]].append(loc)

    out_bed.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with out_bed.open("w") as out:
        for chrom, group in sorted(by_chrom.items()):
            if chrom not in sizes:
                raise SystemExit(f"{chrom} not in fai")
            group = sorted(group, key=lambda x: x["start"])
            i_min = group[0]["start"]
            i_max = max(g["end"] for g in group)
            up = group[0]
            down = max(group, key=lambda x: x["end"])
            wstart = max(0, i_min - 1 - flank)
            wend = min(sizes[chrom], i_max + flank)
            name = (
                f"cluster_{chrom}|n{len(group)}|"
                f"up={chrom}:{up['start']}-{up['end']}({up['strand']})|"
                f"down={chrom}:{down['start']}-{down['end']}({down['strand']})"
            )
            out.write(f"{chrom}\t{wstart}\t{wend}\t{name}\t0\t.\n")
            n += 1
            span = wend - wstart
            print(
                f"{chrom}: {len(group)} I-domains  "
                f"core={i_min}-{i_max} ({i_max - i_min + 1} bp)  "
                f"window={wstart}-{wend} ({span} bp; flank={flank})"
            )
            print(f"  most-upstream I:   {up['chrom']}:{up['start']}-{up['end']}({up['strand']})")
            print(f"  most-downstream I: {down['chrom']}:{down['start']}-{down['end']}({down['strand']})")
    print(f"Wrote {n} cluster window(s) -> {out_bed}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("confirmed_fa", type=Path, help="Confirmed I-domain protein FASTA")
    p.add_argument("-g", "--genome-fai", type=Path, required=True)
    p.add_argument("-o", "--out-bed", type=Path, required=True)
    p.add_argument("--flank", type=int, default=100000)
    p.add_argument(
        "--mode",
        choices=("cluster", "per-locus"),
        default="cluster",
        help="cluster = outermost I ± flank (default); per-locus = one window per I",
    )
    args = p.parse_args()

    loci = parse_loci(args.confirmed_fa)
    if not loci:
        raise SystemExit("No loci parsed from confirmed FASTA")
    sizes = load_fai(args.genome_fai)

    if args.mode == "cluster":
        write_cluster_bed(loci, sizes, args.flank, args.out_bed)
    else:
        write_per_locus_bed(loci, sizes, args.flank, args.out_bed)


if __name__ == "__main__":
    main()
