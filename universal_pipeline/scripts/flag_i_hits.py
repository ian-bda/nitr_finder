#!/usr/bin/env python3
"""Score TBLASTN I-domain hit peptides for manual confirmation.

This does **not** replace eyes on the alignment. It flags cysteine count and
writes a TSV you can sort, then a keeper FASTA from a list of IDs.

Canonical NITR I domains usually have several cysteines (often ~4–6). Allow
~1 missing. Hits with 0–1 C are almost never I domains.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def load_fa(path: Path) -> list[tuple[str, str, str]]:
    recs, name, hdr, chunks = [], None, "", []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, hdr, "".join(chunks)))
                hdr = line[1:].strip()
                name = hdr.split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            recs.append((name, hdr, "".join(chunks)))
    return recs


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("aa_fasta", type=Path, help="roundN_hits_aa.fasta")
    p.add_argument("-o", "--out-tsv", type=Path, required=True)
    p.add_argument("--keep", type=Path, help="text file of sequence IDs to keep")
    p.add_argument("--keep-fasta", type=Path, help="write keeper FASTA (needs --keep)")
    args = p.parse_args()

    recs = load_fa(args.aa_fasta)
    args.out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_tsv.open("w") as fh:
        fh.write("id\tn_cys\tcys_positions\taa_len\tnote\theader\n")
        for name, hdr, aa in recs:
            pos = [str(i + 1) for i, x in enumerate(aa) if x == "C"]
            n = len(pos)
            if n >= 4:
                note = "likely_I"
            elif n >= 2:
                note = "inspect"
            else:
                note = "unlikely"
            fh.write(f"{name}\t{n}\t{','.join(pos)}\t{len(aa)}\t{note}\t{hdr}\n")
    print(f"Wrote {args.out_tsv} ({len(recs)} hits)")

    if args.keep and args.keep_fasta:
        keep = {ln.strip() for ln in args.keep.read_text().splitlines() if ln.strip() and not ln.startswith("#")}
        n = 0
        with args.keep_fasta.open("w") as fh:
            for name, hdr, aa in recs:
                if name in keep:
                    fh.write(f">{hdr}\n")
                    for i in range(0, len(aa), 80):
                        fh.write(aa[i : i + 80] + "\n")
                    n += 1
        print(f"Wrote {args.keep_fasta} ({n} kept)")


if __name__ == "__main__":
    main()
