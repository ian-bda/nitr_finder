#!/usr/bin/env python3
"""Parse hmmscan --domtblout for Ig / NITR-I-like domains.

Same hit filters as smart_hmmer_parser.pl:
  HMM name matches ig|I-set|Ig|V-set|IG
  alignment length > 65 aa
  domain i-Evalue < 1e-5

One row per domain hit (catalog of domain exons, not merged ORFs).
Coords are 1-based inclusive genomic on the cluster window's chromosome.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio.Seq import Seq

IG_HMM_RE = re.compile(r"ig|I-set|Ig|V-set|IG")
# SMART names that are real Ig folds; drop Lig/PIG/LIGAN false positives from the regex
FALSE_IG_RE = re.compile(r"^(Lig|PIG|LIGAN|BRIGHT)", re.I)
MIN_AA = 65
MAX_EVALUE = 1e-5


def is_ig_hmm(name: str) -> bool:
    if FALSE_IG_RE.search(name):
        return False
    return bool(IG_HMM_RE.search(name))


def load_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


def load_orf_map(path: Path) -> dict[str, dict]:
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        return {row["orf_id"]: row for row in r}


def load_windows(fasta: Path) -> dict[str, str]:
    """chrom -> window DNA (genome-forward). One window per chrom here."""
    seqs: dict[str, str] = {}
    for rec_id, seq in load_fasta(fasta).items():
        # bedtools -name+ header: name::chrom:start-end
        m = re.search(r"::([^:]+):(\d+)-(\d+)$", rec_id)
        chrom = m.group(1) if m else rec_id.split("::")[-1].split(":")[0]
        seqs[chrom] = seq.upper()
    # also index by full id
    return seqs


def parse_domtblout(path: Path) -> list[dict]:
    hits = []
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 22:
                continue
            hits.append(
                {
                    "hmm_name": p[0],
                    "hmm_acc": p[1],
                    "orf_id": p[3],
                    "full_evalue": float(p[6]),
                    "full_score": float(p[7]),
                    "domain_n": int(p[9]),
                    "domain_of": int(p[10]),
                    "c_evalue": float(p[11]),
                    "i_evalue": float(p[12]),
                    "score": float(p[13]),
                    "hmm_from": int(p[15]),
                    "hmm_to": int(p[16]),
                    "ali_from": int(p[17]),
                    "ali_to": int(p[18]),
                    "env_from": int(p[19]),
                    "env_to": int(p[20]),
                    "hmm_desc": " ".join(p[22:]) if len(p) > 22 else "",
                }
            )
    return hits


def domain_genomic(orf: dict, ali_from: int, ali_to: int) -> tuple[int, int, int, int]:
    """Map 1-based inclusive AA alignment on the ORF to window/genomic DNA."""
    g0 = int(orf["genomic_start"])
    g1 = int(orf["genomic_end"])
    w0 = int(orf["window_start0"])
    w1 = int(orf["window_end"])
    aa_len = ali_to - ali_from + 1
    nt_len = aa_len * 3
    if orf["strand"] == "+":
        dg0 = g0 + (ali_from - 1) * 3
        dg1 = g0 + ali_to * 3 - 1
        dw0 = w0 + (ali_from - 1) * 3
        dw1 = w0 + ali_to * 3
    else:
        dg1 = g1 - (ali_from - 1) * 3
        dg0 = g1 - ali_to * 3 + 1
        dw1 = w1 - (ali_from - 1) * 3
        dw0 = w1 - ali_to * 3
        nt_len = dg1 - dg0 + 1
        dw1 = dw0 + nt_len
    return dg0, dg1, dw0, dw1


def coding_nt(window_dna: str, strand: str, win0: int, win1: int) -> str:
    nt = window_dna[win0:win1]
    if strand == "-":
        nt = str(Seq(nt).reverse_complement())
    return nt


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--domtblout", type=Path, required=True)
    p.add_argument("--orf-map", type=Path, required=True)
    p.add_argument("--orf-fa", type=Path, required=True)
    p.add_argument("--window-fasta", type=Path, required=True)
    p.add_argument("-o", "--out", type=Path, required=True)
    args = p.parse_args()

    orfs = load_orf_map(args.orf_map)
    peptides = load_fasta(args.orf_fa)
    windows = load_windows(args.window_fasta)

    rows = []
    n_ig = 0
    for hit in parse_domtblout(args.domtblout):
        if not is_ig_hmm(hit["hmm_name"]):
            continue
        n_ig += 1
        ali_len = hit["ali_to"] - hit["ali_from"] + 1
        if ali_len <= MIN_AA or hit["i_evalue"] >= MAX_EVALUE:
            continue
        orf = orfs.get(hit["orf_id"])
        if orf is None:
            raise SystemExit(f"orf {hit['orf_id']} not in {args.orf_map}")
        pep = peptides[hit["orf_id"]]
        ig_aa = pep[hit["ali_from"] - 1 : hit["ali_to"]]
        g0, g1, w0, w1 = domain_genomic(orf, hit["ali_from"], hit["ali_to"])
        chrom = orf["chrom"]
        if chrom not in windows:
            if len(windows) == 1:
                window_dna = next(iter(windows.values()))
            else:
                raise SystemExit(f"no window DNA for {chrom}")
        else:
            window_dna = windows[chrom]
        ig_nt = coding_nt(window_dna, orf["strand"], w0, w1)
        strand_word = "forward" if orf["strand"] == "+" else "reverse"
        rows.append(
            {
                "strand": strand_word,
                "scaffold": chrom,
                "scaf_start": str(w0),
                "scaf_end": str(w1),
                "chr_start": str(g0),
                "chr_end": str(g1),
                "hmm_name": hit["hmm_name"],
                "i_evalue": f"{hit['i_evalue']:.3g}",
                "score": f"{hit['score']:.1f}",
                "ali_from": str(hit["ali_from"]),
                "ali_to": str(hit["ali_to"]),
                "orf_id": hit["orf_id"],
                "frame": orf["frame"],
                "Ig-AA": ig_aa,
                "Ig-nucl": ig_nt,
            }
        )

    rows.sort(key=lambda r: (r["scaffold"], int(r["chr_start"]), r["hmm_name"]))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cols = [
        "strand",
        "scaffold",
        "scaf_start",
        "scaf_end",
        "chr_start",
        "chr_end",
        "hmm_name",
        "i_evalue",
        "score",
        "ali_from",
        "ali_to",
        "orf_id",
        "frame",
        "Ig-AA",
        "Ig-nucl",
    ]
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    fa_out = args.out.with_suffix(".fa")
    with fa_out.open("w") as fh:
        for i, r in enumerate(rows, 1):
            hid = f"{r['scaffold']}:{r['chr_start']}-{r['chr_end']}({'+' if r['strand']=='forward' else '-'})_{r['hmm_name']}_{i}"
            fh.write(f">{hid}\n{r['Ig-AA']}\n")

    print(f"Ig-named HMM rows in domtblout: {n_ig}")
    print(f"Significant Ig domains (len>{MIN_AA}, iE<{MAX_EVALUE}): {len(rows)}")
    print(f"  {args.out}")
    print(f"  {fa_out}")
    if rows:
        by_hmm: dict[str, int] = {}
        for r in rows:
            by_hmm[r["hmm_name"]] = by_hmm.get(r["hmm_name"], 0) + 1
        for k, v in sorted(by_hmm.items(), key=lambda kv: -kv[1]):
            print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
