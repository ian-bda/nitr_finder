#!/usr/bin/env python3
"""NITR exon CSV from step 8b GTF: one row per exon, peptide = translation of that NT.

Domain labels (SP, V, I, TM, cyto) come from the cassette / transcript GTF, not
from SignalP or DeepTMHMM. Peptide length matches floor(nt/3). Architecture
TSV is optional (ITIM/RNA notes).
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio.Seq import Seq
from pyfaidx import Fasta

GID_RE = re.compile(r'gene_id "([^"]+)"')
DOM_RE = re.compile(r'domain "([^"]+)"')

DOM_OUT = {
    "SP": "SP",
    "V": "V/D1",
    "I": "I/D2",
    "TM": "TM",
    "cyto": "cyto",
    "C1": "cyto",
    "C2": "cyto",
    "V/D1": "V/D1",
    "I/D2": "I/D2",
}
CANON_ORDER = ["SP", "V/D1", "I/D2", "TM", "cyto"]


def parse_gtf(path: Path) -> dict[str, dict]:
    genes: dict[str, dict] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "exon":
                continue
            chrom, start, end, strand, attrs = p[0], int(p[3]), int(p[4]), p[6], p[8]
            gm = GID_RE.search(attrs)
            if not gm:
                continue
            gid = gm.group(1)
            dm = DOM_RE.search(attrs)
            rec = genes.setdefault(gid, {"chrom": chrom, "strand": strand, "exons": []})
            rec["exons"].append(
                {
                    "start": start,
                    "end": end,
                    "gtf_domain": (dm.group(1) if dm else ""),
                }
            )
    for rec in genes.values():
        rec["exons"].sort(key=lambda x: (x["start"], x["end"]))
    return genes


def fetch_nt(genome: Fasta, chrom: str, start: int, end: int, strand: str) -> str:
    seq = str(genome[chrom][start - 1 : end]).upper()
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())
    return seq


def pep_of(nt: str) -> str:
    nt = nt[: len(nt) - (len(nt) % 3)]
    if not nt:
        return ""
    return str(Seq(nt).translate(to_stop=False)).rstrip("*")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genome", type=Path, required=True)
    p.add_argument("--gtf", type=Path, required=True, help="8b nitr_full_models.gtf")
    p.add_argument("--architecture", type=Path, help="optional nitr_architecture.tsv")
    p.add_argument("-o", "--out", type=Path, required=True)
    args = p.parse_args()

    genome = Fasta(str(args.genome), as_raw=True, sequence_always_upper=True)
    genes = parse_gtf(args.gtf)
    arch = {}
    if args.architecture and args.architecture.exists():
        arch = {r["gene_id"]: r for r in csv.DictReader(args.architecture.open(), delimiter="\t")}

    rows = []
    for gid, rec in genes.items():
        strand, chrom = rec["strand"], rec["chrom"]
        meta = arch.get(gid, {})
        notes = []
        if meta.get("ITIM") == "yes":
            notes.append("ITIM")
        if meta.get("ITSM") == "yes":
            notes.append("ITSM")
        tx = meta.get("RNA") or meta.get("transcript_id") or ""
        if tx and tx not in {"", "-"}:
            notes.append(f"RNA:{tx}")
        gnotes = meta.get("notes") or ""
        if "from_transcript" in gnotes or "RNA:" in gnotes:
            if not any(n.startswith("RNA:") for n in notes):
                notes.append("RNA")
        if "genome_cassette" in gnotes:
            notes.append("genome_cassette")
        if "no_SP" in gnotes:
            notes.append("missing_SP")
        if "no_TM" in gnotes:
            notes.append("missing_TM")
        note = ";".join(notes) if notes else ""
        gstart = min(e["start"] for e in rec["exons"])
        for ex in rec["exons"]:
            raw = ex["gtf_domain"]
            dname = DOM_OUT.get(raw, raw or "exon")
            nt = fetch_nt(genome, chrom, ex["start"], ex["end"], strand)
            pep = pep_of(nt)
            if len(pep) < 6 and dname not in ("SP",):
                if len(pep) < 4:
                    continue
            if dname == "SP" and len(pep) < 8:
                continue
            rows.append(
                {
                    "Scaffold": chrom,
                    "GeneID": gid,
                    "Notes": note,
                    "Exon Start": ex["start"],
                    "Exon Stop": ex["end"],
                    "Orientation": strand,
                    "Protein Domain": dname,
                    "Nucleotide sequence": nt,
                    "Protein sequence": pep,
                    "_gene_start": gstart,
                }
            )

    def sk(r):
        di = CANON_ORDER.index(r["Protein Domain"]) if r["Protein Domain"] in CANON_ORDER else 99
        return (r["Scaffold"], r["_gene_start"], r["Exon Start"], r["Exon Stop"], di)

    rows.sort(key=sk)
    cols = [
        "Scaffold", "GeneID", "Notes",
        "Exon Start", "Exon Stop", "Orientation", "Protein Domain",
        "Nucleotide sequence", "Protein sequence",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"Rows: {len(rows)}  genes: {len({r['GeneID'] for r in rows})}")
    print(f"  {args.out}")
    shown = []
    for r in rows:
        if r["GeneID"] not in shown:
            shown.append(r["GeneID"])
    for gid in shown:
        parts = [x["Protein Domain"] for x in rows if x["GeneID"] == gid]
        print(f"  {gid:12} {', '.join(parts)}")


if __name__ == "__main__":
    main()
