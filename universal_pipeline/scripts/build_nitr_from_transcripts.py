#!/usr/bin/env python3
"""Reconstruct NITR proteins from Ig-domain evidence + StringTie transcripts.

Gene count is NOT fixed. Every I and V Ig locus (HMM, plus BLAST I's HMM missed)
is a domain exon. Those are grouped into genes:

  V + I     canonical membrane NITR (then look for TM/cyto in the ORF)
  V only    possible secreted / excreted NITR
  I only    partial / I-only gene

Never merge two I domains or two V domains into one gene.
Transcripts that overlap two genes are treated as over-stitched and skipped.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio.Seq import Seq
from pyfaidx import Fasta

COORD_RE = re.compile(
    r"(?<![A-Za-z0-9.])(?P<chrom>[A-Za-z][A-Za-z0-9.]*):(?P<start>\d+)-(?P<end>\d+)\((?P<strand>[+-])\)"
)
TID_RE = re.compile(r'transcript_id "([^"]+)"')
GID_RE = re.compile(r'gene_id "([^"]+)"')
COV_RE = re.compile(r'(?:cov|coverage) "([^"]+)"')

# max genomic span between paired V and I on the same gene
MAX_VI_BP = 25000

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def overlaps(a0: int, a1: int, b0: int, b1: int) -> bool:
    return max(a0, b0) <= min(a1, b1)


def parse_blast_i(fasta: Path) -> list[dict]:
    loci, seen = [], set()
    with fasta.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            matches = list(COORD_RE.finditer(header))
            if not matches:
                continue
            m = matches[-1]
            start, end = int(m.group("start")), int(m.group("end"))
            if start > end:
                start, end = end, start
            key = (m.group("chrom"), start, end, m.group("strand"))
            if key in seen:
                continue
            seen.add(key)
            loci.append(
                {
                    "chrom": m.group("chrom"),
                    "start": start,
                    "end": end,
                    "strand": m.group("strand"),
                    "cls": "I",
                    "source": "BLAST_I",
                    "hmms": "BLAST",
                    "nC": None,
                }
            )
    return loci


def parse_hmm_loci(ig_csv: Path) -> list[dict]:
    rows = list(csv.DictReader(ig_csv.open(), delimiter="\t"))
    merged: list[dict] = []
    for h in sorted(rows, key=lambda r: int(r["chr_start"])):
        s, e = int(h["chr_start"]), int(h["chr_end"])
        strand = "+" if h["strand"] == "forward" else "-"
        placed = False
        for loc in merged:
            if loc["chrom"] == h["scaffold"] and loc["strand"] == strand and overlaps(s, e, loc["start"], loc["end"]):
                loc["start"] = min(loc["start"], s)
                loc["end"] = max(loc["end"], e)
                loc["hmms"].add(h["hmm_name"])
                loc["nC"] = max(loc["nC"], h["Ig-AA"].count("C"))
                placed = True
                break
        if not placed:
            merged.append(
                {
                    "chrom": h["scaffold"],
                    "start": s,
                    "end": e,
                    "strand": strand,
                    "hmms": {h["hmm_name"]},
                    "nC": h["Ig-AA"].count("C"),
                    "source": "HMM",
                }
            )
    for loc in merged:
        hmms = loc["hmms"]
        nC = loc["nC"]
        if "IGv" in hmms or nC <= 2:
            loc["cls"] = "V"
        elif nC >= 4 and "IGv" not in hmms:
            loc["cls"] = "I"
        elif nC == 3 and "IGv" not in hmms:
            loc["cls"] = "I"  # divergent I-like; may be revised
        else:
            loc["cls"] = "V"
        loc["hmms"] = ",".join(sorted(hmms))
    return merged


def catalog_domains(ig_csv: Path, blast_fa: Path) -> list[dict]:
    hmm = parse_hmm_loci(ig_csv)
    blast = parse_blast_i(blast_fa)
    out = list(hmm)
    for b in blast:
        hit = next(
            (
                h
                for h in hmm
                if h["chrom"] == b["chrom"] and overlaps(h["start"], h["end"], b["start"], b["end"])
            ),
            None,
        )
        if hit:
            hit["cls"] = "I"
            hit["source"] = "HMM+BLAST_I"
        else:
            out.append(b)
    out.sort(key=lambda x: (x["chrom"], x["start"], x["strand"]))
    for i, loc in enumerate(out, 1):
        loc["dom_id"] = f"Ig{i}"
    return out


def pair_genes(domains: list[dict], prefix: str = "NITR") -> list[dict]:
    """Group I/V exons into genes. Count is data-driven."""
    genes: list[dict] = []
    by_key: dict[tuple[str, str], list[dict]] = {}
    for d in domains:
        by_key.setdefault((d["chrom"], d["strand"]), []).append(d)

    for (chrom, strand), locs in by_key.items():
        locs = sorted(locs, key=lambda x: x["start"])
        Is = [d for d in locs if d["cls"] == "I"]
        Vs = [d for d in locs if d["cls"] == "V"]
        used_v: set[str] = set()

        def vid(v: dict) -> str:
            return v["dom_id"]

        if strand == "-":
            # transcription 5'→3' is high→low coords; V is genomic-right of I
            for i, I in enumerate(Is):
                next_i = Is[i + 1]["start"] if i + 1 < len(Is) else I["end"] + MAX_VI_BP
                limit = min(next_i, I["end"] + MAX_VI_BP)
                cands = [
                    v
                    for v in Vs
                    if vid(v) not in used_v and I["end"] < v["start"] <= limit
                ]
                if cands:
                    v = min(cands, key=lambda x: x["start"])
                    used_v.add(vid(v))
                    genes.append(_gene(chrom, strand, I, v, "V+I"))
                else:
                    genes.append(_gene(chrom, strand, I, None, "I-only"))
        else:
            # plus: V is genomic-left of I
            for i, I in enumerate(Is):
                prev_end = Is[i - 1]["end"] if i else I["start"] - MAX_VI_BP
                limit = max(prev_end, I["start"] - MAX_VI_BP)
                cands = [
                    v
                    for v in Vs
                    if vid(v) not in used_v and limit <= v["end"] < I["start"]
                ]
                if cands:
                    v = max(cands, key=lambda x: x["end"])
                    used_v.add(vid(v))
                    genes.append(_gene(chrom, strand, I, v, "V+I"))
                else:
                    genes.append(_gene(chrom, strand, I, None, "I-only"))

        for v in Vs:
            if vid(v) not in used_v:
                genes.append(_gene(chrom, strand, None, v, "V-only"))

    genes.sort(key=lambda g: (g["chrom"], g["start"]))
    for i, g in enumerate(genes, 1):
        g["gene_id"] = f"{prefix}{i}"
    return genes


def _gene(chrom: str, strand: str, I: dict | None, V: dict | None, arch: str) -> dict:
    parts = [d for d in (I, V) if d]
    return {
        "chrom": chrom,
        "strand": strand,
        "arch": arch,
        "I": I,
        "V": V,
        "start": min(d["start"] for d in parts),
        "end": max(d["end"] for d in parts),
        "domains": parts,
    }


def parse_gtf(gtf: Path) -> dict[str, dict]:
    tx: dict[str, dict] = {}
    with gtf.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            chrom, _src, feat, start, end, _sc, strand, _fr, attrs = p
            start, end = int(start), int(end)
            tm = TID_RE.search(attrs)
            gm = GID_RE.search(attrs)
            if not tm:
                continue
            tid = tm.group(1)
            rec = tx.setdefault(
                tid,
                {
                    "tid": tid,
                    "gid": gm.group(1) if gm else tid,
                    "chrom": chrom,
                    "strand": strand,
                    "exons": [],
                    "cov": 0.0,
                },
            )
            if feat == "exon":
                rec["exons"].append((start, end))
            if feat == "transcript":
                cm = COV_RE.search(attrs)
                if cm:
                    rec["cov"] = float(cm.group(1))
    for rec in tx.values():
        rec["exons"].sort()
        rec["start"] = rec["exons"][0][0] if rec["exons"] else 0
        rec["end"] = rec["exons"][-1][1] if rec["exons"] else 0
    return tx


def tx_overlaps_interval(rec: dict, chrom: str, start: int, end: int) -> bool:
    if rec["chrom"] != chrom:
        return False
    return any(overlaps(a, b, start, end) for a, b in rec["exons"])


def tx_hits_gene(rec: dict, gene: dict) -> bool:
    if rec["chrom"] != gene["chrom"] or rec["strand"] != gene["strand"]:
        return False
    return any(tx_overlaps_interval(rec, d["chrom"], d["start"], d["end"]) for d in gene["domains"])


def n_genes_hit(rec: dict, genes: list[dict]) -> int:
    return sum(1 for g in genes if tx_hits_gene(rec, g))


def spliced_seq(genome: Fasta, rec: dict) -> str:
    parts = [str(genome[rec["chrom"]][s - 1 : e]) for s, e in rec["exons"]]
    seq = "".join(parts).upper()
    if rec["strand"] == "-":
        seq = str(Seq(seq).reverse_complement())
    return seq


def genomic_to_tx(rec: dict, gstart: int, gend: int) -> tuple[int, int] | None:
    order = rec["exons"] if rec["strand"] == "+" else list(reversed(rec["exons"]))
    tpos = 0
    t0 = t1 = None
    for s, e in order:
        elen = e - s + 1
        if overlaps(s, e, gstart, gend):
            if rec["strand"] == "+":
                local0 = max(0, gstart - s)
                local1 = min(elen, gend - s + 1)
            else:
                local0 = max(0, e - gend)
                local1 = min(elen, e - gstart + 1)
            if t0 is None:
                t0 = tpos + local0
            t1 = tpos + local1
        tpos += elen
    if t0 is None:
        return None
    return t0, t1


def find_orfs(nt: str) -> list[dict]:
    out = []
    for frame in range(3):
        aa = str(Seq(nt[frame:]).translate(to_stop=False))
        i = 0
        while i < len(aa):
            if aa[i] == "*":
                i += 1
                continue
            j = i
            while j < len(aa) and aa[j] != "*":
                j += 1
            pep = aa[i:j]
            if len(pep) >= 60:
                nt0 = frame + i * 3
                nt1 = frame + j * 3
                has_stop = j < len(aa) and aa[j] == "*"
                mets = [k for k, x in enumerate(pep) if x == "M"]
                out.append(
                    {
                        "aa": pep,
                        "nt0": nt0,
                        "nt1": nt1,
                        "has_stop": has_stop,
                        "has_M": bool(mets),
                        "frame": frame,
                    }
                )
                if mets and mets[0] > 0:
                    m = mets[0]
                    out.append(
                        {
                            "aa": pep[m:],
                            "nt0": frame + (i + m) * 3,
                            "nt1": nt1,
                            "has_stop": has_stop,
                            "has_M": True,
                            "frame": frame,
                        }
                    )
            i = j
    return out


def find_tm(aa: str, after: int = 0) -> tuple[int, int] | None:
    window, thresh = 19, 1.6
    best = None
    start = max(0, after)
    for i in range(start, max(start, len(aa) - window + 1)):
        w = aa[i : i + window]
        if "P" in w or "*" in w:
            continue
        score = sum(KD.get(a, 0.0) for a in w) / window
        if score >= thresh and (best is None or i > best[0]):
            best = (i, i + window)
    return best


def choose_orf(orfs: list[dict], must: list[tuple[int, int]]) -> dict | None:
    scored = []
    for o in orfs:
        if not all(overlaps(o["nt0"], o["nt1"] - 1, a, b - 1) for a, b in must if a is not None):
            continue
        score = 100 * int(o["has_M"]) + 50 * int(o["has_stop"]) + len(o["aa"])
        scored.append((score, o))
    if not scored:
        return None
    scored.sort(key=lambda x: -x[0])
    return scored[0][1]


def cds_exons(rec: dict, nt0: int, nt1: int) -> list[tuple[int, int]]:
    exons = rec["exons"] if rec["strand"] == "+" else list(reversed(rec["exons"]))
    tpos = 0
    out = []
    for s, e in exons:
        elen = e - s + 1
        a, b = tpos, tpos + elen
        if b <= nt0 or a >= nt1:
            tpos = b
            continue
        clip0 = max(0, nt0 - a)
        clip1 = min(elen, nt1 - a)
        if rec["strand"] == "+":
            gs, ge = s + clip0, s + clip1 - 1
        else:
            gs, ge = e - clip1 + 1, e - clip0
        out.append((min(gs, ge), max(gs, ge)))
        tpos = b
    out.sort()
    return out


def empty_row(gene: dict, notes: str, rec=None) -> dict:
    I, V = gene["I"], gene["V"]
    return {
        "gene_id": gene["gene_id"],
        "chrom": gene["chrom"],
        "strand": gene["strand"],
        "arch": gene["arch"],
        "v_start": V["start"] if V else "",
        "v_end": V["end"] if V else "",
        "i_start": I["start"] if I else "",
        "i_end": I["end"] if I else "",
        "transcript_id": rec["tid"] if rec else "-",
        "n_exons": len(rec["exons"]) if rec else 0,
        "cov": f"{rec['cov']:.2f}" if rec else "0",
        "protein_len": 0,
        "has_V": "yes" if V else "no",
        "has_I": "yes" if I else "no",
        "has_TM": "no",
        "has_M": "no",
        "has_stop": "no",
        "notes": notes,
        "protein": "",
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genome", type=Path, required=True)
    p.add_argument("--gtf", type=Path, help="StringTie GTF; omit if no RNA-seq")
    p.add_argument("--i-domains", type=Path, required=True, help="BLAST confirmed I FASTA (merged into catalog)")
    p.add_argument("--ig-csv", type=Path, required=True)
    p.add_argument("--gene-prefix", default="NITR", help="GeneID prefix (e.g. PsenNITR)")
    p.add_argument("-o", "--prefix", type=Path, required=True)
    args = p.parse_args()

    domains = catalog_domains(args.ig_csv, args.i_domains)
    genes = pair_genes(domains, args.gene_prefix)
    tx = parse_gtf(args.gtf) if args.gtf and args.gtf.exists() else {}
    genome = Fasta(str(args.genome), as_raw=True, sequence_always_upper=True)

    args.prefix.parent.mkdir(parents=True, exist_ok=True)
    loci_tsv = Path(f"{args.prefix}_ig_loci.tsv")
    prot_fa = Path(f"{args.prefix}_proteins.fa")
    cds_fa = Path(f"{args.prefix}_cds.fa")
    gtf_out = Path(f"{args.prefix}_models.gtf")
    tsv_out = Path(f"{args.prefix}_summary.tsv")
    prot_fa.write_text("")
    cds_fa.write_text("")

    with loci_tsv.open("w", newline="") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["dom_id", "cls", "chrom", "start", "end", "strand", "source", "hmms", "nC"],
            delimiter="\t",
        )
        w.writeheader()
        for d in domains:
            w.writerow({k: d.get(k, "") for k in w.fieldnames})

    rows = []
    gtf_lines = []

    for gene in genes:
        notes = [gene["arch"]]
        cands = [rec for rec in tx.values() if tx_hits_gene(rec, gene)]
        cands = [rec for rec in cands if n_genes_hit(rec, genes) == 1]
        cands.sort(key=lambda r: (-r["cov"], -(r["end"] - r["start"])))
        rec = cands[0] if cands else None
        if rec is None:
            notes.append("no_unique_transcript")
            rows.append(empty_row(gene, ";".join(notes)))
            continue

        nt = spliced_seq(genome, rec)
        must = []
        v_tx = genomic_to_tx(rec, gene["V"]["start"], gene["V"]["end"]) if gene["V"] else None
        i_tx = genomic_to_tx(rec, gene["I"]["start"], gene["I"]["end"]) if gene["I"] else None
        if v_tx:
            must.append(v_tx)
        if i_tx:
            must.append(i_tx)
        if not must:
            notes.append("domains_not_on_transcript")
            rows.append(empty_row(gene, ";".join(notes), rec))
            continue

        orf = choose_orf(find_orfs(nt), must)
        if orf is None:
            notes.append("no_ORF_spanning_domains")
            rows.append(empty_row(gene, ";".join(notes), rec))
            continue

        aa = orf["aa"]
        after = 0
        if i_tx:
            after = max(0, (i_tx[0] - orf["nt0"]) // 3 + 20)
        elif v_tx:
            after = max(0, (v_tx[0] - orf["nt0"]) // 3 + 20)
        tm = find_tm(aa, after=after)
        has_tm = tm is not None

        if not orf["has_M"]:
            notes.append("partial_no_start")
        if not orf["has_stop"]:
            notes.append("partial_no_stop")
        if gene["arch"] == "V-only" and not has_tm:
            notes.append("secreted_candidate")
        elif gene["arch"] == "V-only" and has_tm:
            notes.append("V_only_with_TM")
        if gene["arch"] == "V+I" and has_tm and orf["has_M"] and orf["has_stop"]:
            notes.append("complete_V_I_TM")
        if has_tm:
            notes.append("has_TM")
        else:
            notes.append("no_TM")

        cds = nt[orf["nt0"] : orf["nt1"] + (3 if orf["has_stop"] else 0)]
        exons_cds = cds_exons(rec, orf["nt0"], orf["nt0"] + len(cds))
        g0, g1 = exons_cds[0][0], exons_cds[-1][1]
        gid = gene["gene_id"]
        header = (
            f"{gid} {gene['chrom']}:{g0}-{g1}({gene['strand']}) "
            f"arch={gene['arch']} tx={rec['tid']} aa={len(aa)} notes={';'.join(notes)}"
        )
        with prot_fa.open("a") as fh:
            fh.write(f">{header}\n{wrap(aa)}\n")
        with cds_fa.open("a") as fh:
            fh.write(f">{header}\n{wrap(cds)}\n")
        attr = f'gene_id "{gid}"; transcript_id "{gid}.1"; orig_tx "{rec["tid"]}"; arch "{gene["arch"]}";'
        gtf_lines.append(f"{gene['chrom']}\tnitr\ttranscript\t{g0}\t{g1}\t.\t{gene['strand']}\t.\t{attr}\n")
        for es, ee in exons_cds:
            gtf_lines.append(f"{gene['chrom']}\tnitr\texon\t{es}\t{ee}\t.\t{gene['strand']}\t.\t{attr}\n")

        I, V = gene["I"], gene["V"]
        rows.append(
            {
                "gene_id": gid,
                "chrom": gene["chrom"],
                "strand": gene["strand"],
                "arch": gene["arch"],
                "v_start": V["start"] if V else "",
                "v_end": V["end"] if V else "",
                "i_start": I["start"] if I else "",
                "i_end": I["end"] if I else "",
                "transcript_id": rec["tid"],
                "n_exons": len(rec["exons"]),
                "cov": f"{rec['cov']:.2f}",
                "protein_len": len(aa),
                "has_V": "yes" if V else "no",
                "has_I": "yes" if I else "no",
                "has_TM": "yes" if has_tm else "no",
                "has_M": "yes" if orf["has_M"] else "no",
                "has_stop": "yes" if orf["has_stop"] else "no",
                "notes": ";".join(notes),
                "protein": aa,
            }
        )

    cols = [
        "gene_id", "chrom", "strand", "arch",
        "v_start", "v_end", "i_start", "i_end",
        "transcript_id", "n_exons", "cov", "protein_len",
        "has_V", "has_I", "has_TM", "has_M", "has_stop", "notes", "protein",
    ]
    with tsv_out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    gtf_out.write_text("".join(gtf_lines))

    n_vi = sum(1 for g in genes if g["arch"] == "V+I")
    n_v = sum(1 for g in genes if g["arch"] == "V-only")
    n_i = sum(1 for g in genes if g["arch"] == "I-only")
    print(f"Ig loci: {len(domains)}  (I={sum(1 for d in domains if d['cls']=='I')} V={sum(1 for d in domains if d['cls']=='V')})")
    print(f"Gene models: {len(genes)}  V+I={n_vi}  V-only={n_v}  I-only={n_i}")
    print(f"StringTie transcripts: {len(tx)}")
    print(f"  {loci_tsv}")
    print(f"  {tsv_out}")
    print(f"  {prot_fa}")
    print(f"  {cds_fa}")
    print(f"  {gtf_out}")
    for r in rows:
        print(
            f"  {r['gene_id']} {r['arch']} {r['chrom']}:{r['v_start'] or r['i_start']}-{r['i_end'] or r['v_end']}"
            f"({r['strand']}) aa={r['protein_len']} TM={r['has_TM']} {r['notes']}"
        )


if __name__ == "__main__":
    main()
