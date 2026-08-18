#!/usr/bin/env python3
"""NITR exon CSV: one SP, one V, one I, TM, cyto — labeled from tools, coords from exons.

  SignalP     → SP (only if predicted)
  InterProScan → one V/D1 (5' Ig cluster) and one I/D2 (3' Ig cluster, if distinct)
  DeepTMHMM   → TM helix and cytoplasmic tail

Each GTF exon is assigned to at most one of SP / V / I. TM and cyto are the
DeepTMHMM intervals (split out of the 3' exon). No leftover '-' fragments,
no nested extra Ig hits, no secreted_tail.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq
from pyfaidx import Fasta

GID_RE = re.compile(r'gene_id "([^"]+)"')
DOM_RE = re.compile(r'domain "([^"]+)"')
IG_NAME_RE = re.compile(r"ig|I-set|Ig|V-set|IG|immunoglobulin", re.I)
FALSE_IG_RE = re.compile(r"^(Lig|PIG|LIGAN|BRIGHT)", re.I)
VSET_RE = re.compile(r"V-set|igv|IGV|PF07686|SM00406", re.I)

CANON_ORDER = ["SP", "V/D1", "I/D2", "TM", "cyto"]


def load_fasta(path: Path) -> dict[str, str]:
    seqs, name, chunks = {}, None, []
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
            rec["exons"].append({"start": start, "end": end, "gtf_domain": (dm.group(1) if dm else "")})
    for rec in genes.values():
        rec["exons"].sort(key=lambda x: (x["start"], x["end"]))
    return genes


def tx_positions(exons: list[dict], strand: str) -> list[int]:
    ordered = exons if strand == "+" else list(reversed(exons))
    pos = []
    for ex in ordered:
        s, e = ex["start"], ex["end"]
        if strand == "+":
            pos.extend(range(s, e + 1))
        else:
            pos.extend(range(e, s - 1, -1))
    return pos


def aa_span_of_exon(ex: dict, pos: list[int], strand: str) -> tuple[int, int] | None:
    """1-based AA coordinates covered by this exon."""
    s, e = ex["start"], ex["end"]
    idxs = [i for i, g in enumerate(pos) if s <= g <= e]
    if not idxs:
        return None
    return idxs[0] // 3 + 1, idxs[-1] // 3 + 1


def fetch_nt(genome: Fasta, chrom: str, start: int, end: int, strand: str) -> str:
    seq = str(genome[chrom][start - 1 : end]).upper()
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())
    return seq


def aa_to_genome(pos: list[int], aa0: int, aa1: int) -> list[tuple[int, int]]:
    """Map inclusive 1-based AA span to genomic blocks."""
    i0, i1 = (aa0 - 1) * 3, aa1 * 3
    if i0 >= len(pos) or i1 <= i0:
        return []
    coords = pos[i0 : min(i1, len(pos))]
    if not coords:
        return []
    out = []
    a = b = coords[0]
    for p in coords[1:]:
        if abs(p - b) == 1:
            b = p
        else:
            out.append((min(a, b), max(a, b)))
            a = b = p
    out.append((min(a, b), max(a, b)))
    return out


def parse_signalp(gff: Path, summary: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if summary.exists():
        with summary.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                out[p[0]] = {"is_sp": p[1].upper().startswith("SP"), "score": p[2] if len(p) > 2 else ""}
    if gff.exists():
        with gff.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                if p[2].lower() not in {"signal_peptide", "signal_peptide_region_of_n_terminal"}:
                    continue
                rec = out.setdefault(p[0], {})
                rec["start"], rec["end"] = int(p[3]), int(p[4])
                rec["is_sp"] = True
                rec.setdefault("score", p[5] if len(p) > 5 else "")
    return out


def parse_deeptmhmm(gff: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    if not gff.exists():
        return out
    with gff.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("//"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            gid, kind = p[0], p[1].lower()
            rec = out.setdefault(gid, {})
            span = (int(p[2]), int(p[3]))
            if kind == "tmhelix":
                rec["tm"] = span
            elif kind == "inside":
                rec["cyto"] = span
    return out


def parse_ips_ig(path: Path) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = defaultdict(list)
    if not path.exists():
        return hits
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            gene, acc, desc = p[0], p[4], p[5]
            if FALSE_IG_RE.search(acc):
                continue
            if not IG_NAME_RE.search(f"{acc} {desc}"):
                continue
            try:
                start, end = int(p[6]), int(p[7])
            except ValueError:
                continue
            hits[gene].append(
                {
                    "start": start,
                    "end": end,
                    "acc": acc,
                    "desc": desc,
                    "vset": bool(VSET_RE.search(f"{acc} {desc}")),
                    "source": f"IPS:{p[3]}",
                    "score": p[8],
                }
            )
    return hits


def cluster_ig(hits: list[dict]) -> list[dict]:
    """Merge overlapping Ig hits; return clusters sorted 5'→3'."""
    if not hits:
        return []
    hits = sorted(hits, key=lambda h: (h["start"], -h["end"]))
    clusters = []
    for h in hits:
        placed = False
        for c in clusters:
            ov = max(0, min(h["end"], c["end"]) - max(h["start"], c["start"]) + 1)
            shorter = min(h["end"] - h["start"] + 1, c["end"] - c["start"] + 1)
            if ov >= 20 and shorter and ov / shorter >= 0.25:
                c["start"] = min(c["start"], h["start"])
                c["end"] = max(c["end"], h["end"])
                c["members"].append(h)
                c["vset"] = c["vset"] or h["vset"]
                placed = True
                break
        if not placed:
            clusters.append(
                {"start": h["start"], "end": h["end"], "members": [h], "vset": h["vset"]}
            )
    clusters.sort(key=lambda c: c["start"])
    return clusters


def pick_v_i(clusters: list[dict]) -> tuple[dict | None, dict | None]:
    """At most one V and one I. Nested/overlapping Ig is not a second domain."""
    if not clusters:
        return None, None
    if len(clusters) == 1:
        c = clusters[0]
        if c["vset"]:
            return c, None
        return None, c
    v, i = clusters[0], clusters[-1]
    ov = max(0, min(v["end"], i["end"]) - max(v["start"], i["start"]) + 1)
    shorter = min(v["end"] - v["start"] + 1, i["end"] - i["start"] + 1)
    if shorter and ov / shorter > 0.5:
        return (v, None) if v["vset"] or not i["vset"] else (None, i)
    if i["start"] <= v["end"]:
        i = {**i, "start": v["end"] + 1}
        if i["end"] - i["start"] + 1 < 30:
            return v, None
    return v, i


def overlap_aa(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0) + 1)


def best_hit_seq(cluster: dict, aa: str) -> tuple[str, str, str]:
    """Representative IPS accession/score and peptide for a cluster."""
    best = max(cluster["members"], key=lambda h: h["end"] - h["start"])
    pep = aa[cluster["start"] - 1 : cluster["end"]]
    return best["source"], best["acc"], pep


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genome", type=Path, required=True)
    p.add_argument("--gtf", type=Path, required=True)
    p.add_argument("--cds", type=Path, required=True)
    p.add_argument("--proteins", type=Path, required=True)
    p.add_argument("--ips", type=Path, required=True)
    p.add_argument("--signalp-gff", type=Path, required=True)
    p.add_argument("--signalp-summary", type=Path, required=True)
    p.add_argument("--deeptmhmm-gff", type=Path, required=True)
    p.add_argument("--architecture", type=Path, required=True)
    p.add_argument("-o", "--out", type=Path, required=True)
    args = p.parse_args()

    genome = Fasta(str(args.genome), as_raw=True, sequence_always_upper=True)
    genes = parse_gtf(args.gtf)
    proteins = load_fasta(args.proteins)
    arch = {r["gene_id"]: r for r in csv.DictReader(args.architecture.open(), delimiter="\t")}
    sigp = parse_signalp(args.signalp_gff, args.signalp_summary)
    dtm = parse_deeptmhmm(args.deeptmhmm_gff)
    ips_ig = parse_ips_ig(args.ips)

    rows = []
    for gid, rec in genes.items():
        aa = proteins.get(gid, "")
        strand, chrom = rec["strand"], rec["chrom"]
        pos = tx_positions(rec["exons"], strand)
        meta = arch.get(gid, {})
        notes = []
        if meta.get("ITIM") == "yes":
            notes.append("ITIM")
        if meta.get("ITSM") == "yes":
            notes.append("ITSM")
        tx = meta.get("RNA", "-")
        if tx and tx not in {"", "-"}:
            notes.append(f"RNA:{tx}")
        elif meta.get("notes") and "no_RNA" in meta["notes"]:
            notes.append("no_RNA")

        domains: dict[str, dict] = {}
        srec = sigp.get(gid, {})
        if srec.get("is_sp") and "start" in srec:
            domains["SP"] = {
                "aa0": srec["start"],
                "aa1": srec["end"],
                "pep": aa[srec["start"] - 1 : srec["end"]],
            }
        v_cl, i_cl = pick_v_i(cluster_ig(ips_ig.get(gid, [])))
        if v_cl:
            src, acc, pep = best_hit_seq(v_cl, aa)
            domains["V/D1"] = {"aa0": v_cl["start"], "aa1": v_cl["end"], "pep": pep}
        if i_cl:
            src, acc, pep = best_hit_seq(i_cl, aa)
            domains["I/D2"] = {"aa0": i_cl["start"], "aa1": i_cl["end"], "pep": pep}
        drec = dtm.get(gid, {})
        if drec.get("tm"):
            t0, t1 = drec["tm"]
            domains["TM"] = {"aa0": t0, "aa1": t1, "pep": aa[t0 - 1 : t1]}
        if drec.get("cyto"):
            c0, c1 = drec["cyto"]
            domains["cyto"] = {"aa0": c0, "aa1": c1, "pep": aa[c0 - 1 : c1]}

        # GTF I exon fallback when IPS missed (NITR10)
        if "I/D2" not in domains:
            for ex in rec["exons"]:
                if ex["gtf_domain"] == "I":
                    asp = aa_span_of_exon(ex, pos, strand)
                    if asp:
                        domains["I/D2"] = {
                            "aa0": asp[0],
                            "aa1": asp[1],
                            "pep": aa[asp[0] - 1 : asp[1]],
                        }
                    break

        # Assign SP / V / I to at most one exon each (the exon with most overlapping AA)
        best_exon: dict[str, tuple[int, dict]] = {}
        for dname in ("SP", "V/D1", "I/D2"):
            if dname not in domains:
                continue
            d = domains[dname]
            scored = []
            for ex in rec["exons"]:
                asp = aa_span_of_exon(ex, pos, strand)
                if not asp:
                    continue
                ov = overlap_aa(asp[0], asp[1], d["aa0"], d["aa1"])
                if ov >= 8:
                    scored.append((ov, ex))
            if scored:
                scored.sort(key=lambda x: -x[0])
                best_exon[dname] = (scored[0][0], scored[0][1])

        note = ";".join(notes)
        for dname in ("SP", "V/D1", "I/D2"):
            if dname not in best_exon:
                continue
            ex = best_exon[dname][1]
            rows.append(
                {
                    "Scaffold": chrom,
                    "GeneID": gid,
                    "Notes": note,
                    "Exon Start": ex["start"],
                    "Exon Stop": ex["end"],
                    "Orientation": strand,
                    "Protein Domain": dname,
                    "Nucleotide sequence": fetch_nt(genome, chrom, ex["start"], ex["end"], strand),
                    "Protein sequence": domains[dname]["pep"],
                    "_gene_start": min(e["start"] for e in rec["exons"]),
                }
            )

        for dname in ("TM", "cyto"):
            if dname not in domains:
                continue
            d = domains[dname]
            for gs, ge in aa_to_genome(pos, d["aa0"], d["aa1"]):
                if ge - gs + 1 < 6:
                    continue
                rows.append(
                    {
                        "Scaffold": chrom,
                        "GeneID": gid,
                        "Notes": note,
                        "Exon Start": gs,
                        "Exon Stop": ge,
                        "Orientation": strand,
                        "Protein Domain": dname,
                        "Nucleotide sequence": fetch_nt(genome, chrom, gs, ge, strand),
                        "Protein sequence": d["pep"],
                        "_gene_start": min(e["start"] for e in rec["exons"]),
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
