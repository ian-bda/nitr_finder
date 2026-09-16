#!/usr/bin/env python3
"""Step 10: gene-level table for NITRs + in-cluster non-NITRs + flanks.

Does NOT use NCBI Polypterus annotation. Loci come from:
  - our NITR 8b models
  - miniprot of related-species proteins (zebrafish / gar) onto a window
  - StringTie transcripts in that window that are not NITRs

Overlapping homology hits are collapsed to one locus (best identity).
Hits overlapping a NITR are dropped. ~n_flank genes are kept on each side
of the NITR span; every non-overlapping locus between the first and last
NITR is kept.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio.Seq import Seq
from pyfaidx import Fasta

GID_RE = re.compile(r'gene_id "([^"]+)"')
TID_RE = re.compile(r'transcript_id "([^"]+)"')
ATTR_RE = re.compile(r"(?:^|;)([^=;]+)=([^;]*)")
TARGET_RE = re.compile(r"Target=([^\s;]+)")
ID_RE = re.compile(r"(?:^|;)ID=([^;]+)")
IDENT_RE = re.compile(r"Identity=([0-9.]+)", re.I)
NITR_QUERY_RE = re.compile(r"\bnitr\b|novel immune-type|immune-type receptor", re.I)


def parse_gtf_genes(path: Path, feature: str = "exon") -> list[dict]:
    genes: dict[str, dict] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != feature:
                continue
            chrom, start, end, strand, attrs = p[0], int(p[3]), int(p[4]), p[6], p[8]
            gm = GID_RE.search(attrs)
            if not gm:
                continue
            gid = gm.group(1)
            rec = genes.setdefault(
                gid,
                {
                    "gene_id": gid,
                    "chrom": chrom,
                    "strand": strand,
                    "start": start,
                    "end": end,
                    "exons": [],
                    "transcripts": {},
                },
            )
            rec["start"] = min(rec["start"], start)
            rec["end"] = max(rec["end"], end)
            rec["exons"].append((start, end))
            tm = TID_RE.search(attrs)
            if tm:
                rec["transcripts"].setdefault(tm.group(1), []).append((start, end))
            if strand in "+-":
                rec["strand"] = strand
    return list(genes.values())


def parse_attr(blob: str) -> dict[str, str]:
    out = {}
    for k, v in ATTR_RE.findall(blob):
        out[k.strip()] = v.strip()
    return out


def parse_miniprot_gff(path: Path, bed_start0: int, chrom: str) -> list[dict]:
    """Lift window-local miniprot mRNA/CDS coords onto the scaffold."""
    by_id: dict[str, dict] = {}
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            kind = p[2]
            local_s, local_e = int(p[3]), int(p[4])
            strand = p[6] if p[6] in "+-" else "."
            attrs = p[8]
            if kind in {"mRNA", "gene"}:
                ident_m = IDENT_RE.search(attrs)
                ident = float(ident_m.group(1)) if ident_m else 0.0
                tgt_m = TARGET_RE.search(attrs)
                target = tgt_m.group(1) if tgt_m else parse_attr(attrs).get("Target", "")
                target = target.split()[0] if target else ""
                im = ID_RE.search(attrs)
                mp_id = im.group(1) if im else f"anon{len(by_id)}"
                try:
                    score = float(p[5]) if p[5] not in {".", ""} else 0.0
                except ValueError:
                    score = 0.0
                rec = by_id.setdefault(mp_id, {"cds": []})
                rec.update(
                    {
                        "chrom": chrom,
                        "start": bed_start0 + local_s,
                        "end": bed_start0 + local_e,
                        "strand": strand,
                        "score": score,
                        "identity": ident,
                        "query": target,
                        "mp_id": mp_id,
                    }
                )
            elif kind == "CDS":
                parent = parse_attr(attrs).get("Parent", "")
                try:
                    phase = int(p[7]) if p[7] not in {".", ""} else 0
                except ValueError:
                    phase = 0
                rec = by_id.setdefault(parent, {"cds": []})
                rec["cds"].append((bed_start0 + local_s, bed_start0 + local_e, phase))
    return [h for h in by_id.values() if "query" in h]


def load_faa_names(path: Path) -> dict[str, str]:
    """query_id -> description (protein FASTA headers)."""
    names = {}
    with path.open() as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            hid = line[1:].split()[0]
            desc = line[1:].strip()
            names[hid] = desc
            names[hid.split(".")[0]] = desc
    return names


def nice_name(header: str, query: str) -> str:
    """Pull a gene-like token from an NCBI/Ensembl protein header."""
    if not header:
        return query or "homology"
    h = header
    m = re.search(r"gene[:=]([A-Za-z0-9._-]+)", h, re.I)
    if m and not m.group(1).startswith("ENS"):
        return m.group(1)
    # NCBI: "solute carrier family 22 member 4 [Danio rerio]"
    rest = h.split(None, 1)[1] if " " in h else h
    rest = re.sub(r"\s*\[[^\]]+\]\s*$", "", rest).strip()
    # symbol sometimes at the end after "isoform X"
    m = re.search(r"\(([^)]+)\)\s*$", rest)
    if m and 2 <= len(m.group(1)) <= 20 and " " not in m.group(1):
        return m.group(1)
    slug = re.sub(r"[^A-Za-z0-9]+", "_", rest)
    slug = re.sub(r"_+", "_", slug).strip("_")[:40]
    return slug or (query or "homology")


def span(d: dict) -> tuple[int, int]:
    if "start" in d:
        return int(d["start"]), int(d["end"])
    return int(d["Gene Start"]), int(d["Gene Stop"])


def overlap_bp(a: dict, b: dict, pad: int = 0) -> int:
    a0, a1 = span(a)
    b0, b1 = span(b)
    lo = max(a0 - pad, b0 - pad)
    hi = min(a1 + pad, b1 + pad)
    return max(0, hi - lo + 1)


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


def fetch_nt(genome: Fasta, chrom: str, start: int, end: int, strand: str) -> str:
    seq = str(genome[chrom][start - 1 : end]).upper()
    if strand == "-":
        seq = str(Seq(seq).reverse_complement())
    return seq


STOP_CODONS = {"TAG", "TAA", "TGA"}


def nt_with_coords(
    genome: Fasta, chrom: str, intervals: list[tuple[int, int]], strand: str
) -> tuple[str, list[int]]:
    """Sense-strand bases and the genomic position of each base."""
    bases: list[str] = []
    coords: list[int] = []
    for s, e in intervals:
        if e < s:
            continue
        if strand == "+":
            seq = str(genome[chrom][s - 1 : e]).upper()
            pos = range(s, e + 1)
        else:
            seq = str(Seq(str(genome[chrom][s - 1 : e])).reverse_complement()).upper()
            pos = range(e, s - 1, -1)
        if len(seq) != (e - s + 1):
            continue
        bases.append(seq)
        coords.extend(pos)
    return "".join(bases), coords


def allowed_range(
    start: int, end: int, nitrs: list[dict], chrom_n: int, pad: int = 50
) -> tuple[int, int]:
    """Do not walk start/stop into a NITR model."""
    prev = max((n["end"] for n in nitrs if n["end"] < start), default=0)
    nxt = min((n["start"] for n in nitrs if n["start"] > end), default=chrom_n + 1)
    lo, hi = prev + pad, nxt - pad
    if lo > hi:
        return start, end
    return lo, hi


def codon_at(
    genome: Fasta, chrom: str, strand: str, g0: int, g1: int
) -> str:
    if strand == "+":
        return str(genome[chrom][g0 - 1 : g1]).upper()
    return str(Seq(str(genome[chrom][g0 - 1 : g1])).reverse_complement()).upper()


def extend_terminal(
    genome: Fasta,
    chrom: str,
    strand: str,
    five_g: int,
    three_g: int,
    allow_lo: int,
    allow_hi: int,
    max_bp: int,
) -> tuple[str, list[int], str, list[int]]:
    """In-frame 5'/3' extension of the terminal exons, stopping at in-frame stops."""
    up: list[tuple[str, list[int]]] = []
    down_nt, down_coords = [], []
    if strand == "+":
        pos = five_g
        while pos - 3 >= allow_lo and five_g - (pos - 3) <= max_bp:
            a, b = pos - 3, pos - 1
            codon = codon_at(genome, chrom, strand, a, b)
            if len(codon) < 3 or codon in STOP_CODONS:
                break
            up.append((codon, [a, a + 1, b]))
            pos = a
        pos = three_g
        while pos + 3 <= allow_hi and (pos + 3) - three_g <= max_bp:
            a, b = pos + 1, pos + 3
            codon = codon_at(genome, chrom, strand, a, b)
            if len(codon) < 3 or codon in STOP_CODONS:
                break
            down_nt.append(codon)
            down_coords.extend([a, a + 1, b])
            pos = b
    else:
        pos = five_g
        while pos + 3 <= allow_hi and (pos + 3) - five_g <= max_bp:
            a, b = pos + 1, pos + 3
            codon = codon_at(genome, chrom, strand, a, b)
            if len(codon) < 3 or codon in STOP_CODONS:
                break
            up.append((codon, [b, b - 1, a]))
            pos = b
        pos = three_g
        while pos - 3 >= allow_lo and three_g - (pos - 3) <= max_bp:
            a, b = pos - 3, pos - 1
            codon = codon_at(genome, chrom, strand, a, b)
            if len(codon) < 3 or codon in STOP_CODONS:
                break
            down_nt.append(codon)
            down_coords.extend([b, b - 1, a])
            pos = a
    up.reverse()
    return (
        "".join(c for c, _ in up),
        [p for _, ps in up for p in ps],
        "".join(down_nt),
        down_coords,
    )


def orf_from_frame0(nt: str, coords: list[int]) -> tuple[str, int | None, int | None]:
    """First Met through in-frame stop (or end) in frame 0."""
    nt = nt[: len(nt) - len(nt) % 3]
    if len(nt) < 6:
        return "", None, None
    pep = str(Seq(nt).translate(to_stop=False))
    i = pep.find("M")
    if i < 0:
        return "", None, None
    j = pep.find("*", i)
    piece = pep[i:] if j < 0 else pep[i:j]
    nt_a, nt_b = i * 3, (i + len(piece)) * 3
    gpos = coords[nt_a:nt_b]
    if not gpos:
        return piece, None, None
    return piece, min(gpos), max(gpos)


def longest_orf_span(nt: str, coords: list[int]) -> tuple[str, int, int, int | None, int | None]:
    """Longest Met→stop in three frames; return pep, nt_start, nt_end, g0, g1."""
    best: tuple[str, int, int, int | None, int | None] = ("", 0, 0, None, None)
    for fr in range(3):
        usable = nt[fr : len(nt) - (len(nt) - fr) % 3]
        if len(usable) < 6:
            continue
        pep = str(Seq(usable).translate(to_stop=False))
        i = 0
        while i < len(pep):
            if pep[i] != "M":
                i += 1
                continue
            j = pep.find("*", i)
            piece = pep[i:] if j < 0 else pep[i:j]
            nt_a = fr + i * 3
            nt_b = nt_a + len(piece) * 3
            if len(piece) > len(best[0]):
                gpos = coords[nt_a:nt_b]
                g0 = min(gpos) if gpos else None
                g1 = max(gpos) if gpos else None
                best = (piece, nt_a, nt_b, g0, g1)
            i = j + 1 if j >= 0 else len(pep)
    return best


def finish_met_stop(
    genome: Fasta,
    chrom: str,
    cds: list[tuple[int, int, int]],
    strand: str,
    allowed: tuple[int, int] | None = None,
    max_bp: int = 1500,
) -> tuple[str, int | None, int | None]:
    """Complete a miniprot alignment to Met→stop in the alignment frame.

    Walk the first/last exon codon-by-codon until an in-frame stop, but only
    within max_bp of the aligned CDS so a short RS-repeat exon cannot eat a
    neighboring gene. If the aligned CDS itself is frameshifted, take the
    longest clean ORF inside it and extend that frame the same way.
    """
    if not cds or strand not in "+-":
        return "", None, None
    chrom_n = len(genome[chrom])
    cds_lo = min(s for s, e, _p in cds)
    cds_hi = max(e for s, e, _p in cds)
    allow_lo, allow_hi = (1, chrom_n) if allowed is None else allowed
    allow_lo = max(1, allow_lo, cds_lo - max_bp)
    allow_hi = min(chrom_n, allow_hi, cds_hi + max_bp)

    tx = sorted(cds, key=lambda x: x[0])
    if strand == "-":
        tx = list(reversed(tx))
    first_s, first_e, phase = tx[0]
    phase = phase if phase in (1, 2) else 0

    core_iv: list[tuple[int, int]] = []
    if strand == "+":
        a = first_s + phase
        if a <= first_e:
            core_iv.append((a, first_e))
    else:
        b = first_e - phase
        if first_s <= b:
            core_iv.append((first_s, b))
    for s, e, _ph in tx[1:]:
        if s <= e:
            core_iv.append((s, e))
    if not core_iv:
        return "", None, None

    core_nt, core_coords = nt_with_coords(genome, chrom, core_iv, strand)
    core_nt = core_nt[: len(core_nt) - len(core_nt) % 3]
    core_coords = core_coords[: len(core_nt)]
    if not core_nt:
        return "", None, None

    core_pep = str(Seq(core_nt).translate(to_stop=False))
    if "*" in core_pep:
        piece, nt_a, nt_b, _g0, _g1 = longest_orf_span(core_nt, core_coords)
        if not piece:
            return "", None, None
        core_nt = core_nt[nt_a:nt_b]
        core_coords = core_coords[nt_a:nt_b]

    if strand == "+":
        five_g, three_g = core_coords[0], core_coords[-1]
    else:
        five_g, three_g = core_coords[0], core_coords[-1]

    up_nt, up_coords, down_nt, down_coords = extend_terminal(
        genome, chrom, strand, five_g, three_g, allow_lo, allow_hi, max_bp
    )
    nt = up_nt + core_nt + down_nt
    coords = up_coords + core_coords + down_coords
    pep, g0, g1 = orf_from_frame0(nt, coords)
    if pep:
        return pep, g0, g1
    return longest_orf(core_nt), None, None


def longest_orf(nt: str) -> str:
    """Longest Met→stop peptide in three frames (nt already sense-strand)."""
    best = ""
    for fr in range(3):
        seq = nt[fr : len(nt) - (len(nt) - fr) % 3]
        if len(seq) < 6:
            continue
        pep = str(Seq(seq).translate(to_stop=False))
        i = 0
        while i < len(pep):
            if pep[i] != "M":
                i += 1
                continue
            j = pep.find("*", i)
            piece = pep[i:] if j < 0 else pep[i:j]
            if len(piece) > len(best):
                best = piece
            i = j + 1 if j >= 0 else len(pep)
    return best


def splice_exons(genome: Fasta, chrom: str, exons: list[tuple[int, int]], strand: str) -> str:
    ordered = sorted(exons)
    if strand == "-":
        ordered = list(reversed(ordered))
    return "".join(fetch_nt(genome, chrom, s, e, strand) for s, e in ordered)


def protein_from_stringtie(genome: Fasta, rec: dict) -> str:
    chrom, strand = rec["chrom"], rec["strand"]
    txs = rec.get("transcripts") or {}
    candidates = list(txs.values()) if txs else [rec.get("exons") or []]
    best = ""
    strands = [strand] if strand in "+-" else ["+", "-"]
    for exons in candidates:
        if not exons:
            continue
        for st in strands:
            pep = longest_orf(splice_exons(genome, chrom, exons, st))
            if len(pep) > len(best):
                best = pep
    return best


def merge_hits(hits: list[dict], min_ident: float) -> list[dict]:
    hits = [h for h in hits if h["identity"] >= min_ident - 0.005]
    hits.sort(key=lambda h: (-h["identity"], -h["score"]))
    kept: list[dict] = []
    for h in hits:
        rival = next((k for k in kept if overlap_bp(h, k) >= 50), None)
        if rival is None:
            kept.append(dict(h))
        elif h["identity"] > rival["identity"] or (
            h["identity"] == rival["identity"] and h["score"] > rival["score"]
        ):
            kept.remove(rival)
            kept.append(dict(h))
    kept.sort(key=lambda h: h["start"])
    return kept


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nitr-gtf", type=Path, required=True)
    p.add_argument("--miniprot-gff", type=Path, required=True)
    p.add_argument("--window-bed", type=Path, required=True)
    p.add_argument("--proteins", type=Path, required=True, help="query protein FASTA (headers)")
    p.add_argument("--genome", type=Path, required=True, help="reference genome FASTA")
    p.add_argument("--nitr-proteins", type=Path, required=True, help="NITR proteins from step 8b")
    p.add_argument("--stringtie-gtf", type=Path, help="StringTie GTF in the same window")
    p.add_argument("--n-flank", type=int, default=6)
    p.add_argument("--min-identity", type=float, default=0.30)
    p.add_argument("-o", "--out", type=Path, required=True)
    args = p.parse_args()

    genome = Fasta(str(args.genome))
    nitr_prot = load_fasta(args.nitr_proteins)

    with args.window_bed.open() as fh:
        bed = next(l for l in fh if l.strip() and not l.startswith("#")).split("\t")
    chrom = bed[0]
    bed_start0 = int(bed[1])

    nitrs = parse_gtf_genes(args.nitr_gtf)
    nitrs.sort(key=lambda g: g["start"])
    if not nitrs:
        raise SystemExit("no NITR genes in GTF")
    nitr_min = min(g["start"] for g in nitrs)
    nitr_max = max(g["end"] for g in nitrs)

    names = load_faa_names(args.proteins)
    raw = parse_miniprot_gff(args.miniprot_gff, bed_start0, chrom)
    loci = merge_hits(raw, args.min_identity)
    chrom_n = len(genome[chrom])

    rows = []
    for n in nitrs:
        pep = nitr_prot.get(n["gene_id"], "")
        rows.append(
            {
                "GeneID": n["gene_id"],
                "Scaffold": n["chrom"],
                "Gene Start": n["start"],
                "Gene Stop": n["end"],
                "Strand": n["strand"],
                "Gene Type": "NITR",
                "Notes": f"8b_model;aa={len(pep)};protein=8b_model",
                "Protein": pep,
            }
        )

    def hits_nitr(h: dict) -> bool:
        return any(overlap_bp(h, n, pad=2000) >= 50 for n in nitrs)

    homology = []
    for h in loci:
        if hits_nitr(h):
            continue
        header = names.get(h["query"], names.get(h["query"].split(".")[0], h["query"]))
        if NITR_QUERY_RE.search(header or "") and nitr_min <= h["start"] <= nitr_max:
            continue
        gid = nice_name(header, h["query"])
        pep, g0, g1 = finish_met_stop(
            genome,
            h["chrom"],
            h.get("cds") or [],
            h["strand"],
            allowed=allowed_range(h["start"], h["end"], nitrs, chrom_n),
        )
        gstart, gstop = h["start"], h["end"]
        if g0 is not None and g1 is not None:
            gstart, gstop = min(gstart, g0), max(gstop, g1)
        homology.append(
            {
                "GeneID": gid,
                "Scaffold": h["chrom"],
                "start": gstart,
                "end": gstop,
                "Gene Start": gstart,
                "Gene Stop": gstop,
                "Strand": h["strand"],
                "Gene Type": "",  # fill after
                "Notes": f"miniprot:{h['query']};identity={h['identity']:.2f}",
                "header": header,
                "Protein": pep,
                "protein_src": "miniprot_MetStop" if pep else "",
            }
        )

    stringtie = []
    if args.stringtie_gtf and args.stringtie_gtf.exists():
        for g in parse_gtf_genes(args.stringtie_gtf, feature="exon"):
            if g["strand"] not in "+-":
                g["strand"] = "."
            if hits_nitr(g):
                continue
            stringtie.append(g)

    # merge StringTie into homology when they overlap; else keep as STRG
    used_st = set()
    for h in homology:
        for i, s in enumerate(stringtie):
            if i in used_st:
                continue
            if overlap_bp(h, s) >= 50:
                h["Gene Start"] = h["start"] = min(h["Gene Start"], s["start"])
                h["Gene Stop"] = h["end"] = max(h["Gene Stop"], s["end"])
                h["Notes"] += f";RNA:{s['gene_id']}"
                if not h.get("Protein"):
                    pep = protein_from_stringtie(genome, s)
                    h["Protein"] = pep
                    h["protein_src"] = "StringTie_ORF" if pep else ""
                used_st.add(i)
    for i, s in enumerate(stringtie):
        if i in used_st:
            continue
        pep = protein_from_stringtie(genome, s)
        homology.append(
            {
                "GeneID": s["gene_id"],
                "Scaffold": s["chrom"],
                "Gene Start": s["start"],
                "Gene Stop": s["end"],
                "Strand": s["strand"],
                "Gene Type": "",
                "Notes": "StringTie_only",
                "Protein": pep,
                "protein_src": "StringTie_ORF" if pep else "",
            }
        )

    # unique GeneIDs
    seen = {}
    for h in homology:
        base = h["GeneID"]
        n = seen.get(base, 0)
        seen[base] = n + 1
        if n:
            h["GeneID"] = f"{base}_{n + 1}"

    up = [h for h in homology if h["Gene Stop"] < nitr_min]
    down = [h for h in homology if h["Gene Start"] > nitr_max]
    inside = [h for h in homology if not (h["Gene Stop"] < nitr_min or h["Gene Start"] > nitr_max)]
    up.sort(key=lambda h: h["Gene Stop"], reverse=True)
    up = list(reversed(up[: args.n_flank]))
    down.sort(key=lambda h: h["Gene Start"])
    down = down[: args.n_flank]
    inside.sort(key=lambda h: h["Gene Start"])

    for h in up:
        h["Gene Type"] = "flanking_upstream"
    for h in inside:
        h["Gene Type"] = "nonNITR_in_cluster"
    for h in down:
        h["Gene Type"] = "flanking_downstream"

    for h in (*up, *inside, *down):
        pep = h.get("Protein") or ""
        src = h.get("protein_src") or ("miniprot_CDS" if pep else "none")
        h["Notes"] += f";aa={len(pep)};protein={src}"
        h["Protein"] = pep

    rows.extend(up)
    rows.extend(inside)
    rows.extend(down)
    rows.sort(key=lambda r: (r["Scaffold"], r["Gene Start"], r["Gene Stop"]))

    cols = [
        "GeneID",
        "Scaffold",
        "Gene Start",
        "Gene Stop",
        "Strand",
        "Gene Type",
        "Notes",
        "Protein",
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    print(f"Genes: {len(rows)}")
    print(f"  NITR: {sum(1 for r in rows if r['Gene Type']=='NITR')}")
    print(f"  in-cluster non-NITR: {sum(1 for r in rows if r['Gene Type']=='nonNITR_in_cluster')}")
    print(f"  flanking_upstream: {sum(1 for r in rows if r['Gene Type']=='flanking_upstream')}")
    print(f"  flanking_downstream: {sum(1 for r in rows if r['Gene Type']=='flanking_downstream')}")
    print(f"  {args.out}")
    for r in rows:
        print(
            f"  {r['Gene Start']}-{r['Gene Stop']} {r['Strand']:1} "
            f"{r['Gene Type']:22} {r['GeneID']:24} aa={len(r.get('Protein') or '')} "
            f"{r['Notes'][:90]}"
        )


if __name__ == "__main__":
    main()
