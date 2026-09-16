#!/usr/bin/env python3
"""Step 8b: NITR exons from RNA, or V/I splice sites from the genome.

V/I from step 8 are the anchors. If a same-strand transcript covers this gene,
those spliced exons are the model (fusions are clipped to the gene window).
Without RNA, snap V/I to GT–AG and stop: do not invent SP/TM/cyto from the
genome (that was the old Met→stop stitch). Missing ends stay missing.
"""

from __future__ import annotations

import argparse
import csv
import re
from itertools import product
from pathlib import Path

from Bio.Seq import Seq
from pyfaidx import Fasta

MAX_5P = 16000
MAX_3P = 16000
MIN_INTRON = 50
MAX_INTRON_5P = 14000
MAX_INTRON_3P = 15000
MIN_V = 240
MAX_V = 430
MIN_I = 220
MAX_I = 430
MIN_SP = 36
MAX_SP = 72

TID_RE = re.compile(r'transcript_id "([^"]+)"')
COV_RE = re.compile(r'(?:cov|coverage) "([^"]+)"')
REPEAT_AA = re.compile(r"(IY){5,}|(YI){5,}|I{8,}|Y{8,}|Q{8,}|X{6,}")
ITIM_RE = re.compile(r"[SIVL][A-Z]Y[A-Z]{2}[ILV]")
STOPS = {"TAA", "TAG", "TGA"}
HYD = set("AILMFVWY")

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def wrap(seq: str, n: int = 80) -> str:
    return "\n".join(seq[i : i + n] for i in range(0, len(seq), n))


def rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def translate(nt: str) -> str:
    nt = nt[: len(nt) - (len(nt) % 3)]
    if not nt:
        return ""
    return str(Seq(nt).translate(to_stop=False))


def hyd_frac(aa: str) -> float:
    if not aa:
        return 0.0
    return sum(a in HYD for a in aa) / len(aa)


def hyd_score(aa: str) -> float:
    if not aa:
        return -99.0
    return sum(KD.get(a, 0.0) for a in aa) / len(aa)


def has_tm_helix(aa: str, window: int = 18, thresh: float = 2.20) -> bool:
    if len(aa) < window:
        return False
    return any(hyd_score(aa[i : i + window]) >= thresh for i in range(len(aa) - window + 1))


def is_sp(pep: str) -> bool:
    if not pep.startswith("M") or "*" in pep:
        return False
    if not (12 <= len(pep) <= 24):
        return False
    hyd = set("AILMFVW")
    best = max(sum(1 for a in pep[i : i + 8] if a in hyd) for i in range(0, len(pep) - 7))
    n_charged = sum(1 for a in pep[:20] if a in "DEKR")
    return best >= 5 and n_charged <= 2


def is_tm_pep(pep: str) -> bool:
    pep = pep.replace("*", "")
    if not (15 <= len(pep) <= 55):
        return False
    if ITIM_RE.search(pep):
        return False
    return has_tm_helix(pep)


def low_complexity(aa: str) -> bool:
    return bool(REPEAT_AA.search(aa))


def overlap(a0, a1, b0, b1) -> bool:
    return max(a0, b0) <= min(a1, b1)


class Locus:
    """Genomic interval in 5'→3' transcript coordinates."""

    def __init__(self, genome: Fasta, chrom: str, g0: int, g1: int, strand: str):
        self.chrom = chrom
        self.g0 = g0
        self.g1 = g1
        self.strand = strand
        raw = str(genome[chrom][g0 - 1 : g1]).upper()
        self.seq = rc(raw) if strand == "-" else raw

    def g2t(self, gpos: int) -> int:
        if self.strand == "+":
            return gpos - self.g0
        return self.g1 - gpos

    def t2g(self, tpos: int) -> int:
        if self.strand == "+":
            return self.g0 + tpos
        return self.g1 - tpos

    def interval_tx(self, gs: int, ge: int) -> tuple[int, int]:
        a, b = self.g2t(gs), self.g2t(ge)
        return min(a, b), max(a, b) + 1

    def exons_tx(self, g_exons: list[tuple[int, int]]) -> list[tuple[int, int]]:
        out = []
        for s, e in g_exons:
            a, b = self.interval_tx(s, e)
            out.append((a, b))
        out.sort()
        return out


def find_motif(seq: str, motif: str, t0: int, t1: int) -> list[int]:
    out, i, s = [], 0, seq[max(0, t0) : max(0, t1)]
    while True:
        j = s.find(motif, i)
        if j < 0:
            break
        out.append(max(0, t0) + j)
        i = j + 1
    return out


def acceptors(seq: str, t0: int, t1: int) -> list[int]:
    return [p + 2 for p in find_motif(seq, "AG", t0 - 2, t1) if 0 <= p + 2 < len(seq)]


def donors(seq: str, t0: int, t1: int) -> list[int]:
    return [p for p in find_motif(seq, "GT", t0, t1) if 0 <= p < len(seq)]


def parse_gtf(gtf: Path) -> list[dict]:
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
            if not tm:
                continue
            tid = tm.group(1)
            rec = tx.setdefault(
                tid,
                {"tid": tid, "chrom": chrom, "strand": strand, "exons": [], "cov": 0.0},
            )
            if feat == "exon":
                rec["exons"].append((start, end))
            cm = COV_RE.search(attrs)
            if cm:
                rec["cov"] = float(cm.group(1))
    for rec in tx.values():
        rec["exons"].sort()
        rec["start"] = rec["exons"][0][0]
        rec["end"] = rec["exons"][-1][1]
    return list(tx.values())


def tx_for_gene(txs: list[dict], g: dict) -> dict | None:
    hits = []
    for rec in txs:
        if rec["chrom"] != g["chrom"] or rec["strand"] != g["strand"]:
            continue
        ok = False
        if g["v_start"] and any(overlap(s, e, g["v_start"], g["v_end"]) for s, e in rec["exons"]):
            ok = True
        if g["i_start"] and any(overlap(s, e, g["i_start"], g["i_end"]) for s, e in rec["exons"]):
            ok = True
        if ok:
            hits.append(rec)
    if not hits:
        return None
    hits.sort(key=lambda r: (-r["cov"], -(r["end"] - r["start"])))
    return hits[0]


def ig_exon_cands(seq: str, t0: int, t1: int, min_len: int, max_len: int) -> list[tuple[int, int]]:
    acc = acceptors(seq, t0 - 150, t0 + 80)
    don = donors(seq, t1 - 80, t1 + 150)
    out = []
    seen = set()
    for a in acc:
        for d in don:
            L = d - a
            if min_len <= L <= max_len:
                key = (a, d)
                if key not in seen:
                    seen.add(key)
                    out.append(key)
    raw = (t0, t1)
    if raw not in seen and 80 <= (t1 - t0) <= max_len + 50:
        out.append(raw)
    out.sort(key=lambda x: abs((x[0] + x[1]) / 2 - (t0 + t1) / 2))
    return out[:12]


def match_8a(aa: str, bait: str) -> float:
    if not bait or len(bait) < 20:
        return 0.0
    hits = 0
    n = 0
    for i in range(0, len(bait) - 24, 20):
        n += 1
        if bait[i : i + 25] in aa:
            hits += 1
    return hits / max(n, 1)


def min_aa_for(arch: str) -> int:
    if arch == "V+I":
        return 160
    return 70


def genomic_exons(loc: Locus, pieces: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    out = []
    for a, b, k in pieces:
        if b <= a:
            continue
        gs, ge = loc.t2g(a), loc.t2g(b - 1)
        out.append((min(gs, ge), max(gs, ge), k))
    return out


def try_core_pair(seq: str, v_tx, i_tx, bait: str, arch: str):
    """Yield (v_exon, i_exon, core_nt) for in-frame V±I concatenations."""
    if v_tx and i_tx:
        vs, ve = v_tx
        i0, i1 = i_tx
        for va, vd in ig_exon_cands(seq, vs, ve, MIN_V, MAX_V):
            for ia, idn in ig_exon_cands(seq, i0, i1, MIN_I, MAX_I):
                if ia < vd + MIN_INTRON:
                    continue
                if seq[vd : vd + 2] != "GT" or (ia >= 2 and seq[ia - 2 : ia] != "AG"):
                    continue
                core = seq[va:vd] + seq[ia:idn]
                pep = translate(core)
                if "*" in pep or len(pep) < 120:
                    continue
                if pep.count("C") < 4:
                    continue
                if bait and match_8a(pep, bait) < 0.3 and len(bait) >= 40:
                    continue
                yield (va, vd), (ia, idn), core
    elif v_tx:
        vs, ve = v_tx
        for va, vd in ig_exon_cands(seq, vs, ve, MIN_V, MAX_V):
            core = seq[va:vd]
            pep = translate(core)
            if "*" in pep or len(pep) < 70:
                continue
            yield (va, vd), None, core
    elif i_tx:
        i0, i1 = i_tx
        for ia, idn in ig_exon_cands(seq, i0, i1, MIN_I, MAX_I):
            core = seq[ia:idn]
            pep = translate(core)
            if "*" in pep or len(pep) < 60:
                continue
            yield None, (ia, idn), core


def pick_sp_exon(seq: str, v0: int, core: str) -> tuple[int, int, str] | None:
    """Separate ATG…GT leader exon, intron, then AG at V. No intron Met fused to V."""
    best = None
    if v0 < 2 or seq[v0 - 2 : v0] != "AG":
        return None
    for don in donors(seq, max(0, v0 - MAX_INTRON_5P), v0 - MIN_INTRON):
        intron = v0 - don
        if not (MIN_INTRON <= intron <= MAX_INTRON_5P):
            continue
        for atg in find_motif(seq, "ATG", max(0, don - MAX_SP), don - MIN_SP + 3):
            if not (MIN_SP <= don - atg <= MAX_SP):
                continue
            if (don - atg) % 3 != 0:
                continue
            sp_nt = seq[atg:don]
            aa = translate(sp_nt + core)
            if not aa.startswith("M") or "*" in aa[: max(80, len(core) // 3)]:
                continue
            sp_pep = translate(sp_nt).split("*")[0]
            if not is_sp(sp_pep):
                continue
            sc = 10 if 12 <= len(sp_pep) <= 22 else 0
            if best is None or sc >= best[0]:
                best = (sc, atg, don, sp_nt)
    if best:
        return best[1], best[2], best[3]
    return None


def internal_exons(seq: str, t0: int, t1: int, min_nt: int, max_nt: int) -> list[tuple[int, int]]:
    """GT–AG cassette exons. Length need not be a multiple of 3 (phase 1/2 splices)."""
    accs = acceptors(seq, t0, t1)
    dons = donors(seq, t0, t1)
    out = []
    for a in accs:
        if a < 2 or seq[a - 2 : a] != "AG":
            continue
        for d in dons:
            L = d - a
            if min_nt <= L <= max_nt and seq[d : d + 2] == "GT":
                out.append((a, d))
    return out


def continued_aa(prefix_nt: str, extra_nt: str) -> str:
    """Amino acids encoded when extra_nt is ligated onto prefix_nt (same ORF)."""
    aa = translate(prefix_nt + extra_nt)
    return aa[len(prefix_nt) // 3 :]


def has_internal_stop(aa: str) -> bool:
    core = aa[:-1] if aa.endswith("*") else aa
    return "*" in core


def exon_to_stop_framed(seq: str, acc: int, limit: int, prefix_nt: str) -> tuple[int, int] | None:
    """Last exon: AG acceptor, then read in the running frame until a stop codon."""
    if acc < 2 or seq[acc - 2 : acc] != "AG":
        return None
    end = min(len(seq), limit, acc + 480)
    extra = seq[acc:end]
    if len(extra) < 36:
        return None
    aa = translate(prefix_nt + extra)
    tail = aa[len(prefix_nt) // 3 :]
    if "*" not in tail:
        return None
    stop_i = tail.find("*")
    stop_nt_end = (len(prefix_nt) // 3 + stop_i + 1) * 3
    exon_end = stop_nt_end - len(prefix_nt)
    if not (36 <= exon_end <= 480):
        return None
    return acc, acc + exon_end


def extend_cassette(seq: str, i1: int, limit: int, core_nt: str) -> list[tuple[int, int, str]]:
    """0–3 extra GT–AG exons after I: TM, optional C1, optional C2 (to stop).

    Each exon continues the V/I reading frame. Do not walk unspliced nt to a stop.
    """
    limit = min(len(seq), i1 + MAX_3P, limit)
    search0 = i1 + MIN_INTRON
    running = core_nt
    tm_hits = []
    for a, d in internal_exons(seq, search0, limit, 36, 165):
        pep = continued_aa(running, seq[a:d])
        if has_internal_stop(pep):
            continue
        pep = pep.replace("*", "")
        if is_tm_pep(pep):
            tm_hits.append((a, d, pep))
    tm_hits.sort(key=lambda x: x[0])
    if not tm_hits:
        return []

    ta, td, _ = tm_hits[0]
    pieces = [(ta, td, "TM")]
    running = core_nt + seq[ta:td]
    after = td + MIN_INTRON

    c1 = None
    for a, d in internal_exons(seq, after, limit, 45, 90):
        pep = continued_aa(running, seq[a:d])
        if has_internal_stop(pep):
            continue
        pep = pep.replace("*", "")
        if 12 <= len(pep) <= 30 and not is_tm_pep(pep):
            c1 = (a, d)
            break
    if c1:
        pieces.append((c1[0], c1[1], "cyto"))
        running = running + seq[c1[0] : c1[1]]
        after = c1[1] + MIN_INTRON

    c2_accs = acceptors(seq, after, min(limit, after + MAX_INTRON_3P))
    best_c2 = None
    for a in c2_accs[:40]:
        hit = exon_to_stop_framed(seq, a, limit, running)
        if not hit:
            continue
        pep = continued_aa(running, seq[hit[0] : hit[1]]).split("*")[0]
        if len(pep) < 8:
            continue
        sc = (20 if ITIM_RE.search(pep) else 0) + len(pep)
        if best_c2 is None or sc > best_c2[0]:
            best_c2 = (sc, hit[0], hit[1])
    if best_c2:
        pieces.append((best_c2[1], best_c2[2], "cyto"))
    return pieces


def pack(g, notes, pep, cds, pieces, loc):
    labs = {k for *_, k in pieces}
    has_stop = bool(pep) and cds[-3:] in STOPS if len(cds) >= 3 else "*" in translate(cds + "NN")
    if pep.endswith("*"):
        pep = pep[:-1]
        has_stop = True
    tm = "TM" in labs
    complete = bool(
        pep
        and (g["arch"] != "V+I" or ("V" in labs and "I" in labs))
        and (g["arch"] != "V+I" or "TM" in labs)
        and pep.startswith("M")
        and has_stop
        and not low_complexity(pep)
    )
    return {
        **{k: g[k] for k in ("gene_id", "chrom", "strand", "arch", "v_start", "v_end", "i_start", "i_end")},
        "notes": ";".join(notes),
        "protein": pep,
        "cds": cds,
        "aa_len": len(pep),
        "has_M": "yes" if pep.startswith("M") else "no",
        "has_stop": "yes" if has_stop else "no",
        "has_TM": "yes" if "TM" in labs else "no",
        "has_SP": "yes" if "SP" in labs else "no",
        "n_exons": len(pieces),
        "exons": genomic_exons(loc, pieces),
        "complete": complete,
        "transcript_id": next((n[4:] for n in notes if n.startswith("RNA:")), "-"),
    }


def empty(g, notes):
    return {
        **{k: g[k] for k in ("gene_id", "chrom", "strand", "arch", "v_start", "v_end", "i_start", "i_end")},
        "notes": ";".join(notes),
        "protein": "",
        "cds": "",
        "aa_len": 0,
        "has_M": "no",
        "has_stop": "no",
        "has_TM": "no",
        "has_SP": "no",
        "n_exons": 0,
        "exons": [],
        "complete": False,
        "transcript_id": "-",
    }


def pieces_to_cds(seq: str, pieces: list[tuple[int, int, str]]) -> tuple[str, str]:
    cds = "".join(seq[a:b] for a, b, _ in pieces)
    aa = translate(cds)
    has_stop = "*" in aa
    pep = aa.split("*")[0]
    if has_stop:
        cds = cds[: (len(pep) + 1) * 3]
    return pep, cds


def stitch_cassette(loc: Locus, g: dict) -> dict:
    seq = loc.seq
    bait = g.get("protein_8a") or ""
    notes = [g["arch"], "genome_cassette"]
    v_tx = loc.interval_tx(g["v_start"], g["v_end"]) if g["v_start"] else None
    i_tx = loc.interval_tx(g["i_start"], g["i_end"]) if g["i_start"] else None
    if v_tx is None and i_tx is None:
        notes.append("no_anchor")
        return empty(g, notes)

    best_core = None
    n_core = 0
    for v_ex, i_ex, core in try_core_pair(seq, v_tx, i_tx, bait, g["arch"]):
        n_core += 1
        if n_core > 8:
            break
        score = match_8a(translate(core), bait) * 100 + len(core)
        if best_core is None or score > best_core[0]:
            best_core = (score, v_ex, i_ex, core)

    if best_core is None:
        pieces = []
        if v_tx:
            pieces.append((v_tx[0], v_tx[1], "V"))
        if i_tx:
            pieces.append((i_tx[0], i_tx[1], "I"))
        notes.append("raw_HMM_exons")
        pep, cds = pieces_to_cds(seq, pieces)
        return pack(g, notes, pep, cds, pieces, loc)

    _sc, v_ex, i_ex, core = best_core
    pieces: list[tuple[int, int, str]] = []
    if v_ex:
        pieces.append((v_ex[0], v_ex[1], "V"))
    if i_ex:
        pieces.append((i_ex[0], i_ex[1], "I"))
    # Without RNA, do not search the genome for SP/TM/cyto. Those exons are
    # only taken from spliced transcripts. GT–AG here is only used to snap
    # V/I to true splice sites.
    notes.append("no_SP")
    notes.append("no_TM")
    notes.append("partial")
    pep, cds = pieces_to_cds(seq, pieces)
    return pack(g, notes, pep, cds, pieces, loc)


def label_tx_pieces(
    loc: Locus,
    g: dict,
    exons_tx: list[tuple[int, int]],
    pep: str,
    cds_from: int,
    cds_to: int,
) -> list[tuple[int, int, str]]:
    v_box = loc.interval_tx(g["v_start"], g["v_end"]) if g["v_start"] else None
    i_box = loc.interval_tx(g["i_start"], g["i_end"]) if g["i_start"] else None
    seq = loc.seq
    clipped = []
    tpos = 0
    for a, b in exons_tx:
        elen = b - a
        if tpos + elen <= cds_from or tpos >= cds_to:
            tpos += elen
            continue
        clip0 = max(0, cds_from - tpos)
        clip1 = min(elen, cds_to - tpos)
        clipped.append((a + clip0, a + clip1))
        tpos += elen
    if not clipped:
        return []

    tagged = []
    for a, b in clipped:
        lab = "?"
        if v_box and overlap(a, b - 1, v_box[0], v_box[1] - 1):
            lab = "V"
        elif i_box and overlap(a, b - 1, i_box[0], i_box[1] - 1):
            lab = "I"
        tagged.append([a, b, lab])

    vi = [i for i, t in enumerate(tagged) if t[2] in ("V", "I")]
    if not vi:
        return [(a, b, "V" if v_box else "I") for a, b in clipped]

    for i in range(vi[0]):
        a, b, _ = tagged[i]
        pep_e = translate(seq[a:b]).split("*")[0]
        tagged[i][2] = (
            "SP"
            if is_sp(pep_e) or (pep_e.startswith("M") and 8 <= len(pep_e) <= 35)
            else "?"
        )
    tm_done = False
    for i in range(vi[-1] + 1, len(tagged)):
        if not tm_done:
            tagged[i][2] = "TM"
            tm_done = True
        else:
            tagged[i][2] = "cyto"

    return [(a, b, k) for a, b, k in tagged if k != "?" and (b - a) >= 18]


TRIMS_5 = tuple(range(0, 19))


def fit_rna_cds(
    seq: str,
    exons_tx: list[tuple[int, int]],
    v_box: tuple[int, int] | None,
    i_box: tuple[int, int] | None,
    bait: str,
) -> list[tuple[int, int]] | None:
    """Trim 5' overhangs on StringTie exons so SP–V–I–TM concat is one ORF.

    StringTie often keeps extra acceptor nt that insert stops between SP and V.
    """
    def overlaps_box(a, b, box):
        return box is not None and overlap(a, b - 1, box[0], box[1] - 1)

    v_idx = [i for i, (a, b) in enumerate(exons_tx) if overlaps_box(a, b, v_box)]
    i_idx = [i for i, (a, b) in enumerate(exons_tx) if overlaps_box(a, b, i_box)]
    if not v_idx and not i_idx:
        return None
    lo = min((v_idx or i_idx) + (i_idx or v_idx))
    hi = max((v_idx or i_idx) + (i_idx or v_idx))
    idxs = list(range(lo, hi + 1))
    if len(idxs) > 4:
        idxs = [idxs[0], idxs[-1]] if len(idxs) > 2 else idxs

    best = None
    for trims in product(TRIMS_5, repeat=len(idxs)):
        parts = []
        chunks = []
        ok = True
        for idx, t in zip(idxs, trims):
            a, b = exons_tx[idx]
            if b - a - t < 60 and idx in set(v_idx + i_idx):
                ok = False
                break
            if b - a - t < 18:
                ok = False
                break
            parts.append((a + t, b))
            chunks.append(seq[a + t : b])
        if not ok:
            continue
        core = "".join(chunks)
        pep = translate(core)
        if has_internal_stop(pep):
            continue
        pep = pep.replace("*", "")
        if len(pep) < 60:
            continue
        sc = 100 * match_8a(pep, bait) + len(pep) - 0.1 * sum(trims)
        if best is None or sc > best[0]:
            best = (sc, parts, core)

    if best is None:
        return None

    _sc, pieces, running = best

    for idx in range(lo - 1, -1, -1):
        a, b = exons_tx[idx]
        placed = None
        for t in TRIMS_5:
            if b - a - t < 21:
                continue
            extra = seq[a + t : b]
            pep = translate(extra + running)
            if has_internal_stop(pep):
                continue
            head = translate(extra).split("*")[0]
            sc = (40 if pep.startswith("M") else 0) + (
                20 if is_sp(head) or (head.startswith("M") and 8 <= len(head) <= 35) else 0
            )
            if placed is None or sc >= placed[0]:
                placed = (sc, a + t, b, extra)
        if placed is None:
            break
        pieces.insert(0, (placed[1], placed[2]))
        running = placed[3] + running

    for idx in range(hi + 1, len(exons_tx)):
        a, b = exons_tx[idx]
        last = idx == len(exons_tx) - 1
        placed = None
        for t in TRIMS_5:
            if b - a - t < 18:
                continue
            extra = seq[a + t : b]
            pep = continued_aa(running, extra)
            if last:
                pep0 = pep.split("*")[0]
                if len(pep0) < 8:
                    continue
                if "*" in pep:
                    hit = exon_to_stop_framed(seq, a + t, b + 3, running)
                    if hit:
                        extra = seq[hit[0] : hit[1]]
                        a = hit[0]
                        b = hit[1]
                        t = 0
                sc = len(pep0) + (15 if ITIM_RE.search(pep0) else 0)
            else:
                if has_internal_stop(pep):
                    continue
                sc = len(pep.replace("*", ""))
            if placed is None or sc >= placed[0]:
                placed = (sc, a + t, b, extra)
        if placed is None:
            break
        pieces.append((placed[1], placed[2]))
        running = running + placed[3]
    return pieces


def stitch_from_tx(loc: Locus, g: dict, rec: dict) -> dict | None:
    exons_tx = loc.exons_tx(rec["exons"])
    if not exons_tx:
        return None
    seq = loc.seq
    bait = g.get("protein_8a") or ""
    v_tx = loc.interval_tx(g["v_start"], g["v_end"]) if g["v_start"] else None
    i_tx = loc.interval_tx(g["i_start"], g["i_end"]) if g["i_start"] else None
    if v_tx is None and i_tx is None:
        return None

    clipped = fit_rna_cds(seq, exons_tx, v_tx, i_tx, bait)
    if not clipped:
        return None
    cds = "".join(seq[a:b] for a, b in clipped)
    pep = translate(cds).split("*")[0]
    if len(pep) < 60:
        return None
    tlen = sum(b - a for a, b in clipped)
    pieces = label_tx_pieces(loc, g, clipped, pep, 0, tlen)
    if not pieces:
        return None
    cds = "".join(seq[a:b] for a, b, _ in pieces)
    pep = translate(cds).split("*")[0]
    notes = [g["arch"], "from_transcript", f"RNA:{rec['tid']}"]
    labs = {k for *_, k in pieces}
    if "SP" not in labs:
        notes.append("no_SP")
    if "TM" not in labs:
        notes.append("no_TM")
    if "V" in labs and "I" in labs and "TM" in labs:
        notes.append("RNA_full")
    else:
        notes.append("partial")
    rec_out = pack(g, notes, pep, cds, pieces, loc)
    rec_out["transcript_id"] = rec["tid"]
    return rec_out


def load_genes(summary: Path) -> list[dict]:
    rows = list(csv.DictReader(summary.open(), delimiter="\t"))
    genes = []
    for r in rows:
        def coord(k):
            v = r.get(k, "")
            return int(v) if v else None
        genes.append(
            {
                "gene_id": r["gene_id"],
                "chrom": r["chrom"],
                "strand": r["strand"],
                "arch": r["arch"],
                "v_start": coord("v_start"),
                "v_end": coord("v_end"),
                "i_start": coord("i_start"),
                "i_end": coord("i_end"),
                "protein_8a": (r.get("protein") or "").strip(),
            }
        )
    genes.sort(key=lambda g: (g["chrom"], min(x for x in [g["v_start"], g["i_start"]] if x)))
    return genes


def _span(g: dict) -> tuple[int, int]:
    if g.get("v_start") or g.get("i_start"):
        a = [x for x in [g["v_start"], g["v_end"], g["i_start"], g["i_end"]] if x]
        return min(a), max(a)
    return int(g["start"]), int(g["end"])


def dist_to_ig(ex_s: int, ex_e: int, lo: int, hi: int) -> int:
    """0 if the exon overlaps the Ig span, else the gap in bp."""
    if overlap(ex_s, ex_e, lo, hi):
        return 0
    if ex_e < lo:
        return lo - ex_e
    return ex_s - hi


def _within_flank(ex_s: int, ex_e: int, lo: int, hi: int, strand: str) -> bool:
    d = dist_to_ig(ex_s, ex_e, lo, hi)
    if d == 0:
        return True
    if strand == "+":
        if ex_e < lo:
            return d <= MAX_5P
        if ex_s > hi:
            return d <= MAX_3P
        return True
    if ex_s > hi:
        return d <= MAX_5P
    if ex_e < lo:
        return d <= MAX_3P
    return True


def owner_gene(ex_s: int, ex_e: int, genes_same: list[dict]) -> dict | None:
    """Assign an exon to exactly one same-strand gene.

    Overlap with a gene's V/I span wins. In the gap between two genes (e.g.
    after NITR4's I and before NITR5's V), the nearer Ig edge wins: closer to
    the downstream V → that gene's SP; closer to the upstream I → TM/CT.
    """
    scored: list[tuple[int, int, int, dict]] = []
    for g in genes_same:
        lo, hi = _span(g)
        if not _within_flank(ex_s, ex_e, lo, hi, g["strand"]):
            continue
        d = dist_to_ig(ex_s, ex_e, lo, hi)
        ov = min(ex_e, hi) - max(ex_s, lo) + 1 if overlap(ex_s, ex_e, lo, hi) else 0
        scored.append((d, -ov, lo, g))
    if not scored:
        return None
    scored.sort(key=lambda x: (x[0], x[1], x[2]))
    return scored[0][3]


def exons_owned_by(exons: list[tuple[int, int]], gene: dict, genes: list[dict]) -> list[tuple[int, int]]:
    same = [x for x in genes if x["chrom"] == gene["chrom"] and x["strand"] == gene["strand"]]
    return [(s, e) for s, e in exons if owner_gene(s, e, same) is gene]


def neighbor_bounds(genes: list[dict], i: int) -> tuple[int, int]:
    g = genes[i]
    lo, hi = _span(g)
    if g["strand"] == "+":
        lo, hi = lo - MAX_5P, hi + MAX_3P
    else:
        lo, hi = lo - MAX_3P, hi + MAX_5P
    same = [x for x in genes if x["chrom"] == g["chrom"] and x["strand"] == g["strand"]]
    idx = same.index(g)
    if idx > 0:
        phi = _span(same[idx - 1])[1]
        lo = max(lo, (phi + _span(g)[0]) // 2)
    if idx + 1 < len(same):
        nlo = _span(same[idx + 1])[0]
        hi = min(hi, (_span(g)[1] + nlo) // 2)
    return lo, hi


def stitch_one(genome: Fasta, genes: list[dict], i: int, txs: list[dict]) -> dict:
    g = genes[i]
    w0, w1 = neighbor_bounds(genes, i)
    w0 = max(1, w0)
    rec = tx_for_gene(txs, g) if txs else None
    if rec is not None:
        rec = dict(rec)
        owned = exons_owned_by(rec["exons"], g, genes)
        rec["exons"] = owned
        if owned:
            w0 = min(w0, min(s for s, _e in owned))
            w1 = max(w1, max(e for _s, e in owned))
        hits_ig = False
        if g["v_start"]:
            hits_ig = hits_ig or any(overlap(s, e, g["v_start"], g["v_end"]) for s, e in rec["exons"])
        if g["i_start"]:
            hits_ig = hits_ig or any(overlap(s, e, g["i_start"], g["i_end"]) for s, e in rec["exons"])
        if hits_ig and rec["exons"]:
            loc = Locus(genome, g["chrom"], max(1, w0), w1, g["strand"])
            got = stitch_from_tx(loc, g, rec)
            if got and got["aa_len"] >= 60:
                return got
    loc = Locus(genome, g["chrom"], max(1, w0), w1, g["strand"])
    return stitch_cassette(loc, g)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genome", type=Path, required=True)
    p.add_argument("--summary", type=Path, required=True, help="step 8 nitr_loci_summary.tsv")
    p.add_argument("--gtf", type=Path, help="StringTie cluster transcripts (preferred exons)")
    p.add_argument("-o", "--prefix", type=Path, required=True)
    args = p.parse_args()

    genome = Fasta(str(args.genome), as_raw=True, sequence_always_upper=True)
    genes = load_genes(args.summary)
    txs = parse_gtf(args.gtf) if args.gtf and args.gtf.exists() else []
    args.prefix.parent.mkdir(parents=True, exist_ok=True)

    results = [stitch_one(genome, genes, i, txs) for i in range(len(genes))]

    fa = Path(f"{args.prefix}_proteins.fa")
    cds_fa = Path(f"{args.prefix}_cds.fa")
    tsv = Path(f"{args.prefix}_summary.tsv")
    gtf = Path(f"{args.prefix}_models.gtf")
    fa.write_text("")
    cds_fa.write_text("")
    gtf_lines = []

    cols = [
        "gene_id", "chrom", "strand", "arch", "aa_len",
        "has_SP", "has_M", "has_stop", "has_TM", "n_exons",
        "complete", "transcript_id", "notes", "protein",
    ]
    with tsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow({**r, "complete": "yes" if r["complete"] else "no"})
            if r["protein"]:
                hdr = (
                    f"{r['gene_id']} {r['chrom']}({r['strand']}) arch={r['arch']} "
                    f"aa={r['aa_len']} M={r['has_M']} stop={r['has_stop']} TM={r['has_TM']} "
                    f"notes={r['notes']}"
                )
                with fa.open("a") as out:
                    out.write(f">{hdr}\n{wrap(r['protein'])}\n")
                with cds_fa.open("a") as out:
                    out.write(f">{hdr}\n{wrap(r['cds'])}\n")
            gid = r["gene_id"]
            for s, e, k in r["exons"]:
                attr = f'gene_id "{gid}"; transcript_id "{gid}.1"; domain "{k}";'
                gtf_lines.append(
                    f"{r['chrom']}\tnitr8b\texon\t{s}\t{e}\t.\t{r['strand']}\t.\t{attr}\n"
                )
    gtf.write_text("".join(gtf_lines))

    n_ok = sum(1 for r in results if r["complete"])
    print(f"Genes: {len(results)}  cassette_full: {n_ok}  transcripts_loaded: {len(txs)}")
    print(f"  {tsv}")
    print(f"  {fa}")
    print(f"  {cds_fa}")
    print(f"  {gtf}")
    for r in results:
        flag = "full" if r["complete"] else "partial"
        print(
            f"  {r['gene_id']:12} {r['arch']:7} {flag:8} aa={r['aa_len']:<4} "
            f"SP={r['has_SP']} TM={r['has_TM']}  {r['notes']}"
        )
        if r["protein"]:
            print(f"    {r['protein'][:80]}{'...' if len(r['protein'])>80 else ''}")


if __name__ == "__main__":
    main()
