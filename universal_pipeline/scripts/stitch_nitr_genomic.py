#!/usr/bin/env python3
"""Step 8b: genomic stitch of full-length NITRs (Met → stop).

Ig V/I intervals from 8a are the anchors. StringTie exons are used when a
transcript uniquely overlaps the gene; otherwise exons are placed on GT–AG
around the HMM hits. The ORF is always translated from the ATG, so SP, V, I
and TM stay in one frame. A model is complete only if it starts with Met,
ends at stop, and is long enough to contain the Ig domains.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio.Seq import Seq
from pyfaidx import Fasta

MAX_5P = 8000
MAX_3P = 16000
MIN_INTRON = 50
MAX_INTRON_5P = 4000
MAX_INTRON_3P = 15000
MIN_V = 240
MAX_V = 430
MIN_I = 220
MAX_I = 430
MIN_SP = 36
MAX_SP = 150

TID_RE = re.compile(r'transcript_id "([^"]+)"')
GID_RE = re.compile(r'gene_id "([^"]+)"')
COV_RE = re.compile(r'(?:cov|coverage) "([^"]+)"')
REPEAT_AA = re.compile(r"(IY){5,}|(YI){5,}|I{8,}|Y{8,}|Q{8,}|X{6,}")

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


def hyd_score(aa: str) -> float:
    if not aa:
        return -99.0
    return sum(KD.get(a, 0.0) for a in aa) / len(aa)


def has_tm(aa: str, window: int = 18, thresh: float = 1.55) -> bool:
    if len(aa) < window:
        return False
    return any(hyd_score(aa[i : i + window]) >= thresh for i in range(len(aa) - window + 1))


def is_sp(pep: str) -> bool:
    if not pep.startswith("M") or "*" in pep:
        return False
    if not (12 <= len(pep) <= 48):
        return False
    hyd = set("AILMFVW")
    best = max(sum(1 for a in pep[i : i + 8] if a in hyd) for i in range(0, len(pep) - 7))
    n_charged = sum(1 for a in pep[:20] if a in "DEKR")
    return best >= 5 and n_charged <= 4


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
    """Exon starts: AG immediately upstream."""
    return [p + 2 for p in find_motif(seq, "AG", t0 - 2, t1) if 0 <= p + 2 < len(seq)]


def donors(seq: str, t0: int, t1: int) -> list[int]:
    """Exon ends: GT immediately downstream."""
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
    acc = acceptors(seq, t0 - 90, t0 + 60)
    don = donors(seq, t1 - 60, t1 + 90)
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
        return 200
    return 90


def score_protein(aa: str, arch: str, bait: str, tm: bool, has_stop: bool) -> float:
    if not aa.startswith("M") or low_complexity(aa):
        return -1e9
    if "*" in aa:
        return -1e9
    sc = len(aa)
    sc += 80 * match_8a(aa, bait)
    sc += 40 if has_stop else -30
    sc += 25 if is_sp(aa[: min(28, len(aa))]) else -15
    sc += 20 if tm else 0
    if arch == "V+I" and tm:
        sc += 15
    if arch == "V-only" and not tm:
        sc += 10  # secreted OK
    c = aa.count("C")
    if arch == "V+I" and c >= 5:
        sc += 20
    elif arch != "V+I" and c >= 2:
        sc += 10
    if len(aa) < min_aa_for(arch):
        sc -= 200
    return sc


def pick_sp(seq: str, v0: int, core: str, bait: str, arch: str) -> tuple[int, int, str] | None:
    """Return (atg, donor, sp_nt) such that SP+core is an in-frame Ig ORF."""
    best = None
    core_aa = len(core) // 3
    need = max(80, min(core_aa, 160))

    def consider(atg: int, don: int) -> None:
        nonlocal best
        if not (MIN_SP <= don - atg <= MAX_SP):
            return
        cds = seq[atg:don] + core
        aa = translate(cds)
        if not aa.startswith("M"):
            return
        if "*" in aa[:need]:
            return
        aa_ns = aa.split("*")[0]
        if not is_sp(aa_ns[: min(28, len(aa_ns))]):
            return
        if low_complexity(aa_ns[:40]):
            return
        sc = len(aa_ns) + 80 * match_8a(aa_ns, bait)
        sc += 10 if 15 <= (don - atg) // 3 <= 30 else 0
        if best is None or sc > best[0]:
            best = (sc, atg, don, seq[atg:don])

    # spliced SP exon
    if v0 >= 2 and seq[v0 - 2 : v0] == "AG":
        for don in donors(seq, max(0, v0 - MAX_INTRON_5P), v0 - MIN_INTRON):
            intron = v0 - don
            if not (MIN_INTRON <= intron <= MAX_INTRON_5P):
                continue
            for atg in find_motif(seq, "ATG", max(0, don - MAX_SP), don - MIN_SP + 3):
                consider(atg, don)
    # ATG on the V/I exon itself (SP fused to Ig)
    for atg in find_motif(seq, "ATG", max(0, v0 - 90), v0 + 30):
        consider(atg, v0)
    if best:
        return best[1], best[2], best[3]
    return None


STOPS = {"TAA", "TAG", "TGA"}


def extend_3p(seq: str, cds0: str, i1: int, limit: int, arch: str, bait: str) -> tuple[list[tuple[int, int]], str, bool, bool]:
    """Append 0–1 extra exons after last Ig, in the ATG frame."""
    limit = min(len(seq), i1 + MAX_3P, limit)
    core_aa = len(translate(cds0).split("*")[0])
    cands: list[tuple[float, list[tuple[int, int]], str, bool, bool]] = []

    def consider(extra_nt: str, exons: list[tuple[int, int]]) -> None:
        aa = translate(cds0 + extra_nt)
        if "*" not in aa:
            pep, has_stop = aa, False
        else:
            stop = aa.index("*")
            if stop < max(core_aa - 8, min_aa_for(arch) - 30):
                return
            pep, has_stop = aa[:stop], True
        tail = pep[core_aa:]
        if low_complexity(pep) or low_complexity(tail):
            return
        if has_stop and not (8 <= len(tail) <= 140):
            return
        tm = has_tm(tail) or has_tm(pep[-40:])
        sc = score_protein(pep, arch, bait, tm, has_stop)
        cands.append((sc, exons, pep, has_stop, tm))

    consider(seq[i1 : min(limit, i1 + 300)], [(i1, min(limit, i1 + 300))])
    hits: list[tuple[int, int, str, int, int]] = []
    for don in donors(seq, i1, min(limit, i1 + 180))[:8]:
        pre = seq[i1:don]
        if not (0 <= len(pre) <= 180):
            continue
        P = len(cds0) + len(pre)
        need = (3 - (P % 3)) % 3
        search1 = min(limit - 3, don + MAX_INTRON_3P + 210)
        for s in range(don + MIN_INTRON + 2, search1 + 1):
            if seq[s : s + 3] not in STOPS:
                continue
            L0 = 60 + (need - 60 % 3) % 3
            for L in range(L0, 181, 3):
                acc = s - L
                if acc < don + MIN_INTRON or acc < 2:
                    continue
                if seq[acc - 2 : acc] != "AG":
                    continue
                intron = acc - don
                rank = abs(L - 110) + (0 if 400 <= intron <= 14000 else 800)
                hits.append((rank, don, pre, acc, s))
    hits.sort()
    for _rank, don, pre, acc, s in hits[:30]:
        extra = pre + seq[acc : s + 3]
        exons = [(i1, don), (acc, s + 3)] if pre else [(acc, s + 3)]
        consider(extra, exons)
    if not cands:
        pep = translate(cds0 + seq[i1:limit]).split("*")[0]
        return [(i1, min(limit, i1 + max(3, (len(pep) - core_aa) * 3)))], pep, False, has_tm(pep[core_aa:])
    cands.sort(reverse=True)
    _sc, exons, pep, has_stop, tm = cands[0]
    return exons, pep, has_stop, tm


def genomic_exons(loc: Locus, pieces: list[tuple[int, int, str]]) -> list[tuple[int, int, str]]:
    out = []
    for a, b, k in sorted(pieces):
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


def stitch_from_tx(loc: Locus, g: dict, rec: dict) -> dict | None:
    """Complete an 8a/StringTie transcript: Met in SP, keep Ig, extend to stop."""
    exons_tx = loc.exons_tx(rec["exons"])
    if not exons_tx:
        return None
    seq = loc.seq
    bait = g.get("protein_8a") or ""
    nt = "".join(seq[a:b] for a, b in exons_tx)
    # map tx-exon concat coords back
    best = None
    search_end = min(len(nt), 500)
    for atg in find_motif(nt, "ATG", 0, search_end):
        aa = translate(nt[atg:])
        if not aa.startswith("M"):
            continue
        has_stop = "*" in aa
        pep = aa.split("*")[0]
        if len(pep) < 80:
            continue
        if not is_sp(pep[: min(28, len(pep))]) and match_8a(pep, bait) < 0.5:
            continue
        tm = has_tm(pep)
        sc = score_protein(pep, g["arch"], bait, tm, has_stop)
        if best is None or sc > best[0]:
            best = (sc, atg, pep, has_stop, tm, nt[atg : atg + (len(pep) + int(has_stop)) * 3])
    if best is None:
        return None
    _sc, atg, pep, has_stop, tm, cds = best
    # if transcript is 5'-incomplete (no SP), caller will genomic-extend
    if not pep.startswith("M") or (g["arch"] == "V+I" and len(pep) < 180 and match_8a(pep, bait) < 0.5):
        return None
    pieces = []
    tpos = 0
    cds_from = atg
    cds_to = atg + len(cds)
    for a, b in exons_tx:
        elen = b - a
        if tpos + elen <= cds_from or tpos >= cds_to:
            tpos += elen
            continue
        clip0 = max(0, cds_from - tpos)
        clip1 = min(elen, cds_to - tpos)
        pieces.append((a + clip0, a + clip1, "tx"))
        tpos += elen
    notes = [g["arch"], "from_transcript"]
    if is_sp(pep[: min(28, len(pep))]):
        notes.append("SP_genomic")
    notes.append("stop_genomic" if has_stop else "no_stop")
    notes.append("has_TM" if tm else "no_TM")
    complete = bool(
        pep.startswith("M")
        and has_stop
        and len(pep) >= min_aa_for(g["arch"])
        and not low_complexity(pep)
        and (not bait or match_8a(pep, bait) >= 0.4 or len(bait) < 30)
    )
    notes.append("Met_to_stop" if complete else "partial")
    return pack(g, notes, pep, cds, tm, has_stop, pieces, loc, complete)


def pack(g, notes, pep, cds, tm, has_stop, pieces, loc, complete):
    return {
        **{k: g[k] for k in ("gene_id", "chrom", "strand", "arch", "v_start", "v_end", "i_start", "i_end")},
        "notes": ";".join(notes),
        "protein": pep,
        "cds": cds,
        "aa_len": len(pep),
        "has_M": "yes" if pep.startswith("M") else "no",
        "has_stop": "yes" if has_stop else "no",
        "has_TM": "yes" if tm else "no",
        "has_SP": "yes" if "SP_genomic" in notes or (pep.startswith("M") and is_sp(pep[: min(28, len(pep))])) else "no",
        "n_exons": len(pieces),
        "exons": genomic_exons(loc, pieces),
        "complete": complete,
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
    }


def stitch_genomic(loc: Locus, g: dict) -> dict:
    seq = loc.seq
    bait = g.get("protein_8a") or ""
    notes = [g["arch"], "genomic_stitch"]
    v_tx = loc.interval_tx(g["v_start"], g["v_end"]) if g["v_start"] else None
    i_tx = loc.interval_tx(g["i_start"], g["i_end"]) if g["i_start"] else None
    if v_tx is None and i_tx is None:
        notes.append("no_anchor")
        return empty(g, notes)

    best = None
    n_core = 0
    for v_ex, i_ex, core in try_core_pair(seq, v_tx, i_tx, bait, g["arch"]):
        n_core += 1
        if n_core > 8:
            break
        first0 = (v_ex or i_ex)[0]
        last1 = (i_ex or v_ex)[1]
        sp = pick_sp(seq, first0, core, bait, g["arch"])
        if sp:
            atg, don, sp_nt = sp
            cds0 = sp_nt + core
            notes_sp = "SP_genomic"
        else:
            atg = don = None
            cds0 = core
            notes_sp = "no_SP"
        exons3, pep, has_stop, tm = extend_3p(seq, cds0, last1, len(seq), g["arch"], bait)
        extra = "".join(seq[a:b] for a, b in exons3)
        cds = cds0 + extra
        aa = translate(cds)
        has_stop = "*" in aa
        pep = aa.split("*")[0]
        if "M" in pep:
            m = pep.find("M")
            if m > 0 and is_sp(pep[m : m + 28] if len(pep) - m >= 12 else pep[m:]):
                pep = pep[m:]
                cds = cds[m * 3 :]
        tm = has_tm(pep)
        sc = score_protein(pep, g["arch"], bait, tm, has_stop)
        pieces = []
        if atg is not None:
            pieces.append((atg, don, "SP"))
        if v_ex:
            pieces.append((v_ex[0], v_ex[1], "V"))
        if i_ex:
            pieces.append((i_ex[0], i_ex[1], "I"))
        for a, b in exons3:
            pieces.append((a, b, "3p"))
        if best is None or sc > best[0]:
            best = (sc, pep, cds, tm, has_stop, pieces, notes_sp)

    if best is None:
        notes.append("no_inframe_core")
        return empty(g, notes)

    _sc, pep, cds, tm, has_stop, pieces, notes_sp = best
    notes.append(notes_sp)
    notes.append("stop_genomic" if has_stop else "no_stop")
    notes.append("has_TM" if tm else "no_TM")
    complete = bool(
        pep.startswith("M")
        and has_stop
        and len(pep) >= min_aa_for(g["arch"])
        and not low_complexity(pep)
        and (not bait or match_8a(pep, bait) >= 0.4 or len(bait) < 30)
    )
    notes.append("Met_to_stop" if complete else "partial")
    return pack(g, notes, pep, cds, tm, has_stop, pieces, loc, complete)


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
    a = [x for x in [g["v_start"], g["v_end"], g["i_start"], g["i_end"]] if x]
    return min(a), max(a)


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
    loc = Locus(genome, g["chrom"], w0, w1, g["strand"])
    rec = tx_for_gene(txs, g) if txs else None
    got = None
    if rec is not None:
        rec = dict(rec)
        rec["exons"] = [(s, e) for s, e in rec["exons"] if s <= w1 and e >= w0]
        covers_v = True
        if g["v_start"]:
            covers_v = any(overlap(s, e, g["v_start"], g["v_end"]) for s, e in rec["exons"])
        if covers_v and rec["exons"]:
            got = stitch_from_tx(loc, g, rec)
            if got and got["complete"]:
                return got
    geo = stitch_genomic(loc, g)
    if got and geo:
        if geo["complete"] and not got["complete"]:
            return geo
        if got["complete"] or got["aa_len"] >= geo["aa_len"]:
            return got
        return geo
    return got or geo


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--genome", type=Path, required=True)
    p.add_argument("--summary", type=Path, required=True, help="8a nitr_summary.tsv")
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
        "complete", "notes", "protein",
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
    print(f"Genes: {len(results)}  Met-to-stop: {n_ok}  transcripts_loaded: {len(txs)}")
    print(f"  {tsv}")
    print(f"  {fa}")
    print(f"  {cds_fa}")
    print(f"  {gtf}")
    for r in results:
        flag = "OK" if r["complete"] else "partial"
        print(
            f"  {r['gene_id']:12} {r['arch']:7} {flag:8} aa={r['aa_len']:<4} "
            f"M={r['has_M']} stop={r['has_stop']} TM={r['has_TM']}  {r['notes']}"
        )
        if r["protein"]:
            print(f"    {r['protein'][:80]}{'...' if len(r['protein'])>80 else ''}")


if __name__ == "__main__":
    main()
