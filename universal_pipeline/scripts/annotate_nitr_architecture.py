#!/usr/bin/env python3
"""Step 9: label NITR domain architecture (SP, V/D1, I/D2, TM, cyto).

Ig domains come from InterProScan + SMART hmmscan.
SP comes from SignalP (GFF + protein_type summary) when provided.
TM and cyto come from DeepTMHMM (TMRs.gff3) when provided.
Heuristics are used only if those files are absent.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

ITIM_RE = re.compile(r"[IVL].Y..[LV]")
ITSM_RE = re.compile(r"T.Y..[VIL]")
IG_NAME_RE = re.compile(r"ig|I-set|Ig|V-set|IG|immunoglobulin", re.I)
FALSE_IG_RE = re.compile(r"^(Lig|PIG|LIGAN|BRIGHT)", re.I)

KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


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


def hyd(aa: str) -> float:
    if not aa:
        return -99.0
    return sum(KD.get(a, 0.0) for a in aa) / len(aa)


def is_ig_sig(name: str, desc: str = "") -> bool:
    blob = f"{name} {desc}"
    if FALSE_IG_RE.search(name):
        return False
    return bool(IG_NAME_RE.search(blob))


def call_sp(aa: str) -> tuple[int, int] | None:
    """Heuristic SP: Met + hydrophobic h-region in the first 12–40 aa. 1-based."""
    if not aa.startswith("M") or len(aa) < 16:
        return None
    window = aa[: min(40, len(aa))]
    if "*" in window:
        return None
    hydset = set("AILMFVW")
    best_i, best_n = 0, 0
    for i in range(0, max(1, len(window) - 7)):
        n = sum(1 for a in window[i : i + 8] if a in hydset)
        if n > best_n:
            best_i, best_n = i, n
    n_charged = sum(1 for a in window[:20] if a in "DEKR")
    if best_n < 5 or n_charged > 4 or best_i > 12:
        return None
    end = min(len(window), max(16, best_i + 12))
    while end < min(len(window), 35) and window[end - 1] not in "AGS":
        end += 1
    return 1, end


def call_tm(aa: str, after: int) -> tuple[int, int] | None:
    """Best KD helix after `after` (0-based). Returns 1-based inclusive."""
    w, thresh = 18, 1.55
    start = max(0, after)
    if len(aa) - start < w:
        return None
    best = None
    for i in range(start, len(aa) - w + 1):
        sc = hyd(aa[i : i + w])
        if sc >= thresh and (best is None or sc > best[0]):
            best = (sc, i, i + w)
    if best is None:
        return None
    return best[1] + 1, best[2]


def parse_signalp(gff: Path | None, summary: Path | None) -> dict[str, dict]:
    """SignalP 5/6: SP only when Prediction is SP (not OTHER). Coords from GFF."""
    out: dict[str, dict] = {}
    if summary and summary.exists():
        with summary.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 2:
                    continue
                gid, pred = p[0], p[1]
                score = p[2] if len(p) > 2 else ""
                cs = p[4] if len(p) > 4 else ""
                rec = out.setdefault(gid, {})
                rec["pred"] = pred
                rec["score"] = score
                rec["cs"] = cs
                rec["is_sp"] = pred.upper().startswith("SP")
    if gff and gff.exists():
        with gff.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 5:
                    continue
                if p[2].lower() not in {"signal_peptide", "signal_peptide_region_of_n_terminal"}:
                    continue
                gid, start, end, score = p[0], int(p[3]), int(p[4]), p[5] if len(p) > 5 else ""
                rec = out.setdefault(gid, {})
                rec["start"], rec["end"] = start, end
                rec.setdefault("score", score)
                rec.setdefault("is_sp", True)
                rec.setdefault("pred", "SP")
    return out


def parse_deeptmhmm(gff: Path | None) -> dict[str, dict]:
    """DeepTMHMM TMRs.gff3: signal / outside / TMhelix / inside (1-based)."""
    out: dict[str, dict] = {}
    if not gff or not gff.exists():
        return out
    with gff.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#") or line.startswith("//"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            gid, kind, start, end = p[0], p[1].lower(), int(p[2]), int(p[3])
            rec = out.setdefault(
                gid, {"signal": None, "tm": None, "inside": None, "outside": []}
            )
            span = (start, end)
            if kind == "signal":
                rec["signal"] = span
            elif kind == "tmhelix":
                rec["tm"] = span
            elif kind == "inside":
                rec["inside"] = span
            elif kind == "outside":
                rec["outside"].append(span)
    return out


def parse_ips(path: Path) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = {}
    if not path or not path.exists():
        return hits
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            gene, analysis, acc, desc = p[0], p[3], p[4], p[5]
            try:
                start, end = int(p[6]), int(p[7])
            except ValueError:
                continue
            score = p[8]
            ipr = p[11] if len(p) > 11 else ""
            ipr_desc = p[12] if len(p) > 12 else ""
            hits.setdefault(gene, []).append(
                {
                    "source": f"IPS:{analysis}",
                    "acc": acc,
                    "desc": desc or ipr_desc,
                    "start": start,
                    "end": end,
                    "score": score,
                    "ipr": ipr,
                }
            )
    return hits


def parse_domtbl(path: Path) -> dict[str, list[dict]]:
    hits: dict[str, list[dict]] = {}
    if not path or not path.exists():
        return hits
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.split()
            if len(p) < 23:
                continue
            hmm, gene = p[0], p[3]
            try:
                ie = float(p[12])
                ali_from, ali_to = int(p[17]), int(p[18])
            except ValueError:
                continue
            if ie > 1e-4:
                continue
            hits.setdefault(gene, []).append(
                {
                    "source": "SMART_hmmscan",
                    "acc": hmm,
                    "desc": hmm,
                    "start": ali_from,
                    "end": ali_to,
                    "score": f"{ie:.2e}",
                    "ipr": "",
                }
            )
    return hits


def load_summary(path: Path) -> dict[str, dict]:
    rows = {}
    with path.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows[r["gene_id"]] = r
    return rows


def classify_ig(hits: list[dict], aa: str) -> list[dict]:
    ig = []
    for h in hits:
        if not is_ig_sig(h["acc"], h["desc"]):
            continue
        pep = aa[h["start"] - 1 : h["end"]]
        nC = pep.count("C")
        name = h["acc"].upper()
        if "IGV" in name or nC <= 2:
            cls = "V"
        elif nC >= 4:
            cls = "I"
        else:
            cls = "Ig"
        ig.append({**h, "cls": cls, "nC": nC, "pep": pep})
    ig.sort(key=lambda x: (x["start"], -x["end"]))
    # drop nested weaker hits
    kept = []
    for h in ig:
        if any(k["start"] <= h["start"] and h["end"] <= k["end"] and k is not h for k in kept):
            continue
        kept.append(h)
    # if two Ig along the protein, force 5' V / 3' I (canonical NITR)
    if len(kept) >= 2:
        kept = sorted(kept, key=lambda x: x["start"])
        kept[0]["cls"] = "V"
        kept[-1]["cls"] = "I"
        for mid in kept[1:-1]:
            mid["cls"] = "Ig"
    return kept


def architecture(parts: list[str]) -> str:
    return "-".join(parts) if parts else "none"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--proteins", type=Path, required=True)
    p.add_argument("--summary8b", type=Path, required=True)
    p.add_argument("--summary8a", type=Path, help="optional RNA notes")
    p.add_argument("--ips", type=Path, help="InterProScan TSV")
    p.add_argument("--domtblout", type=Path, help="hmmscan SMART domtblout")
    p.add_argument("--signalp-gff", type=Path, help="SignalP output.gff3")
    p.add_argument("--signalp-summary", type=Path, help="SignalP output_protein_type.txt")
    p.add_argument("--deeptmhmm-gff", type=Path, help="DeepTMHMM TMRs.gff3")
    p.add_argument("-o", "--outdir", type=Path, required=True)
    args = p.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    seqs = load_fasta(args.proteins)
    s8b = load_summary(args.summary8b)
    s8a = load_summary(args.summary8a) if args.summary8a and args.summary8a.exists() else {}
    ips = parse_ips(args.ips) if args.ips else {}
    hmm = parse_domtbl(args.domtblout) if args.domtblout else {}
    sigp = parse_signalp(args.signalp_gff, args.signalp_summary)
    dtm = parse_deeptmhmm(args.deeptmhmm_gff)
    use_sigp = bool(sigp)
    use_dtm = bool(dtm)

    gene_rows = []
    seg_rows = []

    for gid, aa in seqs.items():
        meta = s8b.get(gid, {})
        rna = s8a.get(gid, {})
        hits = ips.get(gid, []) + hmm.get(gid, [])
        igs = classify_ig(hits, aa)
        has_v = any(h["cls"] == "V" for h in igs)
        has_i = any(h["cls"] == "I" for h in igs)
        last_ig = max((h["end"] for h in igs), default=40)

        notes = []
        arch_bits = []

        sp = None
        sp_src, sp_acc, sp_score = "heuristic", "SP_hydrophobic", ""
        if use_sigp:
            srec = sigp.get(gid, {})
            if srec.get("is_sp") and "start" in srec:
                sp = (srec["start"], srec["end"])
                sp_src, sp_acc, sp_score = "SignalP", srec.get("pred", "SP"), srec.get("score", "")
            elif srec.get("is_sp") and srec.get("cs"):
                m = re.search(r"CS pos:\s*(\d+)-", srec["cs"])
                if m:
                    sp = (1, int(m.group(1)))
                    sp_src, sp_acc, sp_score = "SignalP", srec.get("pred", "SP"), srec.get("score", "")
            else:
                notes.append("no_SP_SignalP")
        else:
            sp = call_sp(aa)
            if not sp:
                notes.append("no_SP_heuristic")

        tm = None
        cyto_span = None
        tm_src, tm_acc, tm_score = "heuristic_KD", "KD_helix", ""
        if use_dtm:
            drec = dtm.get(gid, {})
            tm = drec.get("tm")
            cyto_span = drec.get("inside")
            if tm:
                tm_src, tm_acc = "DeepTMHMM", "TMhelix"
                tm_score = ""
            else:
                notes.append("no_TM_DeepTMHMM")
            dsig = drec.get("signal")
            if dsig and not sp:
                notes.append(f"DeepTMHMM_signal_unconfirmed:{dsig[0]}-{dsig[1]}")
        else:
            tm = call_tm(aa, after=last_ig)
            if tm:
                tm_score = f"{hyd(aa[tm[0] - 1 : tm[1]]):.2f}"
            else:
                notes.append("no_TM_heuristic")

        if sp:
            arch_bits.append("SP")
            seg_rows.append(
                {
                    "gene_id": gid, "domain": "SP", "aa_start": sp[0], "aa_end": sp[1],
                    "source": sp_src, "acc": sp_acc, "score": sp_score,
                    "sequence": aa[sp[0] - 1 : sp[1]],
                }
            )

        for h in igs:
            label = "V/D1" if h["cls"] == "V" else ("I/D2" if h["cls"] == "I" else "Ig")
            if label not in arch_bits:
                arch_bits.append(label)
            seg_rows.append(
                {
                    "gene_id": gid, "domain": label, "aa_start": h["start"], "aa_end": h["end"],
                    "source": h["source"], "acc": h["acc"], "score": h["score"],
                    "sequence": h["pep"],
                }
            )

        if tm:
            arch_bits.append("TM")
            tm_pep = aa[tm[0] - 1 : tm[1]]
            if tm_src == "heuristic_KD" and not tm_score:
                tm_score = f"{hyd(tm_pep):.2f}"
            seg_rows.append(
                {
                    "gene_id": gid, "domain": "TM", "aa_start": tm[0], "aa_end": tm[1],
                    "source": tm_src, "acc": tm_acc, "score": tm_score,
                    "sequence": tm_pep,
                }
            )
            if cyto_span:
                c0, c1 = cyto_span
            else:
                c0, c1 = tm[1] + 1, len(aa)
            cyto = aa[c0 - 1 : c1]
            if cyto:
                arch_bits.append("cyto")
                itims = ",".join(ITIM_RE.findall(cyto))
                isms = ",".join(ITSM_RE.findall(cyto))
                note_cy = []
                if itims:
                    note_cy.append(f"ITIM:{itims}")
                if isms:
                    note_cy.append(f"ITSM:{isms}")
                seg_rows.append(
                    {
                        "gene_id": gid, "domain": "cyto", "aa_start": c0, "aa_end": c1,
                        "source": "DeepTMHMM" if use_dtm else "after_TM",
                        "acc": ";".join(note_cy) or "inside",
                        "score": "",
                        "sequence": cyto,
                    }
                )
        else:
            if igs:
                tail = aa[last_ig :]
                if 5 <= len(tail) <= 80:
                    arch_bits.append("secreted_tail")
                    seg_rows.append(
                        {
                            "gene_id": gid, "domain": "secreted_tail",
                            "aa_start": last_ig + 1, "aa_end": len(aa),
                            "source": "after_Ig", "acc": "", "score": "",
                            "sequence": tail,
                        }
                    )

        tx = rna.get("transcript_id", "-") or "-"
        if tx in {"", "-"}:
            notes.append("no_RNA")
        else:
            notes.append(f"RNA:{tx}")

        if not igs:
            notes.append("no_Ig_HMM")
            if (meta.get("arch") or "").startswith("I"):
                notes.append("possible_decaying_I")
        elif not has_v and has_i:
            notes.append("I_only")
        elif has_v and not has_i:
            notes.append("V_only")

        if cyto_span:
            cyto_aa = aa[cyto_span[0] - 1 : cyto_span[1]]
        elif tm:
            cyto_aa = aa[tm[1] :]
        else:
            cyto_aa = ""
        itim = bool(ITIM_RE.search(cyto_aa))
        itsm = bool(ITSM_RE.search(cyto_aa))

        gene_rows.append(
            {
                "gene_id": gid,
                "chrom": meta.get("chrom", ""),
                "strand": meta.get("strand", ""),
                "aa_len": len(aa),
                "arch_8b": meta.get("arch", ""),
                "architecture": architecture(arch_bits),
                "has_SP": "yes" if sp else "no",
                "has_V": "yes" if has_v else "no",
                "has_I": "yes" if has_i else "no",
                "has_TM": "yes" if tm else "no",
                "n_Ig_hits": len(igs),
                "ITIM": "yes" if itim else "no",
                "ITSM": "yes" if itsm else "no",
                "RNA": tx,
                "notes": ";".join(notes) if notes else "",
                "protein": aa,
            }
        )

    gtsv = args.outdir / "nitr_architecture.tsv"
    stsv = args.outdir / "nitr_domain_segments.tsv"
    gcols = [
        "gene_id", "chrom", "strand", "aa_len", "arch_8b", "architecture",
        "has_SP", "has_V", "has_I", "has_TM", "n_Ig_hits", "ITIM", "ITSM",
        "RNA", "notes", "protein",
    ]
    scols = ["gene_id", "domain", "aa_start", "aa_end", "source", "acc", "score", "sequence"]
    with gtsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=gcols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(gene_rows)
    with stsv.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=scols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(seg_rows)

    print(f"Genes: {len(gene_rows)}")
    print(f"  {gtsv}")
    print(f"  {stsv}")
    if use_sigp and use_dtm:
        print("SP from SignalP; TM/cyto from DeepTMHMM; V/I from InterProScan/SMART.")
    elif use_sigp:
        print("SP from SignalP; TM is heuristic (no DeepTMHMM).")
    elif use_dtm:
        print("TM/cyto from DeepTMHMM; SP is heuristic (no SignalP).")
    else:
        print("SP/TM are heuristics (SignalP and DeepTMHMM not provided).")
    print("--- architecture ---")
    for r in gene_rows:
        print(
            f"  {r['gene_id']:12} {r['architecture']:28} "
            f"SP={r['has_SP']} V={r['has_V']} I={r['has_I']} TM={r['has_TM']} "
            f"ITIM={r['ITIM']} {r['notes']}"
        )


if __name__ == "__main__":
    main()
