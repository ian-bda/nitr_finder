#!/usr/bin/env python3
"""Name cluster_genes.tsv loci from BLASTP vs zebrafish + spotted gar."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

UNCHAR_RE = re.compile(
    r"uncharacterized|hypothetical|unnamed|low quality protein|unknown product",
    re.I,
)
LOC_RE = re.compile(r"\bLOC(\d+)\b")
ISOFORM_RE = re.compile(r",?\s*isoform\s+X?\d+.*$", re.I)
PREFIX_RE = re.compile(r"^(PREDICTED:|LOW QUALITY PROTEIN:)\s*", re.I)
SPECIES_RE = re.compile(r"\s*\[[^\]]+\]\s*$")
ACC_RE = re.compile(r"^[A-Z]{2}_\d+\.\d+\s+")

# NCBI product phrases → gene symbol
PRODUCT_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"downstream neighbor of son", re.I), "DONSON"),
    (re.compile(r"\bprotein SON\b|^SON$", re.I), "SON"),
    (re.compile(r"runt-related transcription factor\s*(\d+[a-z]?)", re.I), "RUNX{0}"),
    (
        re.compile(r"serine(?:/| and )arginine-?rich splicing factor\s*(\d+[a-z]?)", re.I),
        "SRSF{0}",
    ),
    (re.compile(r"serine/threonine-protein kinase pim-(\d+)", re.I), "PIM{0}"),
    (re.compile(r"transposon TX1", re.I), "TX1"),
    (re.compile(r"novel immune[- ]type receptor|\bnitr\b", re.I), "NITR"),
    (re.compile(r"axoneme-associated protein mst101", re.I), "MST101"),
    (re.compile(r"heparan-sulfate 6-O-sulfotransferase\s*(\d+)", re.I), "HS6ST{0}"),
]


def locus_key(row: dict) -> str:
    return f"{row['Gene Start']}_{row['Gene Stop']}_{row['Strand']}"


def dump_fasta(rows: list[dict], path: Path) -> None:
    with path.open("w") as fh:
        for row in rows:
            pep = (row.get("Protein") or "").replace("*", "")
            if not pep:
                continue
            fh.write(f">{locus_key(row)}\n")
            for i in range(0, len(pep), 80):
                fh.write(pep[i : i + 80] + "\n")


def clean_product(title: str) -> str:
    t = SPECIES_RE.sub("", title).strip()
    t = ACC_RE.sub("", t).strip()
    t = PREFIX_RE.sub("", t).strip()
    t = ISOFORM_RE.sub("", t).strip()
    if t.lower().startswith("protein "):
        rest = t[8:]
        if rest and rest[0].isupper() and " " not in rest.split()[0]:
            t = rest
    return t


def slug(text: str, n: int = 50) -> str:
    s = re.sub(r"[^A-Za-z0-9]+", "_", text)
    s = re.sub(r"_+", "_", s).strip("_")
    return s[:n] or "uncharacterized"


def symbol_from_product(product: str) -> str | None:
    p = clean_product(product)
    for rx, fmt in PRODUCT_PATTERNS:
        m = rx.search(p)
        if not m:
            continue
        if "{0}" in fmt:
            return fmt.format(m.group(1).upper() if fmt.startswith("SRSF") else m.group(1))
        return fmt
    loc = LOC_RE.search(p)
    if loc and UNCHAR_RE.search(p):
        return f"LOC{loc.group(1)}"
    if UNCHAR_RE.search(p):
        return None
    return slug(p)


def parse_blast(path: Path, species: str) -> dict[str, dict]:
    """Best hit per query (highest bitscore, then lowest evalue)."""
    best: dict[str, dict] = {}
    if not path.exists():
        return best
    with path.open() as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10:
                continue
            qid, sid = p[0], p[1]
            try:
                pident = float(p[2])
                qlen = int(p[4])
                slen = int(p[5])
                qcov = float(p[6])
                evalue = float(p[7])
                bits = float(p[8])
            except ValueError:
                continue
            title = p[9]
            rec = {
                "species": species,
                "accession": sid.split("|")[1] if sid.startswith("ref|") else sid,
                "pident": pident,
                "qlen": qlen,
                "slen": slen,
                "qcov": qcov,
                "evalue": evalue,
                "bits": bits,
                "title": title,
                "product": clean_product(title),
            }
            prev = best.get(qid)
            if prev is None or bits > prev["bits"] or (
                bits == prev["bits"] and evalue < prev["evalue"]
            ):
                best[qid] = rec
    return best


def is_informative(sym: str | None) -> bool:
    if not sym:
        return False
    u = sym.upper()
    if u.startswith("LOC") or u.startswith("UNCHARACTERIZED"):
        return False
    if u.startswith("XP_") or u.startswith("NP_") or u.startswith("YP_"):
        return False
    return True


def choose_name(zf: dict | None, gar: dict | None) -> tuple[str | None, dict | None]:
    """Prefer a real gene symbol; LOC/uncharacterized only if that is all we have."""
    cands = [h for h in (zf, gar) if h and h["evalue"] <= 1e-5 and h["qcov"] >= 20]
    if not cands:
        cands = [h for h in (zf, gar) if h and h["evalue"] <= 1e-3]
    if not cands:
        return None, None
    scored = []
    for h in cands:
        sym = symbol_from_product(h["title"])
        scored.append((sym, h, is_informative(sym)))
    named = [s for s in scored if s[2]]
    if named:
        named.sort(key=lambda x: (-x[1]["bits"], x[1]["evalue"]))
        return named[0][0], named[0][1]
    scored.sort(key=lambda x: (-x[1]["bits"], x[1]["evalue"]))
    return scored[0][0], scored[0][1]


def fmt_hit(h: dict | None) -> str:
    if not h:
        return "none"
    return (
        f"{h['species']}|{h['accession']}|{h['product'][:80]}|"
        f"pident={h['pident']:.1f}|qcov={h['qcov']:.0f}|evalue={h['evalue']:.0e}"
    )


def unique_ids(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for n in names:
        k = seen.get(n, 0) + 1
        seen[n] = k
        out.append(n if k == 1 else f"{n}_{k}")
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--table", type=Path, required=True)
    p.add_argument("--fasta-out", type=Path, help="write query FASTA and exit")
    p.add_argument("--blast-zf", type=Path)
    p.add_argument("--blast-gar", type=Path)
    p.add_argument("-o", "--out", type=Path)
    p.add_argument("--hits-out", type=Path, help="one-row-per-gene BLAST summary")
    args = p.parse_args()

    with args.table.open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit("empty cluster table")

    if args.fasta_out:
        dump_fasta(rows, args.fasta_out)
        print(f"Wrote {args.fasta_out} ({sum(1 for r in rows if r.get('Protein'))} seqs)")
        return

    if not args.blast_zf or not args.blast_gar or not args.out:
        raise SystemExit("need --blast-zf --blast-gar -o (or --fasta-out)")

    zf = parse_blast(args.blast_zf, "zebrafish")
    gar = parse_blast(args.blast_gar, "gar")

    new_ids = []
    summaries = []
    for row in rows:
        key = locus_key(row)
        old = row["GeneID"]
        if row.get("Gene Type") == "NITR":
            new_ids.append(old)
            zhit, ghit = zf.get(key), gar.get(key)
            summaries.append(
                {
                    "Locus": key,
                    "OldGeneID": old,
                    "NewGeneID": old,
                    "Zebrafish": fmt_hit(zhit),
                    "Gar": fmt_hit(ghit),
                }
            )
            zf_note = fmt_hit(zhit)
            gar_note = fmt_hit(ghit)
            notes = row.get("Notes") or ""
            notes = re.sub(r";blastp:[^;]*", "", notes)
            row["Notes"] = notes + f";blastp_zf={zf_note};blastp_gar={gar_note}"
            continue

        sym, used = choose_name(zf.get(key), gar.get(key))
        sym, used = choose_name(zf.get(key), gar.get(key))
        if not is_informative(sym):
            m = re.match(r"(STRG\.\d+|PsenNITR\d+)", old)
            if m:
                sym = m.group(1)
        if not sym:
            loc = LOC_RE.search(old)
            sym = f"LOC{loc.group(1)}" if loc else f"uncharacterized_{row['Gene Start']}"
        if (
            used
            and is_informative(sym)
            and used["pident"] < 50
            and not str(sym).endswith("-like")
            and not str(sym).startswith("STRG.")
            and str(sym) not in {"TX1"}
        ):
            sym = f"{sym}-like"
        new_ids.append(sym)
        summaries.append(
            {
                "Locus": key,
                "OldGeneID": old,
                "NewGeneID": sym,
                "Zebrafish": fmt_hit(zf.get(key)),
                "Gar": fmt_hit(gar.get(key)),
            }
        )
        notes = row.get("Notes") or ""
        notes = re.sub(r";blastp:[^;]*", "", notes)
        notes = re.sub(r";blastp_zf=[^;]*", "", notes)
        notes = re.sub(r";blastp_gar=[^;]*", "", notes)
        row["Notes"] = notes + f";blastp_zf={fmt_hit(zf.get(key))};blastp_gar={fmt_hit(gar.get(key))}"

    final_ids = unique_ids(new_ids)
    for row, gid, s in zip(rows, final_ids, summaries):
        row["GeneID"] = gid
        s["NewGeneID"] = gid

    cols = list(rows[0].keys())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    if args.hits_out:
        with args.hits_out.open("w", newline="") as fh:
            w = csv.DictWriter(
                fh,
                fieldnames=["Locus", "OldGeneID", "NewGeneID", "Zebrafish", "Gar"],
                delimiter="\t",
            )
            w.writeheader()
            w.writerows(summaries)

    print(f"Wrote {args.out}")
    for row, s in zip(rows, summaries):
        print(f"  {row['Gene Start']}-{row['Gene Stop']}  {s['OldGeneID']}  ->  {row['GeneID']}")


if __name__ == "__main__":
    main()
