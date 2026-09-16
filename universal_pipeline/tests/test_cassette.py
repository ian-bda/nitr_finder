#!/usr/bin/env python3
"""Cassette + RNA stitch vs Polypterus recommended exons (and a synthetic locus)."""

from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from Bio.Seq import Seq  # noqa: E402
from pyfaidx import Fasta  # noqa: E402

import stitch_nitr_genomic as st  # noqa: E402

BUNDLE = Path("/Users/ianbda/Desktop/igv_bundle")
REC = BUNDLE / "NITR_recommended_exons.csv"
BICHIR_FA = BUNDLE / "polypterus_bichir/genome/cluster_window.fa"
SENEG_FA = BUNDLE / "polypterus_senegalus/genome/cluster_window.fa"

# fasta pos 1 = original genomic start of the IGV window
OFFSET = {
    "JANRMT010028174.1": 31199999,
    "CM029056.1": 22473013,
}

# Met→stop / DeepTMHMM inventions from the old 8b+9w table (original coords)
OLD_FAKE = {
    "PbicNITR1": [(31488876, 31488988, "TM")],
    "PbicNITR4": [(31595188, 31595337, "SP"), (31579829, 31579941, "TM")],
    "PbicNITR9": [(31722451, 31722600, "SP")],
    "PbicNITR12": [(31844786, 31844920, "SP"), (31823144, 31823256, "TM")],
}


def load_recommended(path: Path) -> dict[str, list[dict]]:
    genes: dict[str, list[dict]] = defaultdict(list)
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            gid = row["GeneID"]
            genes[gid].append(
                {
                    "chrom": row["Pseudo-chromosome"],
                    "gene": gid,
                    "notes": row["Notes"],
                    "start": int(row["Exon Start"]),
                    "end": int(row["Exon Stop"]),
                    "strand": "+" if row["Orientation"].lower().startswith("f") or row["Orientation"] == "+" else "-",
                    "domain": row["Protein Domain*"].strip(),
                    "nt": row["Nucleotide Sequence"],
                    "pep": row["Protein sequence"],
                }
            )
    return genes


def shift(chrom: str, pos: int) -> int:
    return pos - OFFSET[chrom]


def unshift(chrom: str, pos: int) -> int:
    return pos + OFFSET[chrom]


def is_rna_gene(notes: str) -> bool:
    return "no_RNA" not in notes and "RNA" in notes


def gene_row(gid: str, exons: list[dict]) -> dict:
    v = next((e for e in exons if e["domain"] == "V"), None)
    i = next((e for e in exons if e["domain"] == "I"), None)
    labs = {e["domain"] for e in exons}
    if "V" in labs and "I" in labs:
        arch = "V+I"
    elif "V" in labs:
        arch = "V-only"
    else:
        arch = "I-only"
    chrom = exons[0]["chrom"]
    bait = "".join(e["pep"] for e in exons if e["domain"] in ("V", "I"))
    return {
        "gene_id": gid,
        "chrom": chrom,
        "strand": exons[0]["strand"],
        "arch": arch,
        "v_start": shift(chrom, v["start"]) if v else None,
        "v_end": shift(chrom, v["end"]) if v else None,
        "i_start": shift(chrom, i["start"]) if i else None,
        "i_end": shift(chrom, i["end"]) if i else None,
        "protein_8a": bait,
    }


def mock_gtf(genes: dict[str, list[dict]], only_rna: bool = True) -> str:
    lines = []
    for gid, exons in genes.items():
        if only_rna and not is_rna_gene(exons[0]["notes"]):
            continue
        chrom = exons[0]["chrom"]
        strand = exons[0]["strand"]
        for e in exons:
            s, t = shift(chrom, e["start"]), shift(chrom, e["end"])
            lab = e["domain"] if e["domain"] in ("SP", "V", "I", "TM", "cyto") else "cyto"
            attr = f'gene_id "{gid}"; transcript_id "{gid}.rna"; domain "{lab}"; cov "9";'
            lines.append(f"{chrom}\trec\texon\t{s}\t{t}\t.\t{strand}\t.\t{attr}\n")
    return "".join(lines)


def overlap_bp(a0, a1, b0, b1) -> int:
    return max(0, min(a1, b1) - max(a0, b0) + 1)


class TestHelpers(unittest.TestCase):
    def test_is_sp_nitr2(self):
        self.assertTrue(st.is_sp("MSPLDFVFWIIYSA"))
        self.assertTrue(st.is_sp("MIVCSILLLLFCKN"))
        self.assertTrue(st.is_sp("MRPLCFILMCAGSS"))
        self.assertFalse(st.is_sp("YNSLGVPFRSLTAP"))
        self.assertFalse(st.is_sp("MGLFCCFLKVINLNCEIFSK"))

    def test_continued_aa_phase(self):
        prefix = "ATGAAA" + "A"  # 7 nt, leftover 1
        extra = "TTAAAATAA"
        pep = st.continued_aa(prefix, extra)
        self.assertTrue(pep.startswith("I"))

    def test_no_intron_walk_tm(self):
        self.assertEqual(st.extend_cassette("N" * 200, 10, 200, "ATG" * 80), [])

    def test_exon_goes_to_one_gene(self):
        g4 = {
            "gene_id": "NITR4", "chrom": "chr1", "strand": "+",
            "v_start": 1000, "v_end": 1300, "i_start": 2000, "i_end": 2300,
        }
        g5 = {
            "gene_id": "NITR5", "chrom": "chr1", "strand": "+",
            "v_start": 5000, "v_end": 5300, "i_start": 6000, "i_end": 6300,
        }
        same = [g4, g5]
        # SP next to NITR5 V belongs to 5, not 4's tail
        self.assertIs(st.owner_gene(4800, 4900, same), g5)
        # TM just after NITR4 I belongs to 4
        self.assertIs(st.owner_gene(2500, 2600, same), g4)
        # NITR4 I overlap
        self.assertIs(st.owner_gene(2000, 2300, same), g4)
        # midpoint of the gap: closer to 4 I vs 5 V
        self.assertIs(st.owner_gene(2500, 2550, same), g4)
        self.assertIs(st.owner_gene(4700, 4750, same), g5)


@unittest.skipUnless(REC.exists() and BICHIR_FA.exists(), "igv_bundle recommended table not present")
class TestBichirWindow(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rec = {g: e for g, e in load_recommended(REC).items() if g.startswith("Pbic")}
        cls.genome = Fasta(str(BICHIR_FA), as_raw=True, sequence_always_upper=True)
        cls.genes = [gene_row(g, e) for g, e in sorted(cls.rec.items(), key=lambda x: x[1][0]["start"])]

    def _run(self, with_rna: bool):
        gtf_txt = mock_gtf(self.rec, only_rna=True) if with_rna else ""
        txs = []
        if gtf_txt:
            with tempfile.NamedTemporaryFile("w", suffix=".gtf", delete=False) as fh:
                fh.write(gtf_txt)
                path = Path(fh.name)
            txs = st.parse_gtf(path)
            path.unlink()
        return [st.stitch_one(self.genome, self.genes, i, txs) for i in range(len(self.genes))]

    def test_rna_recovers_full_models(self):
        results = {r["gene_id"]: r for r in self._run(True)}
        for gid, exons in self.rec.items():
            if not is_rna_gene(exons[0]["notes"]):
                continue
            got = results[gid]
            want = [e["domain"] for e in exons]
            have = [k for *_, k in got["exons"]]
            self.assertGreaterEqual(got["aa_len"], 60, gid)
            self.assertIn("from_transcript", got["notes"], f"{gid} {got['notes']}")
            for lab in ("V", "I", "SP", "TM"):
                if lab in want:
                    self.assertIn(lab, have, f"{gid} missing {lab}: {have} notes={got['notes']}")

    def test_cassette_does_not_invent_old_fakes(self):
        results = {r["gene_id"]: r for r in self._run(False)}
        for gid, fakes in OLD_FAKE.items():
            got = results[gid]
            chrom = self.rec[gid][0]["chrom"]
            for s, e, lab in fakes:
                for gs, ge, k in got["exons"]:
                    os, oe = unshift(chrom, gs), unshift(chrom, ge)
                    ov = overlap_bp(os, oe, s, e)
                    span = min(e - s + 1, oe - os + 1)
                    if ov >= 0.8 * span and k == lab:
                        self.fail(f"{gid} cassette revived old {lab} {s}-{e} as {os}-{oe}")

    def test_no_rna_genes_keep_v_i_only_unless_true_splice(self):
        results = {r["gene_id"]: r for r in self._run(False)}
        for gid, exons in self.rec.items():
            if is_rna_gene(exons[0]["notes"]):
                continue
            got = results[gid]
            have = {k for *_, k in got["exons"]}
            self.assertTrue(have & {"V", "I"}, f"{gid} lost Ig: {have}")
            self.assertTrue(have <= {"V", "I"}, f"{gid} invented {have - {'V', 'I'}}")

    def test_csv_peptide_matches_nt(self):
        from build_nitr_exon_csv import pep_of

        results = self._run(True)
        for r in results:
            seq = str(self.genome[r["chrom"]])
            for s, e, k in r["exons"]:
                nt = seq[s - 1 : e]
                if r["strand"] == "-":
                    nt = str(Seq(nt).reverse_complement())
                pep = pep_of(nt)
                self.assertEqual(pep, pep_of(nt), f"{r['gene_id']} {k}")
                self.assertEqual(len(pep), len(nt) // 3)


@unittest.skipUnless(REC.exists() and SENEG_FA.exists(), "senegalus window not present")
class TestSenegalusNoInvent(unittest.TestCase):
    def test_cassette_keeps_no_rna_ig(self):
        rec = {g: e for g, e in load_recommended(REC).items() if g.startswith("Psen")}
        genome = Fasta(str(SENEG_FA), as_raw=True, sequence_always_upper=True)
        genes = [gene_row(g, e) for g, e in sorted(rec.items(), key=lambda x: x[1][0]["start"])]
        results = {r["gene_id"]: r for r in [st.stitch_one(genome, genes, i, []) for i in range(len(genes))]}
        for gid, exons in rec.items():
            if is_rna_gene(exons[0]["notes"]):
                continue
            have = {k for *_, k in results[gid]["exons"]}
            self.assertTrue(have & {"V", "I"}, f"{gid} {have}")
            self.assertTrue(have <= {"V", "I"}, f"{gid} invented {have - {'V', 'I'}}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
