#!/usr/bin/env python3
"""Tests for src/bin/rescue_flnc.py (RESCUE_FLNC_FROM_UNPREDICTED_ORFS / RESCUE_FLNC_FROM_ARTIFACTS).

Usage: python3 assets/ci/test_rescue_flnc.py
"""

import doctest
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "src/bin/rescue_flnc.py"


def bed(name, start, end, exons, thick=None, strand="+", chrom="chr1", extra=()):
    """BED12 line from genomic exons [(start, end), ...]."""
    thick_start, thick_end = thick or (start, end)
    sizes = ",".join(str(b - a) for a, b in exons) + ","
    starts = ",".join(str(a - start) for a, _ in exons) + ","
    fields = [chrom, start, end, name, 0, strand, thick_start, thick_end, "0,0,0", len(exons), sizes, starts]
    return "\t".join(str(x) for x in list(fields) + list(extra)) + "\n"


def polish(name, start, end, exons, thick, r=("PASS", "A"), t="PASS", i="PASS"):
    """VEREDICT row: BED12 + R/T/I (status, code, html) + ORF html."""
    extra = [r[0], r[1], "<r>", t, "A" if t == "PASS" else "T", "<t>", i, "A" if i == "PASS" else "I", "<i>", "<o>"]
    return bed(name, start, end, exons, thick, extra=extra)


def run(*args):
    result = subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return result.stderr


def names(path):
    return [line.split("\t")[3] for line in Path(path).read_text().splitlines() if line]


# INFO: reference: tx1 coding 1100-4500 with CDS introns (2000,3000) and (3500,4000);
# INFO: tx2 single-exon CDS; nc1 non-coding
REF = (
    bed("tx1", 1000, 5000, [(1000, 2000), (3000, 3500), (4000, 5000)], (1100, 4500))
    + bed("tx2", 8000, 9000, [(8000, 9000)], (8100, 8900))
    + bed("nc1", 10000, 12000, [(10000, 10500), (11000, 12000)], (10000, 10000))
)
TX1_EXONS = [(1000, 2000), (3000, 3500), (4000, 5000)]


def test_unpredicted(tmp):
    models = tmp / "models.bed"
    models.write_text(
        bed("A__FC0#CN3", 1000, 5000, TX1_EXONS)                                  # got an ORF, in pass
        + bed("B__FC0#CN2", 900, 5200, [(900, 2000), (3000, 3500), (4000, 5200)])  # ORF dropped later; carries tx1 CDS
        + bed("C__FC0#CN1", 1050, 4700, [(1050, 2000), (3000, 3500), (4000, 4700)])  # no ORF; carries tx1 CDS
        + bed("D__FC0#CN1", 1000, 5000, [(1000, 2000), (3000, 3400), (4000, 5000)])  # no ORF; CDS chain broken
        + bed("E__FC0#CN1", 8000, 9000, [(8000, 9000)])                           # no ORF; single-exon CDS only
        + bed("F__FC0#CN1", 1000, 5000, TX1_EXONS)                                 # ORF went to NMD
    )
    predicted = tmp / "xorf.bed"
    predicted.write_text(
        bed("A__FC0#CN3#SC90#OR1@P1", 1000, 5000, TX1_EXONS, (1100, 4500))
        + bed("B__FC0#CN2#SC25#OR1@P1", 900, 5200, [(900, 2000), (3000, 3500), (4000, 5200)], (1100, 4500))
        + bed("F__FC0#CN1#SC90#OR2", 1000, 5000, TX1_EXONS, (1100, 3200))
    )
    passed = tmp / "s.pass.bed"
    passed.write_text(polish("A__FC0#CN3#SC90#OR1@P1", 1000, 5000, TX1_EXONS, (1100, 4500)))
    nmd = tmp / "nmd.bed"
    nmd.write_text(bed("F__FC0#CN1#SC90#OR2", 1000, 5000, TX1_EXONS, (1100, 3200)))
    ref = tmp / "ref.bed"
    ref.write_text(REF)

    out = tmp / "u"
    run("unpredicted", "--models", models, "--predicted", predicted, "--out", passed, nmd,
        "--reference", ref, "--prefix", "s", "--outdir", out)

    assert sorted(names(out / "s.unpredicted.bed")) == ["D__FC0#CN1#UPNO", "E__FC0#CN1#UPNO"], names(out / "s.unpredicted.bed")
    for line in (out / "s.unpredicted.bed").read_text().splitlines():
        f = line.split("\t")
        assert len(f) == 12 and f[6] == f[7] == f[1], line  # INFO: no CDS
    rescued = [line.split("\t") for line in (out / "s.rescued_unpredicted.bed").read_text().splitlines()]
    assert sorted(f[3] for f in rescued) == ["B__FC0#CN2#RSU", "C__FC0#CN1#RSU"], rescued
    for f in rescued:
        assert len(f) == 22 and (f[6], f[7]) == ("1100", "4500"), f
    audit = (out / "s.rescued_unpredicted.tsv").read_text()
    assert "unpredicted:VD\tpass\tcds_transfer\ttx1" in audit and "unpredicted:NO\tpass\tcds_transfer\ttx1" in audit, audit

    # INFO: --no-rescue (fusions): every unpredicted model is kept
    run("unpredicted", "--models", models, "--predicted", predicted, "--out", passed, nmd,
        "--no-rescue", "--prefix", "f", "--outdir", out)
    assert len(names(out / "f.unpredicted.bed")) == 4
    assert "B__FC0#CN2#UPVD" in names(out / "f.unpredicted.bed")
    assert (out / "f.rescued_unpredicted.bed").read_text() == ""


def test_artifacts(tmp):
    tx1_cds = (1100, 4500)
    rt_exons = [(1000, 2000), (3000, 3500), (3600, 3700), (4000, 5000)]
    files = {
        "pass": polish("P#OR1@X", 1000, 5000, TX1_EXONS, tx1_cds),
        "retentions": polish("R1#OR1@X", 1000, 5000, [(1000, 2000), (3000, 5000)], (1100, 4500), r=("RETENTION", "R"))
        + polish("R2#OR1@X", 1000, 5000, TX1_EXONS, (1100, 4500), r=("RETENTION", "R")),  # UTR retention, CDS = tx1
        "truncations": polish("T1#OR1@X", 3000, 5000, [(3000, 3500), (4000, 5000)], (3100, 4500), t="TRUNCATED"),
        "intraprimming": "",
        "rt": polish("X1#OR1@X", 1000, 5000, rt_exons, (1100, 1500), r=("HAS_STRONG_RT", "X"))
        + polish("X2#OR1@X", 1000, 5000, rt_exons, (1100, 1500), r=("HAS_STRONG_RT", "X"), t="TRUNCATED")
        + polish("X3#OR1@X", 1000, 5000, [(1000, 2000), (3000, 3550), (3650, 3700), (4000, 5000)], (1100, 1500), r=("HAS_STRONG_RT", "X")),
        "artifacts": polish("Q1#OR1@X", 8000, 9000, [(8000, 9000)], (8100, 8900), r=("HAS_ARTIFACT", "Q")),
        "trash": polish("F1#OR1@X", 1000, 5000, TX1_EXONS, tx1_cds, t="TRUNCATED", i="INTRAPRIMING"),
    }
    paths = []
    for bucket, content in files.items():
        path = tmp / f"s.{bucket}.bed"
        path.write_text(content)
        paths.append(path)
    # INFO: (3500,3600) annotated STRONG_RT; (3550,3650) unannotated STRONG_RT
    row = "chr1\t{}\t{}\t+" + "\t1" * 11 + "\t{}" + "\t1" * 8 + "\t{}\n"
    introns = tmp / "s.tsv"
    introns.write_text(
        row.format(3500, 3600, "REFERENCE_SUPPORT", "STRONG_RT")
        + row.format(3550, 3650, "NO_REFERENCE_SUPPORT", "STRONG_RT")
        + row.format(2000, 3000, "REFERENCE_SUPPORT", "SPLICED")
    )
    ref = tmp / "ref.bed"
    ref.write_text(REF)

    out = tmp / "a"
    run("artifacts", "--buckets", *paths, "--introns", introns, "--reference", ref, "--annotated-introns",
        "--prefix", "s", "--outdir", out)

    got = {bucket: names(out / f"s.{bucket}.bed") for bucket in files}
    assert sorted(got["pass"]) == ["F1#OR1#RSF@X", "P#OR1@X", "R2#OR1#RSR@X", "X1#OR1#RSX@X"], got
    assert got["retentions"] == ["R1#OR1@X"], got           # CDS has no tx1 intron chain
    assert got["truncations"] == ["T1#OR1@X", "X2#OR1#RSX@X"], got  # T1 CDS differs; X2 re-bucketed
    assert got["rt"] == ["X3#OR1@X"], got                    # flagged intron not annotated
    assert got["artifacts"] == ["Q1#OR1@X"], got             # single-exon CDS never matched
    assert got["trash"] == [], got
    for line in (out / "s.pass.bed").read_text().splitlines():
        assert len(line.split("\t")) == 22, line
    audit = (out / "s.rescued_artifacts.tsv").read_text()
    assert "X2#OR1#RSX@X\tX2#OR1@X\trt\ttruncations\tannotated_introns\t." in audit, audit

    # INFO: without --annotated-introns the rt rows stay
    run("artifacts", "--buckets", *paths, "--reference", ref, "--prefix", "s", "--outdir", tmp / "b")
    assert names(tmp / "b" / "s.rt.bed") == ["X1#OR1@X", "X2#OR1@X", "X3#OR1@X"]


def test_nmd_and_orphans(tmp):
    ref = tmp / "ref.bed"
    ref.write_text(REF)

    # INFO: VEREDICT buckets + NMD (BED12): an NMD row with the tx1 CDS moves to pass, padded
    files = {bucket: "" for bucket in ("pass", "trash", "retentions", "truncations", "intraprimming", "rt", "artifacts")}
    files["pass"] = polish("P#OR1@X", 1000, 5000, TX1_EXONS, (1100, 4500))
    paths = []
    for bucket, content in files.items():
        path = tmp / f"s.{bucket}.bed"
        path.write_text(content)
        paths.append(path)
    nmd = tmp / "in.nmd.bed"  # INFO: staged NMD file; its label is `nmd`
    nmd.write_text(
        bed("N1#OR1@X", 1000, 6000, TX1_EXONS[:2] + [(4000, 5000), (5500, 6000)], (1100, 4500))
        + bed("N2#OR2@X", 1000, 6000, TX1_EXONS[:2] + [(4000, 5000), (5500, 6000)], (1100, 3200))
    )
    out = tmp / "n"
    run("artifacts", "--buckets", *paths, nmd, "--reference", ref, "--prefix", "s", "--outdir", out)
    assert names(out / "s.nmd.bed") == ["N2#OR2@X"], names(out / "s.nmd.bed")
    moved = [l.split("\t") for l in (out / "s.pass.bed").read_text().splitlines() if "#RSN" in l]
    assert [f[3] for f in moved] == ["N1#OR1#RSN@X"] and len(moved[0]) == 22, moved
    assert "nmd\tpass\tcds\ttx1" in (out / "s.rescued_artifacts.tsv").read_text()

    # INFO: NMD only (xORF kept no non-NMD row on this chromosome)
    run("artifacts", "--buckets", nmd, "--reference", ref, "--prefix", "o", "--outdir", out)
    assert sorted(p.name for p in out.glob("o.*.bed")) == ["o.nmd.bed", "o.pass.bed"]
    assert names(out / "o.pass.bed") == ["N1#OR1#RSN@X"]

    # INFO: orphan finder outputs: hq is pass, file labels are kept
    hq = tmp / "smp.pass.hq.bed"
    hq.write_text(files["pass"])
    scraps = tmp / "smp.pass.scraps.bed"
    scraps.write_text(
        polish("S1#OR1@X", 1000, 5000, TX1_EXONS, (1100, 4500))
        + polish("S2#OR1@X", 1000, 5000, TX1_EXONS, (1100, 3200))
    )
    run("artifacts", "--buckets", hq, scraps, "--reference", ref, "--prefix", "smp.pass", "--outdir", out)
    assert sorted(p.name for p in out.glob("smp.pass.*.bed")) == ["smp.pass.hq.bed", "smp.pass.scraps.bed"]
    assert names(out / "smp.pass.hq.bed") == ["P#OR1@X", "S1#OR1#RSS@X"]
    assert names(out / "smp.pass.scraps.bed") == ["S2#OR1@X"]


def test_doctests():
    spec = importlib.util.spec_from_file_location("rescue_flnc", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules["rescue_flnc"] = module  # INFO: dataclasses look the module up while it executes
    sys.dont_write_bytecode = True  # INFO: no __pycache__ in src/bin
    spec.loader.exec_module(module)
    failures, _ = doctest.testmod(module)
    assert failures == 0


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        for test in (test_doctests, test_unpredicted, test_artifacts, test_nmd_and_orphans):
            work = tmp / test.__name__
            work.mkdir()
            test(work) if test is not test_doctests else test()
            print(f"PASS {test.__name__}", flush=True)


if __name__ == "__main__":
    main()
