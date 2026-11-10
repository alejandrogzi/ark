#!/usr/bin/env python3

# Copyright (c) 2026 Alejandro Gonzales-Irribarren <alejandrxgzi@gmail.com>
# Distributed under the terms of the Apache License, Version 2.0.

"""Collect models that lost their ORF and rescue POLISH models backed by the reference annotation.

Two subcommands, one per ARK process (both run per sample@chromosome):

``unpredicted`` (RESCUE_FLNC_FROM_UNPREDICTED_ORFS)
    Models that went into ORF calling (``--models``) minus every model with at least one row in
    any output (``--out``: VEREDICT buckets, NMD). The rest is ``<prefix>.unpredicted.bed``
    (BED12, no CDS, name tagged ``#UPNO`` when xORF kept no ORF for it, ``#UPVD`` when its ORF
    rows were dropped after xORF). With ``--reference``, an unpredicted model whose exon structure
    carries the full CDS of a coding reference transcript takes that CDS and goes to
    ``<prefix>.rescued_unpredicted.bed`` (POLISH schema), tagged ``#RSU``.

``artifacts`` (RESCUE_FLNC_FROM_ARTIFACTS, and RESCUE_FLNC_FROM_ORPHANS on the orphan finder)
    Reads bucket files (``<prefix>.<bucket>.bed``: the seven VEREDICT buckets plus NMD, or the
    orphan finder's ``hq`` (= pass) and ``scraps``) and moves rows out of the rescue categories:
      1. rt / artifacts rows whose flagged introns (RT or ARTIFACT in the iso-classify table)
         are all annotated (REFERENCE_SUPPORT) are re-bucketed by their remaining flaws;
      2. rows of the rescue categories whose CDS (thickStart, thickEnd and the introns between
         them) equals the CDS of a coding reference transcript go to pass.
    Every input bucket is written back as ``<prefix>.<label>.bed`` without its moved rows, and
    moved rows are appended to pass. Moved rows are tagged ``#RS<origin>``: R retentions,
    T truncations, I intraprimming, X rt, Q artifacts, F trash, N nmd, S scraps (U is unpredicted,
    above). Rows without the POLISH columns (NMD is BED12) are padded with NOT_EVALUATED
    descriptors when they join pass.

A CDS match needs at least ``--min-cds-introns`` introns inside the CDS (default 1): a
single-exon CDS matched by coordinates alone would rescue any pre-mRNA lying over it.

Example
-------
    rescue_flnc.py unpredicted --models models.bed --predicted xorf.bed \\
        --out s.pass.bed s.trash.bed nmd.bed --reference ref.bed --prefix s@chr1
    rescue_flnc.py artifacts --buckets s.pass.bed s.trash.bed ... --introns s@chr1.tsv \\
        --reference ref.bed --prefix s@chr1
"""

from __future__ import annotations

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.2"

import argparse
import bisect
import gzip
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Set, Tuple

log = logging.getLogger("rescue_flnc")

# INFO: VEREDICT output layout: BED12 + R/T/I (status, code, html) + ORF html
POLISH_COLUMNS = 22
R_STATUS, R_CODE, T_STATUS, I_STATUS = 12, 13, 15, 18

BUCKETS = (
    "pass",
    "trash",
    "retentions",
    "truncations",
    "intraprimming",
    "rt",
    "artifacts",
    "nmd",
    "scraps",
)
# INFO: file label -> bucket; the orphan finder calls its pass set `hq`
LABELS = {bucket: bucket for bucket in BUCKETS}
LABELS["hq"] = "pass"
ORIGIN_CODES = {
    "retentions": "R",
    "truncations": "T",
    "intraprimming": "I",
    "rt": "X",
    "artifacts": "Q",
    "trash": "F",
    "nmd": "N",
    "scraps": "S",
    "unpredicted": "U",
}
DEFAULT_CATEGORIES = "retentions,truncations,intraprimming,rt,artifacts,trash,nmd,scraps"
RESCUE_TAG = "#RS"
UNPREDICTED_TAG = "#UP"
FLAGGED_SUPPORT = {"STRONG_RT", "WEAK_RT", "ARTIFACT"}

Intron = Tuple[int, int]
CdsKey = Tuple[str, str, int, int, Tuple[Intron, ...]]


@dataclass
class Bed:
    """One BED12 line, extra columns kept verbatim."""

    fields: List[str]

    @property
    def chrom(self) -> str:
        return self.fields[0]

    @property
    def start(self) -> int:
        return int(self.fields[1])

    @property
    def end(self) -> int:
        return int(self.fields[2])

    @property
    def name(self) -> str:
        return self.fields[3]

    @property
    def strand(self) -> str:
        return self.fields[5]

    @property
    def thick_start(self) -> int:
        return int(self.fields[6])

    @property
    def thick_end(self) -> int:
        return int(self.fields[7])

    def exons(self) -> List[Intron]:
        """Return the exon blocks as genomic half-open intervals."""

        sizes = [int(x) for x in self.fields[10].rstrip(",").split(",") if x]
        starts = [int(x) for x in self.fields[11].rstrip(",").split(",") if x]
        return [(self.start + s, self.start + s + z) for s, z in zip(starts, sizes)]

    def introns(self) -> List[Intron]:
        """Return the introns between consecutive exon blocks."""

        exons = self.exons()
        return [(exons[i][1], exons[i + 1][0]) for i in range(len(exons) - 1)]

    def structure(self) -> Tuple[str, str, int, int, str, str]:
        """Return a name-free key of the exon structure."""

        return (
            self.chrom,
            self.strand,
            self.start,
            self.end,
            self.fields[10].rstrip(","),
            self.fields[11].rstrip(","),
        )

    def line(self) -> str:
        return "\t".join(self.fields) + "\n"


def cds_introns(introns: Sequence[Intron], thick_start: int, thick_end: int) -> Tuple[Intron, ...]:
    """Return the introns lying inside [thick_start, thick_end].

    Example
    -------
    >>> cds_introns([(10, 20), (30, 40), (50, 60)], 15, 45)
    ((30, 40),)
    """

    return tuple((a, b) for a, b in introns if a >= thick_start and b <= thick_end)


def in_exon(exons: Sequence[Intron], position: int) -> bool:
    """Return True when the base at `position` is exonic."""

    return any(a <= position < b for a, b in exons)


def open_text(path: Path):
    """Open a plain or gzipped text file."""

    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")


def read_bed(path: Optional[Path]) -> List[Bed]:
    """Read every non-empty, non-header line of a BED file; a missing path reads as empty."""

    if path is None or not Path(path).is_file():
        return []
    rows = []
    with open_text(path) as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            rows.append(Bed(line.split("\t")))
    return rows


def tag_name(name: str, tag: str) -> str:
    """Insert `tag` before a trailing `@PROTEIN` (xORF names) or append it.

    Example
    -------
    >>> tag_name("R1@0_c__FC0#CN1#SC100#OR1@CALM1_HUMAN", "#RSR")
    'R1@0_c__FC0#CN1#SC100#OR1#RSR@CALM1_HUMAN'
    >>> tag_name("R1@0_c__FC0#CN1#SG", "#RSU")
    'R1@0_c__FC0#CN1#SG#RSU'
    """

    if tag in name:
        return name
    last_hash = name.rfind("#")
    at = name.find("@", last_hash + 1) if last_hash >= 0 else -1
    if last_hash >= 0 and at > last_hash:
        return name[:at] + tag + name[at:]
    return name + tag


def is_rescued(name: str) -> bool:
    return RESCUE_TAG in name


class Reference:
    """Coding reference transcripts of the chromosomes in use, indexed by CDS."""

    def __init__(self, path: Path, chroms: Set[str], min_cds_introns: int) -> None:
        self.by_key: Dict[CdsKey, List[str]] = defaultdict(list)
        # INFO: per (chrom, strand): CDS sorted by thickStart for the transfer search
        self.by_locus: Dict[Tuple[str, str], List[Tuple[int, int, Tuple[Intron, ...], int, str]]] = (
            defaultdict(list)
        )
        kept = 0
        for row in read_bed(path):
            if row.chrom not in chroms or row.thick_start >= row.thick_end:
                continue
            introns = cds_introns(row.introns(), row.thick_start, row.thick_end)
            if len(introns) < min_cds_introns:
                continue
            key = (row.chrom, row.strand, row.thick_start, row.thick_end, introns)
            self.by_key[key].append(row.name)
            coding = sum(
                max(0, min(b, row.thick_end) - max(a, row.thick_start)) for a, b in row.exons()
            )
            self.by_locus[(row.chrom, row.strand)].append(
                (row.thick_start, row.thick_end, introns, coding, row.name)
            )
            kept += 1
        for entries in self.by_locus.values():
            entries.sort()
        self.starts = {locus: [e[0] for e in entries] for locus, entries in self.by_locus.items()}
        log.info("Indexed %d coding reference transcripts on %d chromosomes", kept, len(chroms))

    def match(self, row: Bed, min_cds_introns: int) -> List[str]:
        """Return the reference transcripts with exactly the CDS of `row`."""

        if row.thick_start >= row.thick_end:
            return []
        introns = cds_introns(row.introns(), row.thick_start, row.thick_end)
        if len(introns) < min_cds_introns:
            return []
        return self.by_key.get((row.chrom, row.strand, row.thick_start, row.thick_end, introns), [])

    def transfer(self, row: Bed) -> Optional[Tuple[int, int, List[str]]]:
        """Return the longest reference CDS that `row`'s exon structure carries intact.

        The reference CDS fits when its first and last coding bases are exonic in `row` and the
        introns of `row` inside the CDS are exactly the reference CDS introns, so the spliced
        CDS, frame and protein are the reference ones.
        """

        locus = (row.chrom, row.strand)
        entries = self.by_locus.get(locus)
        if not entries:
            return None
        exons, introns = row.exons(), row.introns()
        lo = bisect.bisect_left(self.starts[locus], row.start)
        hi = bisect.bisect_left(self.starts[locus], row.end)
        best: Optional[Tuple[int, int, int, List[str]]] = None
        for thick_start, thick_end, ref_introns, coding, name in entries[lo:hi]:
            if thick_end > row.end:
                continue
            if not (in_exon(exons, thick_start) and in_exon(exons, thick_end - 1)):
                continue
            if cds_introns(introns, thick_start, thick_end) != ref_introns:
                continue
            if best is None or coding > best[2]:
                best = (thick_start, thick_end, coding, [name])
            elif coding == best[2] and (thick_start, thick_end) == best[:2]:
                best[3].append(name)
        if best is None:
            return None
        return best[0], best[1], sorted(best[3])


class Roots:
    """Maps output names (model name + xORF suffix) back to the input model names."""

    def __init__(self, models: Sequence[Bed]) -> None:
        self.names = {row.name for row in models}
        self.structures: Dict[Tuple, List[str]] = defaultdict(list)
        for row in models:
            self.structures[row.structure()].append(row.name)
        self.unresolved = 0

    def resolve(self, row: Bed) -> List[str]:
        """Return the model names behind one output row.

        xORF writes `[PREFIX#...]<model>#SC..#OR..[#DU][#IN..][@PROTEIN]` (`UNMASK#` when no
        selenocysteine masking is run, plus `--custom-prefix`), so the root is the longest
        `#`-bounded substring, starting after a `#` and ending before a `#` or `@`, that is a model
        name. Renamed roots (`--rebase`) fall back to the exon structure, which xORF keeps.

        Example
        -------
        >>> roots = Roots([Bed("c 0 9 R1@0_c__FC0#CN3 0 + 0 9 0 1 9, 0,".split())])
        >>> roots.resolve(Bed("c 0 5 UNMASK#R1@0_c__FC0#CN3#SC90#OR1@P 0 + 0 5 0 1 5, 0,".split()))
        ['R1@0_c__FC0#CN3']
        """

        name = row.name
        if name in self.names:
            return [name]
        lefts = [0] + [i + 1 for i, char in enumerate(name) if char == "#"]
        rights = [i for i, char in enumerate(name) if char in "#@"][::-1]
        for left in lefts:
            for right in rights:
                if right <= left:
                    break
                if name[left:right] in self.names:
                    return [name[left:right]]
        by_structure = self.structures.get(row.structure())
        if by_structure:
            return by_structure
        self.unresolved += 1
        return []


def write_rows(path: Path, rows: Iterable[str]) -> int:
    count = 0
    with open(path, "w", encoding="utf-8") as handle:
        for line in rows:
            handle.write(line)
            count += 1
    return count


def write_audit(path: Path, records: Sequence[Tuple[str, ...]]) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("id\toriginal_id\torigin\tdestination\trule\treference\n")
        for record in records:
            handle.write("\t".join(record) + "\n")


def padded(row: Bed, source: str, references: Sequence[str]) -> Bed:
    """Give a BED12 row the POLISH columns (descriptors NOT_EVALUATED) so it can join pass."""

    if len(row.fields) >= POLISH_COLUMNS:
        return row
    fields = row.fields[:12]
    note = escape(", ".join(references))
    not_evaluated = ["NOT_EVALUATED", ".", f"<h2>Rescue</h2><br>- not evaluated: {escape(source)}<br>"]
    orf_html = (
        f"<h2>ORF</h2><br>- source: {escape(source)}<br>"
        f"- reference: {note}<br>- cds: {fields[6]}-{fields[7]}<br>"
    )
    return Bed(fields + not_evaluated * 3 + [orf_html])


def polish_row(row: Bed, thick_start: int, thick_end: int, references: Sequence[str]) -> Bed:
    """Turn an unpredicted model into a POLISH-schema pass row with the transferred CDS."""

    fields = row.fields[:12]
    fields[3] = tag_name(row.name, RESCUE_TAG + ORIGIN_CODES["unpredicted"])
    fields[6], fields[7] = str(thick_start), str(thick_end)
    return padded(Bed(fields), "reference annotation (rescued, no xORF ORF)", references)


def run_unpredicted(args: argparse.Namespace) -> None:
    """RESCUE_FLNC_FROM_UNPREDICTED_ORFS: models in minus models out, then CDS transfer."""

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    models = read_bed(args.models)
    roots = Roots(models)

    went_out: Set[str] = set()
    for path in args.out or []:
        for row in read_bed(path):
            went_out.update(roots.resolve(row))
    predicted: Set[str] = set()
    for row in read_bed(args.predicted):
        predicted.update(roots.resolve(row))
    if roots.unresolved:
        log.warning("%d output rows did not map back to an input model", roots.unresolved)

    unpredicted = [row for row in models if row.name not in went_out]
    reference = None
    if args.reference and not args.no_rescue:
        reference = Reference(args.reference, {row.chrom for row in unpredicted}, args.min_cds_introns)

    kept, rescued, audit = [], [], []
    for row in unpredicted:
        reason = "VD" if row.name in predicted else "NO"
        hit = reference.transfer(row) if reference else None
        if hit is None:
            fields = row.fields[:12]
            fields[3] = tag_name(row.name, UNPREDICTED_TAG + reason)
            fields[6] = fields[7] = fields[1]  # INFO: no CDS
            kept.append(Bed(fields).line())
            continue
        thick_start, thick_end, references = hit
        new = polish_row(row, thick_start, thick_end, references)
        rescued.append(new.line())
        audit.append(
            (new.name, row.name, f"unpredicted:{reason}", "pass", "cds_transfer", ",".join(references))
        )

    n_kept = write_rows(outdir / f"{args.prefix}.unpredicted.bed", kept)
    write_rows(outdir / f"{args.prefix}.rescued_unpredicted.bed", rescued)
    write_audit(outdir / f"{args.prefix}.rescued_unpredicted.tsv", audit)
    log.info(
        "%d models in, %d with an output row, %d unpredicted: %d rescued with a reference CDS, %d kept",
        len(models),
        len(models) - len(unpredicted),
        len(unpredicted),
        len(rescued),
        n_kept,
    )


def flagged_introns(path: Optional[Path]) -> Dict[Tuple[str, int, int, str], bool]:
    """Read the iso-classify table: RT/ARTIFACT introns -> annotated (REFERENCE_SUPPORT)."""

    flagged: Dict[Tuple[str, int, int, str], bool] = {}
    if path is None or not Path(path).is_file():
        return flagged
    with open_text(path) as handle:
        for raw in handle:
            parts = raw.rstrip("\n").split("\t")
            if len(parts) < 25 or parts[24] not in FLAGGED_SUPPORT:
                continue
            key = (parts[0], int(parts[1]), int(parts[2]), parts[3])
            flagged[key] = parts[15] == "REFERENCE_SUPPORT"
    return flagged


def rebucket(row: Bed, flaws: int) -> str:
    """Bucket an rt/artifact row by its remaining flaws once its flagged introns are cleared."""

    fails = {
        "retentions": "R" in row.fields[R_CODE],
        "truncations": row.fields[T_STATUS] != "PASS",
        "intraprimming": row.fields[I_STATUS] != "PASS",
    }
    count = sum(fails.values())
    if count == 0:
        return "pass"
    if count >= flaws:
        return "trash"
    return next(bucket for bucket, failed in fails.items() if failed)


def run_artifacts(args: argparse.Namespace) -> None:
    """RESCUE_FLNC_FROM_ARTIFACTS: annotated-intron re-bucketing, then reference CDS matches."""

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    buckets: Dict[str, List[Bed]] = {bucket: [] for bucket in BUCKETS}
    labels: Dict[str, str] = {"pass": "pass"}  # INFO: bucket -> output file label
    for path in args.buckets:
        bucket, label = bucket_of(Path(path).name)
        buckets[bucket].extend(read_bed(path))
        labels[bucket] = label

    unknown = [c for c in args.categories.split(",") if c and (c not in ORIGIN_CODES or c == "unpredicted")]
    if unknown:
        raise SystemExit(f"ERROR: unknown rescue categories: {unknown}")
    # INFO: only categories with an input file; the orphan run sees pass + scraps only
    categories = [c for c in args.categories.split(",") if c in labels]

    audit: List[Tuple[str, ...]] = []
    moves: Dict[Tuple[str, str], int] = defaultdict(int)

    def move(row: Bed, origin: str, destination: str, rule: str, references: Sequence[str]) -> None:
        original = row.name
        row.fields[3] = tag_name(original, RESCUE_TAG + ORIGIN_CODES[origin])
        if destination == "pass":
            row = padded(row, f"xORF (rescued from {origin})", references)
        buckets[destination].append(row)
        audit.append((row.name, original, origin, destination, rule, ",".join(references) or "."))
        moves[(origin, destination)] += 1

    # INFO: 1) rt / artifact rows whose RT/ARTIFACT introns are all annotated
    if args.annotated_introns:
        table = flagged_introns(args.introns)
        for origin in ("rt", "artifacts"):
            if origin not in categories:
                continue
            stay = []
            for row in buckets[origin]:
                hits = [table.get((row.chrom, a, b, row.strand)) for a, b in row.introns()]
                hits = [annotated for annotated in hits if annotated is not None]
                if hits and all(hits):
                    move(row, origin, rebucket(row, args.flaws), "annotated_introns", [])
                else:
                    stay.append(row)
            buckets[origin] = stay

    # INFO: 2) reference CDS matches in every rescue category go to pass
    chroms = {row.chrom for bucket in categories for row in buckets[bucket]}
    reference = Reference(args.reference, chroms, args.min_cds_introns)
    for origin in categories:
        stay = []
        for row in buckets[origin]:
            references = reference.match(row, args.min_cds_introns)
            if not references:
                stay.append(row)
            elif is_rescued(row.name):
                # INFO: already moved by rule 1; it keeps that origin tag
                buckets["pass"].append(padded(row, f"xORF (rescued from {origin})", references))
                audit.append((row.name, row.name, origin, "pass", "cds", ",".join(sorted(references))))
                moves[(origin, "pass")] += 1
            else:
                move(row, origin, "pass", "cds", sorted(references))
        buckets[origin] = stay

    for bucket, label in labels.items():
        write_rows(outdir / f"{args.prefix}.{label}.bed", (row.line() for row in buckets[bucket]))
    write_audit(outdir / f"{args.prefix}.rescued_artifacts.tsv", audit)
    for (origin, destination), count in sorted(moves.items()):
        log.info("Rescued %d rows: %s -> %s", count, origin, destination)
    log.info("Rescued %d rows in total", len(audit))


def bucket_of(filename: str) -> Tuple[str, str]:
    """Return the (bucket, file label) of `<prefix>.<label>.bed`.

    Example
    -------
    >>> bucket_of("s@chr1.retentions.bed"), bucket_of("s.pass.hq.bed"), bucket_of("s.pass.scraps.bed")
    (('retentions', 'retentions'), ('pass', 'hq'), ('scraps', 'scraps'))
    """

    for label, bucket in LABELS.items():
        if filename.endswith(f".{label}.bed"):
            return bucket, label
    raise SystemExit(f"ERROR: {filename} is not a bucket file (<prefix>.<{'|'.join(LABELS)}>.bed)")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("-v", "--version", action="version", version=f"rescue_flnc.py {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-p", "--prefix", required=True, help="Prefix of every output file")
    common.add_argument("-o", "--outdir", default="rescue", help="Output directory [rescue]")
    common.add_argument("-r", "--reference", type=Path, help="Reference annotation BED12[.gz]")
    common.add_argument(
        "-m",
        "--min-cds-introns",
        type=int,
        default=1,
        help="Introns a CDS needs to be matched [1; a single-exon CDS is never matched]",
    )

    unpredicted = sub.add_parser("unpredicted", parents=[common], help=run_unpredicted.__doc__)
    unpredicted.add_argument("--models", type=Path, required=True, help="Models that went into ORF calling")
    unpredicted.add_argument("--predicted", type=Path, help="xORF output BED (models that kept an ORF)")
    unpredicted.add_argument("--out", type=Path, nargs="*", help="Every BED a model can leave in")
    unpredicted.add_argument("--no-rescue", action="store_true", help="Only collect unpredicted models")
    unpredicted.set_defaults(func=run_unpredicted)

    artifacts = sub.add_parser("artifacts", parents=[common], help=run_artifacts.__doc__)
    artifacts.add_argument("--buckets", type=Path, nargs="+", required=True, help="Bucket files <prefix>.<bucket>.bed (VEREDICT, nmd, orphan hq/scraps)")
    artifacts.add_argument("--introns", type=Path, help="iso-classify intron table (annotated-intron rule)")
    artifacts.add_argument(
        "--categories",
        default=DEFAULT_CATEGORIES,
        help=f"Buckets to rescue from [{DEFAULT_CATEGORIES}]",
    )
    artifacts.add_argument(
        "--annotated-introns",
        action="store_true",
        help="Re-bucket rt/artifact rows whose RT/ARTIFACT introns are all annotated",
    )
    artifacts.add_argument("--flaws", type=int, default=2, help="Failed descriptors that make trash [2, as VEREDICT]")
    artifacts.set_defaults(func=run_artifacts)

    args = parser.parse_args(argv)
    if args.command == "artifacts" and args.reference is None:
        parser.error("artifacts needs --reference")
    return args


def main(argv: Optional[Sequence[str]] = None) -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s", stream=sys.stderr)
    args = parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
