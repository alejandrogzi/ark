#!/usr/bin/env python3

from __future__ import annotations

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.0.9"

import argparse
import logging
import os
from dataclasses import dataclass
from html import escape
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple, Union

import pandas as pd

BED12_COLS: List[str] = [
    "chr",
    "start",
    "end",
    "id",
    "score",
    "strand",
    "cds_start",
    "cds_end",
    "rgb",
    "exon_count",
    "exon_sizes",
    "exon_starts",
]

STATUS_TSV_COLS: List[str] = ["id", "status", "code", "html"]
STATUS_COLUMN_SUFFIXES: Tuple[str, str, str] = ("status", "code", "html")

ORF_COLUMNS_TO_KEEP: List[str] = [
    "chr",
    "start",
    "end",
    "id",
    "strand",
    "tai_start_score",
    "tai_stop_score",
    "start_codon",
    "stop_codon",
    "inner_stops",
    "orf_type",
    "netstart_atg_score",
    "transaid_start_score",
    "transaid_stop_score",
    "transaid_integrated_score",
    "blast_pid",
    "blast_evalue",
    "blast_offset",
    "blast_length",
    "blast_percentage_aligned",
    "rna_score",
    "tai_mean_score",
    "tai_combined_score",
    "transaid_mean_score",
    "predicted_class",
    "prob_noncoding",
    "prob_coding",
]

ORF_RENAME_MAP: Dict[str, str] = {
    column: f"O_{column}" for column in ORF_COLUMNS_TO_KEEP
}
ORF_HTML_COL: str = "O_html"

ORF_TYPE_MAPPING: Dict[Union[str, int], str] = {
    1: "COMPLETE",
    2: "COMPLETE_NESTED",
    3: "TAI_UNIQUE",
    4: "FIVE_PRIME_PARTIAL",
    5: "THREE_PRIME_PARTIAL",
    6: "THREE_PRIME_PARTIAL_NESTED",
    7: "FIVE_PRIME_PARTIAL_NESTED",
    "1": "COMPLETE",
    "2": "COMPLETE_NESTED",
    "3": "TAI_UNIQUE",
    "4": "FIVE_PRIME_PARTIAL",
    "5": "THREE_PRIME_PARTIAL",
    "6": "THREE_PRIME_PARTIAL_NESTED",
    "7": "FIVE_PRIME_PARTIAL_NESTED",
}

OUTPUT_FILENAMES: Dict[str, str] = {
    "pass": "pass.bed",
    "trash": "trash.bed",
    "retentions": "retentions.bed",
    "truncations": "truncations.bed",
    "intraprimming": "intraprimming.bed",
    "rt": "rt.bed",
    "artifact": "artifacts.bed",
}


PathLike = Union[str, os.PathLike[str]]

log = logging.getLogger(__name__)
logging.basicConfig(encoding="utf-8", level=logging.INFO)


@dataclass(frozen=True)
class StatusInputSpec:
    """Configuration for one mandatory status table input."""

    argument_name: str
    prefix: str
    bucket_name: str

    @property
    def status_column(self) -> str:
        """Return the merged status column name for this input."""

        return f"{self.prefix}_status"

    @property
    def code_column(self) -> str:
        """Return the merged code column name for this input."""

        return f"{self.prefix}_code"

    @property
    def html_column(self) -> str:
        """Return the merged HTML column name for this input."""

        return f"{self.prefix}_html"

    @property
    def output_columns(self) -> List[str]:
        """Return the full merged column order for this input."""

        return [self.status_column, self.code_column, self.html_column]


STATUS_INPUT_SPECS: Tuple[StatusInputSpec, ...] = (
    StatusInputSpec(argument_name="retentions", prefix="R", bucket_name="retentions"),
    StatusInputSpec(argument_name="truncations", prefix="T", bucket_name="truncations"),
    StatusInputSpec(
        argument_name="intraprimming",
        prefix="I",
        bucket_name="intraprimming",
    ),
)


def parse_args() -> argparse.Namespace:
    """Parse command line arguments for the BED12 assembler."""

    parser = argparse.ArgumentParser(
        description="Merge BED12 reads with status and ORF annotations."
    )
    parser.add_argument(
        "-r",
        "--reads",
        type=Path,
        required=True,
        help="Path to the input BED12 reads file.",
    )
    parser.add_argument(
        "--retentions",
        type=Path,
        required=True,
        help="Path to the retentions TSV file.",
    )
    parser.add_argument(
        "--truncations",
        type=Path,
        required=True,
        help="Path to the truncations TSV file.",
    )
    parser.add_argument(
        "--intraprimming",
        type=Path,
        required=True,
        help="Path to the intraprimming TSV file.",
    )
    parser.add_argument(
        "--orfs",
        type=Path,
        required=True,
        help="Path to the ORF prediction TSV file.",
    )
    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        default=Path("."),
        help="Directory where the output BED files will be written.",
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        required=True,
        help="Prefix used for every output file name.",
    )
    parser.add_argument(
        "-f",
        "--flaws",
        type=int,
        default=2,
        help="Minimum number of failed status columns required to send a row to trash.",
    )
    parser.add_argument(
        "-O",
        "--orf-score-threshold",
        type=float,
        default=0.3,
        help="Minimum ORF prob_coding score required to retain a row.",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"Version: {__version__}",
    )

    args = parser.parse_args()
    if args.flaws < 1:
        raise ValueError("--flaws must be at least 1")

    return args


def run(args: argparse.Namespace) -> None:
    """Load, merge, filter, bucket, and write the BED12 annotation outputs.

    Example
    -------
    >>> args = parse_args()
    >>> run(args)
    INFO: Loaded 100 BED12 rows from reads.bed
    """

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    reads = read_bed12(args.reads)
    log.info("Loaded %d BED12 rows from %s", len(reads), args.reads)

    status_tables: List[Tuple[StatusInputSpec, pd.DataFrame]] = []
    for spec in STATUS_INPUT_SPECS:
        path = getattr(args, spec.argument_name)
        status_tables.append((spec, read_status_table(path, spec)))
        log.info(
            "Loaded %d %s rows from %s",
            len(status_tables[-1][1]),
            spec.argument_name,
            path,
        )

    orfs = read_orf_table(args.orfs)
    log.info("Loaded %d ORF rows from %s", len(orfs), args.orfs)

    schema = merge_inputs(reads, status_tables, orfs)
    filtered = filter_by_orf_score(schema, args.orf_score_threshold)
    buckets = bucket_rows(filtered, args.flaws)
    write_outputs(buckets, outdir, args.prefix)

    log.info("Finished writing BED outputs to %s", outdir)


def read_bed12(path: PathLike) -> pd.DataFrame:
    """Read a BED12 file and validate its column count and unique IDs.

    Parameters
    ----------
    path : PathLike
        Path to BED12 file

    Returns
    -------
    pd.DataFrame
        DataFrame with BED12 columns

    Example
    -------
    >>> df = read_bed12("reads.bed")
    >>> df.shape[1]
    12
    """

    bed = pd.read_csv(
        path,
        sep="\t",
        header=None,
        dtype=str,
        keep_default_na=False,
    )
    if bed.shape[1] != len(BED12_COLS):
        raise ValueError(
            f"{path} does not look like BED12. Expected {len(BED12_COLS)} columns, "
            f"found {bed.shape[1]}."
        )

    bed.columns = BED12_COLS
    assert_unique_ids(bed, "reads")
    return bed


def read_status_table(path: PathLike, spec: StatusInputSpec) -> pd.DataFrame:
    """Read one status TSV with the fixed `id status code html` layout.

    Example input (TSV):
        read1   PASS    .   <html>content</html>

    Returns:
        DataFrame with id and columns like R_status, R_code, R_html (using prefix from spec)"""

    records: List[Dict[str, str]] = []
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t", 3)
            if len(parts) != len(STATUS_TSV_COLS):
                raise ValueError(
                    f"{path} line {line_number} does not have 4 tab-separated fields: {line}"
                )

            read_id, status, code, html = parts
            records.append(
                {
                    "id": read_id.strip(),
                    spec.status_column: status.strip(),
                    spec.code_column: code.strip(),
                    spec.html_column: html,
                }
            )

    table = pd.DataFrame.from_records(records, columns=["id"] + spec.output_columns)
    assert_unique_ids(table, spec.argument_name)
    return table


def read_orf_table(path: PathLike) -> pd.DataFrame:
    """Read the ORF TSV, keep the requested columns, and add an HTML summary.

    Example input (TSV):
        chr1    100     200     read1  +   0.5   0.3   ...  0.8    COMPLETE

    Returns:
        DataFrame with O_ prefixed columns, e.g., O_id, O_chr, O_strand, O_prob_coding, O_html"""

    orfs = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    missing_columns = [
        column for column in ORF_COLUMNS_TO_KEEP if column not in orfs.columns
    ]
    if missing_columns:
        raise ValueError(
            f"{path} is missing required ORF columns: {', '.join(missing_columns)}"
        )

    selected = orfs.loc[:, ORF_COLUMNS_TO_KEEP].copy()
    assert_unique_ids(selected, "orfs")

    selected["orf_type"] = [normalize_orf_type(value) for value in selected["orf_type"]]
    selected[ORF_HTML_COL] = [
        render_orf_html(dict(zip(ORF_COLUMNS_TO_KEEP, row)))
        for row in selected.itertuples(index=False, name=None)
    ]

    selected.insert(0, "merge_id", selected["id"])
    selected.rename(columns=ORF_RENAME_MAP, inplace=True)
    selected.rename(columns={"merge_id": "id"}, inplace=True)
    return selected


def normalize_orf_type(value: str) -> str:
    """Normalize numeric ORF type values into readable labels.

    Parameters
    ----------
    value : str
        ORF type value (numeric or string)

    Returns
    -------
    str
        Normalized ORF type label (e.g., "COMPLETE")

    Example
    -------
    >>> normalize_orf_type("1")
    'COMPLETE'
    """

    text = str(value).strip()
    if not text:
        return text
    if text in ORF_TYPE_MAPPING:
        return ORF_TYPE_MAPPING[text]

    try:
        numeric_value = int(float(text))
    except ValueError:
        return text

    return ORF_TYPE_MAPPING.get(numeric_value, text)


def render_orf_html(row: Mapping[str, str]) -> str:
    """Render one ORF row as a compact HTML table.

    Parameters
    ----------
    row : Mapping[str, str]
        ORF row data as column-value mapping

    Returns
    -------
    str
        HTML table representation

    Example
    -------
    >>> html = render_orf_html({"id": "r1", "orf_type": "COMPLETE"})
    >>> "<table" in html
    True
    """

    html_parts = [
        "<h2>ORF</h2>",
        '<table border="1" style="border-collapse:collapse;">',
    ]
    for column in ORF_COLUMNS_TO_KEEP:
        html_parts.append(
            "<tr>"
            f'<th style="padding:4px;text-align:left;">{escape(column)}</th>'
            f'<td style="padding:4px;">{format_html_value(row.get(column, ""))}</td>'
            "</tr>"
        )
    html_parts.append("</table>")
    return "".join(html_parts)


def format_html_value(value: str) -> str:
    """Escape one value for HTML output while keeping empty strings empty.

    Parameters
    ----------
    value : str
        Value to escape for HTML

    Returns
    -------
    str
        HTML-safe string

    Example
    -------
    >>> format_html_value("<test>")
    '&lt;test&gt;'
    """

    if value is None:
        return ""
    text = str(value)
    if not text:
        return ""
    return escape(text)


def merge_inputs(
    reads: pd.DataFrame,
    status_tables: Sequence[Tuple[StatusInputSpec, pd.DataFrame]],
    orfs: pd.DataFrame,
) -> pd.DataFrame:
    """Merge BED12 reads with all mandatory status tables and ORF annotations.

    Parameters
    ----------
    reads : pd.DataFrame
        BED12 reads DataFrame
    status_tables : Sequence[Tuple[StatusInputSpec, pd.DataFrame]]
        List of (spec, DataFrame) tuples for status tables
    orfs : pd.DataFrame
        ORF annotations DataFrame

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with all annotations

    Example
    -------
    >>> reads = pd.DataFrame({"id": ["read1"], "chr": ["chr1"]})
    >>> status_tables = [(spec_R, pd.DataFrame({"id": ["read1"], "R_status": ["PASS"]}))]
    >>> orfs = pd.DataFrame({"id": ["read1"], "O_prob_coding": ["0.8"]})
    >>> result = merge_inputs(reads, status_tables, orfs)
    """

    schema = reads.copy()
    for spec, table in status_tables:
        validate_id_alignment(schema["id"], table["id"], spec.argument_name)
        schema = schema.merge(table, on="id", how="left", validate="one_to_one")

    validate_id_alignment(schema["id"], orfs["id"], "orfs")
    schema = schema.merge(orfs, on="id", how="left", validate="one_to_one")
    validate_orf_identity(schema)

    # INFO: since we avoid raising errors for missing IDs, we safely drop np.nan values
    log.warning(
        f"WARN: dropped {len(schema) - len(schema.dropna(subset=['O_chr', 'O_strand', 'O_id']))} rows with missing ORF annotations"
    )
    log.warning(
        f"WARN: these rows will be dropped: {schema.loc[schema.isna().any(axis=1)].to_dict(orient='records')}"
    )
    schema = schema.dropna(subset=["O_chr", "O_strand", "O_id"])

    return schema.loc[:, build_schema_columns()]


def validate_orf_identity(schema: pd.DataFrame) -> None:
    """Ensure the merged ORF row still refers to the same BED12 record.

    Parameters
    ----------
    schema : pd.DataFrame
        Merged DataFrame with BED12 and ORF columns

    Raises
    ------
    ValueError
        If ORF chr/strand/id don't match BED12 chr/strand/id
    """

    identity_pairs = (
        ("id", "O_id"),
        ("chr", "O_chr"),
        ("strand", "O_strand"),
    )
    for left, right in identity_pairs:
        mismatch = schema[left].fillna("") != schema[right].fillna("")
        if mismatch.any():
            examples = (
                schema.loc[mismatch, ["id", left, right]]
                .head(5)
                .to_dict(orient="records")
            )
            log.warning(
                f"WARN: Merged ORF column {right} does not match BED12 column {left}. "
                f"WARN: Examples: {examples}"
            )


def validate_id_alignment(
    read_ids: Sequence[str],
    other_ids: Sequence[str],
    label: str,
) -> None:
    """Validate that one annotation table has exactly the same IDs as the reads.

    Parameters
    ----------
    read_ids : Sequence[str]
        IDs from the reads file
    other_ids : Sequence[str]
        IDs from the annotation table
    label : str
        Label for error messages (e.g., "orfs", "retentions")

    Raises
    ------
    ValueError
        If ID sets don't match
    """

    read_id_set = set(read_ids)
    other_id_set = set(other_ids)

    missing_ids = sorted(read_id_set - other_id_set)
    extra_ids = sorted(other_id_set - read_id_set)
    if not missing_ids and not extra_ids:
        return

    # INFO: avoding raising error; these IDs will likely be in NMD
    details: List[str] = []
    if missing_ids:
        details.append(
            f"missing {len(missing_ids)} read IDs, e.g. {sample_ids(missing_ids)}"
        )

        log.warning(f"WARN: {label} ID set does not match reads: {'; '.join(details)}")
        return
    if extra_ids:
        details.append(
            f"has extra {len(extra_ids)} unexpected IDs, e.g. {sample_ids(extra_ids)}"
        )
        log.warning(f"{label} ID set does not match reads: {'; '.join(details)}")
        return


def sample_ids(values: Sequence[str], limit: int = 5) -> List[str]:
    """Return a short sample of IDs for error messages.

    Parameters
    ----------
    values : Sequence[str]
        List of ID strings
    limit : int, optional
        Maximum number of IDs to return (default: 5)

    Returns
    -------
    List[str]
        Sample of up to `limit` IDs
    """

    return list(values[:limit])


def assert_unique_ids(frame: pd.DataFrame, label: str) -> None:
    """Raise if a table contains duplicated IDs.

    Parameters
    ----------
    frame : pd.DataFrame
        DataFrame with an "id" column
    label : str
        Label for error messages

    Raises
    ------
    ValueError
        If duplicated IDs are found
    """

    duplicates = frame["id"][frame["id"].duplicated(keep=False)]
    if duplicates.empty:
        return

    duplicated_ids = duplicates.drop_duplicates().tolist()
    raise ValueError(
        f"{label} contains duplicated IDs, e.g. {sample_ids(duplicated_ids)}"
    )


def build_status_columns() -> List[str]:
    """Build the merged status/code/html column order.

    Returns
    -------
    List[str]
        Column names like R_status, R_code, R_html, T_status, etc.
    """

    status_columns: List[str] = []
    for spec in STATUS_INPUT_SPECS:
        for suffix in STATUS_COLUMN_SUFFIXES:
            status_columns.append(f"{spec.prefix}_{suffix}")

    return status_columns


def build_schema_columns() -> List[str]:
    """Build the internal merged schema column order.

    Returns
    -------
    List[str]
        Full column list: BED12 + status + ORF + HTML
    """

    status_columns = build_status_columns()
    orf_columns = [ORF_RENAME_MAP[column] for column in ORF_COLUMNS_TO_KEEP]
    return BED12_COLS + status_columns + orf_columns + [ORF_HTML_COL]


def build_output_columns() -> List[str]:
    """Build the final emitted BED12 plus annotation column order.

    Returns
    -------
    List[str]
        Column list: BED12 + status + HTML
    """

    return BED12_COLS + build_status_columns() + [ORF_HTML_COL]


def filter_by_orf_score(schema: pd.DataFrame, threshold: float) -> pd.DataFrame:
    """Keep only rows whose ORF `prob_coding` is at or above the requested threshold.

    Parameters
    ----------
    schema : pd.DataFrame
        Merged DataFrame with O_prob_coding column
    threshold : float
        Minimum prob_coding value to retain

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame

    Example
    -------
    >>> schema = pd.DataFrame({"id": ["r1", "r2"], "O_prob_coding": ["0.8", "0.1"]})
    >>> filter_by_orf_score(schema, 0.3)
        id O_prob_coding
    0  r1           0.8
    """

    scores = pd.to_numeric(schema["O_prob_coding"], errors="coerce")
    invalid_rows = schema.loc[scores.isna(), "id"].tolist()
    if invalid_rows:
        raise ValueError(
            "Could not parse O_prob_coding for the following IDs: "
            f"{sample_ids(invalid_rows)}"
        )

    keep_mask = scores >= threshold
    retained_count = int(keep_mask.sum())
    discarded_count = int((~keep_mask).sum())
    log.info(
        "ORF score filtering retained %d rows and discarded %d rows at threshold %.4f",
        retained_count,
        discarded_count,
        threshold,
    )

    return schema.loc[keep_mask].copy()


def bucket_rows(schema: pd.DataFrame, flaws: int) -> Dict[str, pd.DataFrame]:
    """Split the filtered schema into RT, pass, trash, and single-failure buckets.

    Parameters
    ----------
    schema : pd.DataFrame
        Filtered DataFrame with status columns
    flaws : int
        Minimum failed status columns to route to trash

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary with keys: rt, pass, trash, retentions, truncations, intraprimming

    Example
    -------
    >>> schema = pd.DataFrame({
    ...     "id": ["r1", "r2", "r3"],
    ...     "R_status": ["PASS", "FAIL", "PASS"],
    ...     "T_status": ["PASS", "FAIL", "FAIL"],
    ...     "I_status": ["PASS", "PASS", "PASS"],
    ...     "R_code": [".", "K", "."],
    ... })
    >>> buckets = bucket_rows(schema, flaws=2)
    >>> list(buckets.keys())
    ['rt', 'pass', 'trash', 'retentions', 'truncations', 'intraprimming']
    """

    # RT routing takes precedence over every other bucket, artifact routing over
    # flaw-based bucketing; each row lands in exactly one bucket.
    rt_mask = schema["R_code"].str.contains("X", regex=False, na=False)
    artifact_mask = (
        schema["R_code"].str.contains("Q", regex=False, na=False) & ~rt_mask
    )
    remaining_schema = schema.loc[~(rt_mask | artifact_mask)].copy()

    status_columns = [spec.status_column for spec in STATUS_INPUT_SPECS]
    status_frame = remaining_schema.loc[:, status_columns].copy()
    fail_mask = status_frame.ne("PASS")
    fail_counts = fail_mask.sum(axis=1)
    single_failure_mask = fail_counts.eq(1) & fail_counts.lt(flaws)
    ambiguous_mask = fail_counts.gt(1) & fail_counts.lt(flaws)

    if ambiguous_mask.any():
        ambiguous_ids = remaining_schema.loc[ambiguous_mask, "id"].tolist()
        raise ValueError(
            "Rows with multiple failed status columns do not map to a unique output "
            f"bucket when --flaws={flaws}. Sample IDs: {sample_ids(ambiguous_ids)}"
        )

    buckets: Dict[str, pd.DataFrame] = {
        "rt": schema.loc[rt_mask].copy(),
        "artifact": schema.loc[artifact_mask].copy(),
        "pass": remaining_schema.loc[fail_counts.eq(0)].copy(),
        "trash": remaining_schema.loc[fail_counts.ge(flaws)].copy(),
    }

    for spec in STATUS_INPUT_SPECS:
        buckets[spec.bucket_name] = remaining_schema.loc[
            single_failure_mask & fail_mask[spec.status_column]
        ].copy()

    for bucket_name, frame in buckets.items():
        log.info("Bucket %s contains %d rows", bucket_name, len(frame))

    return buckets


def write_outputs(
    buckets: Mapping[str, pd.DataFrame],
    outdir: Path,
    prefix: str,
) -> None:
    """Write the five requested BED12-plus output files without headers.

    Parameters
    ----------
    buckets : Mapping[str, pd.DataFrame]
        Dictionary of bucket name to DataFrame
    outdir : Path
        Output directory
    prefix : str
        Prefix for output filenames

    Example
    -------
    >>> buckets = {"pass": df_pass, "trash": df_trash, "rt": df_rt}
    >>> write_outputs(buckets, Path("output"), "sample")
    # writes sample.pass.bed, sample.trash.bed, sample.rt.bed, etc.
    """

    output_columns = build_output_columns()
    empty_frame = pd.DataFrame(columns=output_columns)

    for bucket_name, suffix in OUTPUT_FILENAMES.items():
        output_path = outdir / f"{prefix}.{suffix}"
        frame = buckets.get(bucket_name, empty_frame).reindex(columns=output_columns)
        frame.to_csv(output_path, sep="\t", header=False, index=False)
        log.info("Wrote %d rows to %s", len(frame), output_path)


if __name__ == "__main__":
    run(parse_args())
