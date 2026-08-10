"""Derive read counts from artefacts the run already produced.

hostile counts by inserting two `tee >(samtools view -c)` stages into the
alignment pipeline, which means the stream is duplicated twice purely to be
counted. Here the host count is read out of the flagstat that is generated
anyway for MultiQC, so counting costs nothing extra.
"""

from __future__ import annotations

import re
from pathlib import Path

# flagstat's first line, e.g. "1234 + 0 in total (QC-passed reads + QC-failed reads)"
_TOTAL = re.compile(r"^(\d+)\s*\+\s*(\d+)\s+in total\b")

# Secondary and supplementary records are alignments of reads already counted,
# so they must be excluded to get a count of *reads*.
_SECONDARY = re.compile(r"^(\d+)\s*\+\s*(\d+)\s+secondary\b")
_SUPPLEMENTARY = re.compile(r"^(\d+)\s*\+\s*(\d+)\s+supplementary\b")


def host_reads_from_flagstat(flagstat: Path) -> int:
    """Number of host reads, from `samtools flagstat` output.

    Equivalent to hostile's `samtools view -F 2304 -c`, which excludes the
    same two record classes (256 secondary + 2048 supplementary = 2304), but
    read from a file that was produced for MultiQC regardless.
    """
    text = flagstat.read_text()
    total = _first(_TOTAL, text, "total")
    secondary = _first(_SECONDARY, text, "secondary")
    supplementary = _first(_SUPPLEMENTARY, text, "supplementary")

    reads = total - secondary - supplementary
    if reads < 0:
        raise ValueError(
            f"{flagstat}: secondary+supplementary ({secondary + supplementary}) "
            f"exceeds total ({total})"
        )
    return reads


def _first(pattern: re.Pattern[str], text: str, label: str) -> int:
    """Sum the QC-passed and QC-failed columns of a flagstat line."""
    for line in text.splitlines():
        match = pattern.match(line.strip())
        if match:
            return int(match.group(1)) + int(match.group(2))
    raise ValueError(f"no '{label}' line found in flagstat output")


def host_reads_from_count_file(path: Path) -> int:
    """Read the count written by `samtools view -c` when no BAM was requested."""
    text = path.read_text().strip()
    if not text:
        raise ValueError(f"{path}: empty count file")
    try:
        return int(text.splitlines()[0])
    except ValueError as exc:
        raise ValueError(f"{path}: not a read count: {text[:40]!r}") from exc


def reads_in_fastq(path: Path) -> int:
    """Count reads in a gzipped FASTQ.

    Only used for single-end input, where there is no mate to infer from.
    """
    import gzip

    lines = 0
    with gzip.open(path, "rb") as handle:
        for _ in handle:
            lines += 1
    if lines % 4:
        raise ValueError(f"{path}: {lines} lines is not a multiple of 4")
    return lines // 4
