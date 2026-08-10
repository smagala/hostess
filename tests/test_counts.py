"""Counts come from the flagstat that is generated for MultiQC anyway."""

from __future__ import annotations

from pathlib import Path

import pytest

from hostess import counts

# Real `samtools flagstat` output, trimmed to the lines that matter.
FLAGSTAT = """\
1000 + 0 in total (QC-passed reads + QC-failed reads)
900 + 0 primary
60 + 0 secondary
40 + 0 supplementary
0 + 0 duplicates
1000 + 0 mapped (100.00% : N/A)
900 + 0 primary mapped (100.00% : N/A)
900 + 0 paired in sequencing
450 + 0 read1
450 + 0 read2
"""


def test_excludes_secondary_and_supplementary(tmp_path: Path):
    """Those are extra alignments of reads already counted, not extra reads.

    Equivalent to hostile's `samtools view -F 2304 -c` (256 + 2048).
    """
    path = tmp_path / "f.txt"
    path.write_text(FLAGSTAT)
    assert counts.host_reads_from_flagstat(path) == 1000 - 60 - 40


def test_counts_qc_failed_reads_too(tmp_path: Path):
    """Both columns count: a QC-failed read was still removed."""
    path = tmp_path / "f.txt"
    path.write_text(
        "100 + 5 in total (QC-passed reads + QC-failed reads)\n"
        "0 + 0 secondary\n"
        "0 + 0 supplementary\n"
    )
    assert counts.host_reads_from_flagstat(path) == 105


def test_zero_host_reads(tmp_path: Path):
    path = tmp_path / "f.txt"
    path.write_text(
        "0 + 0 in total (QC-passed reads + QC-failed reads)\n"
        "0 + 0 secondary\n"
        "0 + 0 supplementary\n"
    )
    assert counts.host_reads_from_flagstat(path) == 0


def test_missing_line_is_an_error(tmp_path: Path):
    path = tmp_path / "f.txt"
    path.write_text("0 + 0 in total (QC-passed reads + QC-failed reads)\n")
    with pytest.raises(ValueError, match="secondary"):
        counts.host_reads_from_flagstat(path)


def test_impossible_counts_are_rejected(tmp_path: Path):
    path = tmp_path / "f.txt"
    path.write_text(
        "10 + 0 in total (QC-passed reads + QC-failed reads)\n"
        "20 + 0 secondary\n"
        "0 + 0 supplementary\n"
    )
    with pytest.raises(ValueError, match="exceeds total"):
        counts.host_reads_from_flagstat(path)


def test_fastq_read_count(tmp_path: Path):
    import gzip

    path = tmp_path / "r.fastq.gz"
    with gzip.open(path, "wt") as handle:
        for i in range(3):
            handle.write(f"@r{i}\nACGT\n+\nIIII\n")
    assert counts.reads_in_fastq(path) == 3


def test_truncated_fastq_is_an_error(tmp_path: Path):
    import gzip

    path = tmp_path / "r.fastq.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("@r0\nACGT\n+\n")
    with pytest.raises(ValueError, match="multiple of 4"):
        counts.reads_in_fastq(path)
