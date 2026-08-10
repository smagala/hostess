"""The pipeline must read the alignment stream exactly once."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from hostess import pipeline


def _paths(tmp_path: Path, paired: bool = True, bam: bool = True) -> pipeline.Paths:
    return pipeline.Paths.for_sample(tmp_path, "SAMPLE", paired=paired, bam=bam)


def _cmd(tmp_path: Path, **overrides) -> str:
    kwargs = dict(
        index=Path("/refs/human-t2t-hla"),
        fastq1=Path("/data/SAMPLE_R1.fastq.gz"),
        fastq2=Path("/data/SAMPLE_R2.fastq.gz"),
        paths=_paths(tmp_path),
        aligner=pipeline.BOWTIE2,
        align_threads=4,
        sort_threads=2,
        compression_threads=1,
        sort_memory="1G",
    )
    kwargs.update(overrides)
    return pipeline.build_command(**kwargs)


def test_no_tee(tmp_path: Path):
    """The whole point: the stream is split, not duplicated.

    Checked as a command token and as process substitution, rather than as a
    substring: pytest's tmp_path embeds the test name, so a naive `"tee" not
    in cmd` matches its own directory.
    """
    cmd = _cmd(tmp_path)
    assert not re.search(r"(^|[|;&(]\s*)tee\b", cmd), "stream is duplicated with tee"
    assert ">(" not in cmd, "process substitution implies a duplicated stream"


def test_splits_with_output_unselected(tmp_path: Path):
    """One samtools view writes both the kept and the removed reads."""
    cmd = _cmd(tmp_path)
    assert "-U " in cmd
    assert cmd.count("bowtie2 ") == 1


def test_aligns_once(tmp_path: Path):
    assert _cmd(tmp_path).count("bowtie2 -x") == 1


def test_paired_selects_both_mates_unmapped(tmp_path: Path):
    """flag 12 = read unmapped + mate unmapped, i.e. the clean pairs."""
    assert f"-f {pipeline.FLAG_PAIR_BOTH_UNMAPPED} " in _cmd(tmp_path)


def test_single_end_selects_read_unmapped(tmp_path: Path):
    cmd = _cmd(tmp_path, fastq2=None, paths=_paths(tmp_path, paired=False))
    assert f"-f {pipeline.FLAG_READ_UNMAPPED} " in cmd
    assert "-1 " not in cmd.split("bowtie2")[1].split("|")[0]


def test_reorder_is_off_by_default(tmp_path: Path):
    """--reorder costs time and memory; nothing downstream needs read order."""
    assert "--reorder" not in _cmd(tmp_path)


def test_reorder_can_be_enabled(tmp_path: Path):
    assert "--reorder" in _cmd(tmp_path, reorder=True)


def test_produces_sorted_bam_and_index(tmp_path: Path):
    cmd = _cmd(tmp_path)
    assert "samtools sort" in cmd
    assert "samtools index" in cmd


def test_emits_multiqc_parsable_bam_stats(tmp_path: Path):
    """flagstat and stats are what MultiQC's samtools module reads."""
    cmd = _cmd(tmp_path)
    assert "samtools flagstat" in cmd
    assert "samtools stats" in cmd


def test_pipefail_is_set(tmp_path: Path):
    """Without this a failing bowtie2 is masked by a succeeding samtools."""
    assert _cmd(tmp_path).startswith("set -o pipefail")


def test_unsorted_intermediate_is_removed(tmp_path: Path):
    cmd = _cmd(tmp_path)
    assert "rm -f" in cmd


def test_thread_counts_are_honoured(tmp_path: Path):
    cmd = _cmd(tmp_path, align_threads=8, sort_threads=3)
    assert "-p 8" in cmd
    assert "-@ 3" in cmd


def test_paths_with_spaces_are_quoted(tmp_path: Path):
    """Staged filenames can contain spaces; the pipeline is a shell string."""
    cmd = _cmd(tmp_path, fastq1=Path("/data/my sample_R1.fastq.gz"))
    assert "'/data/my sample_R1.fastq.gz'" in cmd


def test_single_quotes_in_paths_are_escaped(tmp_path: Path):
    cmd = _cmd(tmp_path, fastq1=Path("/data/o'brien_R1.fastq.gz"))
    assert "'/data/o'\\''brien_R1.fastq.gz'" in cmd


def test_paired_paths_name_both_mates(tmp_path: Path):
    paths = _paths(tmp_path)
    assert paths.clean_1.name.endswith("clean_1.fastq.gz")
    assert paths.clean_2 is not None


def test_single_end_paths_have_no_mate(tmp_path: Path):
    paths = _paths(tmp_path, paired=False)
    assert paths.clean_2 is None
    assert paths.clean_1.name.endswith(".clean.fastq.gz")
