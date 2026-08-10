"""Aligner selection, long reads, and optional BAM output."""

from __future__ import annotations

from pathlib import Path

import pytest

from hostess import pipeline
from hostess.cli import build_parser


def _cmd(tmp_path: Path, **overrides) -> str:
    paired = overrides.pop("paired", True)
    bam = overrides.pop("bam", False)
    kwargs = dict(
        aligner=pipeline.BOWTIE2,
        index=Path("/refs/host"),
        fastq1=Path("/data/S_R1.fastq.gz"),
        fastq2=Path("/data/S_R2.fastq.gz") if paired else None,
        paths=pipeline.Paths.for_sample(tmp_path, "S", paired=paired, bam=bam),
        align_threads=4,
        sort_threads=2,
        compression_threads=2,
        sort_memory="768M",
    )
    kwargs.update(overrides)
    return pipeline.build_command(**kwargs)


# --- aligner selection -------------------------------------------------

def test_auto_picks_bowtie2_for_paired():
    """Matches hostile: two FASTQs means short paired reads."""
    assert pipeline.choose_aligner(fastq2=Path("r2.fq")) == pipeline.BOWTIE2


def test_auto_picks_minimap2_for_unpaired():
    """Matches hostile: one FASTQ is assumed to be long reads."""
    assert pipeline.choose_aligner(fastq2=None) == pipeline.MINIMAP2


@pytest.mark.parametrize("requested", ["bowtie2", "minimap2"])
def test_explicit_aligner_wins(requested: str):
    """Single-end short reads need this, since auto would choose minimap2."""
    assert pipeline.choose_aligner(fastq2=None, requested=requested) == requested


def test_unknown_aligner_rejected():
    with pytest.raises(ValueError, match="unknown aligner"):
        pipeline.choose_aligner(fastq2=None, requested="bwa")


# --- long reads --------------------------------------------------------

def test_long_read_uses_map_ont_preset(tmp_path: Path):
    cmd = _cmd(tmp_path, aligner=pipeline.MINIMAP2, paired=False)
    assert "minimap2 -ax map-ont" in cmd


def test_long_read_selects_read_unmapped(tmp_path: Path):
    """No mate exists, so the paired both-unmapped flag would never match."""
    cmd = _cmd(tmp_path, aligner=pipeline.MINIMAP2, paired=False)
    assert f"-f {pipeline.FLAG_READ_UNMAPPED} " in cmd


def test_minimap2_paired_uses_short_read_preset(tmp_path: Path):
    cmd = _cmd(tmp_path, aligner=pipeline.MINIMAP2)
    assert "minimap2 -ax sr" in cmd


def test_minimap2_suppresses_secondary(tmp_path: Path):
    """One alignment per read is all a host/not-host decision needs."""
    assert "--secondary no" in _cmd(tmp_path, aligner=pipeline.MINIMAP2, paired=False)


def test_minimap2_does_not_get_bowtie2_flags(tmp_path: Path):
    cmd = _cmd(tmp_path, aligner=pipeline.MINIMAP2, paired=False)
    assert "--mm" not in cmd
    assert "-k 1" not in cmd


# --- optional BAM ------------------------------------------------------

def test_bam_off_by_default_matches_hostile_output(tmp_path: Path):
    paths = pipeline.Paths.for_sample(tmp_path, "S", paired=True, bam=False)
    assert paths.host_bam is None
    cmd = _cmd(tmp_path, bam=False)
    assert "samtools sort" not in cmd
    assert "samtools index" not in cmd


def test_bam_off_still_counts_host_reads(tmp_path: Path):
    """The report needs the count even when the reads are discarded."""
    cmd = _cmd(tmp_path, bam=False)
    assert f"samtools view -c -F {pipeline.FLAG_NOT_PRIMARY}" in cmd


def test_bam_on_produces_sorted_indexed_bam(tmp_path: Path):
    cmd = _cmd(tmp_path, bam=True)
    assert "samtools sort" in cmd
    assert "samtools index" in cmd
    assert "samtools flagstat" in cmd
    assert "samtools stats" in cmd


def test_bam_is_sorted_straight_from_the_fifo(tmp_path: Path):
    """No unsorted intermediate: the host reads touch disk once."""
    cmd = _cmd(tmp_path, bam=True)
    assert "unsorted" not in cmd
    assert 'samtools sort' in cmd and '< "$fifo"' in cmd


# --- single pass, however configured -----------------------------------

@pytest.mark.parametrize("bam", [True, False])
@pytest.mark.parametrize("aligner", [pipeline.BOWTIE2, pipeline.MINIMAP2])
def test_always_exactly_one_alignment(tmp_path: Path, bam: bool, aligner: str):
    cmd = _cmd(tmp_path, bam=bam, aligner=aligner)
    assert cmd.count("bowtie2 -x") + cmd.count("minimap2 -ax") == 1


@pytest.mark.parametrize("bam", [True, False])
def test_consumer_failure_is_not_swallowed(tmp_path: Path, bam: bool):
    """bash does not propagate a background job's status on its own."""
    cmd = _cmd(tmp_path, bam=bam)
    assert "wait $consumer_pid" in cmd
    assert 'exit "$consumer_status"' in cmd


@pytest.mark.parametrize("bam", [True, False])
def test_fifo_is_created_before_the_producer_runs(tmp_path: Path, bam: bool):
    """Opening a FIFO for writing blocks until a reader exists."""
    cmd = _cmd(tmp_path, bam=bam)
    lines = cmd.splitlines()
    mkfifo = next(i for i, line in enumerate(lines) if "mkfifo" in line)
    consumer = next(i for i, line in enumerate(lines) if line.endswith(") &"))
    producer = next(i for i, line in enumerate(lines) if "samtools view -u" in line)
    assert mkfifo < consumer < producer


# --- hostile compatibility --------------------------------------------

def test_clean_subcommand_mirrors_hostile():
    args = build_parser().parse_args(
        ["clean", "--fastq1", "a.fq.gz", "--index", "idx"]
    )
    assert args.command == "clean"


def test_hostile_short_flags_are_accepted():
    """-o and -t are how hostile spells output and threads."""
    args = build_parser().parse_args(
        ["clean", "--fastq1", "a.fq.gz", "--index", "idx", "-o", "out", "-t", "8"]
    )
    assert args.output == Path("out")
    assert args.threads == 8


def test_airplane_accepted_and_inert():
    """hostess never downloads an index, so this is always effectively on."""
    args = build_parser().parse_args(
        ["clean", "--fastq1", "a.fq.gz", "--index", "idx", "--airplane"]
    )
    assert args.airplane is True
