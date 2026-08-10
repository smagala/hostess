"""Build the shell pipeline that does the dehosting.

The whole point of hostess is that a read is aligned once and the resulting
alignment stream is consumed once. Everything downstream is derived from that
single pass:

    bowtie2 ─┬─ (selected: unaligned pairs) ─→ samtools fastq ─→ clean FASTQs
             └─ (unselected: host reads)    ─→ BAM ─→ sort ─→ BAM + BAI

The split is `samtools view -U`, which writes records matching the filter to
its normal output and everything else to a second file. That is one pass over
one stream, as opposed to `tee`, which duplicates the entire stream so each
consumer can filter it again.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# Both mates unmapped. These reads aligned to nothing in the host index, so
# they are what the caller wants to keep.
FLAG_PAIR_BOTH_UNMAPPED = 12

# Read unmapped, for single-end input where there is no mate to consider.
FLAG_READ_UNMAPPED = 4


@dataclass(frozen=True)
class Paths:
    """Where every artefact of a single sample's run goes."""

    clean_1: Path
    clean_2: Path | None
    host_bam: Path
    host_bai: Path
    host_flagstat: Path
    host_stats: Path
    report: Path

    @classmethod
    def for_sample(cls, output: Path, prefix: str, paired: bool) -> "Paths":
        return cls(
            clean_1=output / (f"{prefix}.clean_1.fastq.gz" if paired else f"{prefix}.clean.fastq.gz"),
            clean_2=output / f"{prefix}.clean_2.fastq.gz" if paired else None,
            host_bam=output / f"{prefix}.host.bam",
            host_bai=output / f"{prefix}.host.bam.bai",
            host_flagstat=output / f"{prefix}.host.flagstat.txt",
            host_stats=output / f"{prefix}.host.stats.txt",
            report=output / f"{prefix}.hostess.json",
        )


def build_command(
    *,
    index: Path,
    fastq1: Path,
    fastq2: Path | None,
    paths: Paths,
    align_threads: int,
    sort_threads: int,
    compression_threads: int,
    sort_memory: str,
    aligner_args: str = "",
    reorder: bool = False,
) -> str:
    """Return the bash pipeline for one sample.

    `reorder` is off by default, unlike hostile. bowtie2's --reorder buffers
    finished alignments so output order matches input order, and it costs both
    time and memory: measured peak RSS rose from under 12 GB at one alignment
    thread to 43.6 GB at four. Nothing downstream here depends on read order,
    because the host reads are coordinate-sorted anyway and the clean FASTQs
    are consumed as a set. Callers that need byte-reproducible FASTQ ordering
    can turn it back on.
    """
    paired = fastq2 is not None
    select_flag = FLAG_PAIR_BOTH_UNMAPPED if paired else FLAG_READ_UNMAPPED

    reads = f"-1 {_q(fastq1)} -2 {_q(fastq2)}" if paired else f"-U {_q(fastq1)}"
    extra = f" {aligner_args}" if aligner_args.strip() else ""
    extra += " --reorder" if reorder else ""

    # -k 1: one alignment per read is all a yes/no host decision needs.
    # --mm: memory-map the index so concurrent tasks on one machine share it
    #       through the page cache instead of each loading a private copy.
    align = (
        f"bowtie2 -x {_q(index)} {reads} -k 1 --mm -p {align_threads}{extra}"
    )

    if paired:
        to_fastq = (
            f"samtools fastq --threads {compression_threads} -c 6 -N "
            f"-1 {_q(paths.clean_1)} -2 {_q(paths.clean_2)} -0 /dev/null -s /dev/null -"
        )
    else:
        to_fastq = (
            f"samtools fastq --threads {compression_threads} -c 6 "
            f"-0 {_q(paths.clean_1)} -"
        )

    # The single split. Selected records (unaligned) go to stdout and on to
    # samtools fastq; unselected records (host) go straight to -U as BAM.
    # Neither stream is duplicated.
    split = (
        f"samtools view -u -f {select_flag} "
        f"-U {_q(paths.host_bam.with_suffix('.unsorted.bam'))} -o - -"
    )

    sort = (
        f"samtools sort -@ {sort_threads} -m {sort_memory} "
        f"-o {_q(paths.host_bam)} {_q(paths.host_bam.with_suffix('.unsorted.bam'))}"
    )

    return "\n".join(
        [
            "set -o pipefail",
            f"{align} | {split} | {to_fastq}",
            sort,
            f"rm -f {_q(paths.host_bam.with_suffix('.unsorted.bam'))}",
            f"samtools index -@ {sort_threads} {_q(paths.host_bam)}",
            # Written verbatim so MultiQC's samtools module parses them with no
            # post-processing; it keys off content signatures, not filenames.
            f"samtools flagstat -@ {sort_threads} {_q(paths.host_bam)} > {_q(paths.host_flagstat)}",
            f"samtools stats -@ {sort_threads} {_q(paths.host_bam)} > {_q(paths.host_stats)}",
        ]
    )


def _q(path: Path | str) -> str:
    """Single-quote a path for bash, escaping any embedded single quotes."""
    text = str(path)
    return "'" + text.replace("'", "'\\''") + "'"
