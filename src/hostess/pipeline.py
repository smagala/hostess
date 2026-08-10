"""Build the shell pipeline that does the dehosting.

A read is aligned once and the resulting alignment stream is consumed once.
Everything downstream is derived from that single pass:

    aligner ─┬─ selected: reads that hit nothing  ─→ samtools fastq ─→ clean FASTQs
             └─ unselected: reads that hit the host ─→ counted, and optionally
                                                       sorted into a BAM + BAI

The split is `samtools view -U`, which writes records matching the filter to
its normal output and everything else elsewhere. That is one pass over one
stream, as opposed to `tee`, which duplicates the entire stream so each
consumer can filter it again.

The unselected side is a FIFO rather than a temporary file, so the host reads
are never written twice. The consumer runs in the background and is waited on
explicitly, because bash does not reliably propagate the exit status of a
process substitution.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# Both mates unmapped: the pair aligned to nothing in the host index.
FLAG_PAIR_BOTH_UNMAPPED = 12

# Read unmapped, for single-end and long-read input where there is no mate.
FLAG_READ_UNMAPPED = 4

# Secondary (256) + supplementary (2048). These are additional alignments of
# reads already counted, so they are excluded when counting reads.
FLAG_NOT_PRIMARY = 2304

BOWTIE2 = "bowtie2"
MINIMAP2 = "minimap2"


@dataclass(frozen=True)
class Paths:
    """Where every artefact of a single sample's run goes."""

    clean_1: Path
    clean_2: Path | None
    report: Path
    host_count: Path
    host_bam: Path | None
    host_bai: Path | None
    host_flagstat: Path | None
    host_stats: Path | None

    @classmethod
    def for_sample(
        cls, output: Path, prefix: str, *, paired: bool, bam: bool
    ) -> "Paths":
        host_bam = output / f"{prefix}.host.bam" if bam else None
        return cls(
            clean_1=output
            / (f"{prefix}.clean_1.fastq.gz" if paired else f"{prefix}.clean.fastq.gz"),
            clean_2=output / f"{prefix}.clean_2.fastq.gz" if paired else None,
            report=output / f"{prefix}.hostess.json",
            host_count=output / f"{prefix}.host.count.txt",
            host_bam=host_bam,
            host_bai=output / f"{prefix}.host.bam.bai" if bam else None,
            host_flagstat=output / f"{prefix}.host.flagstat.txt" if bam else None,
            host_stats=output / f"{prefix}.host.stats.txt" if bam else None,
        )


def choose_aligner(*, fastq2: Path | None, requested: str = "auto") -> str:
    """Resolve `auto` the way hostile does, and for the same reason.

    hostile keys off how many FASTQs were given: two means short paired reads
    and bowtie2, one means long reads and minimap2. That is a guess about read
    length made from file count, and it is wrong for single-end Illumina, which
    is why hostile's own help says to override with bowtie2 for those. The
    behaviour is matched so that an existing hostile invocation resolves the
    same way; pass --aligner explicitly to be sure.
    """
    if requested in (BOWTIE2, MINIMAP2):
        return requested
    if requested != "auto":
        raise ValueError(f"unknown aligner: {requested}")
    return BOWTIE2 if fastq2 is not None else MINIMAP2


def build_command(
    *,
    aligner: str,
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
    """Return the bash pipeline for one sample."""
    paired = fastq2 is not None
    # A long-read alignment has no mate to consider, so the flag that means
    # "this fragment hit nothing" differs.
    select_flag = FLAG_PAIR_BOTH_UNMAPPED if paired else FLAG_READ_UNMAPPED

    align = _align_command(
        aligner=aligner,
        index=index,
        fastq1=fastq1,
        fastq2=fastq2,
        threads=align_threads,
        aligner_args=aligner_args,
        reorder=reorder,
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

    want_bam = paths.host_bam is not None
    if want_bam:
        # Sort straight from the FIFO. The unselected records are never
        # written unsorted first, so the host reads touch disk once.
        consumer = (
            f"samtools sort -@ {sort_threads} -m {sort_memory} "
            f'-o {_q(paths.host_bam)} - < "$fifo"'
        )
    else:
        # Counting only. -F 2304 excludes secondary and supplementary records
        # so the result is a count of reads, matching how hostile counts.
        consumer = (
            f'samtools view -c -F {FLAG_NOT_PRIMARY} - < "$fifo" '
            f"> {_q(paths.host_count)}"
        )

    lines = [
        "set -o pipefail",
        f"fifo={_q(paths.host_count.with_suffix('.fifo'))}",
        'rm -f "$fifo"',
        'mkfifo "$fifo"',
        # The consumer must be running before the producer opens the FIFO for
        # writing, otherwise the open blocks forever.
        f"( {consumer} ) &",
        "consumer_pid=$!",
        f"{align} | samtools view -u -f {select_flag} -U \"$fifo\" -o - - | {to_fastq}",
        "pipeline_status=$?",
        # Waited on explicitly: a failure here means the host reads are wrong,
        # and bash would otherwise let it pass silently.
        "wait $consumer_pid",
        "consumer_status=$?",
        'rm -f "$fifo"',
        'if [ "$pipeline_status" -ne 0 ]; then exit "$pipeline_status"; fi',
        'if [ "$consumer_status" -ne 0 ]; then exit "$consumer_status"; fi',
    ]

    if want_bam:
        lines += [
            f"samtools index -@ {sort_threads} {_q(paths.host_bam)}",
            # Written verbatim so MultiQC's samtools module parses them with no
            # post-processing; it keys off content signatures, not filenames.
            f"samtools flagstat -@ {sort_threads} {_q(paths.host_bam)} > {_q(paths.host_flagstat)}",
            f"samtools stats -@ {sort_threads} {_q(paths.host_bam)} > {_q(paths.host_stats)}",
        ]

    return "\n".join(lines)


def _align_command(
    *,
    aligner: str,
    index: Path,
    fastq1: Path,
    fastq2: Path | None,
    threads: int,
    aligner_args: str,
    reorder: bool,
) -> str:
    extra = f" {aligner_args.strip()}" if aligner_args.strip() else ""

    if aligner == MINIMAP2:
        if fastq2 is not None:
            # Short paired reads through minimap2 is what hostile does when
            # asked; -ax sr is its short-read preset.
            preset = "-ax sr"
            reads = f"{_q(fastq1)} {_q(fastq2)}"
        else:
            preset = "-ax map-ont"
            reads = _q(fastq1)
        # --secondary no keeps one alignment per read, which is all a yes/no
        # host decision needs, and keeps the BAM free of duplicate records.
        return (
            f"minimap2 {preset} --secondary no -t {threads}{extra} "
            f"{_q(index)} {reads}"
        )

    if aligner != BOWTIE2:
        raise ValueError(f"unknown aligner: {aligner}")

    # -k 1: one alignment per read is all a yes/no host decision needs.
    # --mm: memory-map the index so concurrent tasks on one machine share it
    #       through the page cache instead of each loading a private copy.
    reads = (
        f"-1 {_q(fastq1)} -2 {_q(fastq2)}" if fastq2 is not None else f"-U {_q(fastq1)}"
    )
    extra += " --reorder" if reorder else ""
    return f"bowtie2 -x {_q(index)} {reads} -k 1 --mm -p {threads}{extra}"


def _q(path: Path | str) -> str:
    """Single-quote a path for bash, escaping any embedded single quotes."""
    text = str(path)
    return "'" + text.replace("'", "'\\''") + "'"
