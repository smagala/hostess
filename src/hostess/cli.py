"""Command line interface.

The signature follows `hostile clean` where that makes sense, so an existing
invocation mostly transfers. It does not follow it where the behaviour differs:

  --airplane   accepted and ignored. hostile uses it to disable downloading a
               reference index. hostess never downloads one, so it is always
               effectively on. Accepted so existing command lines do not break.
  --invert     not supported. hostile's inversion exists because it discards
               the host reads; --bam keeps them, which is strictly more useful
               than swapping which set is thrown away. Rejected with a message
               pointing at --bam rather than silently ignored.
  --bam        new. Off by default, so the default output matches hostile's.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from hostess import __version__, counts, pipeline, report


def allocate_threads(total: int) -> tuple[int, int, int]:
    """Split the thread budget between alignment, sorting and compression.

    Two measurements shaped this, both on 533k read pairs:

    Compression is not free and starving it dominates the run. Holding
    everything else constant and varying only `samtools fastq --threads`:

        1 thread   2456 ms
        2 threads  1286 ms
        3 threads   936 ms
        4 threads   772 ms
        6 threads   696 ms

    So compression gets a real share rather than whatever is left over, and
    the split targets roughly 50 percent alignment, 35 percent compression,
    15 percent sorting.

    Dropping the gzip level is the wrong lever. At one thread, `-c 1` ran in
    679 ms against `-c 6` at 2456 ms, but four threads at `-c 6` reached
    772 ms while keeping the output 15 percent smaller. Threads buy the same
    speed without the size penalty.

    Alignment never drops below two threads. hostile's 2:1 split yields
    `(2 * 2) // 3 == 1` alignment thread at two CPUs, which is what made a
    62M read-pair sample take half an hour.
    """
    total = max(1, int(total))
    if total == 1:
        return 1, 1, 1
    if total <= 3:
        return 2, 1, 1

    align = max(2, total // 2)
    remaining = total - align
    compress = max(1, (remaining * 7) // 10)
    sort = max(1, remaining - compress)
    return align, sort, compress


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="hostess",
        description=(
            "Dehost reads in a single alignment pass, optionally keeping the "
            "removed reads as a coordinate-sorted BAM."
        ),
    )
    parser.add_argument("--version", action="version", version=f"hostess {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    # Named `clean` so `hostess clean ...` mirrors `hostile clean ...`.
    clean = sub.add_parser("clean", help="remove host reads")
    clean.add_argument("--fastq1", required=True, type=Path)
    clean.add_argument("--fastq2", type=Path, help="omit for single-end or long reads")
    clean.add_argument(
        "--aligner",
        choices=["bowtie2", "minimap2", "auto"],
        default="auto",
        help=(
            "auto resolves as hostile does: bowtie2 given two FASTQs, minimap2 "
            "given one. That infers read length from file count and is wrong "
            "for single-end short reads, so set it explicitly for those."
        ),
    )
    clean.add_argument(
        "--index",
        required=True,
        type=Path,
        help="bowtie2 index prefix, or a reference FASTA/MMI for minimap2",
    )
    clean.add_argument("-o", "--output", type=Path, default=Path("."))
    clean.add_argument("--prefix", help="output basename (default: from --fastq1)")
    clean.add_argument("-t", "--threads", type=int, default=4)
    clean.add_argument("--aligner-args", default="", help="extra aligner arguments")
    clean.add_argument("--force", action="store_true", help="overwrite existing output")
    clean.add_argument(
        "--bam",
        action="store_true",
        help=(
            "also write the removed reads as a coordinate-sorted, indexed BAM "
            "with flagstat and stats. Off by default."
        ),
    )
    clean.add_argument("--sort-memory", default="768M", help="per sorting thread")
    clean.add_argument(
        "--reorder",
        action="store_true",
        help=(
            "make FASTQ output order deterministic (bowtie2 only). Off by "
            "default: --reorder buffers alignments and raised measured peak "
            "RSS from under 12 GB to 43.6 GB at four alignment threads."
        ),
    )
    clean.add_argument(
        "--airplane",
        action="store_true",
        help="accepted for hostile compatibility and ignored; hostess never downloads an index",
    )
    clean.add_argument("--invert", action="store_true", help=argparse.SUPPRESS)
    clean.add_argument("--dry-run", action="store_true", help="print the pipeline and exit")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.invert:
        parser.error(
            "--invert is not supported. hostile inverts because it discards the "
            "host reads; use --bam to keep them as a sorted BAM instead."
        )

    aligner = pipeline.choose_aligner(fastq2=args.fastq2, requested=args.aligner)

    if not args.dry_run:
        for tool in (aligner, "samtools"):
            if shutil.which(tool) is None:
                parser.error(f"{tool} not found on PATH")
        if not args.fastq1.exists():
            parser.error(f"--fastq1 does not exist: {args.fastq1}")
        if args.fastq2 is not None and not args.fastq2.exists():
            parser.error(f"--fastq2 does not exist: {args.fastq2}")

    if args.reorder and aligner != pipeline.BOWTIE2:
        parser.error("--reorder applies to bowtie2 only")

    prefix = args.prefix or _stem(args.fastq1)
    paired = args.fastq2 is not None
    args.output.mkdir(parents=True, exist_ok=True)
    paths = pipeline.Paths.for_sample(args.output, prefix, paired=paired, bam=args.bam)

    if not args.force and not args.dry_run and paths.clean_1.exists():
        parser.error(f"{paths.clean_1} exists. Use --force to overwrite.")

    align, sort, compress = allocate_threads(args.threads)
    command = pipeline.build_command(
        aligner=aligner,
        index=args.index,
        fastq1=args.fastq1,
        fastq2=args.fastq2,
        paths=paths,
        align_threads=align,
        sort_threads=sort,
        compression_threads=compress,
        sort_memory=args.sort_memory,
        aligner_args=args.aligner_args,
        reorder=args.reorder,
    )

    if args.dry_run:
        print(command)
        return 0

    result = subprocess.run(["bash", "-c", command], text=True)
    if result.returncode != 0:
        print(f"hostess: alignment pipeline failed ({result.returncode})", file=sys.stderr)
        return result.returncode

    if args.bam:
        reads_removed = counts.host_reads_from_flagstat(paths.host_flagstat)
    else:
        reads_removed = counts.host_reads_from_count_file(paths.host_count)
        paths.host_count.unlink(missing_ok=True)

    clean = counts.reads_in_fastq(paths.clean_1)
    # For paired input each fragment contributes two reads and the clean FASTQs
    # hold one mate each, so the read total is twice one file's count.
    reads_out = clean * 2 if paired else clean

    options = []
    if args.reorder:
        options.append("reorder")
    if args.bam:
        options.append("bam")

    entry = report.build_entry(
        version=__version__,
        index=args.index,
        fastq1=args.fastq1,
        fastq2=args.fastq2,
        reads_in=reads_out + reads_removed,
        reads_removed=reads_removed,
        aligner=aligner,
        options=options,
        host_bam=paths.host_bam,
    )
    report.write(paths.report, [entry])

    print(f"hostess: {reads_removed} host reads removed, {reads_out} kept")
    if args.bam:
        print(f"hostess: removed reads written to {paths.host_bam}")
    return 0


def _stem(path: Path) -> str:
    """Strip .fastq/.fq and .gz, and any trailing _R1/_1 mate marker."""
    name = path.name
    for suffix in (".gz", ".bz2"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    for suffix in (".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    for marker in ("_R1", "_r1", "_1"):
        if name.endswith(marker):
            return name[: -len(marker)]
    return name


if __name__ == "__main__":
    raise SystemExit(main())
