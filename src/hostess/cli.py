"""Command line interface."""

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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="hostess",
        description=(
            "Dehost short reads in a single alignment pass, keeping the removed "
            "reads as a coordinate-sorted BAM for inspection."
        ),
    )
    parser.add_argument("--version", action="version", version=f"hostess {__version__}")
    parser.add_argument("--fastq1", required=True, type=Path)
    parser.add_argument("--fastq2", type=Path, help="omit for single-end input")
    parser.add_argument(
        "--index",
        required=True,
        type=Path,
        help="bowtie2 index prefix, i.e. the path without the .N.bt2 suffix",
    )
    parser.add_argument("--output", type=Path, default=Path("."))
    parser.add_argument("--prefix", help="output basename (default: derived from --fastq1)")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--sort-memory", default="768M", help="per sorting thread")
    parser.add_argument("--aligner-args", default="", help="extra bowtie2 arguments")
    parser.add_argument(
        "--reorder",
        action="store_true",
        help=(
            "make FASTQ output order deterministic. Off by default: bowtie2's "
            "--reorder buffers alignments and raised measured peak memory from "
            "under 12 GB to 43.6 GB at four alignment threads."
        ),
    )
    parser.add_argument("--dry-run", action="store_true", help="print the pipeline and exit")
    args = parser.parse_args(argv)

    for tool in ("bowtie2", "samtools"):
        if shutil.which(tool) is None and not args.dry_run:
            parser.error(f"{tool} not found on PATH")

    if not args.fastq1.exists() and not args.dry_run:
        parser.error(f"--fastq1 does not exist: {args.fastq1}")
    if args.fastq2 is not None and not args.fastq2.exists() and not args.dry_run:
        parser.error(f"--fastq2 does not exist: {args.fastq2}")

    prefix = args.prefix or _stem(args.fastq1)
    paired = args.fastq2 is not None
    args.output.mkdir(parents=True, exist_ok=True)
    paths = pipeline.Paths.for_sample(args.output, prefix, paired=paired)

    align, sort, compress = allocate_threads(args.threads)
    command = pipeline.build_command(
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

    reads_removed = counts.host_reads_from_flagstat(paths.host_flagstat)
    clean = counts.reads_in_fastq(paths.clean_1)
    # For paired input each fragment contributes two reads, and the clean
    # FASTQs hold one mate each, so the read total is twice one file's count.
    reads_out = clean * 2 if paired else clean

    entry = report.build_entry(
        version=__version__,
        index=args.index,
        fastq1=args.fastq1,
        fastq2=args.fastq2,
        reads_in=reads_out + reads_removed,
        reads_removed=reads_removed,
        options=["reorder"] if args.reorder else [],
        host_bam=paths.host_bam,
    )
    report.write(paths.report, [entry])

    print(f"hostess: {reads_removed} host reads removed, {reads_out} kept")
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
