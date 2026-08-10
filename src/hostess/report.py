"""Emit a report MultiQC parses with no configuration.

MultiQC ships a Hostile module. Its search pattern is any ``*.json`` whose
first 100 lines contain the literal string ``"reads_removed_proportion"``, and
its parser (``multiqc/modules/hostile/hostile.py``) then expects a JSON *list*
of objects with these keys:

    fastq1_in_name            used as the sample name
    index                     basename shown as the database
    reads_in, reads_out, reads_removed
    reads_removed_proportion  also what the search pattern keys off
    version                   optional, reported as the software version

The parser asserts ``reads_out + reads_removed == reads_in`` and warns loudly
if that does not hold, so the counts here are derived rather than measured
independently.

Matching that shape is deliberate: it means hostess output drops into an
existing MultiQC report without a custom module, a plugin, or a config entry.
The BAM statistics are emitted separately as raw ``samtools flagstat`` and
``samtools stats`` output, which MultiQC's samtools module recognises by
content signature.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

SCHEMA_VERSION = 1


def build_entry(
    *,
    version: str,
    index: Path,
    fastq1: Path,
    fastq2: Path | None,
    reads_in: int,
    reads_removed: int,
    aligner: str = "bowtie2",
    options: list[str] | None = None,
    host_bam: Path | None = None,
) -> dict[str, Any]:
    """Build one sample's entry.

    ``reads_out`` is computed rather than counted so the invariant MultiQC
    checks cannot be violated by two counts drifting apart.
    """
    if reads_removed > reads_in:
        raise ValueError(
            f"reads_removed ({reads_removed}) exceeds reads_in ({reads_in}); "
            "counts are inconsistent"
        )

    reads_out = reads_in - reads_removed
    proportion = round(reads_removed / reads_in, 5) if reads_in else 0.0

    entry: dict[str, Any] = {
        "version": version,
        "aligner": aligner,
        "index": str(index),
        "options": options or [],
        "fastq1_in_name": fastq1.name,
        "fastq1_in_path": str(fastq1),
        "reads_in": reads_in,
        "reads_out": reads_out,
        "reads_removed": reads_removed,
        "reads_removed_proportion": proportion,
    }

    if fastq2 is not None:
        entry["fastq2_in_name"] = fastq2.name
        entry["fastq2_in_path"] = str(fastq2)

    # hostess-specific additions. MultiQC's Hostile module reads only the keys
    # it knows and ignores the rest, so these are safe to carry.
    entry["hostess_schema_version"] = SCHEMA_VERSION
    if host_bam is not None:
        entry["host_bam"] = str(host_bam)

    return entry


def write(path: Path, entries: list[dict[str, Any]]) -> None:
    """Write the report as a JSON list, which is the shape MultiQC expects."""
    path.write_text(json.dumps(entries, indent=4) + "\n")
