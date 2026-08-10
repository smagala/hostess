"""The report must stay parsable by MultiQC's stock Hostile module.

These tests encode the contract from multiqc/modules/hostile/hostile.py and
multiqc/search_patterns.yaml, so that a change to our output that would break
an existing MultiQC report fails here first.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from hostess import report


def _entry(**overrides):
    kwargs = dict(
        version="0.1.0",
        index=Path("/refs/human-t2t-hla"),
        fastq1=Path("/data/SAMPLE_R1.fastq.gz"),
        fastq2=Path("/data/SAMPLE_R2.fastq.gz"),
        reads_in=1000,
        reads_removed=250,
    )
    kwargs.update(overrides)
    return report.build_entry(**kwargs)


def test_search_pattern_string_present(tmp_path: Path):
    """MultiQC finds the file by the literal key, within the first 100 lines."""
    path = tmp_path / "sample.hostess.json"
    report.write(path, [_entry()])

    head = "\n".join(path.read_text().splitlines()[:100])
    assert '"reads_removed_proportion"' in head


def test_top_level_is_a_list(tmp_path: Path):
    """The module iterates the parsed JSON; a bare object would break it."""
    path = tmp_path / "sample.hostess.json"
    report.write(path, [_entry()])
    assert isinstance(json.loads(path.read_text()), list)


@pytest.mark.parametrize(
    "key",
    [
        "fastq1_in_name",  # becomes the sample name
        "index",           # basename shown as the database
        "reads_in",
        "reads_out",
        "reads_removed",
        "reads_removed_proportion",
        "version",         # reported as the software version
    ],
)
def test_required_key_present(key: str):
    assert key in _entry()


def test_counts_satisfy_the_invariant_multiqc_checks():
    """hostile.py warns when reads_out + reads_removed != reads_in."""
    entry = _entry(reads_in=1000, reads_removed=250)
    assert entry["reads_out"] + entry["reads_removed"] == entry["reads_in"]


def test_proportion_matches_counts():
    entry = _entry(reads_in=1000, reads_removed=250)
    assert entry["reads_removed_proportion"] == 0.25


def test_zero_reads_does_not_divide_by_zero():
    entry = _entry(reads_in=0, reads_removed=0)
    assert entry["reads_removed_proportion"] == 0.0
    assert entry["reads_out"] == 0


def test_nothing_removed():
    entry = _entry(reads_in=500, reads_removed=0)
    assert entry["reads_out"] == 500
    assert entry["reads_removed_proportion"] == 0.0


def test_everything_removed():
    entry = _entry(reads_in=500, reads_removed=500)
    assert entry["reads_out"] == 0
    assert entry["reads_removed_proportion"] == 1.0


def test_inconsistent_counts_are_rejected():
    """Better to fail here than to emit a report MultiQC warns about."""
    with pytest.raises(ValueError, match="exceeds reads_in"):
        _entry(reads_in=10, reads_removed=11)


def test_single_end_omits_mate_fields():
    entry = _entry(fastq2=None)
    assert "fastq2_in_name" not in entry
    assert "fastq2_in_path" not in entry


def test_hostess_extras_do_not_collide_with_multiqc_keys():
    """Extra keys are ignored by the module, but must not shadow its own."""
    entry = _entry(host_bam=Path("/out/SAMPLE.host.bam"))
    multiqc_keys = {
        "fastq1_in_name",
        "index",
        "reads_in",
        "reads_out",
        "reads_removed",
        "reads_removed_proportion",
        "version",
    }
    extras = {"hostess_schema_version", "host_bam"}
    assert not (extras & multiqc_keys)
    assert entry["hostess_schema_version"] == report.SCHEMA_VERSION
