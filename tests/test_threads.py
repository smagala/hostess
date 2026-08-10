"""Thread allocation, constrained by what was measured.

See allocate_threads' docstring for the numbers these rules come from.
"""

from __future__ import annotations

import pytest

from hostess.cli import allocate_threads


@pytest.mark.parametrize("total", [1, 2, 3, 4, 6, 8, 12, 16, 32, 64])
def test_every_stage_gets_at_least_one_thread(total: int):
    align, sort, compress = allocate_threads(total)
    assert align >= 1 and sort >= 1 and compress >= 1


@pytest.mark.parametrize("total", [2, 3, 4, 6, 8, 12, 16, 32, 64])
def test_alignment_never_drops_to_one_thread(total: int):
    """hostile's split gives bowtie2 a single thread at 2 CPUs.

    That is what made a 62M read-pair sample take half an hour.
    """
    align, _, _ = allocate_threads(total)
    assert align >= 2


@pytest.mark.parametrize("total", [8, 12, 16, 32])
def test_budget_is_not_exceeded(total: int):
    """Oversubscribing a Batch task's CPU allocation just causes contention."""
    align, sort, compress = allocate_threads(total)
    assert align + sort + compress <= total


@pytest.mark.parametrize("total", [8, 12, 16, 32])
def test_compression_is_not_starved(total: int):
    """Compression at one thread was the single biggest avoidable cost.

    1 thread took 2456 ms against 772 ms at four, on identical input.
    """
    _, _, compress = allocate_threads(total)
    assert compress >= 2


@pytest.mark.parametrize("total", [8, 16, 32])
def test_alignment_gets_the_largest_share(total: int):
    align, sort, compress = allocate_threads(total)
    assert align >= compress
    assert align >= sort


def test_more_threads_never_reduces_alignment():
    """Monotonic, so raising a task's CPUs cannot slow alignment down."""
    previous = 0
    for total in range(1, 65):
        align, _, _ = allocate_threads(total)
        assert align >= previous
        previous = align


@pytest.mark.parametrize("bad", [0, -1, -100])
def test_nonsense_input_is_clamped(bad: int):
    align, sort, compress = allocate_threads(bad)
    assert align >= 1 and sort >= 1 and compress >= 1


def test_single_thread_is_honest():
    """With one CPU everything shares it; do not pretend otherwise."""
    assert allocate_threads(1) == (1, 1, 1)
