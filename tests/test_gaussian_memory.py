"""%Mem leaves headroom below the allocation, so the scheduler does not kill
Gaussian for exceeding it."""

import pytest

from gaussian_step.substep import gaussian_memory

GB = 1024**3


@pytest.mark.parametrize(
    "allocation, expected",
    [
        (16 * GB, 13.6 * GB),  # 15 % > 2 GB
        (128 * GB, 108.8 * GB),
        (8 * GB, 6 * GB),  # 2 GB > 15 %
        (2 * GB, 1 * GB),  # never below half
    ],
)
def test_gaussian_memory(allocation, expected):
    assert gaussian_memory(allocation) == pytest.approx(expected)


@pytest.mark.parametrize(
    "step, glob, threads, expected",
    [
        # defaults: the threads' share, below the whole allocation's cap
        ("available", "available", 8, 54.4 * GB),
        ("available", "available", 4, 32 * GB),
        # an explicit amount that fits stays (design's example)
        ("48 GB", "60 GB", 4, 48e9),
        # an explicit amount above the allocation is capped below it
        ("100 GB", "available", 8, 54.4 * GB),
        ("all", "all", 8, 54.4 * GB),
    ],
)
def test_gaussian_memory_for_run(step, glob, threads, expected):
    from gaussian_step.substep import gaussian_memory_for_run

    ce = {"MEM_PER_NODE": 64 * GB, "MEM_PER_CPU": 8 * GB}
    assert gaussian_memory_for_run(step, glob, ce, threads) == pytest.approx(expected)
