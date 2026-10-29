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
