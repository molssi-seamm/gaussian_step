#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `gaussian_step` package."""

import pytest  # noqa: F401
import gaussian_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = gaussian_step.Gaussian()
    assert str(type(result)) == "<class 'gaussian_step.gaussian.Gaussian'>"
