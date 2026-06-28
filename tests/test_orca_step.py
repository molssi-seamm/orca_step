#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `orca_step` package."""

import pytest  # noqa: F401
import orca_step  # noqa: F401


def test_construction():
    """Just create an object and test its type."""
    result = orca_step.Orca()
    assert str(type(result)) == "<class 'orca_step.orca.Orca'>"
