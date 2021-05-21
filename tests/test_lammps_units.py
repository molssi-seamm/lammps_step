#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `lammps_step` package, lammps_units module."""

import lammps_step

# Need pint for formats.
from seamm_util import ureg, Q_, units_class  # noqa: F401


def test_default_system():
    """Testing that we can convert to/from the default unit system."""
    value = "2.0 kJ/mol"
    converted = lammps_step.to_lammps_units(value)
    assert (
        f"{value} in real units is {converted:.4}"
        == "2.0 kJ/mol in real units is 0.478"
    )
    assert (
        f'{converted:.4} --> {lammps_step.from_lammps_units(converted, "kJ/mol"):~.4}'  # noqa: E501
        == "0.478 --> 2.0 kJ / mol"
    )


def test_metal_system():
    """Testing that we can convert to/from the metal unit system."""
    lammps_step.set_lammps_unit_system("metal")
    value = "2.0 kJ/mol"
    converted = lammps_step.to_lammps_units(value)
    assert (
        f"{value} in metal units is {converted:.4}"
        == "2.0 kJ/mol in metal units is 0.02073"
    )
    assert (
        f'{converted:.4} --> {lammps_step.from_lammps_units(converted, "kJ/mol"):~.4}'  # noqa: E501
        == "0.02073 --> 2.0 kJ / mol"
    )
