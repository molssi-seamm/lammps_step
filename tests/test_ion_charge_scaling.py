#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Scaling the charges of ions, for the scaled-charge models of electrolytes."""

from types import SimpleNamespace

import numpy as np
import pytest

import lammps_step


def _configuration():
    """A water, Li+, BF4- and another water, in that order."""
    molecules = [[0, 1, 2], [3], [4, 5, 6, 7, 8], [9, 10, 11]]
    smiles = ["O", "[Li+]", "F[B-](F)(F)F", "O"]
    return SimpleNamespace(
        find_molecules=lambda as_indices=False: molecules,
        get_molecule_smiles=lambda: smiles,
    )


WATER = [-0.82, 0.41, 0.41]
BF4 = [0.9, -0.475, -0.475, -0.475, -0.475]


@pytest.fixture
def initialization():
    step = lammps_step.Initialization()
    step.parent = SimpleNamespace()
    step._id = ("1", "1")
    return step


def test_ions_are_scaled(initialization):
    eex = {"charges": WATER + [1.0] + BF4 + WATER}
    lines = initialization.scale_ion_charges(_configuration(), eex, 0.8)
    q = np.array(eex["charges"])
    np.testing.assert_allclose(q[:3], WATER)
    np.testing.assert_allclose(q[9:], WATER)
    assert q[3] == pytest.approx(0.8)
    np.testing.assert_allclose(q[4:9], 0.8 * np.array(BF4))
    assert q.sum() == pytest.approx(0.0, abs=1e-12)
    assert initialization.parent._ion_charge_scale == 0.8
    assert "[Li+] (+1 -> +0.8)" in lines[0]
    assert "F[B-](F)(F)F (-1 -> -0.8)" in lines[0]


def test_no_ions_leaves_the_charges(initialization):
    configuration = SimpleNamespace(
        find_molecules=lambda as_indices=False: [[0, 1, 2]],
        get_molecule_smiles=lambda: ["O"],
    )
    eex = {"charges": list(WATER)}
    lines = initialization.scale_ion_charges(configuration, eex, 0.8)
    assert eex["charges"] == WATER
    assert lines[0].startswith("Warning")
    assert not hasattr(initialization.parent, "_ion_charge_scale")


def test_the_factor_must_be_positive(initialization):
    with pytest.raises(RuntimeError, match="positive"):
        initialization.scale_ion_charges(_configuration(), {"charges": [0.0]}, 0.0)


def test_only_with_the_forcefields_charges():
    P = lammps_step.InitializationParameters()
    assert P.applies("scale ion charges", {"atomic charges": "default for forcefield"})
    assert not P.applies("scale ion charges", {"atomic charges": "no charges"})


def test_description(initialization):
    initialization.parameters["scale ion charges"].value = 0.8
    assert "scaled by 0.8" in initialization.description_text()
    initialization.parameters["scale ion charges"].value = 1.0
    assert "scaled" not in initialization.description_text()
