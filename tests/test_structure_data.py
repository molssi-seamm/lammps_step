#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the LAMMPS structure file."""

import types

import pytest

import lammps_step

# The charges of the 500 ethylene carbonate molecules of a real job, assigned from
# OPLS-AA via LigParGen. They sum to exactly zero, and carry five decimals: the four
# from the forcefield plus the small adjustment that balances the total.
EC_CHARGES = (
    [-0.36759] * 500
    + [-0.32129] * 500
    + [-0.32079] * 500
    + [-0.04419] * 500
    + [-0.04379] * 500
    + [0.13531] * 1000
    + [0.13541] * 1000
    + [0.55621] * 500
)


def _structure_data(charges):
    """Write a structure file for a system with the given charges."""
    n_atoms = len(charges)
    eex = {
        "n_atoms": n_atoms,
        "n_atom_types": 5,
        "charges": charges,
        "atoms": [(0.0, 0.0, 0.0, 1)] * n_atoms,
        "molecule": [i // 10 for i in range(n_atoms)],
        "periodicity": 0,
        "masses": [(12.011, "C")] * 5,
    }
    configuration = types.SimpleNamespace(
        atoms=types.SimpleNamespace(have_velocities=False)
    )
    node = types.SimpleNamespace(
        eex=eex,
        _data={},
        force_triclinic=False,
        get_system_configuration=lambda: (None, configuration),
    )
    result = lammps_step.LAMMPS.structure_data(node)
    return result[0] if isinstance(result, tuple) else result


def _charges_in(text):
    """The charges as LAMMPS would read them back."""
    lines = text.splitlines()
    start = lines.index("Atoms") + 2
    charges = []
    for line in lines[start:]:
        if line.strip() == "":
            break
        charges.append(float(line.split()[3]))
    return charges


def test_charges_keep_the_system_neutral():
    """The charges must not be rounded until the system carries a net charge.

    They used to be written with three decimals. The forcefield assigns four or
    more, so over a large system the roundings accumulated: these 5000 atoms sum
    to zero but were written as a system with a charge of -1.
    """
    # Summing 5000 floats leaves a little noise, which varies with the platform.
    # The tolerance is far below what this is about: the charges used to add up to
    # a whole electron.
    assert sum(EC_CHARGES) == pytest.approx(0.0, abs=1e-6)

    charges = _charges_in(_structure_data(EC_CHARGES))

    assert len(charges) == len(EC_CHARGES)
    assert sum(charges) == pytest.approx(0.0, abs=1e-6)


def test_charges_are_written_to_six_decimals():
    """Five decimals are needed for this case, six leaves some room."""
    charges = _charges_in(_structure_data(EC_CHARGES))

    assert charges[0] == -0.36759
    assert set(charges) == set(EC_CHARGES)


def test_charges_without_molecules():
    """The same, for the atom style that has no molecule column."""
    n_atoms = len(EC_CHARGES)
    eex = {
        "n_atoms": n_atoms,
        "n_atom_types": 5,
        "charges": EC_CHARGES,
        "atoms": [(0.0, 0.0, 0.0, 1)] * n_atoms,
        "periodicity": 0,
        "masses": [(12.011, "C")] * 5,
    }
    configuration = types.SimpleNamespace(
        atoms=types.SimpleNamespace(have_velocities=False)
    )
    node = types.SimpleNamespace(
        eex=eex,
        _data={},
        force_triclinic=False,
        get_system_configuration=lambda: (None, configuration),
    )
    result = lammps_step.LAMMPS.structure_data(node)
    text = result[0] if isinstance(result, tuple) else result

    lines = text.splitlines()
    start = lines.index("Atoms") + 2
    charges = [float(line.split()[2]) for line in lines[start : start + n_atoms]]

    assert sum(charges) == pytest.approx(0.0, abs=1e-6)
