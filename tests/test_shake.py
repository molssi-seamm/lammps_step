#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of the constraints (fix shake) for rigid bonds and angles and X-H bonds."""

from lammps_step.lammps import LAMMPS


def eex(rigid=True):
    """Two waters and a methanol, CH3-OH.

    Atom types: 1 O(water), 2 H(water), 3 C, 4 H(on C), 5 O(alcohol), 6 H(on O).
    Bond types: 1 O-H water, 2 C-H, 3 C-O, 4 O-H alcohol.
    Angle types: 1 H-O-H water, 2 H-C-H, 3 C-O-H.

    With rigid=True the forcefield makes the water rigid, as e.g. oplsaa+-spc_e.
    """
    types = [1, 2, 2, 1, 2, 2, 3, 4, 4, 4, 5, 6]
    water_bond = "rigid_bond" if rigid else "quadratic_bond"
    water_angle = "rigid_angle" if rigid else "quadratic_angle"
    return {
        "elements": ["O", "H", "C", "H", "O", "H"],
        "atoms": [(0.0, 0.0, 0.0, t) for t in types],
        "n_bonds": 9,
        "bonds": [
            (1, 2, 1),
            (1, 3, 1),
            (4, 5, 1),
            (4, 6, 1),
            (7, 8, 2),
            (7, 9, 2),
            (7, 10, 2),
            (7, 11, 3),
            (11, 12, 4),
        ],
        "bond parameters": [
            (water_bond, {}, ("o", "h"), "explicit", ("o", "h")),
            ("quadratic_bond", {}, ("c", "hc"), "explicit", ("c", "hc")),
            ("quadratic_bond", {}, ("c", "oh"), "explicit", ("c", "oh")),
            ("quadratic_bond", {}, ("oh", "ho"), "explicit", ("oh", "ho")),
        ],
        "n_angles": 4,
        "angles": [(2, 1, 3, 1), (5, 4, 6, 1), (8, 7, 9, 2), (7, 11, 12, 3)],
        "angle parameters": [
            (water_angle, {}, ("h", "o", "h"), "explicit", ("h", "o", "h")),
            ("quadratic_angle", {}, ("hc", "c", "hc"), "explicit", ("hc", "c", "hc")),
            ("quadratic_angle", {}, ("c", "oh", "ho"), "explicit", ("c", "oh", "ho")),
        ],
    }


def fix(rigid=True, xh="none"):
    step = LAMMPS.__new__(LAMMPS)
    P = {"fix_XH_bond_lengths": xh}
    return step.shake_fix(P, eex(rigid))


def test_rigid_from_the_forcefield():
    line = fix()
    assert line.split()[2:] == [
        "all",
        "shake",
        "1.0e-6",
        "50",
        "1000",
        "b",
        "1",
        "a",
        "1",
    ]


def test_flexible_forcefield():
    """A flexible water model stays flexible."""
    assert fix(rigid=False) == ""


def test_ch_bonds():
    line = fix(rigid=False, xh="CH")
    assert line.split()[7:] == ["b", "2"]


def test_all_xh_bonds_and_rigid_waters():
    line = fix(xh="all")
    words = line.split()
    b = words.index("b")
    a = words.index("a")
    assert sorted(words[b + 1 : a]) == ["1", "2", "4"]
    assert words[a + 1 :] == ["1"]
