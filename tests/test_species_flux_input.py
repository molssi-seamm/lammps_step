#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of the sampling of the species' momenta for the thermal conductivity of a
mixture."""

from types import SimpleNamespace

import pytest

import lammps_step
import seamm
from seamm_util import Q_


def _configuration(molecules, smiles, periodicity=3):
    return SimpleNamespace(
        periodicity=periodicity,
        find_molecules=lambda as_indices=False: molecules,
        get_molecule_smiles=lambda: smiles,
    )


@pytest.fixture
def nvt(monkeypatch):
    lammps_step.set_lammps_unit_system("real")
    monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables(), raising=False)
    step = lammps_step.NVT()
    step._id = ("3", "4")
    return step


def _use(monkeypatch, step, configuration):
    monkeypatch.setattr(step, "get_system_configuration", lambda: (None, configuration))


P = {"timestep": Q_(1.0, "fs")}


def test_salt_water(nvt, monkeypatch):
    # 4 waters, then a Na+ and a Cl-, then 2 more waters
    molecules = [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11], [12], [13]]
    molecules += [[14, 15, 16], [17, 18, 19]]
    smiles = ["O"] * 4 + ["[Na+]", "[Cl-]"] + ["O"] * 2
    _use(monkeypatch, nvt, _configuration(molecules, smiles))

    lines, ncomputes, nfixes = nvt.species_flux_input(P, 2, 1000, 3, 5)
    text = "\n".join(lines)
    assert (ncomputes, nfixes) == (6, 6)
    # Water first: it has the most atoms
    assert "group               sflux1 id 1:12 15:20" in text
    assert "group               sflux2 id 13" in text
    assert "group               sflux3 id 14" in text
    assert (
        "compute             4 sflux1 reduce sum v_sflux_px v_sflux_py v_sflux_pz"
        in text
    )
    assert "compute             6 sflux3 reduce" in text
    assert "variable            sflux_px atom mass*vx" in text
    assert "fix                 6 all ave/time 2 1 2" in text
    assert "c_4[1] c_4[2] c_4[3] c_5[1]" in text
    assert '"species":["O","[Na+]","[Cl-]"]' in text
    assert "title2 's1x s1y s1z s2x s2y s2z s3x s3y s3z'" in text
    assert "@4+species_flux.trj" in text

    cleanup = nvt.species_flux_cleanup()
    assert "group               sflux3 delete" in cleanup
    assert "variable            sflux_pz delete" in cleanup
    # Only once
    assert nvt.species_flux_cleanup() == []


def test_heat_flux_cleanup_removes_the_groups(nvt, monkeypatch):
    molecules = [[0, 1, 2], [3]]
    _use(monkeypatch, nvt, _configuration(molecules, ["O", "[Na+]"]))
    nvt.species_flux_input(P, 2, 1000, 0, 0)
    lines = nvt.heat_flux_cleanup({"use centroid stress": True}, ff_form="class1")
    assert "group               sflux2 delete" in lines
    assert lines.index("variable            Jz delete") < lines.index(
        "group               sflux1 delete"
    )


def test_pure_liquid_writes_nothing(nvt, monkeypatch):
    molecules = [[0, 1, 2], [3, 4, 5]]
    _use(monkeypatch, nvt, _configuration(molecules, ["O", "O"]))
    assert nvt.species_flux_input(P, 2, 1000, 3, 5) == ([], 3, 5)
    assert nvt.species_flux_cleanup() == []


def test_not_periodic_writes_nothing(nvt, monkeypatch):
    molecules = [[0, 1, 2], [3]]
    _use(monkeypatch, nvt, _configuration(molecules, ["O", "[Na+]"], periodicity=0))
    assert nvt.species_flux_input(P, 2, 1000, 0, 0) == ([], 0, 0)


def test_long_id_lists_are_split(nvt, monkeypatch):
    # Alternating molecules give one range per molecule
    molecules = [[i] for i in range(50)]
    smiles = ["[Na+]", "[Cl-]"] * 25
    _use(monkeypatch, nvt, _configuration(molecules, smiles))
    lines, _, _ = nvt.species_flux_input(P, 2, 1000, 0, 0)
    groups = [line for line in lines if line.startswith("group               sflux1")]
    assert len(groups) == 3
    assert groups[0].split()[3:5] == ["1", "3"]
