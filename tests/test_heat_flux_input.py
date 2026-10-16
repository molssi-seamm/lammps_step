#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of the LAMMPS input for the Heat Flux sub-step."""

from types import SimpleNamespace

import pytest

import lammps_step
import seamm

SHAKE = "fix                 {} all shake 1.0e-6 50 1000 b  1 a  1"


@pytest.fixture
def heat_flux(monkeypatch):
    """A Heat Flux step for periodic water, without a real flowchart."""
    lammps_step.set_lammps_unit_system("real")
    monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables(), raising=False)
    step = lammps_step.HeatFlux()
    step._id = ("3", "4")
    step.parameters["time"].value = 20
    step.parameters["timestep"].value = 1
    step.parameters["timestep"].units = "fs"
    ff = SimpleNamespace(current_forcefield="nacl_water")
    step.parent = SimpleNamespace(
        forcefield=lambda: ff, ff_form=lambda: "class1", _data={"masses": [16.0, 1.0]}
    )
    configuration = SimpleNamespace(periodicity=3)
    monkeypatch.setattr(step, "get_system_configuration", lambda: (None, configuration))
    monkeypatch.setattr(step, "trajectory_input", lambda P, *args: ([], 0, 0, 0))
    return step


def test_constraints_from_the_start(heat_flux):
    """SHAKE must be on before the runs that settle the system and average it."""
    result = heat_flux.get_input(extras={"shake": SHAKE})
    script = "\n".join(result["script"]).splitlines()
    shake = [i for i, line in enumerate(script) if " shake " in line]
    assert len(shake) == 1
    first_run = next(i for i, line in enumerate(script) if line.startswith("run "))
    assert shake[0] < first_run
    assert script[shake[0]] == SHAKE.format("constraint")
    assert "unfix               constraint" in "\n".join(result["postscript"])


def test_no_constraints(heat_flux):
    result = heat_flux.get_input()
    script = "\n".join(result["script"])
    assert " shake " not in script
    assert "unfix               constraint" not in "\n".join(result["postscript"])
