#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of what the MD sub-steps say they will write, and turning sampling on."""

from types import SimpleNamespace

import pytest

import lammps_step
from seamm_md_analysis import Need
from lammps_step.md_provides import md_provides, request_sampling, sampling_edits


class NVT:
    def __init__(self, parameters_class=lammps_step.NVT_Parameters):
        self.parameters = parameters_class()
        self.title = "NVT dynamics"
        self.parent = SimpleNamespace(title="LAMMPS")


class NPT(NVT):
    def __init__(self):
        super().__init__(lammps_step.NPT_Parameters)


def set_(node, **kwargs):
    for key, value in kwargs.items():
        key = key.replace("_", " ")
        if isinstance(value, tuple):
            node.parameters[key].value = value[0]
            node.parameters[key].units = value[1]
        else:
            node.parameters[key].value = value


def test_state_only_by_default():
    provides = md_provides(NVT())
    assert provides.step == "LAMMPS / NVT dynamics"
    assert provides.data["temperature"].time.m_as("fs") == 50
    assert "pressure_tensor" in provides.data
    assert "volume" not in provides.data
    assert "molecules/position" not in provides.data
    assert provides.timestep is None  # "normal" depends on the masses


def test_npt_has_volume_and_density():
    assert "density" in md_provides(NPT()).data


def test_sampling_by_interval_and_by_number():
    node = NVT()
    set_(node, time=("10", "ps"), timestep=("1", "fs"))
    node.parameters["com positions"].value = "by time interval"
    node.parameters["com positions rate"].value = "200"
    node.parameters["com velocities"].value = "by number of samples"
    node.parameters["com velocities number of samples"].value = "1000"
    provides = md_provides(node)
    assert provides.timestep.m_as("fs") == 1
    assert provides.data["molecules/position"].time.m_as("fs") == 200
    assert provides.data["molecules/velocity"].time.m_as("fs") == pytest.approx(10)


def test_trajectory_goes_to_its_own_group_with_atomic_positions():
    node = NVT()
    node.parameters["trajectory"].value = "by time interval"
    node.parameters["trajectory velocities"].value = "yes"
    assert "all/position" in md_provides(node).data
    assert "all/velocity" in md_provides(node).data
    node.parameters["atomic positions"].value = "by time interval"
    data = md_provides(node).data
    assert "trajectory/position" in data
    assert "trajectory/velocity" in data


def test_variables_are_unknown():
    node = NVT()
    node.parameters["com positions"].value = "$sampling"
    assert md_provides(node).data["molecules/position"].time is None


def test_request_sampling_by_time():
    node = NVT()
    change = request_sampling(
        node, Need(["molecules/position"], max_time="1 ps", sample_time="100 fs")
    )
    assert "molecules/position every 100 fs" in change
    assert md_provides(node).data["molecules/position"].time.m_as("fs") == 100


def test_sampling_edits_leave_the_step_unchanged():
    """The changes as data, for a fix offered before the job runs."""
    node = NVT()
    need = Need(["molecules/position"], max_time="1 ps", sample_time="100 fs")
    edits, description = sampling_edits(node, need)
    assert edits == [
        ("com positions", "by time interval", None),
        ("com positions rate", "100", "fs"),
    ]
    assert description == "sampling of molecules/position every 100 fs in NVT dynamics"
    assert "molecules/position" not in md_provides(node).data
    for parameter, value, units in edits:
        node.parameters[parameter].value = value
        if units is not None:
            node.parameters[parameter].units = units
    assert md_provides(node).data["molecules/position"].time.m_as("fs") == 100


def test_request_sampling_by_steps_needs_the_timestep():
    node = NVT()
    assert request_sampling(node, Need("heat_flux", max_steps=2)) is None
    set_(node, timestep=("0.5", "fs"))
    assert request_sampling(node, Need("heat_flux", max_steps=2)) is not None
    assert md_provides(node).data["heat_flux"].time.m_as("fs") == 1


def test_heat_flux_substep_samples_every_two_steps():
    from lammps_step.heat_flux import HeatFlux

    node = HeatFlux.__new__(HeatFlux)
    node.parameters = lammps_step.HeatFluxParameters()
    node.title = "Heat Flux"
    node.parent = SimpleNamespace(title="LAMMPS")
    assert node.md_provides().data["heat_flux"].steps == 2


def test_pressure_tensor_sampled_on_its_own():
    node = NVT()
    set_(node, time=("10", "ps"), timestep=("1", "fs"))
    # From the state, every 50 fs by default
    assert md_provides(node).data["pressure_tensor"].time.m_as("fs") == 50
    need = Need("pressure_tensor", max_time="5 fs", why="the viscosity")
    change = request_sampling(node, need)
    assert "pressure_tensor" in change
    assert node.parameters["pressure tensor"].value == "by time interval"
    interval = md_provides(node).interval("pressure_tensor")
    assert interval.time.m_as("fs") == 5
    assert interval.steps == 5


def test_old_shear_stress_parameters_are_renamed():
    data = {
        "shear stress": {"value": "by time interval", "units": None},
        "shear stress rate": {"value": "2", "units": "fs"},
        "shear stress number of samples": {"value": "500", "units": None},
    }
    for cls in (
        lammps_step.NVE_Parameters,
        lammps_step.NVT_Parameters,
        lammps_step.NPT_Parameters,
    ):
        P = cls(data={k: dict(v) for k, v in data.items()})
        assert "shear stress" not in P
        assert P["pressure tensor"].value == "by time interval"
        assert P["pressure tensor rate"].value == "2"
        assert P["pressure tensor rate"].units == "fs"
        assert P["pressure tensor number of samples"].value == "500"
