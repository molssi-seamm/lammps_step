#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of what the MD sub-steps say they will write, and turning sampling on."""

from types import SimpleNamespace

import pytest

import lammps_step
from seamm_md_analysis import Need
from lammps_step.md_provides import md_provides, request_sampling


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
