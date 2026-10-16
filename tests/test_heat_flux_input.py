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
    heat_flux.parameters["filter"].value = "per-atom averages"
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
    everything = "\n".join(result["script"] + (result["postscript"] or []))
    assert "unfix               constraint" not in everything


def _defined_and_removed(result):
    """The ids of the computes, named fixes and variables the input defines, and
    of those its postscript removes."""
    import re

    script = "\n".join(result["script"])
    post = "\n".join(result["postscript"] or [])
    defined = set()
    for kind in ("compute", "fix"):
        defined |= {
            (kind, m) for m in re.findall(rf"^{kind}\s+([A-Za-z_]\w*)\s", script, re.M)
        }
    defined |= {
        ("variable", m)
        for m in re.findall(r"^variable\s+(\w+)\s+(?:equal|atom)", script, re.M)
    }
    removed = {("compute", m) for m in re.findall(r"^uncompute\s+(\w+)", post, re.M)}
    removed |= {("fix", m) for m in re.findall(r"^unfix\s+(\w+)", post, re.M)}
    removed |= {
        ("variable", m) for m in re.findall(r"^variable\s+(\w+)\s+delete", post, re.M)
    }
    # and those the script itself removes before the main run
    removed |= {("compute", m) for m in re.findall(r"^uncompute\s+(\w+)", script, re.M)}
    removed |= {("fix", m) for m in re.findall(r"^unfix\s+(\w+)", script, re.M)}
    removed |= {
        ("variable", m) for m in re.findall(r"^variable\s+(\w+)\s+delete", script, re.M)
    }
    return defined, removed


def test_no_filter_is_the_default(heat_flux):
    """Without the filter: no preliminary runs, no Python, the centroid stress
    including the constraint forces, and one copy of the convective term."""
    assert heat_flux.parameters["filter"].value == "none"
    result = heat_flux.get_input(extras={"shake": SHAKE})
    script = "\n".join(result["script"])
    assert result["use python"] is False
    # without Python, input_post.dat is not read: everything is in the script
    assert result["postscript"] is None
    assert "J_filter" not in script and "ave/atom" not in script
    runs = [line for line in script.splitlines() if line.startswith("run ")]
    assert runs == ["run                 20000"]
    assert (
        "compute             S_b all centroid/stress/atom NULL bond angle dihedral "
        "improper fix" in script
    )
    assert "compute             S_p all stress/atom NULL pair kspace\n" in script
    assert "-c_flux_b[4])/vol" in script
    shake = script.index(SHAKE.format("constraint"))
    assert shake < script.index("compute             KE") < script.index(runs[0])
    assert script.index(runs[0]) < script.index("uncompute           KE")


@pytest.mark.parametrize("filter_", ["none", "per-atom averages"])
@pytest.mark.parametrize("centroid", ["yes", "no"])
def test_everything_defined_is_removed(heat_flux, filter_, centroid):
    """A later step must be able to define the same computes, fixes and
    variables again."""
    heat_flux.parameters["filter"].value = filter_
    heat_flux.parameters["use centroid stress"].value = centroid
    result = heat_flux.get_input(extras={"shake": SHAKE})
    defined, removed = _defined_and_removed(result)
    assert defined and defined <= removed, defined - removed


def test_centroid_option_only_without_the_filter():
    P = lammps_step.HeatFluxParameters()
    assert P.applies("use centroid stress", {"filter": "none"})
    assert not P.applies("use centroid stress", {"filter": "per-atom averages"})


def test_filter_is_strict(heat_flux):
    heat_flux.parameters["filter"].value = "per-type averages"
    with pytest.raises(ValueError, match="filter"):
        heat_flux.get_input()


@pytest.mark.parametrize(
    "ff_form, centroid, expected",
    [
        ("class1", True, "centroid/stress/atom NULL bond angle dihedral improper fix"),
        ("class1", False, "stress/atom NULL virial"),
        ("class2", True, "stress/atom NULL virial"),
    ],
)
def test_md_heat_flux_input(ff_form, centroid, expected):
    """The heat flux of NVE, NVT and NPT: the constraint forces in the centroid
    stress, one copy of the convective term, and everything removed afterwards."""
    lammps_step.set_lammps_unit_system("real")
    step = lammps_step.NVE()
    P = {"use centroid stress": centroid}
    text = step.heat_flux_input(P, ff_form)
    assert expected in text
    if centroid and ff_form != "class2":
        assert "-c_flux_b[4])/vol" in text
    cleanup = "\n".join(step.heat_flux_cleanup(P, ff_form))
    defined, removed = _defined_and_removed({"script": [text], "postscript": [cleanup]})
    assert defined and defined <= removed, defined - removed
