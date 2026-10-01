# -*- coding: utf-8 -*-

"""Regression tests for run-time bugs found while adding the shared rules."""

import pytest

import lammps_step
from lammps_step.initialization import explicit_kspace
from seamm_util import Q_

P = {"kspace_accuracy": 1.0e-5, "kspace_smallq": 1.0e-5}


@pytest.fixture(autouse=True)
def real_units():
    """The expected numbers are in LAMMPS 'real' units; other tests change them."""
    from lammps_step.lammps_units import get_lammps_unit_system, set_lammps_unit_system

    previous = get_lammps_unit_system()
    set_lammps_unit_system("real")
    yield
    set_lammps_unit_system(previous)


@pytest.mark.parametrize(
    "method, pair, kspace",
    [
        ("none", "lj/cut", ""),
        ("Ewald summation method", "lj/cut/coul/long", "ewald 1e-05"),
        (
            "PPPM (Particle-particle particle-mesh) method",
            "lj/cut/coul/long",
            "pppm 1e-05",
        ),
        (
            "PPPM method for few charged atoms",
            "lj/cut/coul/long",
            "pppm/cg 1e-05 1e-05",
        ),
        ("MSM (Multilevel summation method)", "lj/cut/coul/msm", "msm 1e-05"),
        ("MSM method for few charged atoms", "lj/cut/coul/msm", "msm/cg 1e-05 1e-05"),
    ],
)
def test_explicit_kspace(method, pair, kspace):
    """An explicitly chosen k-space method (was: KeyError 'kspace_style')."""
    assert explicit_kspace(method, "lj/cut", 10, P) == (pair, kspace)


def test_explicit_kspace_without_charges():
    assert explicit_kspace("Ewald summation method", "lj/cut", 0, P) == ("lj/cut", "")


def test_explicit_msm_needs_an_msm_pair_style():
    with pytest.raises(RuntimeError, match="no MSM version"):
        explicit_kspace("MSM (Multilevel summation method)", "lj/class2", 10, P)


def npt_values(couple, **extra):
    node = lammps_step.NPT()
    values = node.parameters.current_values_to_dict(context={})
    values.update(
        {
            "system type": "solid",
            "use_stress": "general stress",
            "couple": couple,
            "allow shear": False,
            "Panneal": False,
            "Sxx,initial": Q_(1.0, "atm"),
            "Syy,initial": Q_(2.0, "atm"),
            "Szz,initial": Q_(3.0, "atm"),
            "Sxx damp": Q_(1000.0, "fs"),
            "Syy damp": Q_(2000.0, "fs"),
            "Szz damp": Q_(3000.0, "fs"),
        }
    )
    values.update(extra)
    return node, values


def test_npt_damping_times_are_times():
    """The stress damping times were converted as pressures, which failed."""
    node, values = npt_values("none")
    text = node.get_pressure_text(values)
    assert " x -1.000 -1.000 1000.000 " in text
    assert " z -3.000 -3.000 3000.000 " in text


def test_npt_y_and_z_coupling():
    """With y and z coupled, z takes y's stress (was: Szz1 = Szz1)."""
    node, values = npt_values(
        "y and z",
        **{
            "Panneal": True,
            "Syy,final": Q_(5.0, "atm"),
            "Szz,final": Q_(7.0, "atm"),
            "Sxx,final": Q_(1.0, "atm"),
        }
    )
    text = node.get_pressure_text(values)
    assert " y -2.000 -5.000 2000.000 " in text
    assert " z -2.000 -5.000 2000.000 " in text


@pytest.mark.parametrize(
    "mtk, nreset, expected",
    [
        (True, "never", " mtk yes"),
        (False, 100, " mtk no nreset 100"),
        ("yes", "50", " mtk yes nreset 50"),
    ],
)
def test_nose_hoover_barostat_options(mtk, nreset, expected):
    """mtk and nreset were shown in the dialog but never used."""
    assert (
        lammps_step.NPT.nose_hoover_barostat_options({"mtk": mtk, "nreset": nreset})
        == expected
    )


def test_remove_momentum_default_is_a_choice():
    """The default ended in a period, so it was not one of its own choices."""
    P = lammps_step.VelocitiesParameters()
    assert P["remove_momentum"].default in P["remove_momentum"].enumeration
    old = (
        "remove any translational and, for molecular systems, rotational momentum "
        "(default)."
    )
    P = lammps_step.VelocitiesParameters(
        data={"remove_momentum": {"value": old, "units": None}}
    )
    assert P["remove_momentum"].value in P["remove_momentum"].enumeration


@pytest.mark.parametrize(
    "system_type, barostat, expected",
    [
        ("solid", "Nose-Hoover", True),
        ("solid", "Berendsen", False),
        ("fluid", "Nose-Hoover", False),
    ],
)
def test_npt_shear_allowed(system_type, barostat, expected):
    """'allow shear' applies only to a solid with the Nose-Hoover barostat."""
    values = {"allow shear": True, "system type": system_type, "barostat": barostat}
    assert lammps_step.NPT.shear_allowed(values) is expected


def test_berendsen_takes_the_damping_time_but_no_shear():
    """press/berendsen needs the damping time, and cannot control the tilt."""
    node, values = npt_values("none", barostat="Berendsen", **{"allow shear": True})
    values["allow shear"] = node.shear_allowed(values)
    text = node.get_pressure_text(values)
    assert " x -1.000 -1.000 1000.000 " in text
    assert "xy" not in text and "xz" not in text and "yz" not in text


@pytest.mark.parametrize(
    "optimize_cell, system_type, expected",
    [("yes", "solid", True), ("no", "solid", False), ("yes", "fluid", False)],
)
def test_minimization_shear_allowed(optimize_cell, system_type, expected):
    """'allow shear' applies only when optimizing the cell of a solid."""
    values = {
        "allow shear": True,
        "optimize cell": optimize_cell == "yes",
        "system type": system_type,
    }
    assert lammps_step.Minimization.shear_allowed(values) is expected


class _Node:
    def __init__(self, title, uses_berendsen):
        self.title = title
        self.uses_berendsen = uses_berendsen


@pytest.mark.parametrize(
    "forced, triclinic, berendsen, message",
    [
        (True, True, True, "another step in this LAMMPS step allows the cell"),
        (False, True, True, "the cell is not orthorhombic"),
        (False, False, True, None),
        (True, True, False, None),
    ],
)
def test_berendsen_needs_an_orthogonal_box(forced, triclinic, berendsen, message):
    """LAMMPS writes one box for all the steps, and press/berendsen fails with a
    triclinic box; say so before running, not after LAMMPS fails."""
    lammps = lammps_step.LAMMPS()
    lammps.force_triclinic = forced
    lammps._triclinic_box = triclinic
    nodes = [_Node("NPT dynamics", berendsen)]
    if message is None:
        lammps.check_barostats(nodes)
    else:
        with pytest.raises(RuntimeError, match=message):
            lammps.check_barostats(nodes)
