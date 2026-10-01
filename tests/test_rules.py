# -*- coding: utf-8 -*-

"""The rules the dialogs and the flowchart builder share: which settings apply,
given the others (see seamm.Parameters and lammps_step's *_parameters.py)."""

import pytest
import seamm
from seamm.builder import FlowchartBuilder, FlowchartBuildError, set_parameters

import lammps_step


def node_for(substep):
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm.lammps", directory=".")
    node = flowchart.create_node(substep)
    flowchart.add_node(node)
    return node


@pytest.mark.parametrize(
    "substep, settings, refused, reason",
    [
        (
            "Initialization",
            {"kspace_method": "none", "kspace_accuracy": 1.0e-6},
            "kspace_accuracy",
            "it applies when 'kspace_method' is not 'none'",
        ),
        (
            "Initialization",
            {"qeq iterations": 200},
            "qeq iterations",
            "it applies when 'atomic charges' is 'charge equilibration",
        ),
        (
            "Initialization",
            {
                "kspace_method": "PPPM method including dispersion terms",
                "tail_correction": "no",
            },
            "tail_correction",
            "is not 'PPPM method including dispersion terms'",
        ),
        (
            "Velocities",
            {"method": "scaling current velocities", "seed": 1},
            "seed",
            "it applies when 'method' is 'using a random distribution'",
        ),
        (
            "Minimization",
            {"etol": 0.01},
            "etol",
            "it applies when 'convergence' is 'custom'",
        ),
        (
            "Minimization",
            {"timestep": 5.0},
            "timestep",
            "it applies when 'minimizer' is 'Fire' or 'QuickMin'",
        ),
        (
            "Minimization",
            {"optimize cell": "no", "P": 10.0},
            "P",
            "it applies when 'optimize cell' is 'yes'",
        ),
        (
            "Minimization",
            {"Sxx": 10.0},
            "Sxx",
            "it needs 'use_stress', which does not apply (it applies when "
            "'system type' is 'solid')",
        ),
        (
            "NVE",
            {"trajectory extxyz": "yes"},
            "trajectory extxyz",
            "it needs 'trajectory export', which does not apply (it applies when "
            "'trajectory' is not 'never')",
        ),
        (
            "NVE",
            {"trajectory": "by number of samples", "trajectory rate": 10.0},
            "trajectory rate",
            "it applies when 'trajectory' is 'by time interval'",
        ),
        (
            "NVE",
            {"maximum_time": 2.0},
            "maximum_time",
            "it is not implemented yet",
        ),
        (
            "NVT",
            {"thermostat": "Berendsen", "Tchain": 5},
            "Tchain",
            "it applies when 'thermostat' is 'Nose-Hoover'",
        ),
        (
            "NVT",
            {"thermostat": "velocity rescaling", "Tdamp": 50.0},
            "Tdamp",
            "it applies when 'thermostat' is not 'velocity rescaling'",
        ),
        (
            "NPT",
            {"barostat": "Berendsen", "mtk": "no"},
            "mtk",
            "it applies when 'barostat' is 'Nose-Hoover'",
        ),
        (
            "NPT",
            {"allow shear": "yes"},
            "allow shear",
            "it applies when 'system type' is 'solid'",
        ),
        (
            "NPT",
            {"system type": "solid", "use_stress": "general stress", "Pinitial": 5.0},
            "Pinitial",
            "a general stress, not a pressure, is applied",
        ),
        (
            "NPT",
            {
                "system type": "solid",
                "use_stress": "general stress",
                "Syy,initial": 5.0,
            },
            "Syy,initial",
            "'couple' is 'x and z' or 'y and z' or 'none'",
        ),
        (
            "NPT",
            {"Pfinal": 10.0},
            "Pfinal",
            "it applies when 'Panneal' is 'yes'",
        ),
    ],
)
def test_refused_with_reason(substep, settings, refused, reason):
    """A setting that would have no effect is refused, saying why, and the step is
    left as it was."""
    node = node_for(substep)
    before = node.parameters.to_dict()
    with pytest.raises(FlowchartBuildError) as e:
        set_parameters(node, settings)
    message = str(e.value)
    assert f"'{refused}' has no effect" in message
    assert reason in message
    assert node.parameters.to_dict() == before


@pytest.mark.parametrize(
    "substep, settings",
    [
        (
            "Initialization",
            {
                "atomic charges": "charge equilibration (QEq) with point charges",
                "qeq iterations": 200,
                "kspace_method": "MSM method for few charged atoms",
                "kspace_smallq": 1.0e-4,
            },
        ),
        ("Velocities", {"method": "using a random distribution", "seed": 42}),
        (
            "Minimization",
            {
                "convergence": "custom",
                "etol": 0.01,
                "minimizer": "Fire",
                "timestep": 5.0,
                "system type": "solid",
                "use_stress": "general stress",
                "couple": "none",
                "allow shear": "yes",
                "Syy": 10.0,
                "Szz": 10.0,
                "Sxy": 0.0,
            },
        ),
        (
            "NVE",
            {
                "trajectory": "by time interval",
                "trajectory rate": 10.0,
                "trajectory export": "yes",
                "trajectory extxyz": "yes",
                "trajectory extxyz filename": "md.extxyz",
                "heat flux": "by number of samples",
                "heat flux number of samples": 500,
                "use centroid stress": "no",
            },
        ),
        ("NVT", {"thermostat": "Langevin", "seed": 42, "Tdamp": 50.0}),
        (
            "NPT",
            {
                "system type": "solid",
                "barostat": "Nose-Hoover",
                "Panneal": "yes",
                "use_stress": "general stress",
                "couple": "x and y",
                "Szz,initial": 10.0,
                "Szz,final": 20.0,
                "Szz damp": 500.0,
            },
        ),
        ("NPT", {"Panneal": "yes", "Pfinal": 100.0, "barostat": "Berendsen"}),
        ("HeatFlux", {"use centroid stress": "no"}),
    ],
)
def test_accepted_when_they_apply(substep, settings):
    node = node_for(substep)
    set_parameters(node, settings)
    P = node.parameters
    for key, value in settings.items():
        assert P[key].value == value


def test_a_variable_counts_as_met():
    """The value of a variable is known only when the flowchart runs, so the
    settings that depend on it are accepted."""
    node = node_for("NVT")
    set_parameters(node, {"thermostat": "$thermostat", "Tchain": 5, "seed": 3})
    node = node_for("NPT")
    set_parameters(
        node,
        {
            "system type": "$type",
            "use_stress": "$stress",
            "Pinitial": 5.0,
            "Sxx,initial": 2.0,
        },
    )


def test_no_narrowed_choices_or_implied_values():
    """LAMMPS's dialogs neither narrow lists nor set values from other choices, so
    the rules only say what applies."""
    for cls in (
        lammps_step.InitializationParameters,
        lammps_step.MinimizationParameters,
        lammps_step.VelocitiesParameters,
        lammps_step.NPT_Parameters,
        lammps_step.HeatFluxParameters,
    ):
        P = cls()
        assert P.implied() == {}
        assert all(P.choices(key) is None for key in P)
        assert P.problems() == []


def test_keep_orthorhombic_still_translated():
    """Old flowcharts' 'keep orthorhombic' becomes 'allow shear'."""
    data = {"keep orthorhombic": {"value": "no", "units": None}}
    P = lammps_step.NPT_Parameters(data=data)
    assert P["allow shear"].value == "yes"
    values = {**P.current_values(), "system type": "solid"}
    assert P.applies("allow shear", values)


def test_builder():
    """Building a LAMMPS flowchart: the rules are checked as each step is added."""
    fb = FlowchartBuilder("NPT of a solid")
    lammps = fb.add("LAMMPS")
    lammps.add("Initialization")
    lammps.add(
        "NPT",
        {"Syy,initial": 10.0},
        system_type="solid",
        use_stress="general stress",
        couple="none",
    )
    with pytest.raises(FlowchartBuildError, match="'mtk' has no effect"):
        lammps.add("NPT", barostat="Berendsen", mtk="no")
    with pytest.raises(FlowchartBuildError, match="'seed' has no effect"):
        lammps.add("NVT", seed=1)
