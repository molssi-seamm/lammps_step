"""LAMMPS writes its own variables as ${name}: the custom script and a sub-step's
extra commands are passed to LAMMPS as written, not expanded as flowchart
variables (seamm 2026.11.6 expands '${...}' in text parameters)."""

import lammps_step


def test_custom_script_is_not_expanded():
    parameters = lammps_step.CustomParameters()
    parameters["script"].value = "fix 1 all nvt temp ${T} ${T} 100"
    values = parameters.current_values_to_dict(context={})
    assert values["script"] == "fix 1 all nvt temp ${T} ${T} 100"


def test_extra_commands_are_not_expanded():
    for parameters in (
        lammps_step.NVE_Parameters(),
        lammps_step.NVT_Parameters(),
        lammps_step.NPT_Parameters(),
    ):
        parameters["extra commands"].value = "variable x equal ${T}*2"
        values = parameters.current_values_to_dict(context={})
        assert values["extra commands"] == "variable x equal ${T}*2"
