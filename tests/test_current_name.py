"""The structure a step writes keeps its current name, e.g. a title read from a
file, in front of the step's own text (SEAMM 5380, 2026-10-23)."""


def test_the_default_name_keeps_the_current_name():
    from lammps_step import NVT_Parameters

    P = NVT_Parameters()
    choices = list(P["configuration name"].enumeration)
    assert P["configuration name"].default == "{current name} simulated with {model}"
    assert choices[0] == "{current name} simulated with {model}"
    assert "simulated with {model}" in choices
    assert "keep current name" in choices
