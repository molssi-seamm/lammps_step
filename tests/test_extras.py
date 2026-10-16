#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests of the extra commands and extra properties of the MD sub-steps."""

from collections import defaultdict
import json
import re
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import lammps_step
import seamm
from lammps_step import extras
from lammps_step.extras import ExtrasError
from seamm_util import Q_

DATA = Path(__file__).parent / "data" / "md_h5md"

COMMANDS = """\
# The mean squared displacement
compute     msd all msd   # all atoms
fix         msd_out all ave/time 10 1 10 &
            c_msd[4] file msd.txt
variable    twice equal 2*c_msd[4]
"""

PROPERTIES = [
    {"name": "msd", "expression": "c_msd[4]", "units": "distance"},
    {
        "name": "Pxy",
        "expression": "pxy",
        "units": "pressure",
        "treatment": "Green-Kubo integral",
    },
]


# ----------------------------------------------------------------------
# The commands
# ----------------------------------------------------------------------
def test_logical_lines():
    lines = extras.logical_lines(COMMANDS)
    assert lines == [
        "compute msd all msd",
        "fix msd_out all ave/time 10 1 10 c_msd[4] file msd.txt",
        "variable twice equal 2*c_msd[4]",
    ]
    # A '#' inside quotes is not a comment
    assert extras.logical_lines('variable s string "a#b" # c') == [
        'variable s string "a#b"'
    ]


def test_what_the_commands_define_is_removed_in_reverse():
    commands = extras.parse_commands(COMMANDS)
    assert [(c["command"], c["id"]) for c in commands] == [
        ("compute", "msd"),
        ("fix", "msd_out"),
        ("variable", "twice"),
    ]
    assert extras.cleanup_lines(commands) == [
        "variable            twice delete",
        "unfix               msd_out",
        "uncompute           msd",
    ]
    # Redefining an ID removes it once
    twice = extras.parse_commands("compute a all ke\ncompute a all pe")
    assert extras.cleanup_lines(twice) == ["uncompute           a"]


@pytest.mark.parametrize(
    "text, message",
    [
        ("group water type 1 2", "'group' is not allowed"),
        ("run 100", "Use a Custom sub-step"),
        ("fix 3 all nve", "may not be a number"),
        ("compute S_p all stress/atom NULL", "SEAMM uses the compute ID 'S_p'"),
        ("variable seamm_x_a equal 1", "SEAMM uses the variable ID"),
        ("variable a delete", "deleted automatically"),
        ("compute", "has no ID"),
    ],
)
def test_commands_refused(text, message):
    with pytest.raises(ExtrasError, match=message):
        extras.parse_commands(text)


def test_every_problem_is_reported():
    with pytest.raises(ExtrasError) as e:
        extras.parse_commands("group g type 1\nfix 2 all nve\ncompute ok all ke")
    assert len(str(e.value).splitlines()) == 2


# ----------------------------------------------------------------------
# The properties
# ----------------------------------------------------------------------
def test_properties_completed():
    msd, pxy = extras.check_properties(PROPERTIES)
    assert msd["treatment"] == "average" and msd["kind"] == "quantity"
    assert msd["display_units"] == "Å" and msd["column"] == "x_msd"
    assert msd["variable"] == "seamm_x_msd" and msd["store"] is False
    assert pxy["prefactor"] == "V/kT" and pxy["display_units"] == "atm"
    # Explicit units, dimensionless, and empty rows of the table
    rows = [
        {"name": "dip", "expression": "c_d[1]", "units": "debye"},
        {"name": "n", "expression": "count(all)"},
        {"name": "", "expression": ""},
    ]
    dip, n = extras.check_properties(rows)
    assert (dip["kind"], dip["display_units"]) == ("explicit", "debye")
    assert (n["units"], n["display_units"]) == ("dimensionless", "")


@pytest.mark.parametrize(
    "row, message",
    [
        ({"name": "1x", "expression": "pe"}, "must start with a letter"),
        ({"name": "a", "expression": ""}, "has no expression"),
        ({"name": "a", "expression": 'v_"x"'}, "double quote"),
        ({"name": "a", "expression": "pe", "units": "furlongs/fortnight^7x"}, "units"),
        ({"name": "a", "expression": "pe", "treatment": "median"}, "treatment"),
        (
            {
                "name": "a",
                "expression": "pxy",
                "treatment": "Green-Kubo integral",
                "prefactor": "V/kT^3",
            },
            "prefactor",
        ),
    ],
)
def test_properties_refused(row, message):
    with pytest.raises(ExtrasError, match=message):
        extras.check_properties([row])
    with pytest.raises(ExtrasError, match="Two extra properties"):
        extras.check_properties([{"name": "a", "expression": "pe"}] * 2)


def test_extra_input():
    result = extras.extra_input(
        {"extra commands": COMMANDS, "extra properties": PROPERTIES}
    )
    lines = "\n".join(result["lines"])
    assert "compute msd all msd" in lines
    assert 'variable            seamm_x_msd equal "c_msd[4]"' in lines
    assert result["properties"] == "v_seamm_x_msd v_seamm_x_Pxy"
    assert result["titles"] == "x_msd x_Pxy"
    assert result["cleanup"] == [
        "variable            seamm_x_Pxy delete",
        "variable            seamm_x_msd delete",
        "variable            twice delete",
        "unfix               msd_out",
        "uncompute           msd",
    ]
    empty = extras.extra_input({})
    assert empty["lines"] == [] and empty["cleanup"] == []
    assert empty["properties"] == ""


def test_check():
    assert extras.check({}) == ""
    text = extras.check(
        {"extra commands": "group g type 1", "extra properties": [{"name": "T"}]}
    )
    assert "'group' is not allowed" in text and "has no expression" in text


# ----------------------------------------------------------------------
# The input of the MD sub-steps
# ----------------------------------------------------------------------
@pytest.fixture
def md_step(monkeypatch):
    """An NVE, NVT or NPT step for periodic water, without a real flowchart."""

    def make(cls, extra=True):
        lammps_step.set_lammps_unit_system("real")
        monkeypatch.setattr(
            seamm, "flowchart_variables", seamm.Variables(), raising=False
        )
        step = cls()
        step._id = ("3", "4")
        step.parameters["time"].value = 20
        step.parameters["timestep"].value = 1
        step.parameters["timestep"].units = "fs"
        if extra:
            step.parameters["extra commands"].value = COMMANDS
            step.parameters["extra properties"].value = PROPERTIES
        ff = SimpleNamespace(current_forcefield="oplsaa")
        step.parent = SimpleNamespace(
            forcefield=lambda: ff,
            ff_form=lambda: "class1",
            have_dreiding_hbonds=False,
            model_chemistry=lambda task: f"LAMMPS:{task}|FF@oplsaa",
            _data={"masses": [16.0, 1.0]},
        )
        configuration = SimpleNamespace(periodicity=3)
        monkeypatch.setattr(
            step, "get_system_configuration", lambda *a, **k: (None, configuration)
        )
        monkeypatch.setattr(
            step, "trajectory_input", lambda P, dt, n, nc, nd, nf: ([], nc, nd, nf)
        )
        step.references = SimpleNamespace(cite=lambda **kwargs: None)
        step._bibliography = defaultdict(str)
        return step

    return make


@pytest.mark.parametrize("ensemble", ["NVE", "NVT", "NPT"])
def test_md_input(md_step, ensemble):
    step = md_step(getattr(lammps_step, ensemble))
    lines = step.get_input()["script"]
    text = "\n".join(lines)
    ensemble_fix = next(i for i, x in enumerate(lines) if x.startswith("fix "))
    commands = lines.index("compute msd all msd")
    sampling = next(i for i, x in enumerate(lines) if "ave/time" in x and "trj" in x)
    run = next(i for i, x in enumerate(lines) if x.startswith("run "))
    cleanup = lines.index("uncompute           msd")
    assert ensemble_fix < commands < sampling < run < cleanup
    # Both the summary and the state sample the extra properties
    assert text.count("v_seamm_x_msd v_seamm_x_Pxy") == 2
    assert "x_msd x_Pxy" in text
    # SEAMM's own fixes are removed before the variables they sample
    last_unfix = max(i for i, x in enumerate(lines) if re.match(r"unfix +\d+$", x))
    assert last_unfix < lines.index("variable            seamm_x_Pxy delete")


@pytest.mark.parametrize("ensemble", ["NVE", "NVT", "NPT"])
def test_md_input_without_extras_is_unchanged(md_step, ensemble):
    step = md_step(getattr(lammps_step, ensemble), extra=False)
    text = "\n".join(step.get_input()["script"])
    assert "seamm_x_" not in text and "Extra" not in text


def test_md_input_refuses_wrong_extras(md_step):
    step = md_step(lammps_step.NVT)
    step.parameters["extra commands"].value = "group g type 1"
    with pytest.raises(RuntimeError, match="extra commands or properties are wrong"):
        step.get_input()
    problems = step.check_flowchart([])
    assert len(problems) == 1 and "'group' is not allowed" in problems[0]


def test_heat_flux_has_no_extras(md_step):
    step = md_step(lammps_step.HeatFlux)
    assert step.extra_definitions() == []
    assert step.check_flowchart([]) == []


# ----------------------------------------------------------------------
# The treatments
# ----------------------------------------------------------------------
def ar1(n, phi, seed=1):
    """An AR(1) series of unit variance: its ACF is phi**k."""
    from scipy.signal import lfilter

    rng = np.random.default_rng(seed)
    noise = rng.normal(scale=np.sqrt(1 - phi**2), size=n)
    noise[0] = rng.normal()
    return lfilter([1.0], [1.0, -phi], noise)


def test_autocorrelation():
    phi = 0.9
    result = extras.autocorrelation(ar1(100_000, phi), Q_(2.0, "fs"))
    acf = np.array(result["acf"])
    assert acf[0] == pytest.approx(1.0)
    assert acf[1] == pytest.approx(phi, abs=0.01)
    assert acf[10] == pytest.approx(phi**10, abs=0.03)
    # The trapezoid integral of phi**k: 1/(1-phi) - 1/2 samples
    tau = result["tau"].m_as("fs")
    assert tau == pytest.approx(2.0 * (1 / (1 - phi) - 0.5), rel=0.1)
    with pytest.raises(ExtrasError, match="constant"):
        extras.autocorrelation(np.ones(100), Q_(1.0, "fs"))


def test_green_kubo_integral():
    """An AR(1) series correlated over about a picosecond, sampled every 10 fs:
    the integral of phi**k is 1/(1-phi) - 1/2 samples."""
    phi = 0.99
    definition = extras.check_properties(
        [
            {
                "name": "a",
                "expression": "v_a",
                "treatment": "Green-Kubo integral",
                "prefactor": "2 / fs",
            }
        ]
    )[0]
    gk = extras.green_kubo(ar1(1_000_000, phi), Q_(10.0, "fs"), definition)
    expected = 2 * 10 * (1 / (1 - phi) - 0.5)
    value, stderr = gk["value"].m_as(""), gk["stderr"].m_as("")
    # The errors from blocks match the spread over seeds (about 4 %)
    assert gk["error_from"] == "blocks"
    assert 0 < stderr < 0.1 * expected
    assert abs(value - expected) < 3 * stderr
    assert len(gk["ts"]) == len(gk["integral"])


def test_result_units():
    k = Q_("k_B")
    stress = Q_(1000.0, "Å^3") / (k * Q_(300.0, "K")) * Q_(1.0, "atm") ** 2
    assert extras.result_units(stress * Q_(1.0, "fs")) == Q_(1.0, "mPa*s").units
    assert extras.result_units(Q_(2.0, "1/fs") * Q_(1.0, "fs")) == Q_(1.0, "").units
    assert extras.result_units(Q_(1.0, "kg")) == Q_(1.0, "kg").units


def test_prefactors():
    rows = [
        {"name": n, "expression": "x", "treatment": "Green-Kubo integral"}
        for n in ("a", "b", "c")
    ]
    rows[1]["prefactor"] = "V/kT^2"
    rows[2]["prefactor"] = "1/(V kT)"
    a, b, c = extras.check_properties(rows)
    V, T = Q_(1000.0, "Å^3"), Q_(300.0, "K")
    k = Q_("k_B")
    assert extras.prefactor_value(a, V, T).m_as("m*s^2/kg") == pytest.approx(
        (V / (k * T)).m_as("m*s^2/kg")
    )
    assert extras.prefactor_value(b, V, T).check("[length]^3/[energy]/[temperature]")
    assert extras.prefactor_value(c, V, T).check("1/[length]^3/[energy]")
    with pytest.raises(ExtrasError, match="non-periodic"):
        extras.prefactor_value(a, None, T)


# ----------------------------------------------------------------------
# The analysis in the step
# ----------------------------------------------------------------------
def test_analysis_of_the_extra_properties(tmp_path, monkeypatch):
    monkeypatch.setattr(seamm, "flowchart_variables", seamm.Variables(), raising=False)
    step = lammps_step.NVT()
    step.parameters["extra properties"].value = [
        {"name": "msd", "expression": "c_msd[4]", "units": "distance"},
        {
            "name": "acf",
            "expression": "v_a",
            "treatment": "autocorrelation",
        },
        {
            "name": "gk",
            "expression": "v_a",
            "treatment": "Green-Kubo integral",
            "prefactor": "2 / fs",
        },
        {"name": "lost", "expression": "v_b"},
    ]
    monkeypatch.setattr(
        lammps_step.NVT, "directory", property(lambda self: str(tmp_path))
    )
    series = ar1(1_000_000, 0.99)
    properties = {
        "x_msd": {"mean": 1.5, "stderr": 0.1, "timestep": 10.0},
        "x_acf": {"mean": 0.0, "stderr": 0.01, "timestep": 10.0},
        "x_gk": {"mean": 0.0, "stderr": 0.01, "timestep": 10.0},
    }
    step._extra_series = {"x_acf": series, "x_gk": series}
    data = {"T": 300.0}
    text = step._analyze_extra_properties(properties, data, None)
    results = data["extra properties"]
    assert "Extra Properties" in text and "lost" in text
    assert results["msd"]["mean"] == 1.5 and results["msd"]["units"] == "Å"
    assert results["acf"]["correlation time"] == pytest.approx(0.995, rel=0.2)
    gk = results["gk"]
    assert abs(gk["integral"] - 1990.0) < 3 * gk["integral,stderr"]
    assert results["lost"]["error"] == "the property did not equilibrate"
    assert (tmp_path / "acf_acf.graph").exists()
    assert (tmp_path / "gk_gk.graph").exists()
    json.dumps(results)  # stored as JSON


def test_stored_properties():
    store = {}
    names = {}

    class Properties:
        def exists(self, name):
            return name in names

        def add(self, name, _type, units=None, description="", noerror=False):
            names[name] = (units, description)

        def put(self, name, value):
            store[name] = value

    step = lammps_step.NVT()
    configuration = SimpleNamespace(properties=Properties())
    definition = extras.check_properties(
        [{"name": "msd", "expression": "c_msd[4]", "units": "distance"}]
    )[0]
    step._store_extra_property.__func__(
        SimpleNamespace(model="oplsaa"), configuration, definition, 1.5, 0.1, "Å"
    )
    assert store == {"msd#LAMMPS#oplsaa": 1.5, "msd,stderr#LAMMPS#oplsaa": 0.1}
    assert names["msd#LAMMPS#oplsaa"][0] == "Å"


# ----------------------------------------------------------------------
# The H5MD file
# ----------------------------------------------------------------------
def test_extra_properties_in_the_h5md_file(tmp_path):
    from seamm_md_analysis import open_md, read_trajectory
    from .test_md_h5md import NVT, P, configuration

    from lammps_step.md_h5md import write_md_h5md

    lammps_step.set_lammps_unit_system("real")
    for path in (DATA / "water").iterdir():
        shutil.copy(path, tmp_path)
    # Add an extra property to the state file: twice the temperature
    state = tmp_path / "nvt_state.trj"
    lines = state.read_text().splitlines()
    lines[1] += " x_T2"
    for i in range(2, len(lines)):
        if lines[i].strip() != "":
            lines[i] += f" {2 * float(lines[i].split()[2])}"
    state.write_text("\n".join(lines) + "\n")

    node = NVT(tmp_path, [-0.82, 0.41, 0.41] * 64)
    node.extra_definitions = lambda: extras.check_properties(
        [{"name": "T2", "expression": "2*temp", "units": "temperature"}]
    )
    config = configuration(
        64, ["O", "H", "H"], [8, 1, 1], [15.999, 1.008, 1.008], cell=[12.4172] * 3
    )
    config.cell.parameters = [12.417241581015336] * 3 + [90.0] * 3
    write_md_h5md(node, P, config)
    (record,) = node.variables["_md_runs"]
    assert "user/T2" in record["observables"]

    _, columns = read_trajectory(state)
    with open_md(tmp_path / "md.h5md") as md:
        assert "user/T2" in md.observables()
        T2 = md.observable("user/T2")
        assert np.allclose(T2.values, 2 * columns["T"])
        assert T2.units == Q_(1.0, "K").units
        assert T2.attrs["seamm_expression"] == "2*temp"
        assert T2.attrs["seamm_treatment"] == "average"


def test_md_provides_lists_the_extra_properties():
    from lammps_step.md_provides import md_provides

    node = SimpleNamespace(
        parameters=lammps_step.NVT_Parameters(),
        title="NVT dynamics",
        parent=SimpleNamespace(title="LAMMPS"),
        extra_definitions=lambda: extras.check_properties(
            [{"name": "msd", "expression": "c_msd[4]", "units": "distance"}]
        ),
    )
    provides = md_provides(node)
    assert provides.data["user/msd"] == provides.data["temperature"]


@pytest.mark.parametrize("ensemble", ["NVE", "NVT", "NPT"])
def test_summary_averages_each_tenth_of_the_run(md_step, ensemble):
    """Each of the 10 summary lines is the average over its tenth of the run (20000
    steps: every 10 steps, 200 samples, written every 2000), with the time, the
    first value, left as the time of the line."""
    step = md_step(getattr(lammps_step, ensemble), extra=False)
    text = "\n".join(step.get_input()["script"])
    summary = next(block for block in text.split("fix ") if "_summary.trj" in block)
    assert "ave/time 10 200 2000" in summary
    assert "off 1" in summary
    state = next(block for block in text.split("fix ") if "_state.trj" in block)
    assert "off 1" in state
