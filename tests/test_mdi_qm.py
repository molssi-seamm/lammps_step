#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for the QM-MD-over-MDI (Phase C) path in ``lammps_step``.

These exercise the driver-side pieces that read the ``_model_chemistry``
variable published by the Model Chemistry step and launch the QM MDI engine:

* ``_free_tcp_port`` -- the rendezvous port helper.
* ``LAMMPS.ff_form`` -- detecting the QM path before touching ``_forcefield``.
* ``LAMMPS._mdi_engine_launch`` -- MDI-capability gating and resolving the
  owning program step's engine command.
* ``_mdi_launch_script`` -- composing the engine + driver launch script.

The step methods are called as *unbound* methods on a lightweight fake ``self``
so no real flowchart/executor is needed. See
docs/.../campaigns/2026-06-22/NOTES_C.rst.
"""

import shlex
import socket
import types

import pytest

from lammps_step.lammps import LAMMPS, _free_tcp_port, _mdi_launch_script


# ---------------------------------------------------------------------------
# Fakes
# ---------------------------------------------------------------------------
class _FakeStep:
    """Stand-in for a program step exposing get_mdi_engine_command.

    Records the keyword arguments it was called with and echoes the port back
    in the argv, so a test can confirm the port was threaded through.
    """

    def __init__(self, argv):
        self._argv = list(argv)
        self.calls = []

    def get_mdi_engine_command(
        self,
        executor,
        seamm_options,
        *,
        method,
        port,
        hostname="localhost",
        charge=0,
        multiplicity=1,
        n_atoms=None,
    ):
        self.calls.append(
            {
                "executor": executor,
                "seamm_options": seamm_options,
                "method": method,
                "port": port,
                "hostname": hostname,
                "charge": charge,
                "multiplicity": multiplicity,
                "n_atoms": n_atoms,
            }
        )
        return [*self._argv, "--port", str(port)]


class _FakePluginManager:
    def __init__(self, mapping):
        self._mapping = mapping

    def get(self, name):
        return self._mapping[name]


def _fake_self(mc=None, variables=None, step=None):
    """Build a fake LAMMPS node exposing only what the methods under test use."""
    variables = dict(variables or {})
    if mc is not None:
        variables["_model_chemistry"] = mc

    def get_variable(name):
        return variables[name]

    def variable_exists(name):
        return name in variables

    plugin_manager = None
    if step is not None and mc is not None:
        plugin_manager = _FakePluginManager({mc["step"]: step})

    return types.SimpleNamespace(
        get_variable=get_variable,
        variable_exists=variable_exists,
        global_options={"root": "~/SEAMM"},
        flowchart=types.SimpleNamespace(
            plugin_manager=plugin_manager,
            executor="FAKE-EXECUTOR",
        ),
    )


def _mc(mdi_capable=True, periodic_mdi=False):
    # The _model_chemistry wrapper as the Model Chemistry step now stores it:
    # a level spec (owner/type/method/...), not a full model chemistry.
    return {
        "level": "MOPAC:SQM@PM6-ORG",
        "owner": "MOPAC",
        "type": "SQM",
        "method": "PM6-ORG",
        "basis": None,
        "cutoff": None,
        "step": "mopac-step",
        "options": {
            "mdi_capable": mdi_capable,
            "periodic_mdi": periodic_mdi,
            "mdi_method_arg": "PM6-ORG",
        },
    }


# ---------------------------------------------------------------------------
# _free_tcp_port
# ---------------------------------------------------------------------------
def test_free_tcp_port_is_a_usable_port():
    port = _free_tcp_port()
    assert isinstance(port, int)
    assert 1 <= port <= 65535
    # It should be free right now -- we can bind it ourselves.
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        s.bind(("localhost", port))
    finally:
        s.close()


# ---------------------------------------------------------------------------
# ff_form
# ---------------------------------------------------------------------------
def test_ff_form_detects_model_chemistry():
    """_model_chemistry present -> "MDI/QM", without reading _forcefield."""
    me = _fake_self(mc=_mc())
    assert LAMMPS.ff_form(me) == "MDI/QM"


def test_ff_form_does_not_touch_forcefield_on_qm_path():
    """The QM branch must come first: _forcefield does not exist in a QM-MD
    flowchart, so reading it would raise."""
    me = _fake_self(mc=_mc())  # only "_model_chemistry" is defined
    # If ff_form read "_forcefield" it would KeyError; it must not.
    assert LAMMPS.ff_form(me) == "MDI/QM"


def test_ff_form_classical_path_unchanged():
    me = _fake_self(variables={"_forcefield": "OpenKIM"})
    assert LAMMPS.ff_form(me) == "OpenKIM"


# ---------------------------------------------------------------------------
# _mdi_engine_launch
# ---------------------------------------------------------------------------
def test_mdi_engine_launch_happy_path():
    step = _FakeStep(["conda", "run", "-n", "seamm-mopac", "python", "mopac_mdi.py"])
    me = _fake_self(mc=_mc(mdi_capable=True), step=step)
    configuration = types.SimpleNamespace(
        periodicity=0, charge=-1, spin_multiplicity=2, n_atoms=45
    )

    engine_argv, port = LAMMPS._mdi_engine_launch(me, configuration)

    assert isinstance(port, int) and 1 <= port <= 65535
    # The owning step was resolved by its Stevedore handle and called correctly.
    assert len(step.calls) == 1
    call = step.calls[0]
    assert call["method"] == "PM6-ORG"
    assert call["port"] == port
    assert call["hostname"] == "localhost"
    assert call["charge"] == -1  # from the configuration (D6)
    assert call["multiplicity"] == 2  # from the configuration (D6)
    assert call["executor"] == "FAKE-EXECUTOR"
    assert call["seamm_options"] == {"root": "~/SEAMM"}
    # The returned argv is exactly what the step produced.
    assert engine_argv[-2:] == ["--port", str(port)]


def test_mdi_engine_launch_rejects_non_mdi_capable():
    step = _FakeStep(["python", "engine.py"])
    me = _fake_self(mc=_mc(mdi_capable=False), step=step)
    configuration = types.SimpleNamespace(periodicity=0, charge=0, spin_multiplicity=1)
    with pytest.raises(ValueError, match="cannot be driven"):
        LAMMPS._mdi_engine_launch(me, configuration)
    assert step.calls == []  # never reached the engine


def test_mdi_engine_launch_rejects_periodic_without_periodic_mdi():
    step = _FakeStep(["python", "engine.py"])
    me = _fake_self(mc=_mc(mdi_capable=True, periodic_mdi=False), step=step)
    configuration = types.SimpleNamespace(
        periodicity=3, charge=0, spin_multiplicity=1, n_atoms=45
    )
    with pytest.raises(ValueError, match="periodic"):
        LAMMPS._mdi_engine_launch(me, configuration)
    assert step.calls == []


def test_mdi_engine_launch_allows_periodic_when_validated():
    step = _FakeStep(["python", "engine.py"])
    me = _fake_self(mc=_mc(mdi_capable=True, periodic_mdi=True), step=step)
    configuration = types.SimpleNamespace(
        periodicity=3, charge=0, spin_multiplicity=1, n_atoms=45
    )
    engine_argv, port = LAMMPS._mdi_engine_launch(me, configuration)
    assert len(step.calls) == 1
    assert engine_argv[-2:] == ["--port", str(port)]


# ---------------------------------------------------------------------------
# model_chemistry -- the full driver:task|level provenance label
# ---------------------------------------------------------------------------
def test_model_chemistry_label_composes_full_string():
    me = _fake_self(mc=_mc())
    # MDI/QM: LAMMPS drives, MOPAC owns the PES -> owner kept on the level side.
    assert LAMMPS.model_chemistry(me, "MD") == "LAMMPS:MD|MOPAC:SQM@PM6-ORG"
    assert LAMMPS.model_chemistry(me, "OPT") == "LAMMPS:OPT|MOPAC:SQM@PM6-ORG"


def test_model_chemistry_label_classical_falls_back():
    """With no _model_chemistry (classical/MLFF/OpenKIM) the legacy bare label
    is returned unchanged -- no spurious grammar string."""
    me = types.SimpleNamespace(
        variable_exists=lambda name: False,
        model="OPLS-AA",
    )
    assert LAMMPS.model_chemistry(me, "MD") == "OPLS-AA"


# ---------------------------------------------------------------------------
# _mdi_launch_script
# ---------------------------------------------------------------------------
def test_launch_script_shape_and_port_threading():
    port = 54321
    engine_argv = [
        "conda",
        "run",
        "--live-stream",
        "-n",
        "seamm-mopac",
        "python",
        "/abs/mopac_mdi.py",
        "-mdi",
        "-role ENGINE -name MOPAC -method TCP -port 54321 -hostname localhost",
        "--method",
        "PM6-ORG",
    ]
    config = {"code": "mpirun -np {NTASKS} lmp", "installation": "conda"}
    ce = {"NTASKS": 4}

    script = _mdi_launch_script(engine_argv, port, config, ce)
    lines = script.splitlines()

    # Shape: shebang, set -e, engine &, capture pid, driver, wait.
    assert lines[0] == "#!/bin/bash"
    assert lines[1] == "set -e"
    assert lines[2].endswith(" &")  # engine backgrounded first
    assert "mopac_mdi.py" in lines[2]
    assert lines[3] == "ENGINE_PID=$!"
    assert lines[-1] == "wait $ENGINE_PID"

    # The engine line is exactly shlex.join(argv) + " &" -- so the -mdi value
    # (which contains spaces) is quoted as one token.
    assert lines[2] == shlex.join(engine_argv) + " &"
    assert (
        "'-role ENGINE -name MOPAC -method TCP -port 54321 -hostname localhost'"
        in lines[2]
    )

    # The driver line: ini template resolved, MDI DRIVER flag with the SAME
    # port/hostname, and the input file appended.
    driver_line = lines[4]
    assert driver_line.startswith("mpirun -np 4 lmp")  # {NTASKS} -> 4
    assert (
        '-mdi "-role DRIVER -name LAMMPS -method TCP -port 54321 -hostname localhost"'
        in driver_line
    )
    assert driver_line.endswith("-in input.dat")


def test_launch_script_appends_cmd_args():
    script = _mdi_launch_script(
        ["engine"], 7000, {"code": "lmp", "cmd-args": "-sf omp"}, {"NTASKS": 1}
    )
    assert "lmp -sf omp -mdi " in script


def test_launch_script_custom_hostname():
    script = _mdi_launch_script(
        ["engine"], 7000, {"code": "lmp"}, {"NTASKS": 1}, hostname="node07"
    )
    assert "-hostname node07" in script
