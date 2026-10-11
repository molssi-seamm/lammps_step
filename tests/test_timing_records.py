# -*- coding: utf-8 -*-
"""The timing records LAMMPS runs write (seamm_exec campaign 2026-10-05)."""

from types import SimpleNamespace

from lammps_step.lammps import timing_descriptors

LOG = """\
Loop time of 1.234 on 4 procs for 1000 steps with 648 atoms
Performance: 70.0 ns/day
Loop time of 12.5 on 4 procs for 10000 steps with 648 atoms
Total wall time: 0:00:15
"""


def test_descriptors():
    conf = SimpleNamespace(
        atoms=SimpleNamespace(atomic_numbers=[8, 1, 1] * 216),
        charge=0,
        spin_multiplicity=1,
        periodicity=3,
        volume=6450.0,
    )
    d = timing_descriptors([["Initialization", {}], ["NVT dynamics", {}]], LOG, conf)
    assert d["n_calculations"] == 2 and "NVT dynamics" in d["calculations"]
    assert d["n_atoms"] == 648 and d["periodicity"] == 3 and d["volume"] == 6450.0
    assert d["n_runs"] == 2 and d["md_steps"] == 11000
    assert abs(d["code_seconds"] - 13.734) < 1e-9
    assert d["procs"] == 4 and d["total_wall"] == 15
    assert d["terminated_normally"] is True


def test_record_timing_needs_no_run_locals(monkeypatch):
    """2026.10.6 passed a local bound only in the MDI branch to record_timing,
    so every forcefield run raised UnboundLocalError after LAMMPS finished
    (ChemAI jobs 5177-5188). The record now takes only what every run has."""
    import inspect

    import seamm_exec

    from lammps_step.lammps import LAMMPS

    source = inspect.getsource(LAMMPS._execute_single_sim)
    assert "self.record_timing(ce, t, result)" in source
    assert "record_timing(configuration" not in source

    node = LAMMPS.__new__(LAMMPS)
    monkeypatch.setattr(LAMMPS, "directory", property(lambda self: "/nonexistent"))
    node.logger = __import__("logging").getLogger("test")
    node._timing_control = [["NVT dynamics", {}]]
    conf = SimpleNamespace(
        atoms=SimpleNamespace(atomic_numbers=[8, 1, 1]),
        charge=0,
        spin_multiplicity=1,
        periodicity=0,
    )
    node.get_system_configuration = lambda *a, **k: (None, conf)
    recorded = {}
    monkeypatch.setattr(
        seamm_exec,
        "record_timing",
        lambda program, wall, descriptors, **kw: recorded.update(
            program=program, wall=wall, descriptors=descriptors, **kw
        ),
    )
    node.record_timing({"NTASKS": 2, "NGPUS": 0}, 1.5, True)
    assert recorded["program"] == "lammps" and recorded["wall"] == 1.5
    assert recorded["descriptors"]["n_atoms"] == 3 and recorded["ntasks"] == 2


def test_timing_spec():
    from lammps_step import lammps

    assert lammps.TIMING_SPEC["units"] == "md_steps"
    # The forcefield is the class: an MLFF and a classical forcefield differ by
    # orders of magnitude per atom-step
    assert lammps.TIMING_SPEC["klass"] == ["ff_form", "model"]


def _node_and_recorder(monkeypatch, ff_form, model):
    import seamm_exec

    from lammps_step.lammps import LAMMPS

    node = LAMMPS.__new__(LAMMPS)
    monkeypatch.setattr(LAMMPS, "directory", property(lambda self: "/nonexistent"))
    node.logger = __import__("logging").getLogger("test")
    node._timing_control = [["NVT dynamics", {}]]
    conf = SimpleNamespace(
        atoms=SimpleNamespace(atomic_numbers=[8, 1, 1] * 100),
        charge=0,
        spin_multiplicity=1,
        periodicity=3,
        volume=3000.0,
    )
    node.get_system_configuration = lambda *a, **k: (None, conf)
    node.model = model
    if isinstance(ff_form, Exception):

        def raise_it():
            raise ff_form

        node.ff_form = raise_it
    else:
        node.ff_form = lambda: ff_form
    recorded = {}
    monkeypatch.setattr(
        seamm_exec,
        "record_timing",
        lambda program, wall, descriptors, **kw: recorded.update(
            descriptors=descriptors, **kw
        ),
    )
    return node, recorded


def test_record_names_the_forcefield(monkeypatch):
    node, recorded = _node_and_recorder(monkeypatch, "OPLS-AA", "oplsaa+")
    node.record_timing({"NTASKS": 8, "NGPUS": 0}, 100.0, True)
    assert recorded["descriptors"]["ff_form"] == "OPLS-AA"
    assert recorded["descriptors"]["model"] == "oplsaa+"
    assert recorded["gpu_model"] is None  # no GPU, no nvidia-smi


def test_record_names_the_gpu(monkeypatch):
    from lammps_step import lammps

    monkeypatch.setattr(lammps, "_gpu_model", lambda: "NVIDIA A30")
    node, recorded = _node_and_recorder(monkeypatch, "PyTorch", "PyTorch/mace_les_r3")
    node.record_timing({"NTASKS": 1, "NGPUS": 1}, 100.0, True)
    assert recorded["descriptors"]["ff_form"] == "PyTorch"
    assert recorded["descriptors"]["model"] == "PyTorch/mace_les_r3"
    assert recorded["ngpus"] == 1 and recorded["gpu_model"] == "NVIDIA A30"


def test_record_survives_a_missing_forcefield(monkeypatch):
    """No forcefield variable (e.g. a hand-run step): the record is still written."""
    node, recorded = _node_and_recorder(monkeypatch, KeyError("_forcefield"), None)
    node.record_timing({"NTASKS": 2}, 10.0, True)
    assert recorded["descriptors"]["ff_form"] == ""
    assert recorded["descriptors"]["model"] == ""
