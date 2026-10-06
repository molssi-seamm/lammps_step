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
