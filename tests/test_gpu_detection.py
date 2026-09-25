# -*- coding: utf-8 -*-
"""Tests for the nvidia-smi based GPU detection that replaced GPUtil."""

from lammps_step.lammps import _available_gpus, _parse_nvidia_smi

SAMPLE = """0, 3, 512, 40960
1, 95, 39000, 40960
2, 0, 0, 40960
3, 10, 8000, 40960
"""


def test_parse_selects_idle_gpus():
    # load and memory both at or below 20%
    assert _parse_nvidia_smi(SAMPLE, 0.20) == [0, 2, 3]
    # tighter threshold drops GPU 3 (load 10% ok, memory 19.5% ok) only at 5%
    assert _parse_nvidia_smi(SAMPLE, 0.05) == [0, 2]
    # nothing qualifies at 0
    assert _parse_nvidia_smi(SAMPLE, 0.0) == [2]


def test_parse_skips_bad_lines():
    text = "0, [N/A], 100, 1000\n1, 5, 10, 1000\ngarbage\n2, 5, 10, 0\n"
    assert _parse_nvidia_smi(text, 0.5) == [1]


def test_no_nvidia_smi(monkeypatch):
    monkeypatch.setattr("lammps_step.lammps.shutil.which", lambda name: None)
    assert _available_gpus(0.5) == []
