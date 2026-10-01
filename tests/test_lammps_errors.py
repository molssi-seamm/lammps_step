# -*- coding: utf-8 -*-

"""Errors from LAMMPS itself are found in log.lammps and reported clearly, rather
than surfacing later as unrelated errors analyzing missing output."""

import pytest

import lammps_step
from lammps_step.lammps import lammps_error, lammps_error_message

# Trimmed from SEAMM_DEV job 3989
LOST_ATOMS = """\
# Step 3.1: Initialization
units               real
# Step 3.2: NPT dynamics
fix                 2 all press/berendsen  iso 1.000 1.000 500.000 modulus 19738.5
run                 1000
 0   0   704535.7   49266.982
ERROR: Lost atoms: original 2001 current 1999 (src/thermo.cpp:494)
Last command: run                 1000
"""

# Trimmed from SEAMM_DEV job 3987
TRICLINIC = """\
# Step 3.2: NPT dynamics
fix                 2 all press/berendsen  iso 1.000 1.000 500.000 modulus 19738.5
Generated 1 of 1 mixed pair_coeff terms from geometric mixing rule
ERROR: Cannot use fix press/berendsen with triclinic box \
(src/fix_press_berendsen.cpp:253)
Last command: run                 1000
"""


def test_lost_atoms():
    error = lammps_error(LOST_ATOMS)
    assert error["message"] == "Lost atoms: original 2001 current 1999"
    assert error["step"] == "Step 3.2: NPT dynamics"
    assert error["last command"] == "run 1000"
    assert "minimizing" in error["advice"]
    text = lammps_error_message(error, "/jobs/Job_1/3")
    assert text.startswith(
        "LAMMPS failed in Step 3.2: NPT dynamics: Lost atoms: original 2001 "
        "current 1999. The last command was 'run 1000'. Atoms moved too far"
    )
    assert text.endswith("See log.lammps in /jobs/Job_1/3 for details.")


def test_triclinic():
    error = lammps_error(TRICLINIC)
    assert error["message"] == "Cannot use fix press/berendsen with triclinic box"
    assert "triclinic" in error["advice"]


def test_error_on_another_process_is_only_on_the_screen():
    screen = "ERROR on proc 3: Out of range atoms - cannot compute PPPM (src/x.cpp:1)"
    error = lammps_error("run 100\n", screen)
    assert error["message"] == "Out of range atoms - cannot compute PPPM"
    assert error["step"] == ""
    assert "k-space" in error["advice"]
    assert lammps_error_message(error).startswith("LAMMPS failed: Out of range")


def test_no_error():
    assert lammps_error("run 1000\nTotal wall time: 0:00:01\n", "WARNING: x") is None


def test_check_raises_and_writes_no_success(tmp_path, monkeypatch):
    monkeypatch.setattr(
        lammps_step.LAMMPS, "directory", property(lambda self: str(tmp_path))
    )
    lammps = lammps_step.LAMMPS()
    lammps._id = ("3",)
    (tmp_path / "log.lammps").write_text(LOST_ATOMS)
    with pytest.raises(RuntimeError, match="LAMMPS failed in Step 3.2"):
        lammps.check_for_lammps_error({"stdout": "", "stderr": ""})
    assert not (tmp_path / "success.dat").exists()
