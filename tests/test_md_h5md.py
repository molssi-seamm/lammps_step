#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for converting LAMMPS output to SEAMM's H5MD file (md.h5md)."""

from pathlib import Path
import shutil
import types

import numpy as np
import pytest

from seamm_util import Q_
from seamm_md_analysis import open_md, read_trajectory

from lammps_step.md_h5md import MDH5MDConverter, h5md_units, write_md_h5md

DATA = Path(__file__).parent / "data" / "md_h5md"


def without_step0(columns):
    """The native series without the step-0 sample, which md.h5md leaves out."""
    keep = columns["tstep"].astype(int) != 0
    return {k: v[keep] if len(v) == len(keep) else v for k, v in columns.items()}


class NVT:
    """Just what the converter uses of an NVT sub-step."""

    def __init__(self, directory, charges):
        self.directory = str(directory)
        self._id = (3, 3)
        self.title = "NVT dynamics"
        self.model = "oplsaa"
        self.parent = types.SimpleNamespace(eex={"charges": charges})
        self.variables = {}

    def timestep(self, value):
        return 1.0, Q_(1.0, "fs")

    def variable_exists(self, name):
        return name in self.variables

    def get_variable(self, name):
        return self.variables[name]

    def set_variable(self, name, value):
        self.variables[name] = value


def configuration(n_molecules, symbols, numbers, masses, cell=None):
    n = n_molecules * len(symbols)
    per = len(symbols)
    return types.SimpleNamespace(
        periodicity=0 if cell is None else 3,
        atoms=types.SimpleNamespace(
            ids=list(range(1, n + 1)),
            atomic_numbers=numbers * n_molecules,
            atomic_masses=masses * n_molecules,
            symbols=symbols * n_molecules,
        ),
        cell=types.SimpleNamespace(parameters=cell),
        find_molecules=lambda as_indices=True: [
            list(range(i * per, (i + 1) * per)) for i in range(n_molecules)
        ],
        system=types.SimpleNamespace(name="test", id=1),
        name="start",
        id=2,
    )


P = {"time": Q_(2.0, "ps"), "timestep": Q_(1.0, "fs"), "T0": Q_(300.0, "K")}


def test_units_real():
    units = h5md_units("real")
    assert units["time"] == "fs"
    assert units["distance"] == "Angstrom"
    assert units["velocity"] == "Angstrom fs-1"
    assert units["force"] == "kcal mol-1 Angstrom-1"
    assert units["energy"] == "kcal mol-1"
    assert units["pressure"] == "atm"
    assert units["density"] == "g cm-3"


def test_units_metal():
    units = h5md_units("metal")
    assert units["time"] == "ps"
    assert units["velocity"] == "Angstrom ps-1"
    assert units["pressure"] == "bar"


def test_periodic_water(tmp_path):
    for path in (DATA / "water").iterdir():
        shutil.copy(path, tmp_path)
    node = NVT(tmp_path, [-0.82, 0.41, 0.41] * 64)
    config = configuration(
        64, ["O", "H", "H"], [8, 1, 1], [15.999, 1.008, 1.008], cell=[12.4172] * 3
    )
    config.cell.parameters = [12.417241581015336] * 3 + [90.0] * 3
    text = write_md_h5md(node, P, config, wall_time="12.5")
    assert text.startswith("Wrote the MD data to md.h5md")

    (record,) = node.variables["_md_runs"]
    assert record["ensemble"] == "NVT"
    assert record["periodic"]
    assert record["observables"]["heat_flux"] == "2 fs"
    assert record["observables"]["pressure_tensor"] == "5 fs"
    assert record["observables"]["temperature"] == "10 fs"
    assert record["particles"]["all/position"] == "100 fs"
    assert record["particles"]["molecules/position"] == "50 fs"
    assert record["time"] == ["0 fs", "2000 fs"]

    with open_md(tmp_path / "md.h5md") as md:
        assert md.periodic
        assert md.parameters["temperature"] == 300.0
        assert md.provenance["segments"]["wall_time"] == [12.5]

        # The pressure tensor sampled on its own replaces the one in the state
        _, columns = read_trajectory(tmp_path / "pressure_tensor.trj")
        columns = without_step0(columns)
        PT = md.observable("pressure_tensor")
        assert PT.values.shape == (400, 3, 3)
        assert PT.dt == Q_(5.0, "fs")
        assert PT.units == Q_(1, "atm").units
        assert np.allclose(PT.values[:, 0, 0], columns["Pxx"])
        assert np.allclose(PT.values[:, 2, 2], columns["Pzz"])
        assert np.allclose(PT.values[:, 0, 1], columns["Pxy"])
        assert np.allclose(PT.values[:, 1, 0], columns["Pxy"])
        assert np.allclose(PT.values[:, 1, 2], columns["Pyz"])
        assert "pressure_offdiagonal" not in md.observables()

        J = md.observable("heat_flux")
        _, columns = read_trajectory(tmp_path / "heat_flux.trj")
        columns = without_step0(columns)
        assert np.allclose(J.values[:, 2], columns["Jz"])
        assert J.units == Q_(1, "W/m^2").units
        assert J.dt == Q_(2.0, "fs")

        # The trajectory gives positions, velocities, forces and the box
        assert md.n_frames() == 3
        assert md.has_particles("velocity")
        assert md.has_particles("force")
        assert np.allclose(md.charges()[:3], [-0.82, 0.41, 0.41])
        assert np.allclose(np.diag(md.cell()), 12.417241581015336)
        steps, _, frames = next(md.iter_frames(fields=("position", "velocity")))
        assert list(steps) == [0, 100, 200]
        assert np.allclose(frames["position"][0, 0], [10.412, 7.1, 1.3])
        assert np.allclose(frames["velocity"][0, 0], [-0.00718, 0.00657, -0.00326])

        com = md.particle_series("position", group="molecules")
        assert com.values.shape == (3, 64, 3)


def test_pressure_tensor_from_the_state(tmp_path):
    """Without its own sampling, the pressure tensor comes from the state file."""
    for path in (DATA / "water").iterdir():
        if path.name != "pressure_tensor.trj":
            shutil.copy(path, tmp_path)
    node = NVT(tmp_path, [-0.82, 0.41, 0.41] * 64)
    config = configuration(
        64, ["O", "H", "H"], [8, 1, 1], [15.999, 1.008, 1.008], cell=[12.4172] * 3
    )
    config.cell.parameters = [12.417241581015336] * 3 + [90.0] * 3
    write_md_h5md(node, P, config)
    (record,) = node.variables["_md_runs"]
    assert record["observables"]["pressure_tensor"] == "10 fs"

    with open_md(tmp_path / "md.h5md") as md:
        # The pressure tensor is minus the stress, and its trace is the pressure
        _, columns = read_trajectory(tmp_path / "nvt_state.trj")
        columns = without_step0(columns)
        PT = md.observable("pressure_tensor")
        assert np.allclose(PT.values[:, 0, 0], -columns["Sxx"])
        assert np.allclose(PT.values[:, 0, 1], -columns["Sxy"])
        assert np.allclose(PT.values[:, 1, 0], -columns["Sxy"])
        pressure = md.observable("pressure").values
        trace = np.trace(PT.values, axis1=1, axis2=2) / 3
        assert np.allclose(trace, pressure, rtol=1e-4, atol=1)


def test_non_periodic_ethanol(tmp_path):
    for path in (DATA / "ethanol").iterdir():
        shutil.copy(path, tmp_path)
    node = NVT(tmp_path, [-0.18, 0.145, -0.683] + [0.06] * 5 + [0.418])
    config = configuration(
        1,
        ["C", "C", "O"] + ["H"] * 6,
        [6, 6, 8] + [1] * 6,
        [12.011, 12.011, 15.999] + [1.008] * 6,
    )
    P1 = dict(P, time=Q_(1.0, "ps"))
    converter = MDH5MDConverter(node, P1, config)
    record = converter.write()
    assert not record["periodic"]
    with open_md(tmp_path / "md.h5md") as md:
        assert not md.periodic
        assert md.cell() is None
        assert "pressure" not in md.observables()
        assert "temperature" in md.observables()
        assert md.has_particles("velocity")


def test_failure_is_reported_not_raised(tmp_path):
    node = NVT(tmp_path, None)
    text = write_md_h5md(
        node, {"time": Q_(1, "ps")}, configuration(1, ["O"], [8], [16.0])
    )
    assert text.startswith("Warning: could not write")
    assert not (tmp_path / "md.h5md").exists()
    assert not (tmp_path / "md.h5md.tmp").exists()


@pytest.fixture(autouse=True)
def _real_units():
    import lammps_step

    old = lammps_step.get_lammps_unit_system()
    lammps_step.set_lammps_unit_system("real")
    yield
    lammps_step.set_lammps_unit_system(old)


def test_the_kspace_style_is_recorded(tmp_path):
    """The dielectric constant's formula assumes Ewald or PPPM: the run records
    its kspace style."""
    for path in (DATA / "water").iterdir():
        shutil.copy(path, tmp_path)
    (tmp_path / "input.dat").write_text(
        "pair_style          lj/cut/coul/long 12.0\n"
        "kspace_style        pppm 1.0e-5\n"
    )
    node = NVT(tmp_path, [-0.82, 0.41, 0.41] * 64)
    config = configuration(
        64, ["O", "H", "H"], [8, 1, 1], [15.999, 1.008, 1.008], cell=[12.4172] * 3
    )
    config.cell.parameters = [12.417241581015336] * 3 + [90.0] * 3
    write_md_h5md(node, P, config)
    with open_md(tmp_path / "md.h5md") as md:
        assert md.parameters["kspace_style"] == "pppm 1.0e-5"
    (tmp_path / "input.dat").write_text("pair_style lj/cut/coul/cut 12.0\n")
    (tmp_path / "md.h5md").unlink()
    write_md_h5md(node, P, config)
    with open_md(tmp_path / "md.h5md") as md:
        assert md.parameters["kspace_style"] == "none"


def test_the_fix_styles_are_recorded(tmp_path):
    """Whether a run can be treated as isotropic depends on its fixes."""
    from seamm_md_analysis import isotropy

    for path in (DATA / "water").iterdir():
        shutil.copy(path, tmp_path)
    (tmp_path / "input.dat").write_text(
        "fix                 1 all shake 0.0001 20 0 b 1 a 1\n"
        "fix                 2 all nvt temp 298.15 298.15 100.0  # thermostat\n"
        "fix wall all wall/lj93 zlo EDGE 1.0 1.0 2.5\n"
    )
    node = NVT(tmp_path, [-0.82, 0.41, 0.41] * 64)
    config = configuration(
        64, ["O", "H", "H"], [8, 1, 1], [15.999, 1.008, 1.008], cell=[12.4172] * 3
    )
    config.cell.parameters = [12.417241581015336] * 3 + [90.0] * 3
    write_md_h5md(node, P, config)
    with open_md(tmp_path / "md.h5md") as md:
        assert list(md.parameters["fixes"]) == ["nvt", "shake", "wall/lj93"]
        assert isotropy(md) == (False, "the run used wall/lj93")
