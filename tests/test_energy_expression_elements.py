"""The element list in the LAMMPS energy expression is per atom type.

LAMMPS's ``dump_modify ... element`` takes one symbol per atom *type*. The
PyTorch (MDI/ML potential) energy expression used to list one symbol per
*atom*, so any trajectory dump with an ML potential died with
``Unknown dump_modify keyword: H``.
"""

import pytest

from molsystem import SystemDB

import lammps_step


@pytest.fixture
def water_dimer_db():
    db = SystemDB(filename="file:eex_elements?mode=memory&cache=shared")
    system = db.create_system(name="water dimer")
    configuration = system.create_configuration(name="default")
    # O H H O H H -- two types, six atoms
    configuration.atoms.append(
        atno=[8, 1, 1, 8, 1, 1],
        x=[0.0, 0.96, -0.24, 2.9, 3.5, 3.5],
        y=[0.0, 0.0, 0.93, 0.0, 0.8, -0.8],
        z=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    )
    yield db
    db.close()


@pytest.mark.parametrize(
    "method", ["PyTorch_energy_expression", "OpenKIM_energy_expression"]
)
def test_elements_are_per_type(water_dimer_db, method, monkeypatch):
    node = lammps_step.Initialization()
    variables = {
        "_system_db": water_dimer_db,
        "_OpenKIM_Potential": "SW_StillingerWeber_1985_Si__MO_405512056662_006",
    }
    monkeypatch.setattr(node, "get_variable", lambda name: variables[name])

    eex = getattr(node, method)()

    assert eex["n_atoms"] == 6
    assert eex["n_atom_types"] == 2
    assert eex["elements"] == ["O", "H"]
    assert len(eex["elements"]) == len(eex["atom types"])
