# -*- coding: utf-8 -*-

"""Definition of unit systems in LAMMPS."""
from seamm_util import ureg, Q_, units_class  # noqa: F401

lammps_unit_system = "real"

unit_systems = {
    "real": {
        "mass": "grams/mole",
        "distance": "angstroms",
        "time": "femtoseconds",
        "energy": "kcal/mole",
        "velocity": "angstroms/femtosecond",
        "force": "kcal/mole/angstrom",
        "torque": "kcal/mole",
        "temperature": "kelvin",
        "pressure": "atmospheres",
        "dynamic viscosity": "poise",
        "charge": "elementary_charge",
        "dipole": "elementary_charge*angstroms",
        "electric field": "volts/angstrom",
        "density": "gram/cm^3",
    },
    "metal": {
        "mass": "grams/mole",
        "distance": "angstroms",
        "time": "picoseconds",
        "energy": "eV",
        "velocity": "angstroms/picosecond",
        "force": "eV/angstrom",
        "torque": "eV",
        "temperature": "kelvin",
        "pressure": "bars",
        "dynamic viscosity": "poise",
        "charge": "elementary_charge",
        "dipole": "elementary_charge*angstroms",
        "electric field": "volts/angstrom",
        "density": "gram/cm^3",
    },
    "si": {
        "mass": "kilograms",
        "distance": "meters",
        "time": "seconds",
        "energy": "joules",
        "velocity": "meters/second",
        "force": "newtons",
        "torque": "newton*meters",
        "temperature": "kelvin",
        "pressure": "pascals",
        "dynamic viscosity": "pascal*second",
        "charge": "coulombs",
        "dipole": "coulombs*meters",
        "electric field": "volts/meter",
        "density": "kilograms/meter^3",
    },
    "cgs": {
        "mass": "grams",
        "distance": "centimeters",
        "time": "seconds",
        "energy": "ergs",
        "velocity": "centimeters/second",
        "force": "dynes",
        "torque": "dyne*centimeters",
        "temperature": "kelvin",
        "pressure": "dyne/cm^2",
        "dynamic viscosity": "poise",
        "charge": "statcoulombs",
        "dipole": "statcoulombs*cm",
        "electric field": "statvolt/cm",
        "density": "grams/cm^3",
    },
    "electron": {
        "mass": "amu",
        "distance": "bohr",
        "time": "femtoseconds",
        "energy": "hartrees",
        "velocity": "bohr/a_u_time",
        "force": "hartrees/bohr",
        "temperature": "kelvin",
        "pressure": "pascals",
        "charge": "elementary_charge",
        "dipole moment": "debye",
        "electric field": "volts/cm",
    },
    "micro": {
        "mass": "picograms",
        "distance": "micrometers",
        "time": "microseconds",
        "energy": "picogram*micrometer^2/microsecond^2",
        "velocity": "micrometers/microsecond",
        "force": "picogram*micrometer/microsecond^2",
        "torque": "picogram*micrometer^2/microsecond^2",
        "temperature": "kelvin",
        "pressure": "picogram/(micrometer*microsecond^2)",
        "dynamic viscosity": "picogram/(micrometer*microsecond)",
        "charge": "picocoulombs",
        "dipole": "picocoulomb*micrometer",
        "electric field": "volt/micrometer",
        "density": "picograms/micrometer^3",
    },
    "nano": {
        "mass": "attograms",
        "distance": "nanometers",
        "time": "nanoseconds",
        "energy": "attogram*nanometer^2/nanosecond^2",
        "velocity": "nanometers/nanosecond",
        "force": "attogram*nanometer/nanosecond^2",
        "torque": "attogram*nanometer^2/nanosecond^2",
        "temperature": "kelvin",
        "pressure": "attogram/(nanometer*nanosecond^2)",
        "dynamic viscosity": "attogram/(nanometer*nanosecond)",
        "charge": "elementary_charge",
        "dipole": "elementary_charge*nanometer",
        "electric field": "volt/nanometer",
        "density": "attograms/nanometer^3",
    },
}

quantity_to_dimensionality = {
    "mass": ["[mass] / [substance]", "[mass]"],
    "distance": ["[length]"],
    "time": ["[time]"],
    "energy": [
        "[length] ** 2 * [mass] / [substance] / [time] ** 2",
        "[length] ** 2 * [mass] / [time] ** 2",
    ],
    "velocity": ["[length] / [time]"],
    "force": [
        "[mass] * [length] / [time] ** 2 / [substance]",
        "[mass] * [length] / [time] ** 2",
    ],
    "torque": [
        "[length] ** 2 * [mass] / [substance] / [time] ** 2",
        "[length] ** 2 * [mass] / [time] ** 2",
    ],
    "temperature": ["[temperature]"],
    "pressure": ["[mass] / [length] / [time] ** 2"],
    "dynamic viscosity": ["[mass] / [length] / [time]"],
    "charge": ["[current] * [time]", "[length] ** 1.5 * [mass] ** 0.5 / [time]"],
    "dipole": [
        "[current] * [length] * [time]",
        "[length] ** 2.5 * [mass] ** 0.5 / [time]",
    ],
    "electric field": [
        "[length] * [mass] / [current] / [time] ** 3",
        "[mass] ** 0.5 / [length] ** 0.5 / [time]",
    ],
    "density": ["[mass] / [length] ** 3"],
    "dipole moment": ["[current] * [length] * [time]"],
}

dimensionality_to_quantity = {
    "[mass] / [substance]": ["mass"],
    "[mass]": ["mass"],
    "[length]": ["distance"],
    "[time]": ["time"],
    "[length] ** 2 * [mass] / [substance] / [time] ** 2": ["energy", "torque"],
    "[mass] * [length] ** 2 / [time] ** 2 / [substance]": ["energy", "torque"],
    "[length] ** 2 * [mass] / [time] ** 2": ["energy", "torque"],
    "[length] / [time]": ["velocity"],
    "[length] * [mass] / [substance] / [time] ** 2": ["force"],
    "[length] * [mass] / [time] ** 2": ["force"],
    "[mass] * [length] / [time] ** 2 / [substance]": ["force"],
    "[mass] * [length] / [time] ** 2": ["force"],
    "[temperature]": ["temperature"],
    "[mass] / [length] / [time] ** 2": ["pressure"],
    "[mass] / [length] / [time]": ["dynamic viscosity"],
    "[current] * [time]": ["charge"],
    "[length] ** 1.5 * [mass] ** 0.5 / [time]": ["charge"],
    "[current] * [length] * [time]": ["dipole", "dipole moment"],
    "[length] ** 2.5 * [mass] ** 0.5 / [time]": ["dipole"],
    "[length] * [mass] / [current] / [time] ** 3": ["electric field"],
    "[mass] ** 0.5 / [length] ** 0.5 / [time]": ["electric field"],
    "[mass] / [length] ** 3": ["density"],
}


def set_lammps_unit_system(value):
    """Set the default unit system for LAMMPS."""
    global lammps_unit_system
    if value not in unit_systems:
        raise ValueError(f"unit system '{value}' not recognized!")
    lammps_unit_system = value


def get_lammps_unit_system():
    """Get the default unit system for LAMMPS."""
    return lammps_unit_system


def to_lammps_units(value, units=None, quantity=None, unit_system=None):
    """Convert the value with units to the current LAMMPS units."""
    if not isinstance(value, units_class):
        if units is not None:
            value = Q_(value, units)
        else:
            value = Q_(value)

    if quantity is None:
        quantities = dimensionality_to_quantity[str(value.dimensionality)]
        quantity = quantities[0]

    if unit_system is None:
        unit_system = lammps_unit_system

    unit_expr = unit_systems[unit_system][quantity]
    converted = value.to(unit_expr)

    return converted.magnitude


def from_lammps_units(value, units, quantity=None, unit_system=None):
    """Convert the scalar value from LAMMPS units to requested units"""
    if quantity is None:
        dimensions = Q_(1.0, units).dimensionality
        quantities = dimensionality_to_quantity[str(dimensions)]
        quantity = quantities[0]

    if unit_system is None:
        unit_system = lammps_unit_system

    unit_expr = unit_systems[unit_system][quantity]
    converted = Q_(value, unit_expr).to(units)

    return converted


if __name__ == "__main__":  # pragma: no cover
    import pint
    import json

    for unit in ureg:
        try:
            q = ureg(unit)
        except pint.errors.UndefinedUnitError:
            q = "unknown"
        print(f"{unit} --> {q}")
    print("")

    dimensionality = {}
    for system, units in unit_systems.items():
        print(f"{system} system of units")
        for quantity, unit_expr in units.items():
            print(f"\t{quantity} uses {unit_expr}")
            q = ureg(unit_expr)
            print(f"\t\t{q}")
            dimensions = str(q.dimensionality)
            if quantity not in dimensionality:
                dimensionality[quantity] = [dimensions]
            elif dimensions not in dimensionality[quantity]:
                dimensionality[quantity].append(dimensions)
                print(f"{system}: {quantity} <-- {unit_expr}")

    print()
    print(json.dumps(dimensionality, indent=4))

    inverse = {}
    for quantity, unit_exprs in dimensionality.items():
        for unit_expr in unit_exprs:
            if unit_expr not in inverse:
                inverse[unit_expr] = [quantity]
            elif quantity != inverse[unit_expr]:
                inverse[unit_expr].append(quantity)
                print(
                    f"quantity {quantity} != {inverse[unit_expr]} for {unit_expr}"
                )  # noqa: E501

    print()
    print(json.dumps(inverse, indent=4))

    print()
    value = "2.0 kJ/mol"
    converted = to_lammps_units(value)
    print(f"{value} in real units is {converted:.4}")
    print(f'{converted:.4} --> {from_lammps_units(converted, "kJ/mol"):~.4}')

    lammps_unit_system = "metal"
    converted = to_lammps_units(value)
    print(f"{value} in metal units is {converted:.4}")
    print(f'{converted:.4} --> {from_lammps_units(converted, "kJ/mol"):~.4}')
