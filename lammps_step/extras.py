# -*- coding: utf-8 -*-
"""Extra commands and extra properties of the MD sub-steps (NVE, NVT and NPT).

*Extra commands* are LAMMPS ``compute``, ``fix``, ``variable`` and ``dump`` lines
that apply to one sub-step: SEAMM removes what they define when the sub-step ends.
*Extra properties* are named LAMMPS expressions sampled with the built-in
properties, analysed like them, and optionally treated as an autocorrelation or a
Green–Kubo integral. The Custom sub-step remains for commands that change the
system.
"""

import logging
import math
import re

import numpy as np

from seamm_util import Q_

import lammps_step

logger = logging.getLogger(__name__)

#: The commands allowed among the extra commands
ALLOWED_COMMANDS = ("compute", "fix", "variable", "dump")

#: The units in which properties of each LAMMPS quantity are shown, the same as
#: for the built-in properties
DISPLAY_UNITS = {
    "dimensionless": "",
    "energy": "kcal/mol",
    "pressure": "atm",
    "temperature": "K",
    "distance": "Å",
    "volume": "Å^3",
    "time": "fs",
    "density": "g/mL",
    "velocity": "Å/fs",
    "force": "kcal/mol/Å",
    "mass": "g/mol",
    "charge": "e",
    "dipole": "debye",
    "electric field": "V/Å",
    "dynamic viscosity": "mPa*s",
    "torque": "kcal/mol",
}

#: The treatments of an extra property
TREATMENTS = ("average", "autocorrelation", "Green-Kubo integral")

#: The named prefactors of a Green–Kubo integral, V the volume, k Boltzmann's
#: constant and T the temperature
PREFACTORS = ("V/kT", "V/kT^2", "1/(V kT)")

#: Units for the results of Green–Kubo integrals, when they have the
#: dimensions of one of these
RESULT_UNITS = ("mPa*s", "W/(m*K)", "S/m", "m^2/s", "1/ps")

#: The prefix of SEAMM's own variables wrapping the extra properties, and of
#: their columns in the state files
VARIABLE_PREFIX = "seamm_x_"
COLUMN_PREFIX = "x_"

#: IDs and names that SEAMM's own input defines, which the extra commands may not
#: reuse (numbered fixes, computes and dumps are refused separately)
RESERVED_IDS = {
    "S_p",
    "S_b",
    "flux_p",
    "flux_b",
    "flux",
    "hb",
    "ke_atom",
    "pe_atom",
    "stress_atom",
    "dynamics",
    "constraint",
    "time",
    "temp",
    "press",
    "etotal",
    "ke",
    "pe",
    "emol",
    "epair",
    "sxx",
    "syy",
    "szz",
    "syz",
    "sxz",
    "sxy",
    "pxx",
    "pyy",
    "pzz",
    "pyz",
    "pxz",
    "pxy",
    "Jx",
    "Jy",
    "Jz",
    "factor",
    "N_hbond",
    "E_hbond",
}

_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_]*$")
#: An explicit prefactor starts with a number
_NUMBER = re.compile(r"\s*[-+]?(\d|\.\d)")


class ExtrasError(ValueError):
    """A problem with the extra commands or extra properties."""


# ----------------------------------------------------------------------
# The extra commands
# ----------------------------------------------------------------------
def _strip_comment(line):
    """The line without a LAMMPS comment: '#' outside quotes to the end."""
    quote = None
    for i, c in enumerate(line):
        if quote is not None:
            if c == quote:
                quote = None
        elif c in ("'", '"'):
            quote = c
        elif c == "#":
            return line[:i]
    return line


def logical_lines(text):
    """The commands in ``text``, comments removed and continuation lines ('&' at
    the end) joined, as LAMMPS reads them."""
    lines = []
    current = ""
    for raw in (text or "").splitlines():
        line = _strip_comment(raw).rstrip()
        if line.endswith("&"):
            current += line[:-1] + " "
            continue
        current += line
        if current.strip() != "":
            lines.append(" ".join(current.split()))
        current = ""
    if current.strip() != "":
        lines.append(" ".join(current.split()))
    return lines


def parse_commands(text):
    """Check the extra commands and find what they define.

    Parameters
    ----------
    text : str
        The extra commands.

    Returns
    -------
    [dict]
        For each command, in order: ``command`` (compute, fix, variable or dump),
        ``id`` (its ID, or the variable's name) and ``line``.

    Raises
    ------
    ExtrasError
        Listing every problem: a command other than those allowed, a missing ID,
        an ID that SEAMM uses or that is all digits.
    """
    commands = []
    problems = []
    for line in logical_lines(text):
        words = line.split()
        command = words[0]
        if command not in ALLOWED_COMMANDS:
            problems.append(
                f"'{command}' is not allowed in the extra commands, which may only "
                f"be {', '.join(ALLOWED_COMMANDS)}. Use a Custom sub-step for "
                f"commands that change the system: {line}"
            )
            continue
        if len(words) < 3:
            problems.append(f"The {command} has no ID or style: {line}")
            continue
        _id = words[1]
        if command == "variable" and words[2] == "delete":
            problems.append(f"Variables are deleted automatically: {line}")
            continue
        if _id.isdigit():
            problems.append(
                f"The ID of a {command} may not be a number, which SEAMM uses for "
                f"its own: {line}"
            )
            continue
        if _id in RESERVED_IDS or _id.startswith(VARIABLE_PREFIX):
            problems.append(f"SEAMM uses the {command} ID '{_id}' itself: {line}")
            continue
        commands.append({"command": command, "id": _id, "line": line})
    if problems:
        raise ExtrasError("\n".join(problems))
    return commands


def cleanup_lines(commands):
    """The lines removing what the extra commands defined, in reverse order (a fix
    may use a compute, a variable another variable)."""
    lines = []
    seen = set()
    for item in reversed(commands):
        key = (item["command"], item["id"])
        if key in seen:
            continue
        seen.add(key)
        if item["command"] == "compute":
            lines.append(f"uncompute           {item['id']}")
        elif item["command"] == "fix":
            lines.append(f"unfix               {item['id']}")
        elif item["command"] == "dump":
            lines.append(f"undump              {item['id']}")
        else:
            lines.append(f"variable            {item['id']} delete")
    return lines


# ----------------------------------------------------------------------
# The extra properties
# ----------------------------------------------------------------------
def unit_kind(units):
    """Whether ``units`` names a LAMMPS quantity or explicit units.

    Returns
    -------
    (str, str)
        ``("quantity", name)`` or ``("explicit", units)``.

    Raises
    ------
    ExtrasError
        If they are neither a quantity nor units pint understands.
    """
    units = (units or "").strip()
    if units == "" or units == "dimensionless":
        return "quantity", "dimensionless"
    if units in DISPLAY_UNITS:
        return "quantity", units
    try:
        Q_(1.0, units)
    except Exception:
        raise ExtrasError(
            f"'{units}' is neither a LAMMPS quantity ({', '.join(DISPLAY_UNITS)}) "
            "nor units."
        )
    return "explicit", units


def check_properties(rows):
    """Check the extra properties and complete their definitions.

    Parameters
    ----------
    rows : [dict]
        The rows of the table: ``name``, ``expression``, ``units``,
        ``treatment``, ``prefactor`` and ``store``.

    Returns
    -------
    [dict]
        The rows with defaults filled in and ``column``, ``variable``, ``kind``
        (quantity or explicit) and ``display_units`` added.

    Raises
    ------
    ExtrasError
        Listing every problem.
    """
    result = []
    problems = []
    names = set()
    for row in rows or []:
        name = str(row.get("name", "")).strip()
        expression = str(row.get("expression", "")).strip()
        if name == "" and expression == "":
            continue  # an empty row of the table
        if not _NAME.match(name):
            problems.append(
                f"The name '{name}' must start with a letter and contain only "
                "letters, digits and '_'."
            )
            continue
        if name in names:
            problems.append(f"Two extra properties are named '{name}'.")
            continue
        names.add(name)
        if expression == "":
            problems.append(f"The extra property '{name}' has no expression.")
            continue
        if '"' in expression:
            problems.append(
                f"The expression of '{name}' may not contain a double quote."
            )
            continue
        treatment = row.get("treatment") or "average"
        if treatment not in TREATMENTS:
            problems.append(
                f"The treatment '{treatment}' of '{name}' is not one of "
                f"{', '.join(TREATMENTS)}."
            )
            continue
        try:
            kind, units = unit_kind(row.get("units"))
        except ExtrasError as e:
            problems.append(f"The units of '{name}': {e}")
            continue
        prefactor = str(row.get("prefactor") or "").strip()
        if treatment == "Green-Kubo integral":
            if prefactor == "":
                prefactor = "V/kT"
            if prefactor not in PREFACTORS:
                try:
                    if not _NUMBER.match(prefactor):
                        raise ValueError("not a number")
                    Q_(prefactor)
                except Exception:
                    problems.append(
                        f"The prefactor of '{name}' must be one of "
                        f"{', '.join(PREFACTORS)} or a value with units, not "
                        f"'{prefactor}'."
                    )
                    continue
        result.append(
            {
                "name": name,
                "expression": expression,
                "units": units,
                "kind": kind,
                "display_units": DISPLAY_UNITS[units] if kind == "quantity" else units,
                "treatment": treatment,
                "prefactor": prefactor,
                "store": bool(row.get("store", False)),
                "column": COLUMN_PREFIX + name,
                "variable": VARIABLE_PREFIX + name,
            }
        )
    if problems:
        raise ExtrasError("\n".join(problems))
    return result


def extra_input(P):
    """The input for the extra commands and properties of an MD sub-step.

    Parameters
    ----------
    P : dict
        The sub-step's parameters.

    Returns
    -------
    dict
        ``lines``: the commands and the variables wrapping the properties, to put
        after the ensemble's fix and before the sampling; ``properties``: the
        ``v_`` references to add to the sampled properties; ``titles``: their
        column titles; ``cleanup``: the lines removing everything defined, to put
        after the run; ``definitions``: the checked extra properties.
    """
    commands = parse_commands(P.get("extra commands", ""))
    definitions = check_properties(P.get("extra properties", []))

    lines = []
    if commands:
        lines.append("")
        lines.append("# Extra commands for this sub-step")
        lines.extend(item["line"] for item in commands)
    if definitions:
        lines.append("")
        lines.append("# Extra properties")
    for d in definitions:
        lines.append(f'variable            {d["variable"]} equal "{d["expression"]}"')
        commands.append({"command": "variable", "id": d["variable"], "line": ""})
    if lines:
        lines.append("")

    return {
        "lines": lines,
        "properties": " ".join(f"v_{d['variable']}" for d in definitions),
        "titles": " ".join(d["column"] for d in definitions),
        "cleanup": cleanup_lines(commands),
        "definitions": definitions,
    }


def check(P):
    """The problems with the extra commands and properties, as text, or ''."""
    problems = []
    for check_one, key in (
        (parse_commands, "extra commands"),
        (check_properties, "extra properties"),
    ):
        try:
            check_one(P.get(key))
        except ExtrasError as e:
            problems.append(str(e))
    return "\n".join(problems)


def column_metadata(definitions):
    """{column: (title, display units, convert)} for the analysis of the state
    file: ``convert`` is False for explicit units, which are the values' own."""
    return {
        d["column"]: (d["name"], d["display_units"], d["kind"] == "quantity")
        for d in definitions
    }


def h5md_unit(definition, unit_system=None):
    """The units of the values LAMMPS writes for a property, or None if it is
    dimensionless."""
    if definition["kind"] == "explicit":
        return definition["units"]
    if definition["units"] == "dimensionless":
        return None
    return lammps_step.lammps_units(definition["units"], unit_system=unit_system)


# ----------------------------------------------------------------------
# The treatments
# ----------------------------------------------------------------------
def prefactor_value(definition, volume=None, temperature=None):
    """The prefactor of a Green–Kubo integral as a quantity.

    Parameters
    ----------
    definition : dict
        The extra property.
    volume : pint.Quantity
        The (mean) volume of the cell.
    temperature : pint.Quantity
        The mean temperature.
    """
    prefactor = definition["prefactor"]
    k_B = Q_("k_B")
    if prefactor in PREFACTORS:
        if temperature is None:
            raise ExtrasError("The temperature is needed for the prefactor.")
        if volume is None:
            raise ExtrasError(
                "The prefactor needs the volume, which a non-periodic system does "
                "not have."
            )
        if prefactor == "V/kT":
            return volume / (k_B * temperature)
        if prefactor == "V/kT^2":
            return volume / (k_B * temperature**2)
        return 1 / (volume * k_B * temperature)
    return Q_(prefactor)


def result_units(quantity):
    """Readable units for a result: the first of :data:`RESULT_UNITS` with its
    dimensions (the viscosity in mPa*s, ...), else SI base units scaled with a
    prefix."""
    for units in RESULT_UNITS:
        if quantity.dimensionality == Q_(1.0, units).dimensionality:
            return Q_(1.0, units).units
    if quantity.dimensionless:
        return Q_(1.0, "").units
    return quantity.to_base_units().to_compact().units


def autocorrelation(values, dt, m=None):
    """The normalized autocorrelation function of a series, its standard error
    and the integrated correlation time.

    Parameters
    ----------
    values : numpy.ndarray(n)
        The series, sampled every ``dt``.
    dt : pint.Quantity
        The sampling interval.
    m : int
        The lags to keep; a tenth of the series by default.

    Returns
    -------
    dict
        ``acf``, ``stderr`` (lists), ``tau`` (the integral of the ACF up to its
        first zero, a quantity) and ``dt``.
    """
    from seamm_md_analysis.correlate import acf_err, autocorrelation_functions

    x = np.asarray(values, dtype=float)
    if x.std() <= 1.0e-12 * max(1.0, abs(x.mean())):
        raise ExtrasError("The series is constant, so it has no autocorrelation.")
    if m is None:
        m = max(2, x.size // 10)
    acf, _ = autocorrelation_functions(x.reshape(1, -1), m=m)
    acf = acf[0]
    acf = acf / acf[0]
    zero = np.argmax(acf <= 0.0) if np.any(acf <= 0.0) else acf.size
    trapezoid = getattr(np, "trapezoid", None) or np.trapz  # numpy < 2
    tau = (trapezoid(acf[:zero]) if zero > 1 else 0.0) * dt
    return {
        "acf": acf.tolist(),
        "stderr": acf_err(acf).tolist(),
        "tau": tau,
        "dt": dt,
    }


def green_kubo(values, dt, definition, unit="", volume=None, temperature=None, m=None):
    """The Green–Kubo integral of a series: prefactor × ∫⟨δA(0) δA(t)⟩ dt.

    ``values`` are in ``unit`` (empty for a dimensionless property), sampled every
    ``dt``.

    The running integral is fitted with
    :func:`seamm_md_analysis.transport.fit_green_kubo_integral`. Its error is
    the spread of the values fitted to four blocks of the series when they are
    long enough to reach the plateau, else the fit's own (which treats
    correlated points as independent, so is too small).

    The fit smooths the running integral over 0.4 ps, as for the viscosity and
    the thermal conductivity, so the correlations should last at least about a
    picosecond.

    Returns
    -------
    dict
        ``value`` and ``stderr`` (quantities), ``error_from`` ("blocks" or
        "fit"), ``tau`` (the fit's time constants, ps), ``ts`` (ps),
        ``integral`` and ``fit`` (lists, in the value's units) and ``units``.
    """
    from seamm_md_analysis.correlate import autocorrelation_functions
    from seamm_md_analysis.transport import fit_green_kubo_integral

    x = np.asarray(values, dtype=float)
    if m is None:
        m = max(2, x.size // 10)
    factor = prefactor_value(definition, volume, temperature)
    factor = factor * Q_(1.0, unit or "") ** 2 * dt
    units = result_units(factor)
    scale = factor.m_as(units)

    def integral_of(series, lags):
        _, integral = autocorrelation_functions(series.reshape(1, -1), m=lags)
        return integral[0] * scale

    integral = integral_of(x, m)
    ts = np.arange(m) * dt.m_as("ps")
    fit = fit_green_kubo_integral(integral, ts.tolist())

    stderr = fit.stderr
    error_from = "fit"
    # Blocks are used when each still reaches the plateau: 20 of the slowest
    # time constants within a tenth of the block
    blocks = np.array_split(x, 4)
    m_block = min(m, min(b.size for b in blocks) // 10)
    if m_block * dt.m_as("ps") >= 20 * max(fit.tau):
        values = []
        for block in blocks:
            try:
                values.append(
                    fit_green_kubo_integral(
                        integral_of(block, m_block), ts[:m_block].tolist()
                    ).value
                )
            except Exception:
                pass
        values = [v for v in values if math.isfinite(v)]
        if len(values) >= 2:
            stderr = float(np.std(values, ddof=1) / math.sqrt(len(values)))
            error_from = "blocks"
    return {
        "value": Q_(fit.value, units),
        "stderr": Q_(stderr, units),
        "error_from": error_from,
        "tau": list(fit.tau),
        "ts": ts.tolist(),
        "integral": integral.tolist(),
        "fit_ts": list(fit.xs),
        "fit": list(fit.ys),
        "units": f"{units:~P}",
    }
