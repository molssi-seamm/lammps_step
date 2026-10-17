# -*- coding: utf-8 -*-

"""Convert the output of an MD sub-step to SEAMM's H5MD file, md.h5md.

LAMMPS writes its usual text files during the run. When the run ends, the NVE, NVT
and NPT sub-steps (and Heat Flux, an NVE run) stream them into ``md.h5md`` in the
sub-step's directory, with the layout described in seamm_md_analysis ("The H5MD
file"), and add a record of the run to the workspace variable ``_md_runs``. The
native files are kept.
"""

import logging
from pathlib import Path
import platform

import numpy as np

import lammps_step
from lammps_step.lammps_units import unit_systems

logger = logging.getLogger(__name__)

#: The columns of the state files and the observables they become, with the
#: LAMMPS quantity that gives their units.
state_observables = {
    "T": ("temperature", "temperature"),
    "P": ("pressure", "pressure"),
    "V": ("volume", "volume"),
    "density": ("density", "density"),
    "Etot": ("total_energy", "energy"),
    "Eke": ("kinetic_energy", "energy"),
    "Epe": ("potential_energy", "energy"),
    "Emol": ("molecular_energy", "energy"),
    "Epair": ("pair_energy", "energy"),
    "N_hbond": ("hbond_count", None),
    "E_hbond": ("hbond_energy", "energy"),
}

#: The observables that need the volume of a periodic cell
periodic_only = ("pressure", "pressure_tensor", "volume", "density", "heat_flux")


def h5md_units(unit_system=None):
    """The H5MD unit strings for a LAMMPS unit system.

    Parameters
    ----------
    unit_system : str = None
        The LAMMPS unit system, by default the current one.

    Returns
    -------
    {str: str}
        The unit for each LAMMPS quantity ("time", "distance", "energy", ...).
    """
    from seamm_md_analysis.units import to_h5md

    if unit_system is None:
        unit_system = lammps_step.get_lammps_unit_system()
    result = {}
    for quantity, unit in unit_systems[unit_system].items():
        try:
            result[quantity] = to_h5md(unit)
        except Exception:
            logger.debug(f"No H5MD unit for {quantity} = {unit}")
    return result


def _cell_matrix(configuration):
    """The cell vectors (rows) in LAMMPS's orientation."""
    lx, ly, lz, xy, xz, yz = lammps_step.LAMMPS.cell_to_box(
        *configuration.cell.parameters
    )
    return np.array([[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]])


def _cell_from_parameters(a, b, c, alpha=90.0, beta=90.0, gamma=90.0):
    lx, ly, lz, xy, xz, yz = lammps_step.LAMMPS.cell_to_box(a, b, c, alpha, beta, gamma)
    return np.array([[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]])


class MDH5MDConverter:
    """Write md.h5md for one MD sub-step.

    Parameters
    ----------
    node : lammps_step.NVE
        The sub-step (NVE, NVT, NPT or Heat Flux).
    P : dict
        Its parameters, as values.
    configuration : molsystem._Configuration
        The configuration at the start of the run.
    directory : str or pathlib.Path = None
        The sub-step's directory, by default node.directory.
    wall_time : float = None
        The wall time of the run in seconds (LAMMPS's loop time).
    """

    def __init__(self, node, P, configuration, directory=None, wall_time=None):
        self.node = node
        self.wall_time = wall_time
        self.P = P
        self.configuration = configuration
        self.directory = Path(node.directory if directory is None else directory)
        self.units = h5md_units()
        self.time_unit = self.units["time"]
        self.periodic = configuration.periodicity == 3
        self.timestep, self.timestep_q = node.timestep(P["timestep"])
        time = lammps_step.to_lammps_units(P["time"], quantity="time")
        self.nsteps = round(time / self.timestep)
        self.observables = {}  # name -> sampling interval in steps
        self.particles = {}  # "group/name" -> sampling interval in steps
        self.warnings = []

    @property
    def ensemble(self):
        name = type(self.node).__name__
        return "NVE" if name == "HeatFlux" else name

    @property
    def path(self):
        return self.directory / "md.h5md"

    def write(self):
        """Write md.h5md from the native files.

        Returns
        -------
        dict
            The record of the run for ``_md_runs``.
        """
        from seamm_md_analysis import H5MDWriter

        with H5MDWriter(
            self.path,
            creator="SEAMM lammps_step",
            version=lammps_step.__version__,
        ) as writer:
            self._write_particles(writer)
            self._write_state(writer)
            self._write_vector_series(writer)
            self._write_flux_and_stress(writer)
            writer.parameters(self._parameters())
            writer.provenance(self._provenance(), segments=[self._segment()])
        return self.record()

    # Particles

    def _write_particles(self, writer):
        configuration = self.configuration
        eex = self.node.parent.eex
        charges = None
        if isinstance(eex, dict) and "charges" in eex:
            charges = eex["charges"]
        cell = _cell_matrix(configuration) if self.periodic else None
        writer.particles(
            ids=configuration.atoms.ids,
            species=configuration.atoms.atomic_numbers,
            masses=configuration.atoms.atomic_masses,
            charges=charges,
            cell=cell,
            periodic=self.periodic,
            length_unit=self.units["distance"],
        )

        d = self.directory
        positions = d / "atomic_positions.dump_trj"
        velocities = d / "atomic_velocities.dump_trj"
        trajectory = d / "trajectory.dump_trj.gz"

        # The atoms: the explicitly requested positions and velocities if
        # sampled, otherwise those in the trajectory
        if positions.exists():
            self._dump_to_particles(
                writer, positions, "all", {"position": ("xu", "yu", "zu")}
            )
        if velocities.exists():
            self._dump_to_particles(
                writer, velocities, "all", {"velocity": ("vx", "vy", "vz")}
            )
        if trajectory.exists():
            fields = {}
            if not positions.exists():
                fields["position"] = ("xu", "yu", "zu")
            if not velocities.exists():
                fields["velocity"] = ("vx", "vy", "vz")
            fields["force"] = ("fx", "fy", "fz")
            group = "all"
            if positions.exists():
                # Keep the trajectory as its own group; MDAnalysis reads "all",
                # the first group
                group = "trajectory"
                writer.particles(
                    group=group,
                    ids=configuration.atoms.ids,
                    cell=cell,
                    periodic=self.periodic,
                    length_unit=self.units["distance"],
                )
                fields["position"] = ("xu", "yu", "zu")
            self._dump_to_particles(writer, trajectory, group, fields)

        # The molecules
        if (d / "com_positions.trj").exists() or (d / "com_velocities.trj").exists():
            n_molecules = len(configuration.find_molecules(as_indices=True))
            writer.particles(
                group="molecules",
                ids=np.arange(n_molecules),
                cell=cell,
                periodic=self.periodic,
                length_unit=self.units["distance"],
            )

    def _dump_to_particles(self, writer, path, group, fields):
        """Stream a LAMMPS dump file into per-particle elements.

        Parameters
        ----------
        fields : {str: (str, str, str)}
            The H5MD names and the dump columns for them. A field whose columns
            are not in the dump is skipped.
        """
        from seamm_md_analysis import iter_dump_frames

        unit_for = {"position": "distance", "velocity": "velocity", "force": "force"}
        appenders = None
        box = None
        steps = []
        for frame in iter_dump_frames(path):
            columns = frame["columns"]
            if appenders is None:
                appenders = {}
                n = len(columns["id"])
                for name, cols in fields.items():
                    if all(c in columns for c in cols):
                        attrs = {"seamm_unwrapped": 1} if name == "position" else {}
                        appenders[name] = writer.particle_appender(
                            name,
                            n,
                            self.units[unit_for[name]],
                            self.time_unit,
                            group=group,
                            attrs=attrs,
                        )
                if len(appenders) == 0:
                    return
                if self.periodic and "position" in appenders:
                    box = writer.box_appender(
                        group=group,
                        length_unit=self.units["distance"],
                        time_unit=self.time_unit,
                    )
            step = frame["step"]
            time = frame.get("time", step * self.timestep)
            for name, appender in appenders.items():
                cols = fields[name]
                appender.append(step, time, np.column_stack([columns[c] for c in cols]))
            if box is not None:
                box.append(step, time, frame["cell"])
            steps.append(step)
        if appenders:
            interval = int(np.median(np.diff(steps))) if len(steps) > 1 else None
            for name in appenders:
                self.particles[f"{group}/{name}"] = interval

    # Time series

    def _read_trj(self, path):
        """A time series written by the run, without its step-0 sample.

        LAMMPS samples the setup at step 0, before any dynamics. With constraints
        (SHAKE) its virial is not valid there -- pressures and stresses of 1e307
        have been seen -- so the series start at the first step of the dynamics.
        """
        from seamm_md_analysis import read_trajectory

        metadata, columns = read_trajectory(path)
        steps = columns["tstep"].astype(int)
        keep = steps != 0
        if not keep.all():
            columns = {
                key: (value[keep] if len(value) == len(steps) else value)
                for key, value in columns.items()
            }
            steps = steps[keep]
        if "t" in columns:
            times = columns["t"]
        else:
            times = steps * self.timestep
        return steps, times, columns

    def _observable(self, writer, name, steps, times, values, unit, **kwargs):
        if not self.periodic and name in periodic_only:
            return
        writer.observable(name, steps, times, values, unit, self.time_unit, **kwargs)
        self.observables[name] = (
            int(np.median(np.diff(steps))) if len(steps) > 1 else None
        )

    def _write_state(self, writer):
        """The thermodynamic state, from <ensemble>_state.trj or trajectory.trj."""
        paths = sorted(self.directory.glob("*_state.trj"))
        if len(paths) == 0:
            path = self.directory / "trajectory.trj"
            if not path.exists():
                return
        else:
            path = paths[0]
        steps, times, columns = self._read_trj(path)

        for column, (name, quantity) in state_observables.items():
            if column in columns:
                unit = "" if quantity is None else self.units[quantity]
                self._observable(writer, name, steps, times, columns[column], unit)

        # The user's extra properties (lammps_step.extras)
        if hasattr(self.node, "extra_definitions"):
            from seamm_md_analysis.units import to_h5md

            for d in self.node.extra_definitions():
                if d["column"] not in columns:
                    continue
                if d["kind"] == "explicit":
                    unit = to_h5md(d["units"])
                elif d["units"] == "dimensionless":
                    unit = ""
                else:
                    unit = self.units[d["units"]]
                self._observable(
                    writer,
                    f"user/{d['name']}",
                    steps,
                    times,
                    columns[d["column"]],
                    unit,
                    attrs={
                        "seamm_expression": d["expression"],
                        "seamm_treatment": d["treatment"],
                    },
                )

        # The pressure tensor sampled on its own, if it was, replaces this one
        stress = ("Sxx", "Syy", "Szz", "Syz", "Sxz", "Sxy")
        sampled = (self.directory / "pressure_tensor.trj").exists()
        if not sampled and all(s in columns for s in stress):
            xx, yy, zz, yz, xz, xy = (columns[s] for s in stress)
            # The stress is minus the pressure
            tensor = -np.stack(
                [
                    np.stack([xx, xy, xz], axis=-1),
                    np.stack([xy, yy, yz], axis=-1),
                    np.stack([xz, yz, zz], axis=-1),
                ],
                axis=-2,
            )
            self._observable(
                writer, "pressure_tensor", steps, times, tensor, self.units["pressure"]
            )

        # The cell of an NPT run, if there are no frames to take it from
        all_group = writer.file["particles/all"]
        if (
            self.periodic
            and "position" not in all_group
            and all(c in columns for c in ("a", "b", "c"))
        ):
            angles = [
                columns.get(x, np.full(len(steps), 90.0))
                for x in ("alpha", "beta", "gamma")
            ]
            box = writer.box_appender(
                length_unit=self.units["distance"], time_unit=self.time_unit
            )
            for i, step in enumerate(steps):
                cell = _cell_from_parameters(
                    columns["a"][i],
                    columns["b"][i],
                    columns["c"][i],
                    angles[0][i],
                    angles[1][i],
                    angles[2][i],
                )
                box.append(step, times[i], cell)

    def _write_vector_series(self, writer):
        """The centres of mass of the molecules."""
        from seamm_md_analysis import iter_vector_frames

        for filename, name, quantity in (
            ("com_positions.trj", "position", "distance"),
            ("com_velocities.trj", "velocity", "velocity"),
        ):
            path = self.directory / filename
            if not path.exists():
                continue
            frames = iter_vector_frames(path)
            next(frames)
            appender = None
            steps = []
            for step, values in frames:
                if appender is None:
                    attrs = {"seamm_unwrapped": 1} if name == "position" else {}
                    appender = writer.particle_appender(
                        name,
                        values.shape[0],
                        self.units[quantity],
                        self.time_unit,
                        group="molecules",
                        attrs=attrs,
                    )
                appender.append(step, step * self.timestep, values)
                steps.append(step)
            if appender is not None:
                self.particles[f"molecules/{name}"] = (
                    int(np.median(np.diff(steps))) if len(steps) > 1 else None
                )

    def _write_flux_and_stress(self, writer):
        path = self.directory / "heat_flux.trj"
        if path.exists():
            steps, _, columns = self._read_trj(path)
            times = steps * self.timestep
            flux = np.column_stack([columns["Jx"], columns["Jy"], columns["Jz"]])
            # Converted to W/m^2 in the LAMMPS input
            self._observable(writer, "heat_flux", steps, times, flux, "W m-2")

        path = self.directory / "pressure_tensor.trj"
        if path.exists():
            steps, _, columns = self._read_trj(path)
            times = steps * self.timestep
            xx, yy, zz, xy, xz, yz = (
                columns[c] for c in ("Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz")
            )
            tensor = np.stack(
                [
                    np.stack([xx, xy, xz], axis=-1),
                    np.stack([xy, yy, yz], axis=-1),
                    np.stack([xz, yz, zz], axis=-1),
                ],
                axis=-2,
            )
            self._observable(
                writer, "pressure_tensor", steps, times, tensor, self.units["pressure"]
            )

    # Parameters, provenance and the record

    def _quantity(self, key):
        value = self.P.get(key, None)
        if value is None:
            return None, None
        if hasattr(value, "magnitude"):
            return float(value.magnitude), f"{value.units:~P}"
        return value, ""

    def _parameters(self):
        result = {
            "ensemble": self.ensemble,
            "timestep": float(self.timestep),
            "timestep_unit": self.time_unit,
            "n_steps": int(self.nsteps),
            "periodic": self.periodic,
        }
        for key, name in (
            ("T0", "temperature"),
            ("T1", "final_temperature"),
            ("Tdamp", "thermostat_damping"),
            # NPT's pressures are 'Pinitial' and 'Pfinal'
            ("Pinitial", "pressure"),
            ("Pfinal", "final_pressure"),
            ("P0", "pressure"),
            ("P1", "final_pressure"),
            ("Pdamp", "barostat_damping"),
        ):
            value, unit = self._quantity(key)
            if value is not None:
                result[name] = value
                result[f"{name}_unit"] = unit
        for key, name in (("thermostat", "thermostat"), ("barostat", "barostat")):
            if key in self.P:
                result[name] = str(self.P[key])
        # The long-range electrostatics, for analyses whose formulas depend on it
        # (the dielectric constant assumes Ewald or PPPM with conducting
        # boundaries)
        kspace = self._kspace_style()
        if kspace is not None:
            result["kspace_style"] = kspace
        # The styles of the fixes, for analyses that need to know whether anything
        # directional (walls, fields, deformation, flow) acted on the system
        fixes = self._fix_styles()
        if fixes is not None:
            result["fixes"] = fixes
        step = {}
        for key, value in self.P.items():
            step[key] = f"{value:~P}" if hasattr(value, "magnitude") else str(value)
        result["step"] = step
        return result

    def _kspace_style(self):
        """The kspace_style of the run, from the LAMMPS input ("none" if it has
        none), or None if the input cannot be found."""
        for path in (self.directory / "input.dat", self.directory.parent / "input.dat"):
            if path.exists():
                style = "none"
                for line in path.read_text().splitlines():
                    words = line.split()
                    if len(words) > 1 and words[0] == "kspace_style":
                        style = " ".join(words[1:])
                return style
        return None

    def _fix_styles(self):
        """The styles of the fixes in the LAMMPS input, e.g. ["nvt", "shake"], or
        None if the input cannot be found. The input may hold several sub-steps, so
        this is every fix in it: a cautious answer for whether one was directional."""
        for path in (self.directory / "input.dat", self.directory.parent / "input.dat"):
            if path.exists():
                styles = set()
                for line in path.read_text().splitlines():
                    words = line.split("#", 1)[0].split()
                    if len(words) > 3 and words[0] == "fix":
                        styles.add(words[3])
                return sorted(styles)
        return None

    def _provenance(self):
        node = self.node
        result = {
            "code": "LAMMPS",
            "lammps_step": lammps_step.__version__,
            "step_id": ".".join(str(i) for i in node._id),
            "step_title": node.title,
            "periodic": self.periodic,
            "unit_system": lammps_step.get_lammps_unit_system(),
        }
        try:
            import seamm_md_analysis

            result["seamm_md_analysis"] = seamm_md_analysis.__version__
        except Exception:
            pass
        if node.variable_exists("_job_id"):
            result["job_id"] = str(node.get_variable("_job_id"))
        configuration = self.configuration
        result["system"] = configuration.system.name
        result["configuration"] = configuration.name
        result["configuration_id"] = configuration.id
        model = node.model
        if model is not None:
            result["model"] = str(model)
        return result

    def _segment(self):
        try:
            wall = float(self.wall_time)
        except (TypeError, ValueError):
            wall = -1.0
        return {
            "first_step": 0,
            "last_step": int(self.nsteps),
            "wall_time": wall,
            "node": platform.node(),
        }

    def _interval(self, nsteps):
        if nsteps is None:
            return None
        return f"{nsteps * self.timestep_q.magnitude:g} {self.timestep_q.units:~P}"

    def record(self):
        """The record of this run for the workspace variable ``_md_runs``."""
        node = self.node
        configuration = self.configuration
        return {
            "path": str(self.path.resolve()),
            "step": f"{'.'.join(str(i) for i in node._id)} {node.title}",
            "code": "LAMMPS",
            "ensemble": self.ensemble,
            "periodic": self.periodic,
            "configuration": f"{configuration.system.id}/{configuration.id}",
            "particles": {k: self._interval(v) for k, v in self.particles.items()},
            "observables": {k: self._interval(v) for k, v in self.observables.items()},
            "time": [
                f"0 {self.timestep_q.units:~P}",
                self._interval(self.nsteps),
            ],
        }


def write_md_h5md(node, P, configuration, directory=None, wall_time=None):
    """Write md.h5md for an MD sub-step and record it in ``_md_runs``.

    A failure is reported, not raised: the native files are still there.

    Parameters
    ----------
    node : lammps_step.NVE
        The sub-step.
    P : dict
        Its parameters, as values.
    configuration : molsystem._Configuration
        The configuration at the start of the run.
    directory : str or pathlib.Path = None
        The sub-step's directory, by default node.directory.
    wall_time : float = None
        The wall time of the run in seconds.

    Returns
    -------
    str
        Text for the output.
    """
    try:
        converter = MDH5MDConverter(
            node, P, configuration, directory=directory, wall_time=wall_time
        )
        record = converter.write()
    except Exception as e:
        logger.warning(f"Could not write md.h5md: {e}", exc_info=True)
        return f"Warning: could not write the H5MD file md.h5md: {e}"

    runs = []
    if node.variable_exists("_md_runs"):
        runs = list(node.get_variable("_md_runs"))
    runs.append(record)
    node.set_variable("_md_runs", runs)

    text = "Wrote the MD data to md.h5md"
    n = len(record["observables"])
    groups = sorted(record["particles"])
    if n > 0 or groups:
        parts = []
        if n > 0:
            parts.append(f"{n} time series")
        if groups:
            parts.append("the " + ", ".join(groups))
        text += " (" + " and ".join(parts) + ")"
    return text + "."
