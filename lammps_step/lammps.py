# -*- coding: utf-8 -*-

"""A node or step for LAMMPS in a flowchart"""

import configparser
from contextlib import contextmanager
import copy
import csv
from datetime import datetime, timezone
import importlib
import json
import logging
from math import sqrt, exp, degrees, radians, cos, acos
from pathlib import Path
import os
import os.path
import pkg_resources
import platform
import pprint
import shutil
import string
import sys
import time
import traceback
import warnings

import bibtexparser
from cpuinfo import get_cpu_info
import numpy as np
import pandas
import statsmodels.tsa.stattools as stattools

import lammps_step
import molsystem
import seamm
import seamm_exec
from seamm_ff_util import tabulate_angle
import seamm_util
import seamm_util.printing as printing
from seamm_util import CompactJSONEncoder, Configuration, units_class
from seamm_util.printing import FormattedText as __

# from pymbar import timeseries

logger = logging.getLogger("lammps")
job = printing.getPrinter()
printer = printing.getPrinter("lammps")


# Temporarily used here to stop pymbar's annoying warning.
@contextmanager
def logging_disabled(highest_level=logging.CRITICAL):
    """
    A context manager that will prevent any logging messages
    triggered during the body from being processed.

    :param highest_level: the maximum logging level in use.
      This would only need to be changed if a custom level greater than CRITICAL
      is defined.

    From Simon Weber https://gist.github.com/simon-weber/7853144
    """
    # two kind-of hacks here:
    #    * can't get the highest logging level in effect => delegate to the user
    #    * can't get the current module-level override => use an undocumented
    #       (but non-private!) interface

    previous_level = logging.root.manager.disable

    logging.disable(highest_level)

    try:
        yield
    finally:
        logging.disable(previous_level)


with logging_disabled(highest_level=logging.WARNING):
    from pymbar import timeseries


# Add LAMMPS's properties to the standard properties
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
molsystem.add_properties_from_file(csv_file)

bond_style = {
    "quadratic_bond": "harmonic",
    "quartic_bond": "class2",
    "fene": "fene",
    "morse": "morse",
}

angle_style = {
    "quadratic_angle": "harmonic",
    "quartic_angle": "class2",
    "cosine": "cosine",
    "cosine/squared": "cosine/squared",
    "simple_fourier_angle": "fourier/simple",
    "tabulated_angle": "table",
}

dihedral_style = {
    "torsion_1": "harmonic",
    "torsion_3": "class2",
    "torsion_opls": "opls",
    "torsion_charmm": "charmm",
}

improper_style = {
    "wilson_out_of_plane": "class2",
    "improper_opls": "cvff",
    "dreiding_out_of_plane": "umbrella",
}


class LAMMPS(seamm.Node):
    display_units = {
        "T": "K",
        "P": "atm",
        "t": "fs",
        "V": "Å^3",
        "density": "g/mL",
        "a": "Å",
        "b": "Å",
        "c": "Å",
        "Etot": "kcal/mol",
        "DfH0_reax": "kcal/mol",
        "Eke": "kcal/mol",
        "Epe": "kcal/mol",
        "Emol": "kcal/mol",
        "Epair": "kcal/mol",
        "Jx": "W/m^2",
        "Jy": "W/m^2",
        "Jz": "W/m^2",
        "Kappa_x": "W/K/m",
        "Kappa_y": "W/K/m",
        "Kappa_z": "W/K/m",
        "Kappa": "W/K/m",
        "u": "kcal/mol",
        "N_hbond": "",
        "E_hbond": "kcal/mol",
    }
    display_title = {
        "T": "Temperature",
        "P": "Pressure",
        "t": "Time",
        "V": "Volume",
        "density": "Density",
        "a": "a lattice parameter",
        "b": "b lattice parameter",
        "c": "c lattice parameter",
        "Etot": "Total Energy",
        "DfH0_reax": "Reax \N{GREEK CAPITAL LETTER DELTA}fH\N{SUPERSCRIPT ZERO}",
        "Eke": "Kinetic Energy",
        "Epe": "Potential Energy",
        "Emol": "Molecular Energy, Valence Terms",
        "Epair": "Pair (Nonbond) Energy",
        "Jx": "heat flux in x",
        "Jy": "heat flux in y",
        "Jz": "heat flux in z",
        "Kappa_x": "thermal conductivity in x",
        "Kappa_y": "thermal conductivity in y",
        "Kappa_z": "thermal conductivity in z",
        "Kappa": "thermal conductivity",
        "u": "atom PE",
        "N_hbond": "Number of H-bonds",
        "E_hbond": "Energy of H-bonds",
    }

    def __init__(
        self, flowchart=None, namespace="org.molssi.seamm.lammps", extension=None
    ):
        """Setup the main LAMMPS step

        Keyword arguments:
        """
        logger.debug("Creating LAMMPS {}".format(self))

        # The subflowchart
        self.subflowchart = seamm.Flowchart(
            parent=self, name="LAMMPS", namespace=namespace
        )
        self._initialization_node = None
        self._trajectory = []
        self._data = {}
        self._have_dreiding_hbonds = False
        self._atomic_energy_sum = 0.0

        self._results = {}  # Sotrage for computational and timing results

        super().__init__(
            flowchart=flowchart, title="LAMMPS", extension=extension, logger=logger
        )

        # Set up the timing information
        self._timing_data = []
        self._timing_path = Path("~/.seamm.d/timing/lammps.csv").expanduser()
        self._timing_header = [
            "node",  # 0
            "cpu",  # 1
            "cpu_version",  # 2
            "cpu_count",  # 3
            "cpu_speed",  # 4
            "date",  # 5
            "H_SMILES",  # 6
            "ISOMERIC_SMILES",  # 7
            "formula",  # 8
            "net_charge",  # 9
            "spin_multiplicity",  # 10
            "keywords",  # 11
            "nproc",  # 12
            "time",  # 13
        ]
        try:
            self._timing_path.parent.mkdir(parents=True, exist_ok=True)

            self._timing_data = 14 * [""]
            self._timing_data[0] = platform.node()
            tmp = get_cpu_info()
            if "arch" in tmp:
                self._timing_data[1] = tmp["arch"]
            if "cpuinfo_version_string" in tmp:
                self._timing_data[2] = tmp["cpuinfo_version_string"]
            if "count" in tmp:
                self._timing_data[3] = str(tmp["count"])
            if "hz_advertized_friendly" in tmp:
                self._timing_data[4] = tmp["hz_advertized_friendly"]

            if not self._timing_path.exists():
                with self._timing_path.open("w", newline="") as fd:
                    writer = csv.writer(fd)
                    writer.writerow(self._timing_header)
        except Exception:
            self._timing_data = None

    @property
    def version(self):
        """The semantic version of this module."""
        return lammps_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return lammps_step.__git_revision__

    @property
    def have_dreiding_hbonds(self):
        """Whether the system has Dreiding hydrogen bonds."""
        return self._have_dreiding_hbonds

    @have_dreiding_hbonds.setter
    def have_dreiding_hbonds(self, value):
        self._have_dreiding_hbonds = value

    @property
    def results(self):
        """The storage for results."""
        return self._results

    @staticmethod
    def box_to_cell(lx, ly, lz, xy, xz, yz):
        """Convert the LAMMPS box definition to cell parameters."""
        if xy == 0 and xz == 0 and yz == 0:
            a = lx
            b = ly
            c = lz
            alpha = 0.0
            beta = 0.0
            gamma = 0.0
        else:
            a = lx
            b = sqrt(ly**2 + xy**2)
            c = sqrt(lz**2 + xz**2 + yz**2)
            alpha = degrees(acos((xy * xz + lx * yz) / (b * c)))
            beta = degrees(acos(xz / c))
            gamma = degrees(acos(xy / b))

        return (a, b, c, alpha, beta, gamma)

    @staticmethod
    def cell_to_box(a, b, c, alpha, beta, gamma):
        """Convert cell parameters to the LAMMPS box."""
        if alpha == 90 and beta == 90 and gamma == 90:
            lx = a
            ly = b
            lz = c
            xy = xz = yz = 0.0
        else:
            lx = a
            xy = b * cos(radians(gamma))
            xz = c * cos(radians(beta))
            ly = sqrt(b**2 - xy**2)
            yz = (b * c * cos(radians(alpha)) - xy * xz) / ly
            lz = sqrt(c**2 - xz**2 - yz**2)

        return (lx, ly, lz, xy, xz, yz)

    def create_parser(self):
        """Setup the command-line / config file parser"""
        parser_name = self.step_type
        parser = seamm_util.getParser()

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        result = super().create_parser(name=parser_name)

        if parser_exists:
            return result

        # LAMMPS specific options
        parser.add_argument(
            parser_name,
            "--ncores",
            default="available",
            help=(
                "The maximum number of cores to use for LAMMPS. "
                "Default: all available cores."
            ),
        )
        parser.add_argument(
            parser_name,
            "--atoms-per-core",
            type=int,
            default="1000",
            help="the optimal number of atoms per core for LAMMPS",
        )
        parser.add_argument(
            parser_name,
            "--html",
            action="store_true",
            help="whether to write out html files for graphs, etc.",
        )
        if False:
            parser.add_argument(
                parser_name,
                "--modules",
                nargs="*",
                default=None,
                help="the environment modules to load for LAMMPS",
            )
            parser.add_argument(
                parser_name,
                "--gpu-modules",
                nargs="*",
                default=None,
                help="the environment modules to load for the GPU version of LAMMPS",
            )
            parser.add_argument(
                parser_name,
                "--lammps-path",
                default=None,
                help="the path to the LAMMPS executables",
            )
            parser.add_argument(
                parser_name,
                "--lammps-serial",
                default="lmp_serial",
                help="the serial version of LAMMPS",
            )
            parser.add_argument(
                parser_name,
                "--lammps-mpi",
                default="lmp_mpi",
                help="the mpi version of LAMMPS",
            )
            parser.add_argument(
                parser_name,
                "--cmd-args",
                default="",
                help="the command-line arguments for LAMMPS, e.g. '-k on'",
            )
            parser.add_argument(
                parser_name,
                "--gpu-cmd-args",
                default="",
                help=(
                    "the command-line arguments for GPU version of LAMMPS, e.g. '-k on'"
                ),
            )
            parser.add_argument(
                parser_name, "--mpiexec", default="mpiexec", help="the mpi executable"
            )

        return result

    def set_id(self, node_id):
        """Set the id for node to a given tuple"""
        self._id = node_id

        # and set our subnodes
        self.subflowchart.set_ids(self._id)

        return self.next()

    def description_text(self, P=None):
        """Return a short description of this step.

        Return a nicely formatted string describing what this step will
        do.

        Keyword arguments:
            P: a dictionary of parameter values, which may be variables
                or final values. If None, then the parameters values will
                be used as is.
        """

        self.subflowchart.root_directory = self.flowchart.root_directory

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        text = self.header + "\n\n"

        while node is not None:
            try:
                text += __(node.description_text(), indent=4 * " ").__str__()
            except Exception as e:
                print(
                    "Error describing LAMMPS flowchart: {} in {}".format(
                        str(e), str(node)
                    )
                )
                self.logger.critical(
                    "Error describing LAMMPS flowchart: {} in {}".format(
                        str(e), str(node)
                    )
                )
                raise
            except:  # noqa: E722
                print(
                    "Unexpected error describing LAMMPS flowchart: {} in {}".format(
                        sys.exc_info()[0], str(node)
                    )
                )
                self.logger.critical(
                    "Unexpected error describing LAMMPS flowchart: {} in {}".format(
                        sys.exc_info()[0], str(node)
                    )
                )
                raise
            text += "\n"
            node = node.next()

        return text

    def run(self):
        """Run a LAMMPS simulation"""
        # Set the model
        try:
            ff = self.get_variable("_forcefield")
            if ff == "OpenKIM":
                self.model = "OpenKIM/" + self.get_variable("_OpenKIM_Potential")
            else:
                self.model = ff.current_forcefield
        except Exception:
            self.model = None

        system_db = self.get_variable("_system_db")
        configuration = system_db.system.configuration

        n_atoms = configuration.n_atoms
        if n_atoms == 0:
            self.logger.error("LAMMPS run(): there is no structure!")
            raise RuntimeError("LAMMPS run(): there is no structure!")

        # Initialize storage
        self._results = {}

        next_node = super().run(printer)

        # Get the options
        o = self.options
        global_options = self.global_options

        # Whether to run parallel and if so, how many mpi processes
        if global_options["parallelism"] in ("any", "mpi"):
            np = n_atoms // o["atoms_per_core"] + 1
            if o["ncores"] != "available":
                np = min(np, int(o["ncores"]))
            if global_options["ncores"] != "available":
                np = min(np, int(global_options["ncores"]))
        else:
            np = 1

        # Print headers and get to work
        printer.important(self.header)

        self.subflowchart.root_directory = self.flowchart.root_directory

        files = {}

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        extras = {}

        history_nodes = []

        # Create overall directory for the lammps step

        os.makedirs(self.directory, exist_ok=True)

        files = {}

        control = []
        ff = self.get_variable("_forcefield")
        if ff == "OpenKIM":
            potential = self.get_variable("_OpenKIM_Potential")
            control.append(["forcefield", "OpenKIM " + potential])
        else:
            control.append(["forcefield", ff.current_forcefield])
        while node is not None:
            P = node.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )

            # For the timing data, get the parameters used
            tmp = [node.title]
            for key, value in P.items():
                if isinstance(value, units_class):
                    tmp.append((key, f"{value:~P}"))
                else:
                    tmp.append((key, value))
            control.append(tmp)

            if isinstance(node, lammps_step.Initialization):
                initialization_header, eex = self._get_node_input(
                    node=node, extras={"read_data": True}
                )
                (
                    structure_data,
                    pair_table,
                    bond_table,
                    angle_table,
                    dihedral_table,
                ) = self.structure_data(eex)
                files["structure.dat"] = structure_data
                if "forcefield" in eex:
                    files["forcefield.dat"] = eex["forcefield"]

                if bond_table != "":
                    files["tabulated_bonds.dat"] = bond_table
                if angle_table != "":
                    files["tabulated_angles.dat"] = angle_table
                if dihedral_table != "":
                    files["tabulated_dihedrals.dat"] = dihedral_table

                self.logger.debug("structure.dat:\n" + files["structure.dat"])

                self._initialization_node = node

                files["input.dat"] = copy.deepcopy(initialization_header)

                # Find the bond & angle types as needed for shake/rattle
                shake = self.shake_fix(P, eex)
                if shake != "":
                    extras["shake"] = shake
            else:
                history_nodes.append(node)
            node = node.next()

        files = self._prepare_input(
            files,
            nodes=history_nodes,
            read_restart=False,
            write_restart=True,
            extras=extras,
        )

        self._timing_data[11] = json.dumps(control)
        control = []

        files = self._execute_single_sim(files, np=np)

        self.analyze(nodes=history_nodes)

        self._trajectory = []

        printer.normal("")

        return next_node

    def _execute_single_sim(self, files, np=1, return_files=None):
        """
        Step #1: Execute input file
        """

        return_files = [
            "summary_*.txt",
            "trajectory_*.seamm_trj",
            "*.trj",
            "*.restart.*",
            "*.dump",
            "*.dump.*",
            "*.log",
            "*.dat",
            "log.lammps",
            "log.cite",
            "run_lammps",
        ]

        # Check for already having run
        path = Path(self.directory) / "success.dat"
        if path.exists():
            result = {}
            path = Path(self.directory) / "log.cite"
            if path.exists():
                result["log.cite"] = {
                    "data": path.read_text(),
                }
            path = Path(self.directory) / "stdout.txt"
            if path.exists():
                result["stdout"] = path.read_text()
            result["stderr"] = ""
        else:
            # Set up the computational limits and get the computational enviroment
            cl = {"NTASKS": np}
            ce = seamm_exec.computational_environment(cl)

            executor = self.flowchart.executor

            # Read configuration file for LAMMPS if it exists
            executor_type = executor.name
            full_config = configparser.ConfigParser()
            ini_dir = Path(self.global_options["root"]).expanduser()
            path = ini_dir / "lammps.ini"

            # If the config file doesn't exists, get the default
            if not path.exists():
                resources = importlib.resources.files("lammps_step") / "data"
                ini_text = (resources / "lammps.ini").read_text()
                txt_config = Configuration(path)
                txt_config.from_string(ini_text)

                # Work out the conda info needed
                txt_config.set_value("local", "conda", os.environ["CONDA_EXE"])
                txt_config.set_value("local", "conda-environment", "seamm-lammps")
                txt_config.save()
                printer.normal(f"Wrote the LAMMPS configuration file to {path}")
                printer.normal("")

            full_config.read(ini_dir / "lammps.ini")

            if executor_type not in full_config:
                path = shutil.which("lmp")
                mpi_path = shutil.which("mpirun")
                if path is None or mpi_path is None:
                    raise RuntimeError(
                        f"No section for '{executor_type}' in LAMMPS ini file "
                        f"({ini_dir / 'lammps.ini'}), nor in the defaults, nor "
                        "in the path!"
                    )
                else:
                    txt_config = Configuration(path)
                    txt_config.add_section(executor_type)
                    txt_config.set_value(executor_type, "installation", "local")
                    txt_config.set_value(
                        executor_type,
                        "code",
                        f"{mpi_path} -np {{NTASKS}} {path}",
                    )
                    txt_config.set_value(
                        executor_type,
                        "python",
                        f"mpirun -np {{NTASKS}} {shutil.which('python')}",
                    )
                    txt_config.save()
                    printer.normal(f"Wrote the LAMMPS configuration file to {path}")
                    printer.normal("")
                    full_config.read(ini_dir / "lammps.ini")
                    full_config.set(executor_type, "code", str(path))

            config = dict(full_config.items(executor_type))
            # Use the matching version of the seamm-mopac image by default.
            config["version"] = self.version

            executor_type = executor.name
            if executor_type not in full_config:
                raise RuntimeError(
                    f"No section for '{executor_type}' in LAMMPS ini file "
                    f"({ini_dir / 'lammps.ini'})"
                )
            config = dict(full_config.items(executor_type))

            # Setup the command lines
            cmd = []

            if "run_lammps" in files:
                cmd = ["{python}", "run_lammps"]
                if (
                    "GPUS" not in ce
                    and "cmd_args" in config
                    and config["cmd_args"] != ""
                ):
                    cmd.extend(["--cmd-args", config["cmd_args"]])
                if (
                    "GPUS" in ce
                    and "gpu_cmd_args" in config
                    and config["gpu_cmd_args"] != ""
                ):
                    cmd.extend(["--cmd-args", config["gpu_cmd_args"]])
            else:
                cmd = ["{code}"]
                if (
                    "GPUS" not in ce
                    and "cmd_args" in config
                    and config["cmd_args"] != ""
                ):
                    cmd.extend(config["cmd_args"].split())
                if (
                    "GPUS" in ce
                    and "gpu_cmd_args" in config
                    and config["gpu_cmd_args"] != ""
                ):
                    cmd.extend(config["gpu_cmd_args"].split())
                cmd.extend(["-in", "input.dat"])

            if "NGPUS" in ce:
                printer.important(
                    f"   LAMMPS running with {np} processes and {ce['NGPUS']} gpus."
                )
            else:
                printer.important(f"   LAMMPS using MPI with {np} processes.")
            printer.important("")

            cmd.extend([">", "stdout.txt", "2>", "stderr.txt"])

            if self._timing_data is not None:
                _, configuration = self.get_system_configuration()
                try:
                    self._timing_data[6] = configuration.to_smiles(
                        canonical=True, hydrogens=True
                    )
                except Exception:
                    self._timing_data[6] = ""
                try:
                    self._timing_data[7] = configuration.isomeric_smiles
                except Exception:
                    self._timing_data[7] = ""
                try:
                    self._timing_data[8] = configuration.formula[0]
                except Exception:
                    self._timing_data[7] = ""
                try:
                    self._timing_data[9] = str(configuration.charge)
                except Exception:
                    self._timing_data[9] = ""
                try:
                    self._timing_data[10] = str(configuration.spin_multiplicity)
                except Exception:
                    self._timing_data[10] = ""
                self._timing_data[5] = datetime.now(timezone.utc).isoformat()

            t0 = time.time_ns()

            result = executor.run(
                cmd=cmd,
                config=config,
                directory=self.directory,
                files=files,
                return_files=return_files,
                in_situ=True,
                shell=True,
                ce=ce,
            )

            t = (time.time_ns() - t0) / 1.0e9
            if self._timing_data is not None:
                self._timing_data[13] = f"{t:.3f}"
                self._timing_data[12] = str(ce["NTASKS"])
                try:
                    with self._timing_path.open("a", newline="") as fd:
                        writer = csv.writer(fd)
                        writer.writerow(self._timing_data)
                except Exception:
                    pass

            if not result:
                self.logger.error("There was an error running LAMMPS")
                return None

            self.logger.debug("\n" + pprint.pformat(result))

            f = os.path.join(self.directory, "stdout.txt")
            with open(f, mode="w") as fd:
                fd.write(result["stdout"])

        # Add the citations, getting the version from stdout and any citations
        if "log.cite" in result:
            self._add_lammps_citations(
                result["stdout"], cite=result["log.cite"]["data"]
            )
        else:
            self._add_lammps_citations(result["stdout"])

        initialization_header, eex = self._get_node_input(
            node=self._initialization_node, extras={"read_data": False}
        )
        files["input.dat"] = copy.deepcopy(initialization_header)

        # Write a small file to say that LAMMPS ran successfully, so cancel
        # skip if rerunning.
        path = Path(self.directory) / "success.dat"
        path.write_text("success")

        return files

    def _prepare_input(
        self,
        files,
        nodes=None,
        iteration=0,
        read_restart=False,
        write_restart=False,
        extras=None,
    ):
        _, configuration = self.get_system_configuration()
        python_script = None
        postscript = None
        if isinstance(nodes, list) is False:
            node_ids = [nodes._id[1]]
            todo = self._get_node_input(node=nodes, extras=extras)
            new_input_data = todo["script"]
            if todo["postscript"] is not None:
                postscript = todo["postscript"]
            if todo["use python"] and "python script" in todo:
                python_script = todo["python script"]
        else:
            node_ids = []
            new_input_data = []
            for n in nodes:
                node_ids.append(n._id[1])
                todo = self._get_node_input(node=n, extras=extras)
                new_input_data += todo["script"]
                if todo["postscript"] is not None:
                    postscript = todo["postscript"]
                if todo["use python"] and "python script" in todo:
                    python_script = todo["python script"]

        if postscript is None:
            if configuration.periodicity == 0:
                new_input_data.append(
                    "write_dump         all custom  final.dump id xu yu zu vx vy vz"
                    " modify flush yes sort id"
                )
            else:
                new_input_data.append(
                    "write_dump         all custom  final.dump id xsu ysu zsu vx vy vz"
                    " modify flush yes sort id"
                )
            new_input_data.append("")
            new_input_data.append("")
            new_input_data.append("info               computes fixes dumps out log")

        files["input.dat"] += new_input_data

        files["input.dat"] = "\n".join(files["input.dat"])
        self.logger.debug("input.dat:\n" + files["input.dat"])

        if postscript is not None:
            if configuration.periodicity == 0:
                postscript.append(
                    "write_dump         all custom  final.dump id xu yu zu vx vy vz"
                    " modify flush yes sort id"
                )
            else:
                postscript.append(
                    "write_dump         all custom  final.dump id xsu ysu zsu vx vy vz"
                    " modify flush yes sort id"
                )
            postscript.append("")
            postscript.append("info               computes fixes dumps out log")
            files["input_post.dat"] = "\n".join(postscript)
        if python_script is not None:
            files["run_lammps"] = python_script

        return files

    def _get_node_input(self, node=None, extras=None):
        try:
            ret = node.get_input(extras=extras)
        except Exception as e:
            print("Error running LAMMPS flowchart: {} in {}".format(str(e), str(node)))
            self.logger.critical(
                "Error running LAMMPS flowchart: {} in {}".format(str(e), str(node))
            )
            raise
        except:  # noqa: E722
            print(
                "Unexpected error running LAMMPS flowchart: {} in {}".format(
                    sys.exc_info()[0], str(node)
                )
            )
            self.logger.critical(
                "Unexpected error running LAMMPS flowchart: {} in {}".format(
                    sys.exc_info()[0], str(node)
                )
            )
            raise
        return ret

    def structure_data(self, eex, triclinic=False):
        """Create the LAMMPS structure file from the energy expression"""
        lines = []
        pair_table = []
        bond_table = []
        angle_table = []
        dihedral_table = []

        lines.append("Structure file for LAMMPS generated by a MolSSI flowchart")
        lines.append("{:10d} atoms".format(eex["n_atoms"]))
        lines.append("{:10d} atom types".format(eex["n_atom_types"]))
        if "n_bonds" in eex and eex["n_bonds"] > 0:
            lines.append("{:10d} bonds".format(eex["n_bonds"]))
            lines.append("{:10d} bond types".format(eex["n_bond_types"]))
        if "n_angles" in eex and eex["n_angles"] > 0:
            lines.append("{:10d} angles".format(eex["n_angles"]))
            lines.append("{:10d} angle types".format(eex["n_angle_types"]))
        if "n_torsions" in eex and eex["n_torsions"] > 0:
            lines.append("{:10d} dihedrals".format(eex["n_torsions"]))
            lines.append("{:10d} dihedral types".format(eex["n_torsion_types"]))
        if "n_oops" in eex and eex["n_oops"] > 0:
            lines.append("{:10d} impropers".format(eex["n_oops"]))
            lines.append("{:10d} improper types".format(eex["n_oop_types"]))

        # Find the box limits
        periodicity = eex["periodicity"]
        if periodicity == 3:
            a, b, c, alpha, beta, gamma = eex["cell"]
            lx, ly, lz, xy, xz, yz = LAMMPS.cell_to_box(a, b, c, alpha, beta, gamma)

            lines.append("{} {} xlo xhi".format(0.0, lx))
            lines.append("{} {} ylo yhi".format(0.0, ly))
            lines.append("{} {} zlo zhi".format(0.0, lz))

            xy = xy if abs(xy) > 1.0e-06 else 0.0
            xz = xz if abs(xy) > 1.0e-06 else 0.0
            yz = yz if abs(xy) > 1.0e-06 else 0.0

            if triclinic or xy > 0.0 or xz > 0.0 or yz > 0.0:
                lines.append("{} {} {} xy xz yz".format(xy, xz, yz))
        else:
            x, y, z, index = eex["atoms"][0]
            xlo = xhi = x
            ylo = yhi = y
            zlo = zhi = z
            for x, y, z, index in eex["atoms"]:
                xlo = x if x < xlo else xlo
                xhi = x if x > xhi else xhi
                ylo = y if y < ylo else ylo
                yhi = y if y > yhi else yhi
                zlo = z if z < zlo else zlo
                zhi = z if z > zhi else zhi

            # Some extra space....
            xlo -= 10.0
            xhi += 10.0
            ylo -= 10.0
            yhi += 10.0
            zlo -= 10.0
            zhi += 10.0

            lines.append("{} {} xlo xhi".format(xlo, xhi))
            lines.append("{} {} ylo yhi".format(ylo, yhi))
            lines.append("{} {} zlo zhi".format(zlo, zhi))

        # the atoms and their masses, etc.
        lines.append("")
        lines.append("Atoms")
        lines.append("")

        if "charges" in eex:
            if "molecule" in eex:
                for i, mol, xyz_index, q in zip(
                    range(1, eex["n_atoms"] + 1),
                    eex["molecule"],
                    eex["atoms"],
                    eex["charges"],
                ):
                    x, y, z, index = xyz_index
                    lines.append(
                        f"{i:6d} {mol + 1:6d} {index:6d} {q:6.3f} {x:12.7f} {y:12.7f} "
                        f"{z:12.7f}"
                    )
            else:
                for i, xyz_index, q in zip(
                    range(1, eex["n_atoms"] + 1), eex["atoms"], eex["charges"]
                ):
                    x, y, z, index = xyz_index
                    lines.append(
                        f"{i:6d} {index:6d} {q:6.3f} {x:12.7f} {y:12.7f} {z:12.7f}"
                    )
        else:
            for i, xyz_index in enumerate(eex["atoms"]):
                x, y, z, index = xyz_index
                lines.append(f"{i+1:6d} {index:6d} {x:12.7f} {y:12.7f} {z:12.7f}")

        _, configuration = self.get_system_configuration()
        if configuration.atoms.have_velocities:
            lines.append("")
            lines.append("Velocities")
            lines.append("")
            for i, vxyz in enumerate(
                configuration.atoms.get_velocities(fractionals=False), start=1
            ):
                vx, vy, vz = vxyz
                lines.append(f"{i:6d} {vx:12.7f} {vy:12.7f} {vz:12.7f}")

        lines.append("")
        lines.append("Masses")
        lines.append("")
        self._data["masses"] = []
        for i, parameters in zip(range(1, eex["n_atom_types"] + 1), eex["masses"]):
            mass, itype = parameters
            lines.append("{:6d} {} # {}".format(i, mass, itype))
            self._data["masses"].append(float(mass))

        # nonbonds
        if "nonbond parameters" in eex:
            # If using a hybrid/overlay form for e.g. Dreiding nonbonds, cannot use the
            # section in the structure file because it must have N or N*(N+1) lines,
            # which is not what we want
            forms = set([v[0] for v in eex["nonbond parameters"]])
            use_hybrid = len(forms) > 1

            if not use_hybrid:
                lines.append("")
                if len(eex["nonbond parameters"]) > eex["n_atom_types"]:
                    lines.append("PairIJ Coeffs")
                else:
                    lines.append("Pair Coeffs")
                lines.append("")
                i = 1
                j = 1
                for parameters in eex["nonbond parameters"]:
                    form, values, types, parameters_type, real_types = parameters
                    if form == "nonbond(9-6)":
                        lines.append(
                            f"{i:6d} {values['eps']} {values['rmin']} "
                            f"# {types[0]} --> {real_types[0]}"
                        )
                        i += 1
                    elif form == "nonbond(12-6)":
                        lines.append(
                            f"{i:6d} {values['eps']} {values['sigma']} "
                            f"# {types[0]} --> {real_types[0]}"
                        )
                        i += 1
                    elif form == "buckingham":
                        lines.append(
                            f"{j:6d} {i} {values['A']} {values['rho']} {values['C']}"
                            f" # {types[1]}-{types[0]} --> "
                            f"{real_types[1]}_{real_types[0]}"
                        )
                        if j == i:
                            i += 1
                            j = 1
                        else:
                            j += 1
        # bonds
        if "n_bonds" in eex and eex["n_bonds"] > 0:
            lines.append("")
            lines.append("Bonds")
            lines.append("")
            for counter, tmp in zip(range(1, eex["n_bonds"] + 1), eex["bonds"]):
                i, j, index = tmp
                lines.append("{:6d} {:6d} {:6d} {:6d}".format(counter, index, i, j))

            lines.append("")
            lines.append("Bond Coeffs")
            lines.append("")

            forms = set([v[0] for v in eex["bond parameters"]])
            use_hybrid = len(forms) > 1

            for counter, parameters in zip(
                range(1, eex["n_bond_types"] + 1), eex["bond parameters"]
            ):
                form, values, types, parameters_type, real_types = parameters
                if form == "quadratic_bond":
                    function = "harmonic" if use_hybrid else ""
                    line = f"{counter:6d} {function} {values['K2']} {values['R0']}"
                elif form == "quartic_bond":
                    function = "class2" if use_hybrid else ""
                    line = (
                        f"{counter:6d} {function} {values['R0']} "
                        f"{values['K2']} {values['K3']} {values['K4']}"
                    )
                line += (
                    f" # {types[0]}-{types[1]}" f" --> {real_types[0]}-{real_types[1]}"
                )
                lines.append(line)

        # angles
        if "n_angles" in eex and eex["n_angles"] > 0:
            lines.append("")
            lines.append("Angles")
            lines.append("")
            for counter, tmp in enumerate(eex["angles"], start=1):
                i, j, k, index = tmp
                lines.append(
                    "{:6d} {:6d} {:6d} {:6d} {:6d}".format(counter, index, i, j, k)
                )

            lines.append("")
            lines.append("Angle Coeffs")
            lines.append("")

            quartic_function = "class2" if "n_bond-bond_types" in eex else "quartic"
            forms = set([v[0] for v in eex["angle parameters"]])
            use_hybrid = len(forms) > 1

            for counter, parameters in zip(
                range(1, eex["n_angle_types"] + 1), eex["angle parameters"]
            ):
                form, values, types, parameters_type, real_types = parameters
                if form == "quadratic_angle":
                    function = "harmonic" if use_hybrid else ""
                    line = (
                        f"{counter:6d} {function} " f"{values['K2']} {values['Theta0']}"
                    )
                elif form == "quartic_angle":
                    function = quartic_function if use_hybrid else ""
                    line = (
                        f"{counter:6d} {function} {values['Theta0']} "
                        f"{values['K2']} {values['K3']} {values['K4']}"
                    )
                elif form == "cosine":
                    function = "cosine" if use_hybrid else ""
                    line = f"{counter:6d} {function} {values['K2']}"
                elif form == "cosine/squared":
                    function = "cosine/squared" if use_hybrid else ""
                    line = f"{counter:6d} {function} {values['K2']} {values['Theta0']}"
                elif form == "simple_fourier_angle":
                    function = "fourier/simple" if use_hybrid else ""
                    line = f"{counter:6d} {function} {values['K']} -1 {values['n']}"
                elif form == "tabulated_angle":
                    function = "table" if use_hybrid else ""
                    key = f"{types[0]}-{types[1]}-{types[2]}"
                    line = f"{counter:6d} {function} tabulated_angles.dat {key}"
                    angle_table.extend(self.angle_table(key, values))
                line += (
                    f" # {types[0]} {types[1]} {types[2]}"
                    f" --> {real_types[0]} {real_types[1]} {real_types[2]}"
                )
                lines.append(line)

            # bond-bond coefficients, which must match angles in order & number
            if "n_bond-bond_types" in eex:
                lines.append("")
                lines.append("BondBond Coeffs")
                lines.append("")
                for counter, parameters, angles in zip(
                    range(1, eex["n_bond-bond_types"] + 1),
                    eex["bond-bond parameters"],
                    eex["angle parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    angle_form = angles[0]
                    if angle_form == "quartic_angle":
                        function = "class2" if use_hybrid else ""
                        lines.append(
                            "{:6d} {} {} {} {}".format(
                                counter,
                                function,
                                values["K"],
                                values["R10"],
                                values["R20"],
                            )
                            + " # {}-{}-{} --> {}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                            )
                        )
                    else:
                        lines.append(
                            "{:6d} skip".format(counter)
                            + " # {}-{}-{} --> {}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                            )
                        )

                # bond-angles coefficients, which must match angles in order &
                # number
                lines.append("")
                lines.append("BondAngle Coeffs")
                lines.append("")
                for counter, parameters, angles in zip(
                    range(1, eex["n_bond-angle_types"] + 1),
                    eex["bond-angle parameters"],
                    eex["angle parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    angle_form = angles[0]
                    if angle_form == "quartic_angle":
                        function = "class2" if use_hybrid else ""
                        lines.append(
                            "{:6d} {} {} {} {} {}".format(
                                counter,
                                function,
                                values["K12"],
                                values["K23"],
                                values["R10"],
                                values["R20"],
                            )
                            + " # {}-{}-{} --> {}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                            )
                        )
                    else:
                        lines.append(
                            "{:6d} skip".format(counter)
                            + " # {}-{}-{} --> {}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                            )
                        )

        # torsions
        if "n_torsions" in eex and eex["n_torsions"] > 0:
            lines.append("")
            lines.append("Dihedrals")
            lines.append("")
            for counter, tmp in zip(range(1, eex["n_torsions"] + 1), eex["torsions"]):
                i, j, k, l, index = tmp
                lines.append(
                    "{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}".format(
                        counter, index, i, j, k, l
                    )
                )

            lines.append("")
            lines.append("Dihedral Coeffs")
            lines.append("")

            forms = set([v[0] for v in eex["torsion parameters"]])
            use_hybrid = len(forms) > 1

            for counter, parameters in zip(
                range(1, eex["n_torsion_types"] + 1), eex["torsion parameters"]
            ):
                form, values, types, parameters_type, real_types = parameters
                if form == "torsion_1":
                    KPhi = values["KPhi"]
                    n = values["n"]
                    Phi0 = values["Phi0"]

                    # Discover form is
                    #  KPhi * [1 + cos(n*Phi - Phi0)]
                    #  with trans = 180
                    #
                    #  For ethane, Phi0 = 0 so at Phi=180 E is min. Correct

                    # Lammps for is
                    #  KPhi * [1 + d*cos(n*Phi)]
                    #  with trans = 180
                    #
                    # Again for ethane, d=+1 and at Phi=180, E is min.
                    #
                    # Phi0 = 0   ==> d=+1
                    # Phi0 = 180 ==> d=-1

                    if float(Phi0) == 0.0:
                        d = "-1"
                    elif float(Phi0) == 180.0:
                        d = "+1"
                    else:
                        raise RuntimeError(
                            "LAMMPS cannot handle Phi0 = {}".format(Phi0)
                        )
                    function = "harmonic" if use_hybrid else ""
                    line = f"{counter:6d} {function} {KPhi} {d} {n}"
                elif form == "torsion_3":
                    function = "class2" if use_hybrid else ""
                    line = (
                        f"{counter:6d} {function} "
                        f"{values['V1']} {values['Phi0_1']} "
                        f"{values['V2']} {values['Phi0_2']} "
                        f"{values['V3']} {values['Phi0_3']} "
                    )
                elif form == "torsion_opls":
                    function = "opls" if use_hybrid else ""
                    line = (
                        f"{counter:6d} {function} "
                        f"{values['V1']} "
                        f"{values['V2']} "
                        f"{values['V3']} "
                        f"{values['V4']} "
                    )
                elif form == "torsion_charmm":
                    function = "charmm" if use_hybrid else ""
                    line = (
                        f"{counter:6d} {function} "
                        f"{values['K']} "
                        f"{values['n']} "
                        f"{values['Phi0']} "
                        f"{values['weight']} "
                    )
                line += (
                    f" # {types[0]}-{types[1]}-{types[2]}-{types[3]} "
                    f"--> {real_types[0]}-{real_types[1]}-"
                    f"{real_types[2]}-{real_types[3]}"
                )
                lines.append(line)

            # middle bond-torsion_3 coefficients, which must match torsions
            # in order & number
            if "n_middle_bond-torsion_3_types" in eex:
                lines.append("")
                lines.append("MiddleBondTorsion Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_middle_bond-torsion_3_types"] + 1),
                    eex["middle_bond-torsion_3 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        lines.append(
                            "{:6d} {} {} {} {} {}".format(
                                counter,
                                function,
                                values["V1"],
                                values["V2"],
                                values["V3"],
                                values["R0"],
                            )
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )
                    else:
                        lines.append(
                            "{:6d} skip".format(counter)
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )

                # end bond-torsion_3 coefficients, which must match torsions
                # in order & number
                lines.append("")
                lines.append("EndBondTorsion Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_end_bond-torsion_3_types"] + 1),
                    eex["end_bond-torsion_3 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        lines.append(
                            "{:6d} {} {} {} {} {} {} {} {} {}".format(
                                counter,
                                function,
                                values["V1_L"],
                                values["V2_L"],
                                values["V3_L"],
                                values["V1_R"],
                                values["V2_R"],
                                values["V3_R"],
                                values["R0_L"],
                                values["R0_R"],
                            )
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )
                    else:
                        lines.append(
                            "{:6d} skip".format(counter)
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )

                # angle-torsion_3 coefficients, which must match torsions
                # in order & number
                lines.append("")
                lines.append("AngleTorsion Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_angle-torsion_3_types"] + 1),
                    eex["angle-torsion_3 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        lines.append(
                            "{:6d} {} {} {} {} {} {} {} {} {}".format(
                                counter,
                                function,
                                values["V1_L"],
                                values["V2_L"],
                                values["V3_L"],
                                values["V1_R"],
                                values["V2_R"],
                                values["V3_R"],
                                values["Theta0_L"],
                                values["Theta0_R"],
                            )
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )
                    else:
                        lines.append(
                            "{:6d} skip".format(counter)
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )

                # angle-angle-torsion_1 coefficients, which must match torsions
                # in order & number
                lines.append("")
                lines.append("AngleAngleTorsion Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_angle-angle-torsion_1_types"] + 1),
                    eex["angle-angle-torsion_1 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        lines.append(
                            "{:6d} {} {} {} {}".format(
                                counter,
                                function,
                                values["K"],
                                values["Theta0_L"],
                                values["Theta0_R"],
                            )
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )
                    else:
                        lines.append(
                            "{:6d} skip".format(counter)
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )

                # bond-bond_1_3 coefficients, which must match torsions
                # in order & number
                lines.append("")
                lines.append("BondBond13 Coeffs")
                lines.append("")
                for counter, parameters, torsions in zip(
                    range(1, eex["n_bond-bond_1_3_types"] + 1),
                    eex["bond-bond_1_3 parameters"],
                    eex["torsion parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    torsion_form = torsions[0]
                    if torsion_form == "torsion_3":
                        function = "class2" if use_hybrid else ""
                        lines.append(
                            "{:6d} {} {} {} {}".format(
                                counter,
                                function,
                                values["K"],
                                values["R10"],
                                values["R30"],
                            )
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )
                    else:
                        lines.append(
                            "{:6d} skip".format(counter)
                            + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                                types[0],
                                types[1],
                                types[2],
                                types[3],
                                real_types[0],
                                real_types[1],
                                real_types[2],
                                real_types[3],
                            )
                        )

        # out-of-planes
        if "n_oops" in eex and eex["n_oops"] > 0:
            lines.append("")
            lines.append("Impropers")
            lines.append("")

            # Need forms to reorder impropers....
            form = {i: f[0] for i, f in enumerate(eex["oop parameters"], start=1)}

            for counter, tmp in enumerate(eex["oops"], start=1):
                i, j, k, l, index = tmp
                if form[index] in ("improper_opls", "dreiding_out_of_plane"):
                    lines.append(
                        "{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}".format(
                            counter, index, j, i, k, l
                        )
                    )
                else:
                    lines.append(
                        "{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}".format(
                            counter, index, i, j, k, l
                        )
                    )

            lines.append("")
            lines.append("Improper Coeffs")
            lines.append("")
            for counter, parameters in zip(
                range(1, eex["n_oop_types"] + 1), eex["oop parameters"]
            ):
                form, values, types, parameters_type, real_types = parameters
                if form == "wilson_out_of_plane":
                    lines.append(
                        "{:6d} {} {}".format(counter, values["K"], values["Chi0"])
                        + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                            types[0],
                            types[1],
                            types[2],
                            types[3],
                            real_types[0],
                            real_types[1],
                            real_types[2],
                            real_types[3],
                        )
                    )
                elif form == "dreiding_out_of_plane":
                    # divide by number of oops from atom (3 for all at the moment!)
                    lines.append(
                        f"{counter:6d} {float(values['K2'])/3:.4f} {values['Psi0']} "
                        f"# {types[1]}-{types[0]}-{types[2]}-{types[3]} --> "
                        f"{real_types[1]}-{real_types[0]}-{real_types[2]}-"
                        f"{real_types[3]}"
                    )
                elif form == "improper_opls":
                    # divide by two because OPLS uses V2/2 and CVS V.
                    lines.append(
                        f"{counter:6d} {float(values['V2'])/2:.4f} -1 2 "
                        f"# {types[0]}-{types[1]}-{types[2]}-{types[3]} --> "
                        f"{real_types[0]}-{real_types[1]}-{real_types[2]}-"
                        f"{real_types[3]}"
                    )
                else:
                    raise RuntimeError(f"Can't handle oop form '{form}'")

            # angle-angle
            if "n_angle-angle_types" in eex:
                lines.append("")
                lines.append("AngleAngle Coeffs")
                lines.append("")
                for counter, parameters in zip(
                    range(1, eex["n_angle-angle_types"] + 1),
                    eex["angle-angle parameters"],
                ):
                    form, values, types, parameters_type, real_types = parameters
                    lines.append(
                        "{:6d} {} {} {} {} {} {}".format(
                            counter,
                            values["K1"],
                            values["K2"],
                            values["K3"],
                            values["Theta10"],
                            values["Theta20"],
                            values["Theta30"],
                        )
                        + " # {}-{}-{}-{} --> {}-{}-{}-{}".format(
                            types[0],
                            types[1],
                            types[2],
                            types[3],
                            real_types[0],
                            real_types[1],
                            real_types[2],
                            real_types[3],
                        )
                    )

        lines.append("")
        lines.append("")
        return (
            "\n".join(lines),
            "\n".join(pair_table),
            "\n".join(bond_table),
            "\n".join(angle_table),
            "\n".join(dihedral_table),
        )

    def analyze(self, indent="", nodes=None, **kwargs):
        """Analyze the output of the calculation"""
        if isinstance(nodes, list) is False:
            nodes = [nodes]

        ret = {node._id[1]: None for node in nodes}

        run_dir = Path(self.directory)

        # Divide the log file into sections for the steps
        log_file = run_dir / "log.lammps"
        if not log_file.is_file():
            raise RuntimeError(f"Log file {log_file} is missing!")
        lines = log_file.read_text().split("\n")
        sections = {}
        section = None
        for line in lines:
            if line.startswith("# Step "):
                section = line[2:]
                sections[section] = []
            elif section is not None:
                sections[section].append(line)
            if line.startswith("Total wall time:"):
                try:
                    h, m, s = line.split()[3].split(":")
                    _time = 3600 * float(h) + 60 * float(m) + float(s)
                    self.results["t_lammps_wall"] = _time
                except Exception as _e:
                    print(f"Wall time exception {_e}")

        for node in nodes:
            for value in node.description:
                printer.important(value)
                printer.important(" ")

            subdir = run_dir / str(node._id[-1])

            # Move any old style files
            id_str = "_".join(str(e) for e in node._id)
            paths = sorted(run_dir.glob(f"*trajectory*{id_str}*.seamm_trj"))
            for path in paths:
                ensemble = str(path).split("_")[1]
                path.rename(subdir / f"{ensemble}_state.trj")

            paths = sorted(run_dir.glob(f"*summary*{id_str}*.seamm_trj"))
            for path in paths:
                ensemble = str(path).split("_")[1]
                path.rename(subdir / f"{ensemble}_summary.trj")

            # Find the state trajectory, if any
            paths = sorted(subdir.glob("*state.trj"))
            if len(paths) > 1:
                raise RuntimeError(f"More than one state.trj file for {id_str}.")
            elif len(paths) == 1:
                control_properties = lambda x: x not in ["tstep"]  # noqa: E731
                node_data, table = self.analyze_trajectory(
                    paths[0], control_properties=control_properties, node=node
                )
                # Save the state data as JSON
                data_file = subdir / "state.json"
                with data_file.open("w") as fd:
                    json.dump(
                        node_data, fd, indent=4, cls=CompactJSONEncoder, sort_keys=True
                    )

                # Get just the values from the node data
                values = {k: v["mean"] for k, v in node_data.items()}
                # And the other key values
                for k, v in node_data.items():
                    for key in ("stderr", "tau", "inefficiency", "n_samples"):
                        if key in v:
                            values[f"{k},{key}"] = v[key]
            else:
                node_data = None
                values = {}
                table = None

            values["model"] = self.model
            node.analyze(
                data=values,
                properties=node_data,
                table=table,
                output=sections[node.header],
            )

            ret[node._id[1]] = node_data

        return ret

    def analyze_trajectory(
        self,
        path,
        sampling_rate=20,
        control_properties=None,
        node=None,
    ):
        """Read a trajectory file and do the statistical analysis"""
        ff = self.get_variable("_forcefield")

        write_html = "html" in self.options and self.options["html"]

        results = {}

        table = {
            "Property": [],
            "Value": [],
            " ": [],
            "StdErr": [],
            "Units": [],
            "convergence": [],
            "tau": [],
            "inefficiency": [],
        }

        # Process the trajectory data
        # temporary until we sort out multiple runs
        self._trajectory = []
        with path.open() as fd:
            file_data = pandas.read_csv(
                fd,
                sep=" ",
                header=0,
                comment="!",
                usecols=control_properties,
                index_col="t",
            )
            self._trajectory.append(file_data.iloc[:-1])

        dt = lammps_step.from_lammps_units(
            file_data.index[1] - file_data.index[0], "fs"
        )
        dt_fs = dt.m_as("fs")
        data = pandas.concat(self._trajectory)
        data = data.reset_index(drop=True)
        data.index *= dt_fs

        self.logger.debug("Columns: {}".format(data.columns))
        self.logger.debug("  Types:\n{}".format(data.dtypes))

        printer.normal(f"       Analysis of {path.name}\n")

        # Work out the time step, rather than give the whole vector
        t = data.index
        t_units = "fs"
        len_trj = (len(t) - 1) * dt_fs
        if len_trj >= 4000000000:
            t_units = "ms"
        elif len_trj >= 4000000:
            t_units = "ns"
        elif len_trj >= 4000:
            t_units = "ps"
        t_max = float((len(t) - 1) * dt.m_as(t_units))

        for column in data.columns:
            if "Unnamed:" in column:
                continue
            meta_column = column.rstrip("0123456789")

            if meta_column in LAMMPS.display_title:
                meta_title = LAMMPS.display_title[meta_column]
            elif column in LAMMPS.display_title:
                meta_title = LAMMPS.display_title[column]
            else:
                meta_title = f"unknown: {column}"
            if meta_column in LAMMPS.display_units:
                meta_units = LAMMPS.display_units[meta_column]
            elif column in LAMMPS.display_units:
                meta_units = LAMMPS.display_units[column]
            else:
                meta_units = "???"

            have_warning = False
            have_acf_warning = False
            # Ignore first point, t=0, cause it might not be right.
            yy = data[column].to_numpy()[1:]

            self.logger.info("Analyzing {}, nsamples = {}".format(column, len(yy)))

            # compute indices of uncorrelated timeseries using pymbar
            # Their algorithm is quadratic in length of Y unless you
            # use 'nskip'. I set it so there are about 100 time origins, so
            # the convergence time is accurate to about 1%.
            nskip = yy.size // 100 + 1
            conv, inefficiency, Neff_max = timeseries.detect_equilibration(
                yy, nskip=nskip
            )

            self.logger.info(
                "  converged in {} steps, inefficiency = {}, Neff_max = {}".format(
                    conv, inefficiency, Neff_max
                )
            )

            if np.isnan(inefficiency) or np.isnan(Neff_max):
                # Apparently didn't converge!
                table["Property"].append(column)
                table["Value"].append("")
                table[" "].append("")
                table["StdErr"].append("")
                table["Units"].append("")
                table["convergence"].append("unconverged")
                table["tau"].append("")
                table["inefficiency"].append("")

                have_acf = False
                is_converged = False
            else:
                is_converged = True
                tau = dt_fs * (inefficiency - 1) / 2
                if tau < dt_fs / 2:
                    tau = dt_fs / 2
                t0 = conv * dt_fs
                y_t_equil = yy[conv:]
                indices = timeseries.subsample_correlated_data(
                    y_t_equil, g=inefficiency
                )

                if len(indices) == 0:
                    print("Problem with column " + column)
                    print("yy")
                    print(yy)
                    print("y_t_equil")
                    print(y_t_equil)
                    print("indices")
                    print(indices)
                    continue

                y_n = y_t_equil[indices]
                n_samples = len(y_n)
                mean = float(y_n.mean())
                std = float(y_n.std())
                sem = std / sqrt(n_samples)

                # Get the autocorrelation function
                if len(y_t_equil) < 8:
                    have_acf = False
                    have_acf_warning = True
                    acf_warning = "^"
                elif all(y_t_equil == y_t_equil[0]):
                    # A constant value, so no ACF
                    have_acf = False
                    acf_warning = ""
                else:
                    have_acf = True
                    acf_warning = ""
                    nlags = 4 * int(round(inefficiency + 0.5))

                    if nlags > int(len(y_t_equil) / 2):
                        nlags = int(len(y_t_equil) / 2)
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        acf, confidence = stattools.acf(
                            y_t_equil,
                            nlags=nlags,
                            alpha=0.05,
                            fft=True,
                            adjusted=False,
                        )

                results[column] = {}
                results[column]["mean"] = mean
                results[column]["stderr"] = sem
                results[column]["n_sample"] = n_samples
                results[column]["short_production"] = have_acf_warning
                results[column]["tau"] = float(tau)
                results[column]["inefficiency"] = float(inefficiency)
                results[column]["timestep"] = float(dt_fs)
                results[column]["rootname"] = path.stem
                if have_acf:
                    results[column]["acf"] = acf.tolist()
                    results[column]["acf_confidence"] = confidence.tolist()

                # Don't print or graph some properties, like heat flux
                if column in ("Jx", "Jy", "Jz"):
                    results[column]["values"] = data[column].to_numpy()

                # Work out units on convergence time
                conv_units = "fs"
                t_conv = t0
                if t0 >= 1000000000:
                    conv_units = "ms"
                    t_conv = t0 / 1000000000
                elif t0 >= 1000000:
                    conv_units = "ns"
                    t_conv = t0 / 1000000
                elif t0 >= 1000:
                    conv_units = "ps"
                    t_conv = t0 / 1000

                # Work out units on autocorrelation time
                tau_units = "fs"
                t_tau = tau
                if tau >= 1000000000:
                    tau_units = "ms"
                    t_tau = tau / 1000000000
                elif tau >= 1000000:
                    tau_units = "ns"
                    t_tau = tau / 1000000
                elif tau >= 1000:
                    tau_units = "ps"
                    t_tau = tau / 1000

                if n_samples < 100:
                    have_warning = True
                    warn = "*"
                else:
                    warn = " "

                results[column]["few_neff"] = have_warning

                table["Property"].append(f"{column}{warn}")
                table["Value"].append(f"{mean:.3f}")
                table[" "].append("±")
                table["StdErr"].append(f"{sem:.3f}")
                table["Units"].append(meta_units)
                table["convergence"].append(f"{t_conv:.2f} {conv_units}")
                table["tau"].append(f"{t_tau:.1f} {tau_units}{acf_warning}")
                table["inefficiency"].append(f"{inefficiency:.1f}")

                if column == "Etot" and ff.ff_form == "reaxff":
                    Eat = self._atomic_energy_sum
                    if Eat != 0.0:
                        results["DfH0_reax"] = {}
                        results["DfH0_reax"]["mean"] = mean + Eat
                        results["DfH0_reax"]["stderr"] = sem
                        results["DfH0_reax"]["n_sample"] = n_samples
                        results["DfH0_reax"]["short_production"] = have_acf_warning
                        results["DfH0_reax"]["tau"] = float(tau)
                        results["DfH0_reax"]["inefficiency"] = float(inefficiency)
                        results["DfH0_reax"]["timestep"] = float(dt_fs)
                        results["DfH0_reax"]["rootname"] = path.stem
                        if have_acf:
                            results["DfH0_reax"]["acf"] = acf.tolist()
                            results["DfH0_reax"]["acf_confidence"] = confidence.tolist()
                        results["DfH0_reax"]["few_neff"] = have_warning

                        table["Property"].append(f"DfH0_reax{warn}")
                        table["Value"].append(f"{mean + Eat:.3f}")
                        table[" "].append("±")
                        table["StdErr"].append(f"{sem:.3f}")
                        table["Units"].append(meta_units)
                        table["convergence"].append(f"{t_conv:.2f} {conv_units}")
                        table["tau"].append(f"{t_tau:.1f} {tau_units}{acf_warning}")
                        table["inefficiency"].append(f"{inefficiency:.1f}")

            # Create graphs of the property
            figure = self.create_figure(
                module_path=(self.__module__.split(".")[0], "seamm"),
                template="line.graph_template",
                title=meta_title,
            )

            # The autocorrelation function
            if have_acf:
                plot_acf = figure.add_plot("acf")

                t_acf_units = "fs"
                len_acf = (len(acf) - 1) * dt_fs
                if len_acf >= 2000000000:
                    t_acf_units = "ms"
                elif len_acf >= 2000000:
                    t_acf_units = "ns"
                elif len_acf >= 2000:
                    t_acf_units = "ps"

                x_acf_axis = plot_acf.add_axis(
                    "x", label="Time ({})".format(t_acf_units)
                )
                y_acf_axis = plot_acf.add_axis("y", label="acf", anchor=x_acf_axis)
                x_acf_axis.anchor = y_acf_axis

                # Put the fit to the autocorrelation time in first so the
                # subsequent trajectory trace sits in top
                ts = 0.0
                fit = [1.0]
                for step in range(len(acf) - 1):
                    ts += dt_fs
                    fit.append(exp(-ts / tau))

                plot_acf.add_trace(
                    x_axis=x_acf_axis,
                    y_axis=y_acf_axis,
                    name="fit",
                    x0=0,
                    dx=dt.m_as(t_acf_units),
                    xlabel="t",
                    xunits=t_acf_units,
                    y=fit,
                    ylabel="fit",
                    yunits="",
                    color="gray",
                )

                # the partly transparent error band
                yplus = []
                yminus = []
                t_acf = []
                tmp = 0.0
                for lower, upper in confidence:
                    t_acf.append(tmp)
                    yplus.append(upper)
                    yminus.append(lower)
                    tmp += dt.m_as(t_acf_units)

                plot_acf.add_trace(
                    x_axis=x_acf_axis,
                    y_axis=y_acf_axis,
                    name="stderr",
                    x=t_acf + t_acf[::-1],
                    xlabel="t",
                    xunits=t_acf_units,
                    y=yplus + yminus[::-1],
                    ylabel="stderr",
                    yunits=meta_units,
                    color="rgba(211,211,211,0.5)",
                    fill="toself",
                )

                # And the acf plot last
                plot_acf.add_trace(
                    x_axis=x_acf_axis,
                    y_axis=y_acf_axis,
                    name="acf",
                    x0=0,
                    dx=dt.m_as(t_acf_units),
                    xlabel="t",
                    xunits=t_acf_units,
                    y=list(acf),
                    ylabel="acf",
                    yunits="",
                    color="red",
                )

            # The property data over the trajectory
            y = list(data[column])

            plot = figure.add_plot("trj")

            ylabel = meta_title
            if meta_units != "":
                ylabel += f" ({meta_units})"

            x_axis = plot.add_axis("x", label="Time ({})".format(t_units))
            y_axis = plot.add_axis("y", label=ylabel, anchor=x_axis)
            x_axis.anchor = y_axis

            # Add the trajectory, error band and median value in that order so
            # stack in a nice order.

            # Add the trajectory
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=column,
                x0=0,
                dx=dt.m_as(t_units),
                xlabel="t",
                xunits=t_units,
                y=list(y),
                ylabel=column,
                yunits=meta_units,
                color="#4dbd74",
            )

            if is_converged:
                # the partly transparent error band
                t_min = t0 / dt_fs * dt.m_as(t_units)
                plot.add_trace(
                    x_axis=x_axis,
                    y_axis=y_axis,
                    name="sem",
                    x=[t_min, t_max, t_max, t_min],
                    xlabel="t",
                    xunits=t_units,
                    y=[mean + sem, mean + sem, mean - sem, mean - sem],
                    ylabel="sem",
                    yunits=meta_units,
                    color="rgba(211,211,211,0.5)",
                    fill="toself",
                )

                # and finally the median value so it is on top
                plot.add_trace(
                    x_axis=x_axis,
                    y_axis=y_axis,
                    name="average",
                    x=[t_min, t_max],
                    xlabel="t",
                    xunits=t_units,
                    y=[mean, mean],
                    ylabel="average",
                    yunits=meta_units,
                    color="black",
                )

            if have_acf:
                figure.grid_plots("trj - acf")
            else:
                figure.grid_plots("trj")
            if node is None:
                path = Path(f"{path.stem}_{column}.graph")
            else:
                node_path = Path(node.directory)
                node_path.mkdir(parents=True, exist_ok=True)
                path = node_path / f"{column}.graph"

            figure.dump(path)

            if write_html:
                figure.template = "line.html_template"
                figure.dump(path.with_suffix(".html"))

        # Add citations for pymbar
        self.references.cite(
            raw=self._bibliography["mbar"],
            alias="pymbar-1",
            module="lammps_step",
            level=1,
            note="The main reference for pymbar.",
        )
        self.references.cite(
            raw=self._bibliography["histogram"],
            alias="pymbar-2",
            module="lammps_step",
            level=2,
            note="The second citation for pymbar",
        )
        self.references.cite(
            raw=self._bibliography["convergence"],
            alias="pymbar-3",
            module="lammps_step",
            level=2,
            note="The third citation for pymbar",
        )

        return results, table

    def shake_fix(self, P, eex):
        """Create the 'fix shake' line needed for handling waters and X-H.

        Parameters
        ----------
        P : dict
            The parameters for the initialization step as a dict.
        eex : dict
            The energy expression for this calculation

        Returns
        -------
        line : str
            The correct fix line for LAMMPS
        """

        bond_types = {}
        angle_types = {}

        # Water models
        if P["rigid_waters"]:
            # waters = seamm_util.water_models.Water.find_waters(data.structure)  # noqa: E501

            waters = []

            if len(waters) > 0:
                atoms = []
                for i, j, k in waters:
                    atoms.append(i)
                    atoms.append(j)
                    atoms.append(k)
                if "n_bonds" in eex and eex["n_bonds"] > 0:
                    for i, j, index in eex["bonds"]:
                        if i in atoms and j in atoms:
                            bond_types[index] = 1
                if "n_angles" in eex and eex["n_angles"] > 0:
                    for i, j, k, index in eex["angles"]:
                        if i in atoms and j in atoms and k in atoms:
                            angle_types[index] = 1

        # Fixing bond lengths of X-H bonds...
        if "n_bonds" in eex and eex["n_bonds"] > 0:
            fix_bonds = P["fix_XH_bond_lengths"]
            elements = eex["elements"]
            if fix_bonds == "CH":
                for i, j, index in eex["bonds"]:
                    if (elements[i] == "C" and elements[j] == "H") or (
                        elements[i] == "H" and elements[j] == "C"
                    ):
                        bond_types[index] = 1
            elif fix_bonds == "all":
                for i, j, index in eex["bonds"]:
                    if elements[i] == "H" or elements[j] == "H":
                        bond_types[index] = 1

        # And the result is ....
        if len(bond_types) > 0:
            result = "fix                 {} all rattle 0.001 20 1000 b "
            for bond_type in bond_types.keys():
                result += " " + str(bond_type)
            if len(angle_types) > 0:
                result += " a "
                for angle_type in angle_types.keys():
                    result += " " + str(angle_type)
        else:
            result = ""

        return result

    def read_dump(self, dumpfile):
        """Read the LAMMPS dumpfile and update the system.

        Parameters
        ----------
        dumpfile : str
            The filename (or path) to the dumpfile.
        """
        self.logger.debug("Reading dump file '{}'".format(dumpfile))

        _, configuration = self.get_system_configuration()
        n_atoms = configuration.n_atoms

        cell = None
        xyz = None
        vxyz = None

        section = ""
        section_lines = []
        with open(dumpfile, "r") as fd:
            lineno = 0
            for line in fd:
                line = line.strip()
                lineno += 1
                if lineno == 1:
                    if line[0:5] != "ITEM:":
                        raise RuntimeError(
                            "Error reading dump file '" + dumpfile + "': The "
                            "first line is incorrect! (" + line + ")"
                        )
                    section = line[6:].strip()
                    section_lines = []
                    self.logger.debug("   section = " + section)
                    continue

                if line[0:5] == "ITEM:":
                    # end a section
                    self.logger.debug("  processing section '{}'".format(section))
                    if "BOX BOUNDS" in section:
                        if len(section.split()) == 8:
                            xlo_bound, xhi_bound, xy = section_lines[0].split()
                            ylo_bound, yhi_bound, xz = section_lines[1].split()
                            zlo, zhi, yz = section_lines[2].split()

                            xlo_bound = float(xlo_bound)
                            xhi_bound = float(xhi_bound)
                            ylo_bound = float(ylo_bound)
                            yhi_bound = float(yhi_bound)
                            zlo = float(zlo)
                            zhi = float(zhi)
                            xy = float(xy)
                            xz = float(xz)
                            yz = float(yz)

                            xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
                            xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
                            ylo = ylo_bound - min(0.0, yz)
                            yhi = yhi_bound - max(0.0, yz)
                            cell = LAMMPS.box_to_cell(
                                xhi - xlo, yhi - ylo, zhi - zlo, xy, xz, yz
                            )
                        else:
                            xlo, xhi = section_lines[0].split()
                            ylo, yhi = section_lines[1].split()
                            zlo, zhi = section_lines[2].split()

                            xlo = float(xlo)
                            xhi = float(xhi)
                            ylo = float(ylo)
                            yhi = float(yhi)
                            zlo = float(zlo)
                            zhi = float(zhi)

                            cell = (xhi - xlo, yhi - ylo, zhi - zlo, 90, 90, 90)
                    elif section == "NUMBER OF ATOMS":
                        if int(section_lines[0]) != n_atoms:
                            raise RuntimeError(
                                "Number of atoms has changed! {} to {}".format(
                                    n_atoms, section_lines[0]
                                )
                            )
                    elif "ATOMS" in section:
                        xyz = []
                        vxyz = []
                        keys = section.split()[1:]
                        if keys[1:4] == ["x", "y", "z"] or keys[1:4] == [
                            "xu",
                            "yu",
                            "zu",
                        ]:
                            fractional = False
                        elif keys[1:4] == ["xs", "ys", "zs"] or keys[1:4] == [
                            "xsu",
                            "ysu",
                            "zsu",
                        ]:
                            fractional = True
                        else:
                            logger.error(f"Can't handle dump file, {keys=}")
                        if len(keys) >= 7 and keys[4:7] == ["vx", "vy", "vz"]:
                            have_velocities = True
                            factor = lammps_step.from_lammps_units(1, "fs").magnitude
                            factor = 1 / factor
                            for tmp in section_lines:
                                x, y, z, vx, vy, vz = tmp.split()[1:7]
                                xyz.append([float(x), float(y), float(z)])
                                vxyz.append(
                                    [
                                        factor * float(vx),
                                        factor * float(vy),
                                        factor * float(vz),
                                    ]
                                )
                        else:
                            have_velocities = False
                            for tmp in section_lines:
                                x, y, z = tmp.split()[1:4]
                                xyz.append([float(x), float(y), float(z)])
                    section = line[6:].strip()
                    section_lines = []
                else:
                    section_lines.append(line)

        # Clean up the last section
        xyz = []
        vxyz = []
        if "ATOMS" in section:
            self.logger.debug("  processing section '{}'".format(section))
            self.logger.debug("  handling the atoms")
            keys = section.split()[1:]
            if keys[1:4] == ["x", "y", "z"] or keys[1:4] == ["xu", "yu", "zu"]:
                fractional = False
            elif keys[1:4] == ["xs", "ys", "zs"] or keys[1:4] == ["xsu", "ysu", "zsu"]:
                fractional = True
            else:
                logger.error(f"Can't handle dump file, {keys=}")
            if len(keys) >= 7 and keys[4:7] == ["vx", "vy", "vz"]:
                have_velocities = True
                factor = 1.0 / lammps_step.from_lammps_units(1, "fs").magnitude
                for tmp in section_lines:
                    x, y, z, vx, vy, vz = tmp.split()[1:7]
                    xyz.append([float(x), float(y), float(z)])
                    vxyz.append(
                        [factor * float(vx), factor * float(vy), factor * float(vz)]
                    )
            else:
                have_velocities = False
                for tmp in section_lines:
                    x, y, z = tmp.split()[1:4]
                    xyz.append([float(x), float(y), float(z)])

        if not have_velocities:
            vxyz = None

        if configuration.periodicity == 0:
            cell = None

        return xyz, fractional, cell, vxyz

    def _add_lammps_citations(self, text, cite=None):
        """Add the two main citations for LAMMPS, getting the version from stdout
        text.

        Parameters
        ----------
        text : str
            The standard output from LAMMPS

        Returns
        -------
        None
        """
        # Add the JCP paper
        self.references.cite(
            raw=self._bibliography["PLIMPTON19951"],
            alias="lammps-jcp",
            module="lammps_step",
            level=1,
            note="The principle LAMMPS citation.",
        )

        # And the citation to the LAMMPS code itself
        lines = text.splitlines()

        if len(lines) == 0:
            return

        line = lines[0]
        tmp = line.split(" ")
        if len(tmp) != 4:
            self.logger.info(f"Cannot get LAMMPS version: '{line}'")
            return

        month = tmp[2]
        year = tmp[3].rstrip(")")
        version = " ".join(tmp[1:])

        try:
            template = string.Template(self._bibliography["lammps"])

            citation = template.substitute(month=month, version=version, year=year)

            self.references.cite(
                raw=citation,
                alias="lammps-exe",
                module="lammps_step",
                level=1,
                note="The principle citation for the LAMMPS executable.",
            )

        except Exception as e:
            printer.important(f"Exception in citation {type(e)}: {e}")
            printer.important(traceback.format_exc())

        # If there is a log.cite file, process it
        if cite is not None:
            self.logger.debug("log.cite\n" + cite + "\n")
            bibliography = {}
            tmp = bibtexparser.loads(cite).entries_dict
            writer = bibtexparser.bwriter.BibTexWriter()
            for key, data in tmp.items():
                self.logger.info(f"      {key}")
                bibliography[key] = writer._entry_to_bibtex(data)
            self.logger.debug("Bibliography\n" + pprint.pformat(bibliography))

            for entry in bibliography:
                if entry.lower() in ("commment",):
                    continue
                self.references.cite(
                    raw=bibliography[entry],
                    alias=entry,
                    module="lammps_step",
                    level=1,
                    note="LAMMPS citations from log.cite.",
                )

    def angle_table(self, key, values):
        """Create a section of the tabulated angle file.

        Parameters
        ----------
        key : str
            The name for this section in the file.

        values : {str: int, float, or str}
            The dictionary of the constants for this angle term.

        Returns
        -------
        [str]
            A list of lines of the tabulated angle file.

        values looks like this::

            {
                'reference': '5',
                'Eqn': 'K/8*(1-cos(n*Theta)) + A/(2*Rb*sin(Theta/2))**12',
                'K': 278.4416826003824,
                'n': 4,
                'Rb': 1.606,
                'A': 11.950286806883364,
                'zero-shift': 0.001
            }
        """
        lines = []

        data = {**values}
        if "reference" in data:
            del data["reference"]
        eqn = data.pop("Eqn")
        thetas, Es, dEs = tabulate_angle(eqn, data)

        lines.append(key)
        lines.append(f"N {len(thetas)}")
        lines.append("")
        for i, data in enumerate(zip(thetas, Es, dEs), start=1):
            theta, E, dE = data
            lines.append(f"{i:5d} {theta:9.4f} {E:40}{-dE:40}")
        lines.append("")

        return lines

    def get_dump(self, dumpfile):
        """Read the LAMMPS dumpfile and return the data.

        Parameters
        ----------
        dumpfile : str
            The filename (or path) to the dumpfile.
        """
        self.logger.debug("Reading dump file '{}'".format(dumpfile))

        # system_db = self.get_variable("_system_db")
        # configuration = system_db.system.configuration
        # periodicity = configuration.periodicity
        # n_atoms = configuration.n_atoms

        result = {
            "timestep": [],
            "n_atoms": [],
            "cell": [],
            "data": [],
        }
        section = ""
        section_lines = []
        with open(dumpfile, "r") as fd:
            lineno = 0
            for line in fd:
                line = line.strip()
                lineno += 1
                if lineno == 1:
                    if line[0:5] != "ITEM:":
                        raise RuntimeError(
                            "Error reading dump file '" + dumpfile + "': The "
                            "first line is incorrect! (" + line + ")"
                        )
                    section = line[6:].strip()
                    section_lines = []
                    self.logger.debug("   section = " + section)
                    continue

                if line[0:5] == "ITEM:":
                    # end a section
                    self.logger.debug("  processing section '{}'".format(section))
                    if section == "TIMESTEP":
                        result["timestep"].append(int(section_lines[0]))
                    elif section == "NUMBER OF ATOMS":
                        n_atoms = int(section_lines[0])
                        result["n_atoms"].append(n_atoms)
                    elif "BOX BOUNDS" in section:
                        if len(section.split()) == 8:
                            xlo_bound, xhi_bound, xy = section_lines[0].split()
                            ylo_bound, yhi_bound, xz = section_lines[1].split()
                            zlo, zhi, yz = section_lines[2].split()

                            xlo_bound = float(xlo_bound)
                            xhi_bound = float(xhi_bound)
                            ylo_bound = float(ylo_bound)
                            yhi_bound = float(yhi_bound)
                            zlo = float(zlo)
                            zhi = float(zhi)
                            xy = float(xy)
                            xz = float(xz)
                            yz = float(yz)

                            xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
                            xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
                            ylo = ylo_bound - min(0.0, yz)
                            yhi = yhi_bound - max(0.0, yz)
                            cell = LAMMPS.box_to_cell(
                                xhi - xlo, yhi - ylo, zhi - zlo, xy, xz, yz
                            )
                        else:
                            xlo, xhi = section_lines[0].split()
                            ylo, yhi = section_lines[1].split()
                            zlo, zhi = section_lines[2].split()

                            xlo = float(xlo)
                            xhi = float(xhi)
                            ylo = float(ylo)
                            yhi = float(yhi)
                            zlo = float(zlo)
                            zhi = float(zhi)

                            cell = (xhi - xlo, yhi - ylo, zhi - zlo, 90, 90, 90)
                        result["cell"].append(cell)
                    elif "ATOMS" in section:
                        keys = section.split()[2:]
                        if "fields" not in result:
                            result["fields"] = keys
                        elif keys != result["fields"]:
                            raise ValueError(
                                f"Error reading dump file '{dumpfile}': The fields in "
                                "the ATOMS section do not match the previous fields!\n"
                                f"  timestep: {result['timestep'][-1]}\n"
                                f"  previous: {result['fields']=}\n"
                                f"   current: {keys=}"
                            )
                        tmp = [[None] * n_atoms for key in keys]
                        for line in section_lines:
                            values = line.split()
                            _id = int(values[0]) - 1
                            for index, value in enumerate(values[1:]):
                                key = keys[index]
                                if key in ("mol", "proc", "procp1", "type"):
                                    value = int(value)
                                elif key in ("element",):
                                    value = value.strip()
                                else:
                                    value = float(value)
                                tmp[index][_id] = value
                        result["data"].append(tmp)
                    section = line[6:].strip()
                    section_lines = []
                else:
                    section_lines.append(line)

        # Clean up the last section
        if "ATOMS" in section:
            keys = section.split()[2:]
            if "fields" not in result:
                result["fields"] = keys
            elif keys != result["fields"]:
                raise ValueError(
                    f"Error reading dump file '{dumpfile}': The fields in "
                    "the ATOMS section do not match the previous fields!\n"
                    f"  timestep: {result['timestep'][-1]}\n"
                    f"  previous: {result['fields']=}\n"
                    f"   current: {keys=}"
                )
            tmp = [[None] * n_atoms for key in keys]
            for line in section_lines:
                values = line.split()
                _id = int(values[0]) - 1
                for index, value in enumerate(values[1:]):
                    key = keys[index]
                    if key in ("mol", "proc", "procp1", "type"):
                        value = int(value)
                    elif key in ("element",):
                        value = value.strip()
                    else:
                        value = float(value)
                    tmp[index][_id] = value
            result["data"].append(tmp)

        return result
