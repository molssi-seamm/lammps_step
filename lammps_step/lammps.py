# -*- coding: utf-8 -*-

"""A node or step for LAMMPS in a flowchart"""

import calendar
import copy
import glob
import logging
from math import sqrt, exp, degrees, radians, cos, acos
import os
import os.path
import pprint
import re
import string
import sys
import traceback

import bibtexparser
import numpy
import pandas
import psutil
from scipy.stats import t as t_student
import statsmodels.tsa.stattools as stattools

import lammps_step
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

from pymbar import timeseries

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter('lammps')

bond_style = {
    'quadratic_bond': 'harmonic',
    'quartic_bond': 'class2',
    'fene': 'fene',
    'morse': 'morse',
}

angle_style = {
    'quadratic_angle': 'harmonic',
    'quartic_angle': 'class2',
}

dihedral_style = {
    'torsion_1': 'harmonic',
    'torsion_3': 'class2',
}

improper_style = {
    'wilson_out_of_plane': 'class2',
}


class LAMMPS(seamm.Node):
    display_units = {
        "T": "K",
        "P": "atm",
        "t": "fs",
        "density": "g/mL",
        "a": "Å",
        "b": "Å",
        "c": "Å",
        "Etot": "kcal/mol",
        "Eke": "kcal/mol",
        "Epe": "kcal/mol",
        "Emol": "kcal/mol",
        "Epair": "kcal/mol",
    }
    display_title = {
        "T": "Temperature",
        "P": "Pressure",
        "t": "Time",
        "density": "Density",
        "a": "a lattice parameter",
        "b": "b lattice parameter",
        "c": "c lattice parameter",
        "Etot": "Total Energy",
        "Eke": "Kinetic Energy",
        "Epe": "Potential Energy",
        "Emol": "Molecular Energy, Valence Terms",
        "Epair": "Pair (Nonbond) Energy",
    }

    def __init__(
        self,
        flowchart=None,
        namespace='org.molssi.seamm.lammps',
        extension=None
    ):
        """Setup the main LAMMPS step

        Keyword arguments:
        """
        logger.debug('Creating LAMMPS {}'.format(self))

        # The subflowchart
        self.subflowchart = seamm.Flowchart(
            parent=self, name='LAMMPS', namespace=namespace
        )
        self._initialization_node = None
        self._trajectory = []
        self._data = {}

        super().__init__(
            flowchart=flowchart,
            title='LAMMPS',
            extension=extension,
            logger=logger
        )

    @property
    def version(self):
        """The semantic version of this module.
        """
        return lammps_step.__version__

    @property
    def git_revision(self):
        """The git version of this module.
        """
        return lammps_step.__git_revision__

    @staticmethod
    def box_to_cell(lx, ly, lz, xy, xz, yz):
        """Convert the LAMMPS box definition to cell parameters.
        """
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
            lx = 0
            xy = b * cos(radians(gamma))
            xz = c * cos(radians(beta))
            ly = sqrt(b**2 - xy**2)
            yz = (b * c * cos(radians(alpha)) - xy * xz) / ly
            lz = sqrt(c**2 - xz**2 - yz**2)

        return (lx, ly, lz, xy, xz, yz)

    def create_parser(self):
        """Setup the command-line / config file parser
        """
        parser_name = self.step_type
        parser = seamm.getParser()

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
            '--lammps-serial',
            default='lmp_serial',
            help='the serial version of LAMMPS'
        )
        parser.add_argument(
            parser_name,
            '--mpi-exe',
            default='mpiexec',
            help='the mpi executable'
        )
        parser.add_argument(
            parser_name,
            '--lammps-mpi',
            default='lmp_mpi',
            help='the mpi version of LAMMPS'
        )
        parser.add_argument(
            parser_name,
            '--ncores',
            default='available',
            help=(
                'The maximum number of cores to use for LAMMPS. '
                'Default: all available cores.'
            )
        )
        parser.add_argument(
            parser_name,
            '--atoms-per-core',
            type=int,
            default='1000',
            help='the optimal number of atoms per core for LAMMPS'
        )
        parser.add_argument(
            parser_name,
            '--html',
            action='store_true',
            help='whether to write out html files for graphs, etc.'
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
        node = self.subflowchart.get_node('1').next()

        text = self.header + '\n\n'

        while node is not None:
            try:
                text += __(node.description_text(), indent=3 * ' ').__str__()
            except Exception as e:
                print(
                    'Error describing LAMMPS flowchart: {} in {}'.format(
                        str(e), str(node)
                    )
                )
                self.logger.critical(
                    'Error describing LAMMPS flowchart: {} in {}'.format(
                        str(e), str(node)
                    )
                )
                raise
            except:  # noqa: E722
                print(
                    "Unexpected error describing LAMMPS flowchart: {} in {}"
                    .format(sys.exc_info()[0], str(node))
                )
                self.logger.critical(
                    "Unexpected error describing LAMMPS flowchart: {} in {}"
                    .format(sys.exc_info()[0], str(node))
                )
                raise
            text += '\n'
            node = node.next()

        return text

    def run(self):
        """Run a LAMMPS simulation
        """
        # Add a citation for this plug-in
        try:
            template = string.Template(self._bibliography['lammps_step'])

            version = seamm.__version__
            year, month = version.split('.')[0:2]
            month = calendar.month_abbr[int(month)].lower()
            citation = template.substitute(
                month=month, version=version, year=year
            )

            self.references.cite(
                raw=citation,
                alias='lammps_step',
                module='lammps_step',
                level=2,
                note='The principle citation for the LAMMPS step in SEAMM.'
            )

        except Exception as e:
            printer.important(f'Exception in citation {type(e)}: {e}')
            printer.important(traceback.format_exc())

        system = self.get_variable('_system')

        n_atoms = system.n_atoms()
        if n_atoms == 0:
            self.logger.error('LAMMPS run(): there is no structure!')
            raise RuntimeError('LAMMPS run(): there is no structure!')

        next_node = super().run(printer)

        # Get the options
        o = self.options
        global_options = self.global_options

        # Whether to run parallel and if so, how many mpi processes
        if global_options['parallelism'] in ('any', 'mpi'):
            np = o['ncores']
            if np == 'available':
                np = global_options['ncores']

            if np == 'available':
                np = int(round(n_atoms / o['atoms_per_core']))
                if np < 1:
                    np = 1

                # How many processors does this node have?
                n_cores = psutil.cpu_count(logical=False)
                self.logger.info('The number of cores is {}'.format(n_cores))

                if np > n_cores:
                    self.logger.info(
                        f'LAMMPS could use {np} cores, but only {n_cores} are '
                        'available'
                    )
                    np = n_cores
            else:
                np = int(np)
        else:
            np = 1

        # Print headers and get to work
        printer.important(self.header)
        if np > 1:
            printer.important(
                '    LAMMPS using MPI with {} processes.\n'.format(np)
            )
        else:
            printer.important('   LAMMPS using the serial version.\n')

        self.subflowchart.root_directory = self.flowchart.root_directory

        files = {}

        # Get the first real node
        node = self.subflowchart.get_node('1').next()

        extras = {}

        history_nodes = []

        # Create overall directory for the lammps step

        os.makedirs(self.directory, exist_ok=True)

        files = {}

        while node is not None:

            if isinstance(node, lammps_step.Initialization):
                initialization_header, eex = self._get_node_input(
                    node=node, extras={'read_data': True}
                )
                files['structure'] = {}
                files['structure']['filename'] = 'structure.dat'
                files['structure']['data'] = '\n'.join(
                    self.structure_data(eex)
                )

                f = os.path.join(
                    self.directory, files['structure']['filename']
                )
                with open(f, mode='w') as fd:
                    fd.write(files['structure']['data'])

                self.logger.debug(
                    files['structure']['filename'] + ": " +
                    files['structure']['data']
                )

                self._initialization_node = node

                files['input'] = {}
                files['input']['filename'] = None
                files['input']['data'] = copy.deepcopy(initialization_header)

                # Find the bond & angle types as needed for shake/rattle
                P = node.parameters.current_values_to_dict(
                    context=seamm.flowchart_variables._data
                )

                shake = self.shake_fix(P, eex)
                if shake != '':
                    extras['shake'] = shake

            else:

                P = node.parameters.current_values_to_dict(
                    context=seamm.flowchart_variables._data
                )

                if 'run_control' not in P or 'Until' not in P['run_control']:

                    history_nodes.append(node)

                else:

                    if len(history_nodes) > 0:  # if imcccc

                        files = self._prepare_input(
                            files,
                            nodes=history_nodes,
                            read_restart=False,
                            write_restart=True,
                            extras=extras
                        )

                        files = self._execute_single_sim(files, np=np)

                        self.analyze(nodes=history_nodes)

                        self._trajectory = []

                    iteration = 0

                    extras['nsteps'] = 666

                    while True:

                        extras['nsteps'] = round(1.5 * extras['nsteps'])

                        control_properties = {}

                        for prp in P['control_properties']:
                            k = prp[0]
                            control_properties[k] = {
                                'accuracy': float(prp[1][0].strip('%')),
                                'units': prp[1][1],
                                'enough_accuracy': False
                            }

                        P = node.parameters.current_values_to_dict(
                            context=seamm.flowchart_variables._data
                        )

                        files = self._prepare_input(
                            files,
                            nodes=node,
                            iteration=iteration,
                            read_restart=True,
                            write_restart=True,
                            extras=extras
                        )

                        files = self._execute_single_sim(files, np=np)

                        # Analyze the results
                        analysis = self.analyze(nodes=node)
                        node_id = node._id[1]
                        for prp, v in analysis[node_id].items():
                            if v['short_production'] is False:
                                if v['few_neff'] is False:

                                    accuracy = control_properties[prp][
                                        'accuracy'] / 100
                                    dof = v['n_sample']
                                    mean = v['mean']
                                    ci = t_student.interval(
                                        0.95, dof - 1, loc=0, scale=1
                                    )
                                    interval = (ci[1] - ci[0]) * v['sem']
                                    print(abs(interval / mean))
                                    if abs(interval / mean) < accuracy:
                                        control_properties[prp][
                                            'enough_accuracy'] = True

                        enough_acc = [
                            v['enough_accuracy']
                            for prp, v in control_properties.items()
                        ]
                        if all(enough_acc) is True:
                            history_nodes = []

                            initialization_header, eex = self._get_node_input(
                                node=self._initialization_node,
                                extras={'read_data': False}
                            )
                            files['input']['data'] = copy.deepcopy(
                                initialization_header
                            )
                            files['input']['data'].append(
                                "read_restart       %s" %
                                (os.path.join(files['restart']['filename']))
                            )
                            self._trajectory = []
                            break

                        iteration = iteration + 1

            node = node.next()

        if len(history_nodes) == 0:

            node_initialization = self.subflowchart.get_node('1').next()

            if node_initialization is None:
                raise TypeError(
                    "The initial node in a LAMMPS workflow should be an ",
                    "initialization node"
                )

            if isinstance(
                node_initialization, lammps_step.Initialization
            ) is False:
                raise TypeError(
                    "The initial node in a LAMMPS workflow should be an ",
                    "initialization node"
                )

            base = 'lammps_substep_%s_iter_%d' % (
                node_initialization._id[1], 0
            )

            restart = base + '.restart.*'
            dump = base + '.dump.*'
            input_file = base + '.dat'
            new_input_data = []
            new_input_data.append('run          0')
            new_input_data.append(f'write_restart          {restart}')

            new_input_data.append(
                f'write_dump all custom  {dump} id '
                'xu yu zu modify flush yes sort id'
            )

            files['input']['filename'] = input_file
            files['input']['data'] += new_input_data
            files['input']['data'] = '\n'.join(files['input']['data'])

            self.logger.debug(
                files['input']['filename'] + ':\n' + files['input']['data']
            )

            files = self._execute_single_sim(files, np=np)

        if len(history_nodes) > 0:

            files = self._prepare_input(
                files,
                nodes=history_nodes,
                read_restart=False,
                write_restart=True,
                extras=extras
            )

            files = self._execute_single_sim(files, np=np)

            self.analyze(nodes=history_nodes)

            self._trajectory = []

        self.read_dump(os.path.join(self.directory, files['dump']['filename']))

        return next_node

    def _execute_single_sim(self, files, np=1, return_files=None):
        """
        Step #1: Dump input file
        Step #2: Execute input file
        Step #3: Dump stderr
        Step #4: Dump output files
        """
        tmpdict = {}
        for k, v in files.items():

            if v['filename'] is None or v['data'] is None:
                continue

            filename = os.path.join(self.directory, v['filename'])
            mode = "w" if type(v['data']) is str else "wb"
            with open(filename, mode=mode) as fd:
                fd.write(v['data'])

            tmpdict[v['filename']] = v['data']

        return_files = [
            'summary_*.txt',
            'trajectory_*.seamm_trj',
            '*.restart.*',
            '*.dump.*',
            '*.log',
            'log.cite'
        ]  # yapf: disable

        local = seamm.ExecLocal()

        if np > 1:
            cmd = [
                self.options['mpi_exe'], '-np',
                str(np), self.options['lammps_mpi'], '-in',
                files['input']['filename']
            ]
        else:
            cmd = [
                self.options['lammps_serial'], '-in',
                files['input']['filename']
            ]

        result = local.run(cmd=cmd, files=tmpdict, return_files=return_files)

        if result is None:
            self.logger.error('There was an error running LAMMPS')
            return None

        self.logger.debug('\n' + pprint.pformat(result))

        self.logger.debug('stdout:\n' + result['stdout'])

        f = os.path.join(self.directory, 'stdout.txt')
        with open(f, mode='w') as fd:
            fd.write(result['stdout'])

        # Add the citations, getting the version from stdout and any citations
        if 'log.cite' in result:
            self._add_lammps_citations(
                result['stdout'], cite=result['log.cite']['data']
            )
        else:
            self._add_lammps_citations(result['stdout'])

        if result['stderr'] != '':
            self.logger.warning('stderr:\n' + result['stderr'])
            f = os.path.join(self.directory, 'sstderr.txt')
            with open(f, mode='w') as fd:
                fd.write(result['stderr'])

        for filename in result['files']:
            f = os.path.join(self.directory, filename)
            mode = "wb" if type(result[filename]['data']) is bytes else "w"
            with open(f, mode=mode) as fd:
                if result[filename]['data'] is not None:
                    fd.write(result[filename]['data'])
                else:
                    fd.write(result[filename]['exception'])

        base = os.path.basename(files['input']['filename']).split('.')[0]
        restart_file = base + '.restart.*'
        dump_file = base + '.dump.*'
        filename = os.path.join(self.directory, restart_file)
        restart_filenames = glob.glob(filename)

        # Probably the step didn't run
        if len(restart_filenames) == 0:
            raise FileNotFoundError(
                'Lammps_step: could not find any file with the pattern %s' %
                (filename)
            )

        run_lengths = []

        for f in restart_filenames:
            try:
                pre, ext = os.path.splitext(f)
                ext = int(ext.strip('.'))
            except ValueError:
                raise Exception(
                    'Lammps_step: could not extract run length from %s' % f
                )
            run_lengths.append(ext)

            last_snapshot = str(max(run_lengths))

        restart_file = restart_file.replace('*', last_snapshot)
        dump_file = dump_file.replace('*', last_snapshot)
        files['restart'] = {}
        files['restart']['filename'] = restart_file
        files['dump'] = {}
        files['dump']['filename'] = dump_file
        files['dump']['data'] = None

        filename = os.path.join(self.directory, restart_file)
        with open(filename, mode='rb') as fd:
            files['restart']['data'] = fd.read()

        files['input'] = {}
        files['input']['filename'] = None
        initialization_header, eex = self._get_node_input(
            node=self._initialization_node, extras={'read_data': False}
        )
        files['input']['data'] = copy.deepcopy(initialization_header)

        return files

    def _prepare_input(
        self,
        files,
        nodes=None,
        iteration=0,
        read_restart=False,
        write_restart=False,
        extras=None
    ):

        if isinstance(nodes, list) is False:
            node_ids = [nodes._id[1]]
            new_input_data = self._get_node_input(node=nodes, extras=extras)

        else:
            node_ids = []
            new_input_data = []
            for n in nodes:
                node_ids.append(n._id[1])
                new_input_data += self._get_node_input(node=n, extras=extras)

        base = 'lammps_substep_%s_iter_%d' % ('_'.join(node_ids), iteration)

        input_file = base + '.dat'
        restart = base + '.restart.*'
        dump = base + '.dump.*'

        if read_restart:
            new_input_data.insert(
                0, 'read_restart          %s' % (files['restart']['filename'])
            )

        if write_restart:
            new_input_data.append(f'write_restart          {restart}')

        new_input_data.append(
            f'write_dump all custom  {dump} id '
            'xu yu zu modify flush yes sort id'
        )

        files['input']['data'] += new_input_data

        files['input']['filename'] = input_file
        files['input']['data'] = '\n'.join(files['input']['data'])
        self.logger.debug(
            files['input']['filename'] + ':\n' + files['input']['data']
        )

        return files

    def _get_node_input(self, node=None, extras=None):

        try:
            ret = node.get_input(extras=extras)
        except Exception as e:
            print(
                'Error running LAMMPS flowchart: {} in {}'.format(
                    str(e), str(node)
                )
            )
            self.logger.critical(
                'Error running LAMMPS flowchart: {} in {}'.format(
                    str(e), str(node)
                )
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
        lines.append(
            'Structure file for LAMMPS generated by a MolSSI flowchart'
        )
        lines.append('{:10d} atoms'.format(eex['n_atoms']))
        lines.append('{:10d} atom types'.format(eex['n_atom_types']))
        if 'n_bonds' in eex and eex['n_bonds'] > 0:
            lines.append('{:10d} bonds'.format(eex['n_bonds']))
            lines.append('{:10d} bond types'.format(eex['n_bond_types']))
        if 'n_angles' in eex and eex['n_angles'] > 0:
            lines.append('{:10d} angles'.format(eex['n_angles']))
            lines.append('{:10d} angle types'.format(eex['n_angle_types']))
        if 'n_torsions' in eex and eex['n_torsions'] > 0:
            lines.append('{:10d} dihedrals'.format(eex['n_torsions']))
            lines.append(
                '{:10d} dihedral types'.format(eex['n_torsion_types'])
            )
        if 'n_oops' in eex and eex['n_oops'] > 0:
            lines.append('{:10d} impropers'.format(eex['n_oops']))
            lines.append('{:10d} improper types'.format(eex['n_oop_types']))

        # Find the box limits
        periodicity = eex['periodicity']
        if periodicity == 3:
            a, b, c, alpha, beta, gamma = eex['cell']
            lx, ly, lz, xy, xz, yz = LAMMPS.cell_to_box(
                a, b, c, alpha, beta, gamma
            )

            lines.append('{} {} xlo xhi'.format(0.0, lx))
            lines.append('{} {} ylo yhi'.format(0.0, ly))
            lines.append('{} {} zlo zhi'.format(0.0, lz))

            xy = xy if abs(xy) > 1.0e-06 else 0.0
            xz = xz if abs(xy) > 1.0e-06 else 0.0
            yz = yz if abs(xy) > 1.0e-06 else 0.0

            if triclinic or xy > 0.0 or xz > 0.0 or yz > 0.0:
                lines.append('{} {} {} xy xz yz'.format(xy, xz, yz))
        else:
            x, y, z, index = eex['atoms'][0]
            xlo = xhi = x
            ylo = yhi = y
            zlo = zhi = z
            for x, y, z, index in eex['atoms']:
                xlo = x if x < xlo else xlo
                xhi = x if x > xhi else xlo
                ylo = y if y < ylo else ylo
                yhi = y if y > yhi else ylo
                zlo = z if z < zlo else zlo
                zhi = z if z > zhi else zlo

            # Some extra space....
            xlo -= 10.0
            xhi += 10.0
            ylo -= 10.0
            yhi += 10.0
            zlo -= 10.0
            zhi += 10.0

            lines.append('{} {} xlo xhi'.format(xlo, xhi))
            lines.append('{} {} ylo yhi'.format(ylo, yhi))
            lines.append('{} {} zlo zhi'.format(zlo, zhi))

        # the atoms and their masses, etc.
        lines.append('')
        lines.append('Atoms')
        lines.append('')

        if 'charges' in eex:
            for i, xyz_index, q in zip(
                range(1, eex['n_atoms'] + 1), eex['atoms'], eex['charges']
            ):
                x, y, z, index = xyz_index
                # The '1' is molecule ID ... should correct at some point!
                lines.append(
                    f'{i:6d}      1 {index:6d} {q:6.3f} {x:12.7f} {y:12.7f} '
                    f'{z:12.7f}'
                )
        else:
            for i, xyz_index in enumerate(eex['atoms']):
                x, y, z, index = xyz_index
                lines.append(
                    f'{i+1:6d} {index:6d} {x:12.7f} {y:12.7f} {z:12.7f}'
                )
            pass
        lines.append('')

        lines.append('Masses')
        lines.append('')
        self._data['masses'] = []
        for i, parameters in zip(
            range(1, eex['n_atom_types'] + 1), eex['masses']
        ):
            mass, itype = parameters
            lines.append('{:6d} {} # {}'.format(i, mass, itype))
            self._data['masses'].append(float(mass))

        # nonbonds
        if 'nonbond parameters' in eex:
            lines.append('')
            lines.append('Pair Coeffs')
            lines.append('')
            for i, parameters in zip(
                range(1, eex['n_atom_types'] + 1), eex['nonbond parameters']
            ):
                form, values, types, parameters_type, real_types = \
                    parameters
                if form == 'nonbond(9-6)':
                    lines.append(
                        '{:6d} {} {} # {} --> {}'.format(
                            i, values['eps'], values['rmin'], types[0],
                            real_types[0]
                        )
                    )
                else:
                    lines.append(
                        '{:6d} {} {} # {} --> {}'.format(
                            i, values['eps'], values['sigma'], types[0],
                            real_types[0]
                        )
                    )

        # bonds
        if 'n_bonds' in eex and eex['n_bonds'] > 0:
            lines.append('')
            lines.append('Bonds')
            lines.append('')
            for counter, tmp in zip(
                range(1, eex['n_bonds'] + 1), eex['bonds']
            ):
                i, j, index = tmp
                lines.append(
                    '{:6d} {:6d} {:6d} {:6d}'.format(counter, index, i, j)
                )

            lines.append('')
            lines.append('Bond Coeffs')
            lines.append('')
            for counter, parameters in zip(
                range(1, eex['n_bond_types'] + 1), eex['bond parameters']
            ):
                form, values, types, parameters_type, real_types = \
                    parameters
                if form == 'quadratic_bond':
                    # '{:6d} harmonic {} {}'
                    lines.append(
                        '{:6d} {} {}'
                        .format(counter, values['K2'], values['R0']) +
                        ' # {}-{} --> {}-{}'.format(
                            types[0], types[1], real_types[0], real_types[1]
                        )
                    )
                elif form == 'quartic_bond':
                    lines.append(
                        '{:6d} class2 {} {} {} {}'.format(
                            counter, values['R0'], values['K2'], values['K3'],
                            values['K4']
                        ) + ' # {}-{} --> {}-{}'.format(
                            types[0], types[1], real_types[0], real_types[1]
                        )
                    )

        # angles
        if 'n_angles' in eex and eex['n_angles'] > 0:
            lines.append('')
            lines.append('Angles')
            lines.append('')
            for counter, tmp in zip(
                range(1, eex['n_angles'] + 1), eex['angles']
            ):
                i, j, k, index = tmp
                lines.append(
                    '{:6d} {:6d} {:6d} {:6d} {:6d}'.format(
                        counter, index, i, j, k
                    )
                )

            lines.append('')
            lines.append('Angle Coeffs')
            lines.append('')
            for counter, parameters in zip(
                range(1, eex['n_angle_types'] + 1), eex['angle parameters']
            ):
                form, values, types, parameters_type, real_types = \
                    parameters
                if form == 'quadratic_angle':
                    # '{:6d} harmonic {} {}'
                    lines.append(
                        '{:6d} {} {}'
                        .format(counter, values['K2'], values['Theta0']) +
                        ' # {}-{}-{} --> {}-{}-{}'.format(
                            types[0], types[1], types[2], real_types[0],
                            real_types[1], real_types[2]
                        )
                    )
                elif form == 'quartic_angle':
                    lines.append(
                        '{:6d} class2 {} {} {} {}'.format(
                            counter, values['Theta0'], values['K2'],
                            values['K3'], values['K4']
                        ) + ' # {}-{}-{} --> {}-{}-{}'.format(
                            types[0], types[1], types[2], real_types[0],
                            real_types[1], real_types[2]
                        )
                    )

            # bond-bond coefficients, which must match angles in order & number
            if 'n_bond-bond_types' in eex:
                lines.append('')
                lines.append('BondBond Coeffs')
                lines.append('')
                for counter, parameters, angles in zip(
                    range(1, eex['n_bond-bond_types'] + 1),
                    eex['bond-bond parameters'], eex['angle parameters']
                ):
                    form, values, types, parameters_type, real_types = \
                        parameters
                    angle_form = angles[0]
                    if angle_form == 'quartic_angle':
                        lines.append(
                            '{:6d} class2 {} {} {}'.format(
                                counter, values['K'], values['R10'],
                                values['R20']
                            ) + ' # {}-{}-{} --> {}-{}-{}'.format(
                                types[0], types[1], types[2], real_types[0],
                                real_types[1], real_types[2]
                            )
                        )
                    else:
                        lines.append(
                            '{:6d} skip'.format(counter) +
                            ' # {}-{}-{} --> {}-{}-{}'.format(
                                types[0], types[1], types[2], real_types[0],
                                real_types[1], real_types[2]
                            )
                        )

                # bond-angles coefficients, which must match angles in order &
                # number
                lines.append('')
                lines.append('BondAngle Coeffs')
                lines.append('')
                for counter, parameters, angles in zip(
                    range(1, eex['n_bond-angle_types'] + 1),
                    eex['bond-angle parameters'], eex['angle parameters']
                ):
                    form, values, types, parameters_type, real_types = \
                        parameters
                    angle_form = angles[0]
                    if angle_form == 'quartic_angle':
                        lines.append(
                            '{:6d} class2 {} {} {} {}'.format(
                                counter, values['K12'], values['K23'],
                                values['R10'], values['R20']
                            ) + ' # {}-{}-{} --> {}-{}-{}'.format(
                                types[0], types[1], types[2], real_types[0],
                                real_types[1], real_types[2]
                            )
                        )
                    else:
                        lines.append(
                            '{:6d} skip'.format(counter) +
                            ' # {}-{}-{} --> {}-{}-{}'.format(
                                types[0], types[1], types[2], real_types[0],
                                real_types[1], real_types[2]
                            )
                        )

        # torsions
        if 'n_torsions' in eex and eex['n_torsions'] > 0:
            lines.append('')
            lines.append('Dihedrals')
            lines.append('')
            for counter, tmp in zip(
                range(1, eex['n_torsions'] + 1), eex['torsions']
            ):
                i, j, k, l, index = tmp
                lines.append(
                    '{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}'.format(
                        counter, index, i, j, k, l
                    )
                )

            lines.append('')
            lines.append('Dihedral Coeffs')
            lines.append('')
            for counter, parameters in zip(
                range(1, eex['n_torsion_types'] + 1), eex['torsion parameters']
            ):
                form, values, types, parameters_type, real_types = \
                    parameters
                if form == 'torsion_1':
                    KPhi = values['KPhi']
                    n = values['n']
                    Phi0 = values['Phi0']

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
                        d = '-1'
                    elif float(Phi0) == 180.0:
                        d = '+1'
                    else:
                        raise RuntimeError(
                            'LAMMPS cannot handle Phi0 = {}'.format(Phi0)
                        )

                    lines.append(
                        '{:6d} harmonic {} {} {}'.format(counter, KPhi, d, n) +
                        ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1], real_types[2],
                            real_types[3]
                        )
                    )
                elif form == 'torsion_3':
                    lines.append(
                        '{:6d} class2 {} {} {} {} {} {}'.format(
                            counter, values['V1'], values['Phi0_1'],
                            values['V2'], values['Phi0_2'], values['V3'],
                            values['Phi0_3']
                        ) + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1], real_types[2],
                            real_types[3]
                        )
                    )

            # middle bond-torsion_3 coefficients, which must match torsions
            # in order & number
            if 'n_middle_bond-torsion_3_types' in eex:
                lines.append('')
                lines.append('MiddleBondTorsion Coeffs')
                lines.append('')
                for counter, parameters, torsions in zip(
                    range(1, eex['n_middle_bond-torsion_3_types'] + 1),
                    eex['middle_bond-torsion_3 parameters'],
                    eex['torsion parameters']
                ):
                    form, values, types, parameters_type, real_types = \
                        parameters
                    torsion_form = torsions[0]
                    if torsion_form == 'torsion_3':
                        lines.append(
                            '{:6d} class2 {} {} {} {}'.format(
                                counter, values['V1'], values['V2'],
                                values['V3'], values['R0']
                            ) + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )
                    else:
                        lines.append(
                            '{:6d} skip'.format(counter) +
                            ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )

                # end bond-torsion_3 coefficients, which must match torsions
                # in order & number
                lines.append('')
                lines.append('EndBondTorsion Coeffs')
                lines.append('')
                for counter, parameters, torsions in zip(
                    range(1, eex['n_end_bond-torsion_3_types'] + 1),
                    eex['end_bond-torsion_3 parameters'],
                    eex['torsion parameters']
                ):
                    form, values, types, parameters_type, real_types = \
                        parameters
                    torsion_form = torsions[0]
                    if torsion_form == 'torsion_3':
                        lines.append(
                            '{:6d} class2 {} {} {} {} {} {} {} {}'.format(
                                counter, values['V1_L'], values['V2_L'],
                                values['V3_L'], values['V1_R'], values['V2_R'],
                                values['V3_R'], values['R0_L'], values['R0_R']
                            ) + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )
                    else:
                        lines.append(
                            '{:6d} skip'.format(counter) +
                            ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )

                # angle-torsion_3 coefficients, which must match torsions
                # in order & number
                lines.append('')
                lines.append('AngleTorsion Coeffs')
                lines.append('')
                for counter, parameters, torsions in zip(
                    range(1, eex['n_angle-torsion_3_types'] + 1),
                    eex['angle-torsion_3 parameters'],
                    eex['torsion parameters']
                ):
                    form, values, types, parameters_type, real_types = \
                        parameters
                    torsion_form = torsions[0]
                    if torsion_form == 'torsion_3':
                        lines.append(
                            '{:6d} class2 {} {} {} {} {} {} {} {}'.format(
                                counter, values['V1_L'], values['V2_L'],
                                values['V3_L'], values['V1_R'], values['V2_R'],
                                values['V3_R'], values['Theta0_L'],
                                values['Theta0_R']
                            ) + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )
                    else:
                        lines.append(
                            '{:6d} skip'.format(counter) +
                            ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )

                # angle-angle-torsion_1 coefficients, which must match torsions
                # in order & number
                lines.append('')
                lines.append('AngleAngleTorsion Coeffs')
                lines.append('')
                for counter, parameters, torsions in zip(
                    range(1, eex['n_angle-angle-torsion_1_types'] + 1),
                    eex['angle-angle-torsion_1 parameters'],
                    eex['torsion parameters']
                ):
                    form, values, types, parameters_type, real_types = \
                        parameters
                    torsion_form = torsions[0]
                    if torsion_form == 'torsion_3':
                        lines.append(
                            '{:6d} class2 {} {} {}'.format(
                                counter, values['K'], values['Theta0_L'],
                                values['Theta0_R']
                            ) + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )
                    else:
                        lines.append(
                            '{:6d} skip'.format(counter) +
                            ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )

                # bond-bond_1_3 coefficients, which must match torsions
                # in order & number
                lines.append('')
                lines.append('BondBond13 Coeffs')
                lines.append('')
                for counter, parameters, torsions in zip(
                    range(1, eex['n_bond-bond_1_3_types'] + 1),
                    eex['bond-bond_1_3 parameters'], eex['torsion parameters']
                ):
                    form, values, types, parameters_type, real_types = \
                        parameters
                    torsion_form = torsions[0]
                    if torsion_form == 'torsion_3':
                        lines.append(
                            '{:6d} class2 {} {} {}'.format(
                                counter, values['K'], values['R10'],
                                values['R30']
                            ) + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )
                    else:
                        lines.append(
                            '{:6d} skip'.format(counter) +
                            ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                types[0], types[1], types[2], types[3],
                                real_types[0], real_types[1], real_types[2],
                                real_types[3]
                            )
                        )

        # out-of-planes
        if 'n_oops' in eex and eex['n_oops'] > 0:
            lines.append('')
            lines.append('Impropers')
            lines.append('')
            for counter, tmp in zip(range(1, eex['n_oops'] + 1), eex['oops']):
                i, j, k, l, index = tmp
                lines.append(
                    '{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}'.format(
                        counter, index, i, j, k, l
                    )
                )

            lines.append('')
            lines.append('Improper Coeffs')
            lines.append('')
            for counter, parameters in zip(
                range(1, eex['n_oop_types'] + 1), eex['oop parameters']
            ):
                form, values, types, parameters_type, real_types = \
                    parameters
                lines.append(
                    '{:6d} {} {}'.format(counter, values['K'], values['Chi0'])
                    + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                        types[0], types[1], types[2], types[3], real_types[0],
                        real_types[1], real_types[2], real_types[3]
                    )
                )

            # angle-angle
            if 'n_angle-angle_types' in eex:
                lines.append('')
                lines.append('AngleAngle Coeffs')
                lines.append('')
                for counter, parameters in zip(
                    range(1, eex['n_angle-angle_types'] + 1),
                    eex['angle-angle parameters']
                ):
                    form, values, types, parameters_type, real_types = \
                        parameters
                    lines.append(
                        '{:6d} {} {} {} {} {} {}'.format(
                            counter, values['K1'], values['K2'], values['K3'],
                            values['Theta10'], values['Theta20'],
                            values['Theta30']
                        ) + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1], real_types[2],
                            real_types[3]
                        )
                    )

        lines.append('')
        return lines

    def analyze(self, indent='', nodes=None, **kwargs):
        """Analyze the output of the calculation
        """
        if isinstance(nodes, list) is False:
            nodes = [nodes]

        ret = {node._id[1]: None for node in nodes}

        for node in nodes:

            for value in node.description:
                printer.important(value)
                printer.important(' ')

            id_str = '_'.join(str(e) for e in node._id)

            filenames = glob.glob(
                os.path.join(
                    self.directory, '*trajectory*' + id_str + '*.seamm_trj'
                )
            )

            if len(filenames) > 1:
                filenames.sort(
                    key=lambda x:
                    int(re.search(r'iter_\d+', x).group().split('_')[1])
                )

            if len(filenames) > 0:
                filename = filenames[-1]

            P = node.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )

            node_data = None

            if 'run_control' in P:

                if P['run_control'] == 'For a fixed length of simulated time.':
                    control_properties = lambda x: x not in [  # noqa: E731
                        'tstep'
                    ]
                    # Reset the trajectory data so doesn't carry over
                    self._trajectory = []
                else:

                    if len(P['control_properties']) == 0:
                        raise KeyError(
                            'No physical property selected for automatic',
                            'equilibration detection'
                        )

                    control_properties = [
                        prp[0] for prp in P['control_properties']
                    ]

                node_data = self.analyze_trajectory(
                    filename, control_properties=control_properties
                )
                node.analyze(data=node_data)

            ret[node._id[1]] = node_data

        return ret

    def analyze_trajectory(
        self, filename, sampling_rate=20, control_properties=None
    ):
        """Read a trajectory file and do the statistical analysis
        """

        write_html = ('html' in self.options and self.options['html'])

        rootname = os.path.splitext(filename)[0]

        results = {}

        if isinstance(
            control_properties, list
        ) and 't' not in control_properties:
            control_properties.append('t')

        # Process the trajectory data

        with open(filename, 'r') as fd:
            file_data = pandas.read_csv(
                fd,
                sep=' ',
                header=0,
                comment='!',
                usecols=control_properties,
                index_col='t'
            )
            self._trajectory.append(file_data.iloc[:-1])

        dt_fs = file_data.index[1] - file_data.index[0]
        dt = dt_fs
        data = pandas.concat(self._trajectory)
        data = data.reset_index(drop=True)
        data.index *= dt

        self.logger.debug('Columns: {}'.format(data.columns))
        self.logger.debug('  Types:\n{}'.format(data.dtypes))

        printer.normal(
            '       Analysis of ' + os.path.basename(filename) + '\n'
        )

        printer.normal(
            '                                             Std Error  '
            'Time to\n'
            '               Property           Value       of mean   '
            'convergence     tau    inefficiency\n'
            '          --------------------   ---------  ---------   '
            '-----------  --------  ------------'
        )

        # Work out the time step, rather than give the whole vector
        t = data.index
        t_units = 'fs'
        len_trj = (len(t) - 1) * dt_fs
        divisor = 1
        if len_trj >= 4000000000:
            t_units = 'ms'
            divisor = 1000000000
        elif len_trj >= 4000000:
            t_units = 'ns'
            divisor = 1000000
        elif len_trj >= 4000:
            t_units = 'ps'
            divisor = 1000
        dt /= divisor
        t_max = float((len(t) - 1) * dt)

        for column in data.columns:
            have_warning = False
            have_acf_warning = False
            y = data[column]

            self.logger.info(
                'Analyzing {}, nsamples = {}'.format(column, len(y))
            )

            # compute indices of uncorrelated timeseries using pymbar
            yy = y.to_numpy()
            conv, inefficiency, Neff_max = timeseries.detectEquilibration(yy)

            self.logger.info(
                '  converged in {} steps, inefficiency = {}, Neff_max = {}'
                .format(conv, inefficiency, Neff_max)
            )

            if numpy.isnan(inefficiency) or numpy.isnan(Neff_max):
                # Apparently didn't converge!
                printer.normal(
                    __(
                        '{column:>23s} did not converge to a stationary state',
                        column=column,
                        indent=7 * ' ',
                        wrap=False,
                        dedent=False
                    )
                )

                have_acf = False
                is_converged = False
            else:
                is_converged = True
                tau = dt_fs * (inefficiency - 1) / 2
                if tau < dt_fs / 2:
                    tau = dt_fs / 2
                t0 = conv * dt_fs
                y_t_equil = yy[conv:]
                indices = timeseries.subsampleCorrelatedData(
                    y_t_equil, g=inefficiency
                )

                if len(indices) == 0:
                    print('Problem with column ' + column)
                    print('yy')
                    print(yy)
                    print('y_t_equil')
                    print(y_t_equil)
                    print('indices')
                    print(indices)
                    continue

                y_n = y_t_equil[indices]
                n_samples = len(y_n)
                mean = y_n.mean()
                std = y_n.std()
                sem = std / sqrt(n_samples)

                # Get the autocorrelation function
                if len(y_t_equil) < 8:
                    have_acf = False
                    have_acf_warning = True
                    acf_warning = '^'
                else:
                    have_acf = True
                    acf_warning = ' '
                    nlags = 4 * int(round(inefficiency + 0.5))
                    if nlags > int(len(y_t_equil) / 2):
                        nlags = int(len(y_t_equil) / 2)
                    acf, confidence = stattools.acf(
                        y_t_equil,
                        nlags=nlags,
                        alpha=0.05,
                        fft=nlags > 16,
                        adjusted=False
                    )

                results[column] = {}
                results[column]['mean'] = mean
                results[column]['sem'] = sem
                results[column]['n_sample'] = n_samples
                results[column]['short_production'] = have_acf_warning

                # Work out units on convergence time
                conv_units = 'fs'
                t_conv = t0
                if t0 >= 1000000000:
                    conv_units = 'ms'
                    t_conv = t0 / 1000000000
                elif t0 >= 1000000:
                    conv_units = 'ns'
                    t_conv = t0 / 1000000
                elif t0 >= 1000:
                    conv_units = 'ps'
                    t_conv = t0 / 1000

                # Work out units on autocorrelation time
                tau_units = 'fs'
                t_tau = tau
                if tau >= 1000000000:
                    tau_units = 'ms'
                    t_tau = tau / 1000000000
                elif tau >= 1000000:
                    tau_units = 'ns'
                    t_tau = tau / 1000000
                elif tau >= 1000:
                    tau_units = 'ps'
                    t_tau = tau / 1000

                if n_samples < 100:
                    have_warning = True
                    warn = '*'
                else:
                    warn = ' '

                results[column]['few_neff'] = have_warning

                printer.normal(
                    __(
                        '{column:>23s} = {value:9.3f} ± {stderr:7.3f}{warn}'
                        ' {t0:8.2f} {conv_units} {tau:8.1f} {tau_units}{acf} '
                        '{inefficiency:9.1f}',
                        column=column,
                        value=mean,
                        stderr=sem,
                        warn=warn,
                        t0=t_conv,
                        conv_units=conv_units,
                        tau=t_tau,
                        tau_units=tau_units,
                        acf=acf_warning,
                        inefficiency=inefficiency,
                        indent=7 * ' ',
                        wrap=False,
                        dedent=False
                    )
                )

            # Create graphs of the property
            figure = self.create_figure(
                module_path=(self.__module__.split('.')[0], 'seamm'),
                template='line.graph_template',
                title=LAMMPS.display_title[column]
            )

            # The autocorrelation function
            if have_acf:
                plot_acf = figure.add_plot('acf')

                dt_acf = float(dt_fs)
                t_acf_units = 'fs'
                len_acf = (len(acf) - 1) * dt_fs
                if len_acf >= 2000000000:
                    t_acf_units = 'ms'
                    dt_acf /= 1000000000
                elif len_acf >= 2000000:
                    t_acf_units = 'ns'
                    dt_acf /= 1000000
                elif len_acf >= 2000:
                    t_acf_units = 'ps'
                    dt_acf /= 1000

                x_acf_axis = plot_acf.add_axis(
                    'x', label='Time ({})'.format(t_acf_units)
                )
                y_acf_axis = plot_acf.add_axis(
                    'y', label='acf', anchor=x_acf_axis
                )
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
                    name='fit',
                    x0=0,
                    dx=dt_acf,
                    xlabel='t',
                    xunits=t_acf_units,
                    y=fit,
                    ylabel='fit',
                    yunits='',
                    color='gray'
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
                    tmp += dt_acf

                plot_acf.add_trace(
                    x_axis=x_acf_axis,
                    y_axis=y_acf_axis,
                    name='stderr',
                    x=t_acf + t_acf[::-1],
                    xlabel='t',
                    xunits=t_acf_units,
                    y=yplus + yminus[::-1],
                    ylabel='stderr',
                    yunits=LAMMPS.display_units[column],
                    showlegend='false',
                    color='rgba(211,211,211,0.5)',
                    fill='toself',
                )

                # And the acf plot last
                plot_acf.add_trace(
                    x_axis=x_acf_axis,
                    y_axis=y_acf_axis,
                    name='acf',
                    x0=0,
                    dx=dt_acf,
                    xlabel='t',
                    xunits=t_acf_units,
                    y=list(acf),
                    ylabel='acf',
                    yunits='',
                    color='red'
                )

            # The property data over the trajectory
            y = list(data[column])

            plot = figure.add_plot('trj')

            ylabel = LAMMPS.display_title[column]
            if LAMMPS.display_units[column] != '':
                ylabel += ' ({})'.format(LAMMPS.display_units[column])

            x_axis = plot.add_axis('x', label='Time ({})'.format(t_units))
            y_axis = plot.add_axis('y', label=ylabel, anchor=x_axis)
            x_axis.anchor = y_axis

            # Add the trajectory, error band and median value in that order so
            # stack in a nice order.

            # Add the trajectory
            plot.add_trace(
                x_axis=x_axis,
                y_axis=y_axis,
                name=column,
                x0=0,
                dx=dt,
                xlabel='t',
                xunits=t_units,
                y=list(y),
                ylabel=column,
                yunits=LAMMPS.display_units[column],
                color='#4dbd74'
            )

            if is_converged:
                # the partly transparent error band
                t_min = t0 / divisor
                plot.add_trace(
                    x_axis=x_axis,
                    y_axis=y_axis,
                    name='sem',
                    x=[t_min, t_max, t_max, t_min],
                    xlabel='t',
                    xunits=t_units,
                    y=[mean + sem, mean + sem, mean - sem, mean - sem],
                    ylabel='sem',
                    yunits=LAMMPS.display_units[column],
                    showlegend='false',
                    color='rgba(211,211,211,0.5)',
                    fill='toself',
                )

                # and finally the median value so it is on top
                plot.add_trace(
                    x_axis=x_axis,
                    y_axis=y_axis,
                    name='average',
                    x=[t_min, t_max],
                    xlabel='t',
                    xunits=t_units,
                    y=[mean, mean],
                    ylabel='average',
                    yunits=LAMMPS.display_units[column],
                    color='black'
                )

            if have_acf:
                figure.grid_plots('trj - acf')
            else:
                figure.grid_plots('trj')
            figure.dump('{}_{}.graph'.format(rootname, column))

            if write_html:
                figure.template = 'line.html_template'
                figure.dump('{}_{}.html'.format(rootname, column))

        if have_warning or have_acf_warning:
            printer.normal('\n')
        if have_warning:
            printer.normal(
                __(
                    '          * this property has less than 100 independent '
                    'samples, so may not be accurate.',
                    wrap=False,
                    dedent=False
                )
            )

        if have_acf_warning:
            printer.normal(
                __(
                    '          ^ there are not enough samples after '
                    'equilibration to plot the ACF.',
                    wrap=False,
                    dedent=False
                )
            )

        return results

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
        if P['rigid_waters']:
            # waters = seamm_util.water_models.Water.find_waters(data.structure)  # noqa: E501

            waters = []

            if len(waters) > 0:
                atoms = []
                for i, j, k in waters:
                    atoms.append(i)
                    atoms.append(j)
                    atoms.append(k)
                    if 'n_bonds' in eex and eex['n_bonds'] > 0:
                        for i, j, index in eex['bonds']:
                            if i in atoms and j in atoms:
                                bond_types[index] = 1
                    if 'n_angles' in eex and eex['n_angles'] > 0:
                        for i, j, k, index in eex['angles']:
                            if i in atoms and j in atoms and k in atoms:
                                angle_types[index] = 1

        # Fixing bond lengths of X-H bonds...
        if 'n_bonds' in eex and eex['n_bonds'] > 0:
            fix_bonds = P['fix_XH_bond_lengths']
            elements = eex['elements']
            if fix_bonds == 'CH':
                for i, j, index in eex['bonds']:
                    if (
                        (elements[i] == 'C' and elements[j] == 'H') or
                        (elements[i] == 'H' and elements[j] == 'C')
                    ):
                        bond_types[index] = 1
            elif fix_bonds == 'all':
                for i, j, index in eex['bonds']:
                    if elements[i] == 'H' or elements[j] == 'H':
                        bond_types[index] = 1

        # And the result is ....
        if len(bond_types) > 0:
            result = 'fix                 {} all rattle 0.001 20 1000 b '
            for bond_type in bond_types.keys():
                result += ' ' + str(bond_type)
            if len(angle_types) > 0:
                result += ' a '
                for angle_type in angle_types.keys():
                    result += ' ' + str(angle_type)
        else:
            result = ''

        return result

    def read_dump(self, dumpfile):
        """Read the LAMMPS dumpfile and update the system.

        Parameters
        ----------
        dumpfile : str
            The filename (or path) to the dumpfile.
        """
        self.logger.info("Reading dump file '{}'".format(dumpfile))

        system = self.get_variable('_system')
        periodicity = system.periodicity
        n_atoms = system.n_atoms()

        section = ''
        section_lines = []
        xs = []
        ys = []
        zs = []
        with open(dumpfile, 'r') as fd:
            lineno = 0
            for line in fd:
                line = line.strip()
                lineno += 1
                if lineno == 1:
                    if line[0:5] != 'ITEM:':
                        raise RuntimeError(
                            "Error reading dump file '" + dumpfile + "': The "
                            "first line is incorrect! (" + line + ")"
                        )
                    section = line[6:].strip()
                    section_lines = []
                    self.logger.debug('   section = ' + section)
                    continue

                if line[0:5] == 'ITEM:':
                    # end a section
                    self.logger.debug(
                        "  processing section '{}'".format(section)
                    )
                    if 'BOX BOUNDS' in section:
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

                            cell = (
                                xhi - xlo, yhi - ylo, zhi - zlo, 90, 90, 90
                            )
                    elif section == 'NUMBER OF ATOMS':
                        if int(section_lines[0]) != n_atoms:
                            raise RuntimeError(
                                'Number of atoms has changed! {} to {}'.format(
                                    n_atoms, section_lines[0]
                                )
                            )
                    elif 'ATOMS' in section:
                        for tmp in section_lines:
                            id, x, y, z = tmp.split()
                            xs.append(float(x))
                            ys.append(float(y))
                            zs.append(float(z))
                    section = line[6:].strip()
                    section_lines = []
                else:
                    section_lines.append(line)

        # Clean up the last section
        if 'ATOMS' in section:
            self.logger.debug("  processing section '{}'".format(section))
            self.logger.debug('  handling the atoms')
            for tmp in section_lines:
                id, x, y, z = tmp.split()
                xs.append(float(x))
                ys.append(float(y))
                zs.append(float(z))

        if periodicity == 3:
            system.cell.set_cell(cell)
        system.atoms['x'][0:] = xs
        system.atoms['y'][0:] = ys
        system.atoms['z'][0:] = zs

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
            raw=self._bibliography['PLIMPTON19951'],
            alias='lammps-jcp',
            module='lammps_step',
            level=1,
            note='The principle LAMMPS citation.'
        )

        # And the citation to the LAMMPS code itself
        lines = text.splitlines()

        if len(lines) == 0:
            return

        line = lines[0]
        tmp = line.split(' ')
        if len(tmp) != 4:
            self.logger.info(f"Cannot get LAMMPS version: '{line}'")
            return

        month = tmp[2]
        year = tmp[3].rstrip(')')
        version = ' '.join(tmp[1:])

        try:
            template = string.Template(self._bibliography['lammps'])

            citation = template.substitute(
                month=month, version=version, year=year
            )

            self.references.cite(
                raw=citation,
                alias='lammps-exe',
                module='lammps_step',
                level=1,
                note='The principle citation for the LAMMPS executable.'
            )

        except Exception as e:
            printer.important(f'Exception in citation {type(e)}: {e}')
            printer.important(traceback.format_exc())

        # If there is a log.cite file, process it
        if cite is not None:
            self.logger.warning('log.cite\n' + cite + '\n')
            bibliography = {}
            tmp = bibtexparser.loads(cite).entries_dict
            writer = bibtexparser.bwriter.BibTexWriter()
            for key, data in tmp.items():
                self.logger.warning(f'      {key}')
                bibliography[key] = writer._entry_to_bibtex(data)
            self.logger.warning(
                'Bibliography\n' + pprint.pformat(bibliography)
            )

            for entry in bibliography:
                if entry.lower() in ('commment',):
                    continue
                self.references.cite(
                    raw=bibliography[entry],
                    alias=entry,
                    module='lammps_step',
                    level=2,
                    note='LAMMPS citations from log.cite.'
                )
