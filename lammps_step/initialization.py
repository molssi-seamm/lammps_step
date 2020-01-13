# -*- coding: utf-8 -*-

"""A single-point initialization in LAMMPS"""

import seamm_ff_util
import lammps_step
import logging
import seamm
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __
import seamm_util.smiles
import pprint

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter('lammps')

msm_pair_styles = ['born', 'buck', '', 'lj/charmm', 'lj/cut']

thermo_variables = [
    'step', 'elapsed', 'elaplong', 'dt', 'time', 'cpu', 'tpcpu', 'spcpu',
    'cpuremain', 'part', 'timeremain', 'atoms', 'temp', 'press', 'pe', 'ke',
    'etotal', 'enthalpy', 'evdwl', 'ecoul', 'epair', 'ebond', 'eangle',
    'edihed', 'eimp', 'emol', 'elong', 'etail', 'vol', 'density', 'lx', 'ly',
    'lz', 'xlo', 'xhi', 'ylo', 'yhi', 'zlo', 'zhi', 'xy', 'xz', 'yz', 'xlat',
    'ylat', 'zlat', 'bonds', 'angles', 'dihedrals', 'impropers', 'pxx', 'pyy',
    'pzz', 'pxy', 'pxz', 'pyz', 'fmax', 'fnorm', 'nbuild', 'ndanger', 'cella',
    'cellb', 'cellc', 'cellalpha', 'cellbeta', 'cellgamma'
]


class Initialization(seamm.Node):

    def __init__(self, flowchart=None, title='Initialization', extension=None):
        """Initialize the node"""

        logger.debug('Creating Initialization {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = []
        self.parameters = lammps_step.InitializationParameters()

    @property
    def header(self):
        """A printable header for this section of output"""
        return (
            'Step {}: {}'.format(
                '.'.join(str(e) for e in self._id), self.title
            )
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

    @property
    def kspace_methods(self):
        """The list of avilable methods"""
        return list(lammps_step.kspace_methods)

    def description_text(self, P=None):
        """Return a short description of this step.

        Return a nicely formatted string describing what this step will
        do.

        Keyword arguments:
            P: a dictionary of parameter values, which may be variables
                or final values. If None, then the parameters values will
                be used as is.
        """

        if not P:
            P = self.parameters.values_to_dict()

        text = 'Initialize the calculation with a cutoff of {cutoff} Å'
        if P['shift_nonbond']:
            text += ', shifting the nonbond energies to 0 at the cutoff'
        text += '. If the system is periodic'
        if P['kspace_method'][0] == '$':
            text += " use the variable '{method}' to determine whether "
            text += "and how to accelerate the k-space summations."
        elif P['kspace_method'] == 'none':
            text += ' no k-space acceleration method will be used.'
        elif P['kspace_method'] == 'automatic':
            text += ' the best k-space acceleration method for the '
            text += ' molecular system will be chosen.'
        else:
            text += ' the {method} for k-space acceleration will be used.'

        if P['kspace_method'] != 'none':
            text += ' The accuracy goal is {kspace_accuracy}.'

        return self.header + '\n' + __(text, **P, indent=4 * ' ').__str__()

    def get_input(self):
        """Get the input for the initialization of LAMMPS"""

        self.description = []
        self.description.append('   ' + self.header)

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Fix some things
        P['cutoff'] = P['cutoff'].to('angstrom').magnitude

        # Get the structure
        structure = seamm.data.structure
        logger.debug(
            'Structure in LAMMPS initialization:\n' +
            pprint.pformat(structure)
        )

        # And atom-type if necessary
        ff = seamm.data.forcefield
        ff_name = ff.current_forcefield
        atoms = structure['atoms']
        n_atoms = len(atoms['elements'])
        if 'atom_types' in atoms and ff_name in atoms['atom_types']:
            atom_types = atoms['atom_types'][ff_name]
        else:
            smiles = seamm_util.smiles.from_seamm(structure)
            logger.debug('Atom typing -- smiles = ' + smiles)
            ff_assigner = seamm_ff_util.FFAssigner(ff)
            atom_types = ff_assigner.assign(smiles, add_hydrogens=False)
            logger.info('Atom types: ' + ', '.join(atom_types))
            if 'atom_types' not in atoms:
                atoms['atom_types'] = {}
            atoms['atom_types'][ff_name] = atom_types

        # Get the energy expression. This creates the charges on the
        # atoms as a side-effect.
        eex = ff.energy_expression(structure, style='LAMMPS')
        logger.debug('energy expression:\n' + pprint.pformat(eex))

        # Determine if we have any charges, and if so, if they are sparse
        charges = atoms['charges'][ff_name]
        n_charged_atoms = 0
        smallq = float(P['kspace_smallq'])
        for charge in charges:
            if abs(charge) > smallq:
                n_charged_atoms += 1
        fraction_charged_atoms = n_charged_atoms / n_atoms

        lines = []
        lines.append('')
        lines.append('#     initialization of LAMMPS')
        lines.append('')
        lines.append('units               real')

        periodicity = structure['periodicity']
        if periodicity == 0:
            lines.append('boundary            s s s')
            tail_correction = 'no'
            string = 'Setup for a molecular (non-periodic) system.'
        elif periodicity == 3:
            lines.append('boundary            p p p')
            tail_correction = 'yes' if P['tail_correction'] and \
                              not P['shift_nonbond'] else 'no'
            string = 'Setup for a periodic (crystalline or fluid) system.'
        else:
            raise RuntimeError(
                'The LAMMPS step can only handle 0-'
                ' or 3-D periodicity at the moment!'
            )

        lines.append('atom_style          full')
        lines.append('atom_modify         sort 0 0.0')
        lines.append('newton              on')
        lines.append('')
        lines.append('#    define the style of forcefield')
        lines.append('')

        terms = ff.ff['terms']

        logging.debug(
            'LAMMPS initialization, terms = \n' + pprint.pformat(terms)
        )

        # control of nonbonds
        nonbond_term = None
        if 'pair' in terms:
            if len(terms['pair']) != 1:
                raise RuntimeError('Cannot handle multiple nonbond terms yet!')
            nonbond_term = terms['pair'][0]

        if nonbond_term is None:
            pprint.pprint(terms)
            raise RuntimeError("Cannot find nonbond term in forcefield!")

        if nonbond_term == 'nonbond(9-6)':
            pair_style_base = 'lj/class2'
            mixing = 'sixthpower'
        else:
            raise RuntimeError(
                "Can't handle nonbond term {} yet!".format(nonbond_term)
            )

        shift = 'yes' if P['shift_nonbond'] else 'no'
        if P['kspace_method'] == 'automatic':
            if periodicity == 3:
                kspace_style = ''
                if n_charged_atoms == 0:
                    pair_style = pair_style_base
                    string += (
                        ' The nonbonded interactions will be evaluated using '
                        'a cutoff of {cutoff} Å. Since there are no charges '
                        'on the atoms, no long-range coulomb method will be '
                        'used.'
                    )
                else:
                    string += (
                        ' The nonbonded interactions will be evaluated using '
                        'a cutoff of {cutoff} Å, with the long-range terms '
                    )
                    pair_style = pair_style_base + '/coul/long'
                    if n_atoms < P['ewald_atom_cutoff']:
                        kspace_style = 'ewald {}'.format(P['kspace_accuracy'])
                        string += (
                            'using the Ewald summation method with '
                            'an accuracy of {kspace_accuracy}.'
                        )
                    elif fraction_charged_atoms < \
                            P['charged_atom_fraction_cutoff']:
                        kspace_style = 'pppm/cg {} {}'.format(
                            P['kspace_accuracy'],
                            P['charged_atom_fraction_cutoff']
                        )
                        string += (
                            'using the PPPM method optimized for few '
                            'atoms with charges, with '
                            'an accuracy of {kspace_accuracy}.'
                        )
                    else:
                        kspace_style = 'pppm {}'.format(P['kspace_accuracy'])
                        string += (
                            'using the PPPM method with '
                            'an accuracy of {kspace_accuracy}.'
                        )
                lines.append(
                    'pair_style          {} {}'.format(
                        pair_style, P['cutoff']
                    )
                )
                lines.append(
                    'pair_modify         mix ' + mixing +
                    ' tail {} shift {}'.format(tail_correction, shift)
                )
                if shift:
                    string += (
                        ' The van der Waals terms will be shifted '
                        'to zero energy at the cutoff distance.'
                    )
                if tail_correction:
                    string += (
                        ' A long-range correction for the '
                        'van der Waals terms will be added.'
                    )
                if kspace_style != '':
                    lines.append('kspace_style        ' + kspace_style)
            else:
                kspace_style = ''
                if n_charged_atoms == 0:
                    pair_style = pair_style_base
                    string += (
                        ' The nonbonded interactions will be evaluated using '
                        'a simple cutoff of {cutoff} Å. Since there are no '
                        'charges on the atoms, no long-range coulomb method '
                        'will be used.'
                    )
                elif (
                    n_atoms > P['msm_atom_cutoff'] and
                    pair_style_base in msm_pair_styles
                ):
                    pair_style = pair_style_base + '/coul/msm'
                    string += (
                        'The nonbonded interactions will be handled with '
                        ' a cutoff of {cutoff} Å.'
                    )
                    if fraction_charged_atoms < \
                       P['charged_atom_fraction_cutoff']:
                        kspace_style = (
                            'msm/cg {kspace_accuracy} {kspace_smallq}'
                        )
                        string += (
                            ' The MSM method will be used to handle '
                            'the longer range coulombic interactions, using '
                            'the approach tuned for systems with few charges.'
                            'The accuracy goal is {kspace_accuracy}.'
                        )
                    else:
                        kspace_style = 'msm {kspace_accuracy}'
                        string += (
                            ' The MSM method will be used to handle '
                            'the longer range coulombic interactions.'
                            'The accuracy goal is {kspace_accuracy}.'
                        )
                else:
                    pair_style = pair_style_base + '/coul/cut'
                    string += (
                        'The nonbonded interactions will be handled with '
                        ' a simple cutoff of {cutoff} Å.'
                    )
                if shift:
                    string += (
                        ' The van der Waals terms will be shifted '
                        'to zero energy at the cutoff distance.'
                    )
                lines.append(
                    'pair_style          {} {}'.format(
                        pair_style, P['cutoff']
                    )
                )
                lines.append(
                    'pair_modify         mix ' + mixing +
                    ' shift {}'.format(shift)
                )
                if kspace_style != '':
                    lines.append(
                        'kspace_style        ' + kspace_style.format(**P)
                    )
            self.description.append(__(string, indent=7 * ' ', **P))
        else:
            if periodicity == 3:
                kspace_style = ''
                if n_charged_atoms == 0 or P['kspace_style'] == 'none':
                    pair_style = pair_style_base
                elif fraction_charged_atoms < \
                        P['charged_atom_fraction_cutoff']:
                    pair_style = pair_style_base + '/coul/long'
                    kspace_style = (
                        lammps_step.kspace_methods[P['kspace_method']].format(
                            **P
                        )
                    )
                lines.append(
                    'pair_style          {} {}'.format(
                        pair_style, P['cutoff']
                    )
                )
                lines.append(
                    'pair_modify         mix ' + mixing +
                    ' tail {} shift {}'.format(tail_correction, shift)
                )
                if kspace_style != '':
                    lines.append('kspace_style        ' + kspace_style)
            else:
                if n_charged_atoms == 0:
                    pair_style = pair_style_base
                else:
                    pair_style = pair_style_base + '/coul/cut'
                lines.append(
                    'pair_style          {} {}'.format(
                        pair_style, P['cutoff']
                    )
                )
                lines.append(
                    'pair_modify         mix ' + mixing +
                    ' shift {}'.format(shift)
                )
                if 'msm' in lammps_step.kspace_methods[P['kspace_method']]:
                    kspace_style = (
                        lammps_step.kspace_methods[P['kspace_method']].format(
                            **P
                        )
                    )
                    lines.append('kspace_style        ' + kspace_style)

        if 'bond' in terms and eex['n_bonds'] > 0:
            if len(terms['bond']) == 1:
                bond_style = lammps_step.bond_style[terms['bond'][0]]
                lines.append('bond_style          ' + bond_style)
            else:
                line = 'bond_style          hybrid'
                for term in terms['bond']:
                    line += ' ' + lammps_step.bond_style[term]
                lines.append(line)
        if 'angle' in terms and eex['n_angles'] > 0:
            if len(terms['angle']) == 1:
                angle_style = lammps_step.angle_style[terms['angle'][0]]
                lines.append('angle_style         ' + angle_style)
            else:
                line = 'angle_style         hybrid'
                for term in terms['angle']:
                    line += ' ' + lammps_step.angle_style[term]
                lines.append(line)
        if 'torsion' in terms and eex['n_torsions'] > 0:
            if len(terms['torsion']) == 1:
                dihedral_style = lammps_step.dihedral_style[terms['torsion'][0]
                                                           ]  # noqa: E501
                lines.append('dihedral_style      ' + dihedral_style)
            else:
                line = 'dihedral_style      hybrid'
                for term in terms['torsion']:
                    line += ' ' + lammps_step.dihedral_style[term]
                lines.append(line)
        if 'out-of-plane' in terms and eex['n_oops'] > 0:
            if len(terms['out-of-plane']) == 1:
                improper_style = lammps_step.improper_style[
                    terms['out-of-plane'][0]]
                lines.append('improper_style      ' + improper_style)
            else:
                line = 'improper_style      hybrid'
                for term in terms['out-of-plane']:
                    line += ' ' + lammps_step.improper_style[term]
                lines.append(line)

        lines.append('')
        lines.append('read_data           structure.dat')

        # Set up standard variables
        for variable in thermo_variables:
            lines.append(
                'variable            {var} equal {var}'.format(var=variable)
            )

        return (lines, eex)
