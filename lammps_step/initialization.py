# -*- coding: utf-8 -*-
"""A single-point initialization in LAMMPS"""

import molssi_workflow
import forcefield
import lammps_step
import logging
import molssi_util.smiles
import pprint

logger = logging.getLogger(__name__)

kspace_methods = {
    'automatic': '',
    'none': 'none',
    'Ewald summation': 'ewald {accuracy}',
    'Particle-particle particle-mesh (PPPM)': 'pppm {accuracy}',
    'PPPM for few charged atoms': 'pppm/cg {accuracy} {smallq}',
    'PPPM staggered mesh': 'pppm/stagger {accuracy}',
    'PPPM plus dispersion term': 'pppm/disp {accuracy}',
    'Multilevel summation method (MSM)': 'msm {accuracy}',
    'MSM for few charged atoms': 'msm/cg {accuracy} {smallq}'
}

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


class Initialization(molssi_workflow.Node):
    structures = {
        'current': '',
        'initial': '',
        'other': '',
    }

    def __init__(self,
                 workflow=None,
                 title='Initialization',
                 extension=None):
        """Initialize the node"""

        logger.debug('Creating Initialization {}'.format(self))

        super().__init__(
            workflow=workflow,
            title=title,
            extension=extension)

        self.description = 'Initialization of a LAMMPS calculation'

        self.structure = None
        self.cutoff = 10.0
        self.kspace_method = 'automatic'
        self.kspace_accuracy = 1.0e-05
        self.kspace_smallq = 1.0e-05
        self.charged_atom_fraction_cutoff = 0.1
        self.ewald_atom_cutoff = 3000
        self.use_tail_correction = True
        self.shift_nonbond = False

    @property
    def kspace_methods(self):
        """The list of avilable methods"""
        return list(kspace_methods)

    def get_input(self):
        """Get the input for the initialization of LAMMPS"""

        structure = molssi_workflow.data.structure
        logger.debug('Structure in LAMMPS initialization:\n' +
                     pprint.pformat(structure))

        # Atom-type if necessary
        ff = molssi_workflow.data.forcefield
        ff_name = ff.current_forcefield
        atoms = structure['atoms']
        n_atoms = len(atoms['elements'])
        if 'atom_types' in atoms and ff_name in atoms['atom_types']:
            atom_types = atoms['atom_types'][ff_name]
        else:
            smiles = molssi_util.smiles.from_molssi(structure)
            logger.debug('Atom typing -- smiles = ' + smiles)
            ff_assigner = forcefield.FFAssigner(ff)
            atom_types = ff_assigner.assign(smiles, add_hydrogens=False)
            logger.info('Atom types: ' + ', '.join(atom_types))
            if 'atom_types' not in atoms:
                atoms['atom_types'] = {}
            atoms['atom_types'][ff_name] = atom_types

        # and get the energy expression. This creates the charges on the
        # atoms as a side-effect.
        eex = ff.energy_expression(structure, style='LAMMPS')
        logger.debug('energy expression:\n' + pprint.pformat(eex))

        # Determine if we have any charges, and if so, if they are sparse
        charges = atoms['charges'][ff_name]
        n_charged_atoms = 0
        for charge in charges:
            if abs(charge) > self.kspace_smallq:
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
        elif periodicity == 3:
            lines.append('boundary            p p p')
            tail_correction = 'yes' if self.use_tail_correction and \
                              not self.shift_nonbond else 'no'
        else:
            raise RuntimeError('The LAMMPS step can only handle 0-'
                               ' or 3-D periodicity at the moment!')

        lines.append('atom_style          full')
        lines.append('atom_modify         sort 0 0.0')
        lines.append('newton              on')
        lines.append('')
        lines.append('#    define the style of forcefield')
        lines.append('')

        terms = ff.ff['terms']

        logging.debug(
            'LAMMPS initialization, terms = \n' + pprint.pformat(terms))

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
                "Can't handle nonbond term {} yet!".format(nonbond_term))

        shift = 'yes' if self.shift_nonbond else 'no'
        if self.kspace_method == 'automatic':
            if periodicity == 3:
                kspace_style = ''
                if n_charged_atoms == 0:
                    pair_style = pair_style_base
                else:
                    pair_style = pair_style_base + '/coul/long'
                    if n_atoms < self.ewald_atom_cutoff:
                        kspace_style = 'ewald {}'.format(self.kspace_accuracy)
                    elif fraction_charged_atoms < self.charged_atom_fraction_cutoff:  # nopep8
                        kspace_style = 'pppm/cg {} {}'.format(
                            self.kspace_accuracy,
                            self.charged_atom_fraction_cutoff)
                    else:
                        kspace_style = 'pppm {}'.format(self.kspace_accuracy)
                lines.append('pair_style          {} {}'.format(
                    pair_style, self.cutoff))
                lines.append(
                    'pair_modify         mix ' + mixing +
                    ' tail {} shift {}'.format(tail_correction, shift))
                if kspace_style != '':
                    lines.append('kspace_style        ' + kspace_style)
            else:
                if n_charged_atoms == 0:
                    pair_style = pair_style_base
                elif fraction_charged_atoms < \
                    self.charged_atom_fraction_cutoff \
                    and pair_style_base in msm_pair_styles:  # nopep8
                    pair_style = pair_style_base + '/coul/msm'
                else:
                    pair_style = pair_style_base + '/coul/cut'
                lines.append('pair_style          {} {}'.format(
                    pair_style, self.cutoff))
                lines.append('pair_modify         mix ' + mixing +
                             ' shift {}'.format(shift))
        else:
            if periodicity == 3:
                kspace_style = ''
                if n_charged_atoms == 0 or self.kspace_style == 'none':
                    pair_style = pair_style_base
                elif fraction_charged_atoms < self.charged_atom_fraction_cutoff:  # nopep8
                    pair_style = pair_style_base + '/coul/long'
                    kspace_style = kspace_methods[self.kspace_method].\
                        format(accuracy=self.accuracy,
                               smallq=self.kspace_smallq)
                lines.append('pair_style          {} {}'.format(
                    pair_style, self.cutoff))
                lines.append(
                    'pair_modify         mix ' + mixing +
                    ' tail {} shift {}'.format(tail_correction, shift))
                if kspace_style != '':
                    lines.append('kspace_style        ' + kspace_style)
            else:
                if n_charged_atoms == 0:
                    pair_style = pair_style_base
                else:
                    pair_style = pair_style_base + '/coul/cut'
                lines.append('pair_style          {} {}'.format(
                    pair_style, self.cutoff))
                lines.append('pair_modify         mix ' + mixing +
                             ' shift {}'.format(shift))
                if 'msm' in kspace_methods[self.kspace_method]:
                    kspace_style = kspace_methods[self.kspace_method].\
                                   format(accuracy=self.accuracy,
                                          smallq=self.kspace_smallq)
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
                dihedral_style = lammps_step.dihedral_style[terms['torsion'][
                    0]]
                lines.append('dihedral_style      ' + dihedral_style)
            else:
                line = 'dihedral_style      hybrid'
                for term in terms['torsion']:
                    line += ' ' + lammps_step.dihedral_style[term]
                lines.append(line)
        if 'out-of-plane' in terms and eex['n_oops'] > 0:
            if len(terms['out-of-plane']) == 1:
                improper_style = lammps_step.improper_style[terms[
                    'out-of-plane'][0]]
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
            lines.append('variable            {var} equal {var}'.format(
                var=variable))

        return (lines, eex)
