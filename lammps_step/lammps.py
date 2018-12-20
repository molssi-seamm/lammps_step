# -*- coding: utf-8 -*-
"""A node or step for LAMMPS in a workflow"""

import molssi_workflow
from molssi_workflow import units, Q_, units_class, data  # nopep8
import lammps_step
import logging
import math
import os
import os.path
import pprint

logger = logging.getLogger(__name__)

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

lammps_units = {
    'real': {
        '[mass]': 'g/mol',
        '[distance]': 'Å',
        '[time]': 'fs',
        '[length] ** 2 * [mass] / [substance] / [time] ** 2': 'kcal/mol',
        '[length] ** 2 * [mass] / [time] ** 2': 'kcal/mol',
        '[length] / [time]': 'Å/fs',
        '[length] * [mass] / [substance] / [time] ** 2': 'kcal/mol/Å',
        '[length] * [mass] / [time] ** 2': 'kcal/mol/Å',
        '[temperature]': 'K',
        '[mass] / [length] / [time] ** 2': 'bar',
        '[mass] / [length] / [time]': 'poise',
        '[current] * [time]': 'e',
        '[current] * [length] * [time]': 'e*Å',
        '[length] * [mass] / [current] / [time] ** 3': 'V/Å',
        '[mass] / [length] ** 3': 'g/mL'
    },
    'metal': {
        '[mass]': 'g/mol',
        '[distance]': 'Å',
        '[time]': 'ps',
        '[length] ** 2 * [mass] / [substance] / [time] ** 2': 'eV',
        '[length] ** 2 * [mass] / [time] ** 2': 'eV',
        '[length] / [time]': 'Å/ps',
        '[length] * [mass] / [substance] / [time] ** 2': 'eV/Å',
        '[length] * [mass] / [time] ** 2': 'eV/Å',
        '[temperature]': 'K',
        '[mass] / [length] / [time] ** 2': 'atm',
        '[mass] / [length] / [time]': 'poise',
        '[current] * [time]': 'e',
        '[current] * [length] * [time]': 'e*Å',
        '[length] * [mass] / [current] / [time] ** 3': 'V/Å',
        '[mass] / [length] ** 3': 'g/mL'
    }
}


def cosine(degrees):
    return math.cos(math.radians(degrees))


class LAMMPS(molssi_workflow.Node):
    def __init__(self,
                 workflow=None,
                 namespace='org.molssi.workflow.lammps',
                 extension=None):
        '''Setup the main LAMMPS step

        Keyword arguments:
        '''
        logger.debug('Creating LAMMPS {}'.format(self))

        self.lammps_workflow = molssi_workflow.Workflow(
            parent=self, name='LAMMPS',
            namespace=namespace)
        self.lammps_units = 'real'
        self._data = {}

        super().__init__(
            workflow=workflow,
            title='LAMMPS',
            extension=extension)

    def set_id(self, node_id):
        """Set the id for node to a given tuple"""
        self._id = node_id

        # and set our subnodes
        self.lammps_workflow.set_ids(self._id)

        return self.next()

    def describe(self):
        """Run a LAMMPS simulation
        """

        super().describe()

        self.lammps_workflow.root_directory = self.workflow.root_directory

        # Get the first real node
        node = self.lammps_workflow.get_node('1').next()

        while node is not None:
            node.describe()
            node = node.next()

        return self.next()

    def run(self):
        """Run a LAMMPS simulation
        """

        if data.structure is None:
            logger.error('LAMMPS run(): there is no structure!')
            raise RuntimeError('LAMMPS run(): there is no structure!')

        self.lammps_workflow.root_directory = self.workflow.root_directory

        # Get the first real node
        node = self.lammps_workflow.get_node('1').next()

        input_data = []
        while node is not None:
            if isinstance(node, lammps_step.Initialization):
                lines, eex = node.get_input()
                input_data += lines
            else:
                input_data += node.get_input()
            node = node.next()

        files = {'molssi.dat': '\n'.join(input_data)}
        logger.info('molssi.dat:\n' + files['molssi.dat'])

        # Get the structure file from the eex
        files['structure.dat'] = '\n'.join(self.structure_data(eex))
        # logger.log(0, 'structure.dat:\n' + files['structure.dat'])
        logger.info('structure.dat:\n' + files['structure.dat'])

        os.makedirs(self.directory, exist_ok=True)
        for filename in files:
            with open(os.path.join(self.directory, filename), mode='w') as fd:
                fd.write(files[filename])

        local = molssi_workflow.ExecLocal()
        return_files = ['summary_*.txt', 'trajectory_*.txt']
        result = local.run(
            cmd=['lammps_omp', '-sf', 'omp', '-in', 'molssi.dat'],  # nopep8
            files=files,
            return_files=return_files)

        if result is None:
            logger.error('There was an error running LAMMPS')
            return None

        logger.debug('\n' + pprint.pformat(result))

        logger.debug('stdout:\n' + result['stdout'])
        with open(os.path.join(self.directory, 'stdout.txt'), mode='w') as fd:
            fd.write(result['stdout'])

        if result['stderr'] != '':
            logger.warning('stderr:\n' + result['stderr'])
            with open(os.path.join(self.directory, 'stderr.txt'),
                      mode='w') as fd:
                fd.write(result['stderr'])

        for filename in result['files']:
            with open(os.path.join(self.directory, filename), mode='w') as fd:
                if result[filename]['data'] is not None:
                    fd.write(result[filename]['data'])
                else:
                    fd.write(result[filename]['exception'])

        return super().run()

    def structure_data(self, eex, triclinic=False):
        """Create the LAMMPS structure file from the energy expression"""
        lines = []
        lines.append(
            'Structure file for LAMMPS generated by a MolSSI workflow')
        lines.append('{:10d} atoms'.format(eex['n_atoms']))
        lines.append('{:10d} atom types'.format(eex['n_atom_types']))
        if eex['n_bonds'] > 0:
            lines.append('{:10d} bonds'.format(eex['n_bonds']))
            lines.append('{:10d} bond types'.format(eex['n_bond_types']))
        if eex['n_angles'] > 0:
            lines.append('{:10d} angles'.format(eex['n_angles']))
            lines.append('{:10d} angle types'.format(eex['n_angle_types']))
        if eex['n_torsions'] > 0:
            lines.append('{:10d} dihedrals'.format(eex['n_torsions']))
            lines.append('{:10d} dihedral types'.format(
                eex['n_torsion_types']))
        if eex['n_oops'] > 0:
            lines.append('{:10d} impropers'.format(eex['n_oops']))
            lines.append('{:10d} improper types'.format(eex['n_oop_types']))

        # Find the box limits
        periodicity = eex['periodicity']
        if periodicity == 3:
            a, b, c, alpha, beta, gamma = eex['cell']
            lx = a
            xy = b * cosine(gamma)
            xz = c * cosine(beta)
            ly = math.sqrt(b**2 - xy**2)
            yz = (b * c * cosine(alpha) - xy * xz) / ly
            lz = math.sqrt(c**2 - xz**2 - yz**2)

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

        for i, xyz_index, q in zip(range(1, eex['n_atoms'] + 1), eex['atoms'],
                                   eex['charges']):
            x, y, z, index = xyz_index
            lines.append(
                '{:6d} {:6d} {:6d} {:6.3f} {:12.7f} {:12.7f} {:12.7f}'.format(
                    i, 1, index, q, x, y, z))
        lines.append('')

        lines.append('Masses')
        lines.append('')
        for i, parameters in zip(range(1, eex['n_atom_types'] + 1),
                                 eex['masses']):
            mass, itype = parameters
            lines.append('{:6d} {} # {}'.format(i, mass, itype))

        # nonbonds
        lines.append('')
        lines.append('Pair Coeffs')
        lines.append('')
        for i, parameters in zip(range(1, eex['n_atom_types'] + 1),
                                 eex['nonbond parameters']):
            form, values, types, parameters_type, real_types = \
                parameters
            lines.append('{:6d} {} {} # {} --> {}'.format(
                i, values['eps'], values['r'], types[0], real_types[0]))

        # bonds
        if eex['n_bonds'] > 0:
            lines.append('')
            lines.append('Bonds')
            lines.append('')
            for counter, tmp in zip(
                    range(1, eex['n_bonds'] + 1), eex['bonds']):
                i, j, index = tmp
                lines.append('{:6d} {:6d} {:6d} {:6d}'.format(
                    counter, index, i, j))

            lines.append('')
            lines.append('Bond Coeffs')
            lines.append('')
            for counter, parameters in zip(
                    range(1, eex['n_bond_types'] + 1), eex['bond parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                if form == 'quadratic_bond':
                    lines.append('{:6d} harmonic {} {}'.format(
                        counter, values['K2'],
                        values['R0']) + ' # {}-{} --> {}-{}'.format(
                            types[0], types[1], real_types[0], real_types[1]))
                elif form == 'quartic_bond':
                    lines.append('{:6d} class2 {} {} {} {}'.format(
                        counter, values['R0'], values['K2'], values['K3'],
                        values['K4']) + ' # {}-{} --> {}-{}'.format(
                            types[0], types[1], real_types[0], real_types[1]))

        # angles
        if eex['n_angles'] > 0:
            lines.append('')
            lines.append('Angles')
            lines.append('')
            for counter, tmp in zip(
                    range(1, eex['n_angles'] + 1), eex['angles']):
                i, j, k, index = tmp
                lines.append('{:6d} {:6d} {:6d} {:6d} {:6d}'.format(
                    counter, index, i, j, k))

            lines.append('')
            lines.append('Angle Coeffs')
            lines.append('')
            for counter, parameters in zip(
                    range(1, eex['n_angle_types'] + 1), eex['angle parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                if form == 'quadratic_angle':
                    lines.append('{:6d} harmonic {} {}'.format(
                        counter, values['K2'],
                        values['Theta0']) + ' # {}-{}-{} --> {}-{}-{}'.format(
                            types[0], types[1], types[2],
                            real_types[0], real_types[1], real_types[2]))
                elif form == 'quartic_angle':
                    lines.append('{:6d} class2 {} {} {} {}'.format(
                        counter, values['Theta0'], values['K2'], values['K3'],
                        values['K4']) + ' # {}-{}-{} --> {}-{}-{}'.format(
                            types[0], types[1], types[2],
                            real_types[0], real_types[1], real_types[2]))

            # bond-bond coefficients, which must match angles in order & number
            lines.append('')
            lines.append('BondBond Coeffs')
            lines.append('')
            for counter, parameters, angles in zip(
                    range(1, eex['n_bond-bond_types'] + 1),
                    eex['bond-bond parameters'],
                    eex['angle parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                angle_form = angles[0]
                if angle_form == 'quartic_angle':
                    lines.append('{:6d} class2 {} {} {}'.format(
                        counter, values['K'], values['R10'], values['R20'])
                                 + ' # {}-{}-{} --> {}-{}-{}'.format(
                                     types[0], types[1], types[2],
                                     real_types[0], real_types[1], real_types[2]))
                else:
                    lines.append('{:6d} skip'.format(counter)
                                 + ' # {}-{}-{} --> {}-{}-{}'.format(
                                     types[0], types[1], types[2],
                                     real_types[0], real_types[1], real_types[2]))

            # bond-angles coefficients, which must match angles in order & number
            lines.append('')
            lines.append('BondAngle Coeffs')
            lines.append('')
            for counter, parameters, angles in zip(
                    range(1, eex['n_bond-angle_types'] + 1),
                    eex['bond-angle parameters'],
                    eex['angle parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                angle_form = angles[0]
                if angle_form == 'quartic_angle':
                    lines.append('{:6d} class2 {} {} {} {}'.format(
                        counter, values['K12'], values['K23'],
                        values['R10'], values['R20'])
                                 + ' # {}-{}-{} --> {}-{}-{}'.format(
                                     types[0], types[1], types[2],
                                     real_types[0], real_types[1], real_types[2]))
                else:
                    lines.append('{:6d} skip'.format(counter)
                                 + ' # {}-{}-{} --> {}-{}-{}'.format(
                                     types[0], types[1], types[2],
                                 real_types[0], real_types[1], real_types[2]))

        # torsions
        if eex['n_torsions'] > 0:
            lines.append('')
            lines.append('Dihedrals')
            lines.append('')
            for counter, tmp in zip(
                    range(1, eex['n_torsions'] + 1), eex['torsions']):
                i, j, k, l, index = tmp
                lines.append('{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}'.format(
                    counter, index, i, j, k, l))

            lines.append('')
            lines.append('Dihedral Coeffs')
            lines.append('')
            for counter, parameters in zip(
                    range(1, eex['n_torsion_types'] + 1),
                    eex['torsion parameters']):
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
                            'LAMMPS cannot handle Phi0 = {}'.format(Phi0))

                    lines.append('{:6d} harmonic {} {} {}'.format(
                        counter, KPhi, d, n) +
                        ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1], real_types[2],
                            real_types[3]))
                elif form == 'torsion_3':
                    lines.append('{:6d} class2 {} {} {} {} {} {}'.format(
                        counter, values['V1'], values['Phi0_1'],
                        values['V2'], values['Phi0_2'],
                        values['V3'], values['Phi0_3']) +
                        ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1], real_types[2],
                            real_types[3]))

            # middle bond-torsion_3 coefficients, which must match torsions
            # in order & number
            lines.append('')
            lines.append('MiddleBondTorsion Coeffs')
            lines.append('')
            for counter, parameters, torsions in zip(
                    range(1, eex['n_middle_bond-torsion_3_types'] + 1),
                    eex['middle_bond-torsion_3 parameters'],
                    eex['torsion parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                torsion_form = torsions[0]
                if torsion_form == 'torsion_3':
                    lines.append('{:6d} class2 {} {} {} {}'.format(
                        counter, values['V1'], values['V2'],
                        values['V3'], values['R0'])
                                 + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                     types[0], types[1], types[2], types[3],
                                     real_types[0], real_types[1],
                                     real_types[2], real_types[3]))
                else:
                    lines.append('{:6d} skip'.format(counter)
                                 + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                                     types[0], types[1], types[2], types[3],
                                     real_types[0], real_types[1],
                                     real_types[2], real_types[3]))

            # end bond-torsion_3 coefficients, which must match torsions
            # in order & number
            lines.append('')
            lines.append('EndBondTorsion Coeffs')
            lines.append('')
            for counter, parameters, torsions in zip(
                    range(1, eex['n_end_bond-torsion_3_types'] + 1),
                    eex['end_bond-torsion_3 parameters'],
                    eex['torsion parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                torsion_form = torsions[0]
                if torsion_form == 'torsion_3':
                    lines.append(
                        '{:6d} class2 {} {} {} {} {} {} {} {}'.format(
                            counter,
                            values['V1_L'], values['V2_L'], values['V3_L'],
                            values['V1_R'], values['V2_R'], values['V3_R'],
                            values['R0_L'], values['R0_R'])
                        + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1],
                            real_types[2], real_types[3]))
                else:
                    lines.append(
                        '{:6d} skip'.format(counter)
                        + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1],
                            real_types[2], real_types[3]))

            # angle-torsion_3 coefficients, which must match torsions
            # in order & number
            lines.append('')
            lines.append('AngleTorsion Coeffs')
            lines.append('')
            for counter, parameters, torsions in zip(
                    range(1, eex['n_angle-torsion_3_types'] + 1),
                    eex['angle-torsion_3 parameters'],
                    eex['torsion parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                torsion_form = torsions[0]
                if torsion_form == 'torsion_3':
                    lines.append(
                        '{:6d} class2 {} {} {} {} {} {} {} {}'.format(
                            counter,
                            values['V1_L'], values['V2_L'], values['V3_L'],
                            values['V1_R'], values['V2_R'], values['V3_R'],
                            values['Theta0_L'], values['Theta0_R'])
                        + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1],
                            real_types[2], real_types[3]))
                else:
                    lines.append(
                        '{:6d} skip'.format(counter)
                        + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1],
                            real_types[2], real_types[3]))

            # angle-angle-torsion_1 coefficients, which must match torsions
            # in order & number
            lines.append('')
            lines.append('AngleAngleTorsion Coeffs')
            lines.append('')
            for counter, parameters, torsions in zip(
                    range(1, eex['n_angle-angle-torsion_1_types'] + 1),
                    eex['angle-angle-torsion_1 parameters'],
                    eex['torsion parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                torsion_form = torsions[0]
                if torsion_form == 'torsion_3':
                    lines.append(
                        '{:6d} class2 {} {} {}'.format(
                            counter,
                            values['K'],
                            values['Theta0_L'], values['Theta0_R'])
                        + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1],
                            real_types[2], real_types[3]))
                else:
                    lines.append(
                        '{:6d} skip'.format(counter)
                        + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1],
                            real_types[2], real_types[3]))

            # bond-bond_1_3 coefficients, which must match torsions
            # in order & number
            lines.append('')
            lines.append('BondBond13 Coeffs')
            lines.append('')
            for counter, parameters, torsions in zip(
                    range(1, eex['n_bond-bond_1_3_types'] + 1),
                    eex['bond-bond_1_3 parameters'],
                    eex['torsion parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                torsion_form = torsions[0]
                if torsion_form == 'torsion_3':
                    lines.append(
                        '{:6d} class2 {} {} {}'.format(
                            counter,
                            values['K'],
                            values['R10'], values['R30'])
                        + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1],
                            real_types[2], real_types[3]))
                else:
                    lines.append(
                        '{:6d} skip'.format(counter)
                        + ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                            types[0], types[1], types[2], types[3],
                            real_types[0], real_types[1],
                            real_types[2], real_types[3]))

        # out-of-planes
        if eex['n_oops'] > 0:
            lines.append('')
            lines.append('Impropers')
            lines.append('')
            for counter, tmp in zip(
                    range(1, eex['n_oops'] + 1), eex['oops']):
                i, j, k, l, index = tmp
                lines.append('{:6d} {:6d} {:6d} {:6d} {:6d} {:6d}'.format(
                    counter, index, i, j, k, l))

            lines.append('')
            lines.append('Improper Coeffs')
            lines.append('')
            for counter, parameters in zip(
                    range(1, eex['n_oop_types'] + 1),
                    eex['oop parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                lines.append(
                    '{:6d} {} {}'.format(
                        counter, values['K'], values['Chi0']) +
                    ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                        types[0], types[1], types[2], types[3],
                        real_types[0], real_types[1], real_types[2],
                        real_types[3]))

            # angle-angle
            lines.append('')
            lines.append('AngleAngle Coeffs')
            lines.append('')
            for counter, parameters in zip(
                    range(1, eex['n_angle-angle_types'] + 1),
                    eex['angle-angle parameters']):
                form, values, types, parameters_type, real_types = \
                    parameters
                lines.append(
                    '{:6d} {} {} {} {} {} {}'.format(
                        counter, values['K1'], values['K2'], values['K3'],
                        values['Theta10'], values['Theta20'], values['Theta30'],
                    ) +
                    ' # {}-{}-{}-{} --> {}-{}-{}-{}'.format(
                        types[0], types[1], types[2], types[3],
                        real_types[0], real_types[1], real_types[2],
                        real_types[3]))

        lines.append('')
        return lines

    def to_lammps_units(self, value):
        dimensionality = str(value.dimensionality)
        return value.to(lammps_units[self.lammps_units][dimensionality])

    def magnitude_in_lammps_units(self, value):
        if isinstance(value, units_class):
            return self.to_lammps_units(value).magnitude
        else:
            return value
