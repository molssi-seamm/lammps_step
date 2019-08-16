# -*- coding: utf-8 -*-

"""Minimization step in LAMMPS"""

from seamm import data
from seamm_util import ureg, Q_, units_class  # noqa: F401
import lammps_step
import logging

logger = logging.getLogger(__name__)


class Minimization(lammps_step.Energy):

    def __init__(self, flowchart=None, title='Minimization', extension=None):
        """Initialize the node"""

        logger.debug('Creating Minimization {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = 'Minimization step in LAMMPS'

        self.convergence = 'normal'
        self.etol_method = 'is'
        self.etol = 1.0e-06
        self.etol_variable = ''
        self.ftol_method = 'is'
        self.ftol = Q_(0.1, 'kcal/mol/Å')
        self.ftol_variable = ''
        self.maxiters_method = 'is'
        self.maxiters = 10000
        self.maxiters_variable = ''
        self.maxevals_method = 'is'
        self.maxevals = 30000
        self.maxevals_variable = ''

    def get_input(self):
        """Get the input for a minimization in LAMMPS"""

        n_atoms = len(data.structure['atoms']['elements'])
        nDOF = 3 * n_atoms

        lines = []

        etol = 0.0
        ftol = None
        maxiters = None
        maxevals = None
        if self.convergence == 'normal':
            etol = 0.0
            ftol = 0.1
            maxiters = 3 * nDOF
            maxevals = 5 * maxiters
        elif self.convergence == 'tight':
            etol = 0.0
            ftol = 0.01
            maxiters = 2 * 3 * nDOF
            maxevals = 5 * maxiters
        elif self.convergence == 'loose':
            etol = 0.0
            ftol = 1
            maxiters = nDOF
            maxevals = 5 * maxiters
        elif self.convergence == 'crude':
            etol = 0.0
            ftol = 10
            maxiters = int(nDOF / 10)
            maxevals = 2 * maxiters
        elif 'energy' in self.convergence:
            # energy tolerance
            if self.etol_method == 'is':
                etol = self.etol
            else:
                raise RuntimeError(
                    'Variable handling not implemented for etol'
                )
        elif 'forces' in self.convergence:
            # force tolerance
            if self.ftol_method == 'is':
                ftol = self.ftol.to('kcal/mol/Å').magnitude
            else:
                raise RuntimeError(
                    'Variable handling not implemented for ftol'
                )
        else:
            raise RuntimeError(
                'Variable handling not implemented for convergence'
            )

        # maximum number of iterations
        if maxiters is None:
            if self.maxiters_method == 'default':
                maxiters = 3 * nDOF
            elif self.maxiters_method == 'is':
                maxiters = self.maxiters
            else:
                raise RuntimeError(
                    'Variable handling not implemented for maxiters'
                )

        # maximum number of energy evaluations
        if maxevals is None:
            if self.maxevals_method == 'default':
                maxevals = 3 * nDOF
            elif self.maxevals_method == 'is':
                maxevals = self.maxevals
            else:
                raise RuntimeError(
                    'Variable handling not implemented for maxevals'
                )

        thermo_properties = (
            'press etotal ke pe ebond '
            'eangle edihed eimp evdwl etail ecoul elong'
        )

        lines.append('')
        lines.append('#     Minimization')
        lines.append('')
        lines.append('thermo_style        custom {}'.format(thermo_properties))
        lines.append('thermo              {}'.format(100))
        lines.append(
            'minimize            {} {} {} {}'.format(
                etol, ftol, maxiters, maxevals
            )
        )

        return lines
