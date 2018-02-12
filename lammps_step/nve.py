# -*- coding: utf-8 -*-
"""NVE (microcanonical) dynamics in LAMMPS"""

from molssi_workflow import units, Q_, data  # nopep8
import lammps_step
import logging

logger = logging.getLogger(__name__)


class NVE(lammps_step.Energy):
    def __init__(self,
                 workflow=None,
                 title='NVE dynamics',
                 extension=None):
        """Initialize the node"""

        logger.debug('Creating NVE {}'.format(self))

        super().__init__(
            workflow=workflow,
            title=title,
            extension=extension)

        self.description = 'NVE dynamics step in LAMMPS'

        self.sampling_method = 'is'
        self.sampling_variable = ''
        self.sampling = Q_(20, 'fs')
        self.timestep_method = 'normal'
        self.timestep_variable = ''
        self.timestep = Q_(1, 'fs')
        self.time = Q_(100, 'ps')

    def get_input(self):
        """Get the input for an NVE dynamics run in LAMMPS"""

        lines = []

        if self.timestep == 'automatic':
            timestep = 1.0
        else:
            timestep = self.timestep.to('fs').magnitude

        time = self.time.to('fs').magnitude
        nsteps = round(time / timestep)

        thermo_properties = ('time temp press etotal ke pe ebond '
                             'eangle edihed eimp evdwl etail ecoul elong')
        properties = 'v_time v_temp v_press v_etotal v_ke v_pe v_emol v_epair'

        lines.append('')
        lines.append('#     NVE dynamics')
        lines.append('')
        lines.append('timestep            {}'.format(timestep))
        lines.append('thermo_style        custom {}'.format(thermo_properties))
        lines.append('thermo              {}'.format(int(nsteps/100)))
        lines.append('fix                 1 all nve')
        # summary output written 100 times during run so we can see progress
        nevery = 10
        nfreq = int(nsteps / 100)
        nrepeat = int(nfreq / nevery)
        nfreq = nevery * nrepeat
        lines.append('fix                 2 all  ave/time '
                     '{} {} {} {} file nve_summary.txt'.format(
                         nevery, nrepeat, nfreq, properties))
        # instantaneous output written every 10 steps for averaging
        nevery = 10
        nfreq = int(nsteps / nevery)
        nrepeat = 1
        nfreq = nevery * nrepeat
        lines.append('fix                 3 all  ave/time '
                     '{} {} {} {} file nve.txt'.format(
                         nevery, nrepeat, nfreq, properties))
        lines.append('run                 {}'.format(nsteps))
        lines.append('')
        lines.append('unfix               1')
        lines.append('unfix               2')
        lines.append('unfix               3')

        return lines
