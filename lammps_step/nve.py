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
        titles = 'tstep t T P Etot Eke Epe Emol Epair'

        nfixes = 0
        lines.append('')
        lines.append('#     NVE dynamics')
        lines.append('')
        lines.append('reset_timestep      0')
        lines.append('timestep            {}'.format(timestep))
        lines.append('thermo_style        custom {}'.format(thermo_properties))
        lines.append('thermo              {}'.format(int(nsteps/100)))
        nfixes += 1
        lines.append('fix                 {} all nve'.format(nfixes))
        # summary output written 10 times during run so we can see progress
        nevery = 10
        nfreq = int(nsteps / 10)
        nrepeat = int(nfreq / nevery)
        nfreq = nevery * nrepeat
        nfixes += 1
        lines.append('fix                 {} '.format(nfixes) +
                     'all  ave/time '
                     '{} {} {} {} file summary_nve_{}.txt'.format(
                         nevery, nrepeat, nfreq, properties,
                         '_'.join(str(e) for e in self._id)))
        # instantaneous output written for averaging
        if self.sampling_method != 'none':
            sampling = self.sampling.to('fs').magnitude
            nevery = round(sampling / timestep)
            nfreq = int(nsteps / nevery)
            nrepeat = 1
            nfreq = nevery * nrepeat
            nfixes += 1
            lines.append(
                'fix                 {} '.format(nfixes) +
                'all ave/time ' +
                "{} {} {} {} off 2 title2 '{}' file trajectory_nve_{}.txt"
                .format(
                    nevery, nrepeat, nfreq, properties, titles,
                    '_'.join(str(e) for e in self._id))
            )
                         
        lines.append('run                 {}'.format(nsteps))
        lines.append('')
        for fix in range(1, nfixes+1):
            lines.append('unfix               {}'.format(fix))

        return lines
