# -*- coding: utf-8 -*-

"""NVE (microcanonical) dynamics in LAMMPS"""

import lammps_step
import logging
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter('lammps')


class NVE(lammps_step.Energy):

    def __init__(self, flowchart=None, title='NVE dynamics', extension=None):
        """Initialize the node"""

        logger.debug('Creating NVE {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = 'NVE dynamics step in LAMMPS'

        logger.debug("NVE.init() creating NVE_Parameters object")

        self.parameters = lammps_step.NVE_Parameters()

    def description_text(self):
        """Create the text description of what this step will do.
        """

        text = (
            "{time} of microcanonical (NVE) dynamics using a "
            "timestep of {timestep}. The trajectory will be "
            "sampled every {sampling}."
        )

        return text

    def get_input(self, extras=None):
        """Get the input for an NVE dynamics run in LAMMPS"""

        self.description = []
        self.description.append(__(self.header, indent=3 * ' '))

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        timestep, P['timestep'] = self.timestep(P['timestep'])

        if extras is not None and 'nsteps' in extras:
            nsteps = extras['nsteps']
        else:
            time = P['time'].to('fs').magnitude
            nsteps = round(time / timestep)

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = '{:~P}'.format(PP[key])

        self.description.append(
            __(self.description_text(), **PP, indent=7 * ' ')
        )

        # time = lammps_step.to_lammps_units(P['time'], quantity='time')
        # nsteps = round(time / timestep)

        thermo_properties = (
            'time temp press etotal ke pe ebond '
            'eangle edihed eimp evdwl etail ecoul elong'
        )
        properties = 'v_time v_temp v_press v_etotal v_ke v_pe v_emol v_epair'
        title2 = 'tstep t T P Etot Eke Epe Emol Epair'

        lines = []
        nfixes = 0
        lines.append('')
        lines.append('#     NVE dynamics')
        lines.append('')
        lines.append('reset_timestep      0')
        lines.append('timestep            {}'.format(timestep))
        lines.append('thermo_style        custom {}'.format(thermo_properties))
        lines.append('thermo              {}'.format(int(nsteps / 100)))
        nfixes += 1
        lines.append('fix                 {} all nve'.format(nfixes))
        # summary output written 10 times during run so we can see progress
        nevery = 10
        nfreq = int(nsteps / 10)
        nrepeat = int(nfreq / nevery)
        nfreq = nevery * nrepeat
        nfixes += 1
        lines.append(
            'fix                 {} '.format(nfixes) + 'all  ave/time '
            '{} {} {} {} file summary_nve_{}_iter_0.txt'.format(
                nevery, nrepeat, nfreq, properties,
                '_'.join(str(e) for e in self._id)
            )
        )
        # instantaneous output written for averaging
        if P['sampling'] == 'none':
            self.description.append(
                __(
                    "The run will be {nsteps:n} steps of dynamics.",
                    nsteps=nsteps,
                    indent=7 * ' '
                )
            )
        else:
            sampling = lammps_step.to_lammps_units(
                P['sampling'], quantity='time'
            )
            nevery = round(sampling / timestep)
            nfreq = int(nsteps / nevery)
            nrepeat = 1
            nfreq = nevery * nrepeat
            nfixes += 1
            title1 = (
                '!MolSSI trajectory 1.0 LAMMPS, NVE {} steps of {} fs'.format(
                    int(nsteps / nevery), timestep * nevery
                )
            )
            lines.append(
                (
                    "fix                 {} all ave/time {} {} {} {} off 2 "
                    "title1 '{}' title2 '{}' file "
                    "trajectory_nve_{}_iter_0.seamm_trj"
                ).format(
                    nfixes, nevery, nrepeat, nfreq, properties, title1, title2,
                    '_'.join(str(e) for e in self._id)
                )
            )
            self.description.append(
                __(
                    (
                        "The run will be {nsteps:,d} steps of dynamics "
                        "sampled every {nevery:n} steps."
                    ),
                    nsteps=nsteps,
                    nevery=nevery,
                    indent=7 * ' '
                )
            )

        if extras is not None and 'shake' in extras:
            nfixes += 1
            lines.append(extras['shake'].format(nfixes))

        lines.append('')
        lines.append('run                 {}'.format(nsteps))
        lines.append('')

        for fix in range(1, nfixes + 1):
            lines.append('unfix               {}'.format(fix))
        lines.append('')

        return lines

    def timestep(self, value):
        """Get the timestep in the correct units.

        This handles the 'normal', 'accurate' and 'coarse' values,
        which depend on the mass in an empirical fashion.

        Parameters
        ----------
        value : str or Pint quantity
            The desired timestep, which may be a Pint quantity with units
            or one of 'notmal', 'accurate but slow', or 'coarse but fast'

        Returns
        -------
        timestep : float
            The magnitude of the time step in the appropriate LAMMPS units
        """
        masses = self.parent._data['masses']
        min_mass = min(masses)

        # These are based on masses as a proxy for vibrational frequencies
        if min_mass < 10:
            factor = 1
        elif min_mass < 50:
            factor = 2
        else:
            factor = 4

        if value == 'normal':
            timestep = 1.0 * factor
            value = Q_(timestep, ureg.fs)
        elif value == 'accurate but slow':
            timestep = 0.5 * factor
            value = Q_(timestep, ureg.fs)
        elif value == 'coarse but fast':
            timestep = 2.0 * factor
            value = Q_(timestep, ureg.fs)

        timestep = lammps_step.to_lammps_units(value, quantity='time')

        return (timestep, value)
