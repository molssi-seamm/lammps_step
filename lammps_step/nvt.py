# -*- coding: utf-8 -*-

"""NVT (canonical) dynamics in LAMMPS"""

import lammps_step
import logging
import seamm
from seamm_util import ureg, Q_, units_class  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __
import random

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter('lammps')


class NVT(lammps_step.NVE):

    methods = {
        'Nose-Hoover':
            {
                'documentation': 'http://lammps.sandia.gov/doc/fix_nvt.html',
                'references': ['Shinoda', 'Tuckerman'],
            },
        'Berendsen':
            {
                'documentation':
                    'http://lammps.sandia.gov/doc/fix_temp_berendsen.html',  # noqa: E501
                'references': ['Berendsen'],
            },
        'canonical sampling, velocity rescaling (csvr)':
            {
                'documentation':
                    'http://lammps.sandia.gov/doc/fix_temp_csvr.html',  # noqa: E501
                'references': ['Bussi1'],
            },
        'canonical sampling, langevin dynamics (csld)':
            {
                'documentation':
                    'http://lammps.sandia.gov/doc/fix_temp_csvr.html',  # noqa: E501
                'references': ['Bussi2'],
            },
        'velocity rescaling':
            {
                'documentation':
                    'http://lammps.sandia.gov/doc/fix_temp_rescale.html',  # noqa: E501
                'references': [],
            },
        'Langevin':
            {
                'documentation':
                    'http://lammps.sandia.gov/doc/fix_langevin.html',
                'references': ['Schneider', 'Dunweg']
            },
    }

    references = {
        'Shinoda':
            {
                'bibtex':
                    """
                @article{PhysRevB.69.134103,
                  title = {Rapid estimation of elastic constants by molecular dynamics simulation under constant stress},
                  author = {Shinoda, Wataru and Shiga, Motoyuki and Mikami, Masuhiro},
                  journal = {Phys. Rev. B},
                  volume = {69},
                  issue = {13},
                  pages = {134103},
                  numpages = {8},
                  year = {2004},
                  month = {Apr},
                  publisher = {American Physical Society},
                  doi = {10.1103/PhysRevB.69.134103},
                  url = {https://link.aps.org/doi/10.1103/PhysRevB.69.134103}
            }"""  # noqa: E501
            },
        'Tuckerman':
            {
                'bibtex':
                    """
                @article{0305-4470-39-19-S18,
                  author={Mark E Tuckerman and José Alejandre and Roberto López-Rendón and Andrea L Jochim and Glenn J Martyna},
                  title={A Liouville-operator derived measure-preserving integrator for molecular dynamics simulations in the isothermal–isobaric ensemble},
                  journal={Journal of Physics A: Mathematical and General},
                  volume={39},
                  number={19},
                  pages={5629},
                  url={http://stacks.iop.org/0305-4470/39/i=19/a=S18},
                  year={2006},
                  abstract={The constant-pressure,
                  constant-temperature ( NPT ) molecular dynamics
                  approach is re-examined from the viewpoint of
                  deriving a new measure-preserving reversible
                  geometric integrator for the equations of
                  motion. The underlying concepts of non-Hamiltonian
                  phase-space analysis, measure-preserving integrators
                  and the symplectic property for Hamiltonian systems
                  are briefly reviewed. In addition, current
                  measure-preserving schemes for the constant-volume,
                  constant-temperature ensemble are also reviewed. A
                  new geometric integrator for the NPT method is
                  presented, is shown to preserve the correct
                  phase-space volume element and is demonstrated to
                  perform well in realistic examples. Finally, a
                  multiple time-step version of the integrator is
                  presented for treating systems with motion on
                  several time scales.}
            }"""  # noqa: E501
            },
        'Berendsen':
            {
                'bibtex':
                    """
                @article{doi:10.1063/1.448118,
                author = {H. J. C. Berendsen and J. P. M. Postma and W. F. van Gunsteren and A. DiNola and J. R. Haak},
                title = {Molecular dynamics with coupling to an external bath},
                journal = {The Journal of Chemical Physics},
                volume = {81},
                number = {8},
                pages = {3684-3690},
                year = {1984},
                doi = {10.1063/1.448118},
                URL = {https://doi.org/10.1063/1.448118},
                eprint = {https://doi.org/10.1063/1.448118}
            }"""  # noqa: E501
            },
        'Bussi1':
            {
                'bibtex':
                    """
                @article{doi:10.1063/1.2408420,
                author = {Giovanni Bussi and Davide Donadio and Michele Parrinello},
                title = {Canonical sampling through velocity rescaling},
                journal = {The Journal of Chemical Physics},
                volume = {126},
                number = {1},
                pages = {014101},
                year = {2007},
                doi = {10.1063/1.2408420},
                URL = {https://doi.org/10.1063/1.2408420},
                eprint = {https://doi.org/10.1063/1.2408420}
            }"""  # noqa: E501
            },
        'Bussi2':
            {
                'bibtex':
                    """
                @article{PhysRevE.75.056707,
                title = {Accurate sampling using Langevin dynamics},
                author = {Bussi, Giovanni and Parrinello, Michele},
                journal = {Phys. Rev. E},
                volume = {75},
                issue = {5},
                pages = {056707},
                numpages = {7},
                year = {2007},
                month = {May},
                publisher = {American Physical Society},
                doi = {10.1103/PhysRevE.75.056707},
                url = {https://link.aps.org/doi/10.1103/PhysRevE.75.056707}
            }"""
            },
        'Schneider':
            {
                'bibtex':
                    """
                @article{PhysRevB.17.1302,
                  title = {Molecular-dynamics study of a three-dimensional one-component model for distortive phase transitions},
                  author = {Schneider, T. and Stoll, E.},
                  journal = {Phys. Rev. B},
                  volume = {17},
                  issue = {3},
                  pages = {1302--1322},
                  numpages = {0},
                  year = {1978},
                  month = {Feb},
                  publisher = {American Physical Society},
                  doi = {10.1103/PhysRevB.17.1302},
                  url = {https://link.aps.org/doi/10.1103/PhysRevB.17.1302}
            }"""  # noqa: E501
            },
        'Dunweg':
            {
                'bibtex':
                    """
                @article{doi:10.1142/S0129183191001037,
                author = {DÜNWEG, BURKHARD and PAUL, WOLFGANG},
                title = {BROWNIAN DYNAMICS SIMULATIONS WITHOUT GAUSSIAN RANDOM NUMBERS},
                journal = {International Journal of Modern Physics C},
                volume = {02},
                number = {03},
                pages = {817-827},
                year = {1991},
                doi = {10.1142/S0129183191001037},

                URL = {http://www.worldscientific.com/doi/abs/10.1142/S0129183191001037},
                eprint = {http://www.worldscientific.com/doi/pdf/10.1142/S0129183191001037}
            }"""  # noqa: E501
            },
    }

    def __init__(self, flowchart=None, title='NVT dynamics', extension=None):
        """Initialize the node"""

        logger.debug('Creating NVT {}'.format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        logger.debug('NVT after super init, {}'.format(self))

        self.description = 'NVT dynamics step in LAMMPS'

        logger.debug("NVT.init() creating NVT_Parameters object")

        self.parameters = lammps_step.NVT_Parameters()

        logger.debug("NVT.init() completed")

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

        # What will we do?

        if P['T0'] == P['T1']:
            text = "{time} of canonical (NVT) dynamics at {T0} "
        else:
            text = (
                "{time} of canonical (NVT) dynamics starting "
                " at {T0}, going to {T1}, "
            )
        if P['thermostat'] == 'Nose-Hoover':
            text += "using a Nose-Hoover thermostat."
            if P['Tchain'] != '3':
                if P['Tloop'] != '1':
                    text += (
                        " The thermostat will use a chain of {Tchain} "
                        "thermostats with {Tloop} subcycles and a "
                    )
                else:
                    text += (
                        " The thermostat will use a chain of {Tchain} "
                        "thermostats and a "
                    )
            elif P['Tloop'] != '1':
                text += " The thermostat will use {Tloop} subcycles and a "
            else:
                text += " The thermostat will use a "
            text += "drag factor of {drag}."
        elif P['thermostat'] == 'Berendsen':
            text += (
                "using a Berendsen thermostat with a damping time "
                "of {Tdamp}"
            )
        elif 'csvr' in P['thermostat']:
            text += (
                "using a canonical sampling thermostat using velocity "
                "rescaling (CSVR) with a damping time of {Tdamp} and "
                "a {random_seed}."
            )
        elif 'csld' in P['thermostat']:
            text += (
                "using a canonical sampling thermostat using Langevin "
                "dynamics (CSLD) with a damping time of {Tdamp} and "
                "a {random_seed}."
            )
        elif P['thermostat'] == 'velocity rescaling':
            text += (
                "using velocity rescaling every {frequency} with a "
                "temperature window of {window}."
            )
            if P['fraction'] != 1.0:
                text += (
                    " The velocities will only be scaled a fraction "
                    "({fraction}) of the amount needed to fully correct "
                    "the temperature."
                )
        elif P['thermostat'] == 'Langevin':
            text += (
                "using a Langevin thermostat with a damping time "
                "of {Tdamp} and a {random_seed}"
            )
        else:
            text += ("using the thermostat given by {thermostat}")

        return self.header + '\n' + __(text, **P, indent=4 * ' ').__str__()

    def describe(self, indent='', json_dict=None):
        """Write out information about what this node will do
        If json_dict is passed in, add information to that dictionary
        so that it can be written out by the controller as appropriate.
        """

        # Can't call super() because it will print too much
        self.visited = True
        job.job('\n' + self.indent + self.header)
        next_node = self.next()

        # Local copies of variables in a dictionary

        P = self.parameters.values_to_dict()
        text = self.description_text(P)
        job.job(__(text, indent=self.indent + '    ', **P))

        return next_node

    def get_input(self):
        """Get the input for an NVT dynamics run in LAMMPS"""

        self.description = []

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Fix variables with special cases

        # These need to be based on masses...
        if P['timestep'] == 'normal':
            timestep = 1.0
            P['timestep'] = Q_(timestep, ureg.fs)
        elif P['timestep'] == 'accurate but slow':
            timestep = 0.5
            P['timestep'] = Q_(timestep, ureg.fs)
        elif P['timestep'] == 'coarse but fast':
            timestep = 2.0
            P['timestep'] = Q_(timestep, ureg.fs)
        else:
            timestep = P['timestep'].to('fs').magnitude

        if P['seed'] == 'random':
            P['seed'] = int(random.random() * 2**31)

        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = '{:~P}'.format(PP[key])

        self.description.append(
            __(self.description_text(PP), **PP, indent=3 * ' ').__str__()
        )

        time = P['time'].to('fs').magnitude
        nsteps = round(time / timestep)

        T0 = P['T0'].to('K').magnitude
        T1 = P['T1'].to('K').magnitude
        Tdamp = P['Tdamp'].to('fs').magnitude

        thermo_properties = (
            'time temp press etotal ke pe ebond '
            'eangle edihed eimp evdwl etail ecoul elong'
        )
        properties = 'v_time v_temp v_press v_etotal v_ke v_pe v_epair'
        title2 = 'tstep t T P Etot Eke Epe Epair'

        lines = []
        lines.append('')
        lines.append('#     NVT dynamics')
        lines.append('')
        lines.append('reset_timestep      0')
        lines.append('timestep            {}'.format(timestep))
        lines.append('thermo_style        custom {}'.format(thermo_properties))
        lines.append('thermo              {}'.format(int(nsteps / 100)))

        nfixes = 0
        if P['thermostat'] == 'Nose-Hoover':
            Tchain = P['Tchain']
            Tloop = P['Tloop']
            drag = P['drag']
            nfixes += 1
            lines.append(
                'fix                 {} all nvt '.format(nfixes) +
                'temp {} {} {} '.format(T0, T1, Tdamp) +
                'tchain {} '.format(Tchain) + 'tloop {} '.format(Tloop) +
                'drag {}'.format(drag)
            )
        elif P['thermostat'] == 'Berendsen':
            nfixes += 1
            lines.append(
                'fix                 {} '.format(nfixes) +
                'all temp/berendsen ' + ' {} {} {}'.format(T0, T1, Tdamp)
            )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) + 'all nve')
        elif 'csvr' in P['thermostat']:
            seed = P['seed']
            nfixes += 1
            lines.append(
                'fix                 {} '.format(nfixes) + 'all temp/csvr ' +
                ' {} {} {} {}'.format(T0, T1, Tdamp, seed)
            )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) + 'all nve')
        elif 'csld' in P['thermostat']:
            seed = P['seed']
            nfixes += 1
            lines.append(
                'fix                 {} '.format(nfixes) + 'all temp/csld ' +
                ' {} {} {} {}'.format(T0, T1, Tdamp, seed)
            )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) + 'all nve')
        elif P['thermostat'] == 'velocity rescaling':
            frequency = P['frequency'].to('fs').magnitude
            nevery = round(nsteps / (frequency / timestep))
            window = P['window'].to('K').magnitude
            fraction = P['fraction']
            nfixes += 1
            lines.append(
                'fix                 {} '.format(nfixes) +
                'all temp/rescale ' +
                '{} {} {} {} {}'.format(nevery, T0, T1, window, fraction)
            )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) + 'all nve')
        elif P['thermostat'] == 'Langevin':
            seed = P['seed']
            nfixes += 1
            lines.append(
                'fix                 {} '.format(nfixes) + 'all langevin ' +
                '{} {} {} {} '.format(T0, T1, Tdamp, seed)
            )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) + 'all nve')
        else:
            raise RuntimeError(
                "Don't recognize temperature control " +
                "'{}'".format(P['thermostat'])
            )

        # summary output written 10 times during run so we can see progress
        nevery = 10
        nfreq = int(nsteps / 10)
        nrepeat = int(nfreq / nevery)
        nfreq = nevery * nrepeat
        nfixes += 1
        lines.append(
            (
                "fix                 {} all ave/time {} {} {} {} off 2 "
                "title2 '{}' file summary_nvt_{}.txt"
            ).format(
                nfixes, nevery, nrepeat, nfreq, properties, title2,
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
            sampling = P['sampling'].to('fs').magnitude
            nevery = round(sampling / timestep)
            nfreq = int(nsteps / nevery)
            nrepeat = 1
            nfreq = nevery * nrepeat
            nfixes += 1

            if T0 == T1:
                title1 = (
                    '!MolSSI trajectory 1.0 LAMMPS, NVT {} steps of {} fs, '
                    'T={} K'
                ).format(int(nsteps / nevery), timestep * nevery, T0)
            else:
                title1 = (
                    '!MolSSI trajectory 1.0 LAMMPS, NVT {} steps of {} fs, '
                    'T={}-{} K'
                ).format(nsteps, timestep, T0, T1)
            lines.append(
                (
                    "fix                 {} all ave/time {} {} {} {} off 2"
                    " title1 '{}' title2 '{}' file trajectory_nvt_{}.seamm_trj"
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

        lines.append('run                 {}'.format(nsteps))
        lines.append('')
        for fix in range(1, nfixes + 1):
            lines.append('unfix               {}'.format(fix))

        return lines
