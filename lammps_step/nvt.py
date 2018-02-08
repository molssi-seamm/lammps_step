# -*- coding: utf-8 -*-
"""NVT (canonical) dynamics in LAMMPS"""

from molssi_workflow import units, Q_, units_class, data  # nopep8
import lammps_step
import logging

logger = logging.getLogger(__name__)


class NVT(lammps_step.NVE):
    methods = {
        'Nose-Hoover': {
            'documentation': 'http://lammps.sandia.gov/doc/fix_nvt.html',
            'references': ['Shinoda', 'Tuckerman'],
        },
        'Berendsen': {
            'documentation': 'http://lammps.sandia.gov/doc/fix_temp_berendsen.html',  # nopep8
            'references': ['Berendsen'],
        },
        'canonical sampling, velocity rescaling (csvr)': {
            'documentation': 'http://lammps.sandia.gov/doc/fix_temp_csvr.html',  # nopep8
            'references': ['Bussi1'],
        },
        'canonical sampling, langevin dynamics (csld)': {
            'documentation': 'http://lammps.sandia.gov/doc/fix_temp_csvr.html',  # nopep8
            'references': ['Bussi2'],
        },
        'velocity rescaling': {
            'documentation': 'http://lammps.sandia.gov/doc/fix_temp_rescale.html',  # nopep8
            'references': [],
        },
        'Langevin': {
            'documentation': 'http://lammps.sandia.gov/doc/fix_langevin.html',
            'references': ['Schneider', 'Dunweg']
        },
    }

    references = {
        'Shinoda': {
            'bibtex': """
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
            }"""  # nopep8
        },
        'Tuckerman': {
            'bibtex': """
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
            }"""  # nopep8
        },
        'Berendsen': {
            'bibtex': """
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
            }"""  # nopep8
        },
        'Bussi1': {
            'bibtex': """
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
            }"""  # nopep8
        },
        'Bussi2': {
            'bibtex': """
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
        'Schneider': {
            'bibtex': """
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
            }"""  # nopep8
        },
        'Dunweg': {
            'bibtex': """
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
            }"""  # nopep8
        },
    }

    def __init__(self,
                 workflow=None,
                 gui_object=None,
                 title='NVT dynamics',
                 extension=None):
        """Initialize the node"""

        logger.debug('Creating NVT {}'.format(self))

        super().__init__(
            workflow=workflow,
            title=title,
            gui_object=gui_object,
            extension=extension)

        self.description = 'NVT dynamics step in LAMMPS'

        self.Tcontrol_method = list(lammps_step.NVT.methods)[0]
        self.T0_method = 'is'
        self.T0 = Q_(25, units.degC)
        self.T0_variable = 'T0'
        self.T1_method = 'is'
        self.T1 = self.T0
        self.T1_variable = 'T1'
        self.Tdamp_method = 'is'
        self.Tdamp = Q_(100, 'fs')
        self.Tdamp_variable = 'Tdamp'
        self.Tchain_method = 'is'
        self.Tchain = 3
        self.Tchain_variable = 'Tchain'
        self.Tloop_method = 'is'
        self.Tloop = 1
        self.Tloop_variable = 'Tloop'
        self.drag_method = 'is'
        self.drag = 0.0
        self.drag_variable = 'drag'
        self.seed_method = 'random'
        self.seed = 53
        self.seed_variable = 'seed'
        self.frequency_method = 'is'
        self.frequency = Q_(100, 'fs')
        self.frequency_variable = 'frequency'
        self.window_method = 'is'
        self.window = Q_(20, units.delta_degC)
        self.window_variable = 'window'
        self.fraction_method = 'is'
        self.fraction = 1.0
        self.fraction_variable = 'fraction'

    def get_input(self):
        """Get the input for an NVT dynamics run in LAMMPS"""

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

        T0 = self.get_value('T0')
        T1 = self.get_value('T1')
        Tdamp = self.get_value('Tdamp')

        lines.append('')
        lines.append('#     NVT dynamics')
        lines.append('')
        lines.append('timestep            {}'.format(timestep))
        lines.append('thermo_style        custom {}'.format(thermo_properties))
        lines.append('thermo              {}'.format(int(nsteps/100)))

        nfixes = 0
        if self.Tcontrol_method == 'Nose-Hoover':
            Tchain = self.get_value('Tchain')
            Tloop = self.get_value('Tloop')
            drag = self.get_value('drag')
            nfixes += 1
            lines.append('fix                 {} all nvt '.format(nfixes) +
                         'temp {} {} {} '.format(T0, T1, Tdamp) +
                         'tchain {} '.format(Tchain) +
                         'tloop {} '.format(Tloop) +
                         'drag {}'.format(drag)
                         )
        elif self.Tcontrol_method == 'Berendsen':
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) +
                         'all temp/berendsen ' +
                         ' {} {} {}'.format(T0, T1, Tdamp)
                         )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) +
                         'all nve')
        elif 'csvr' in self.Tcontrol_method:
            seed = self.get_value('seed')
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) +
                         'all temp/csvr ' +
                         ' {} {} {} {}'.format(T0, T1, Tdamp, seed)
                         )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) +
                         'all nve')
        elif 'csld' in self.Tcontrol_method:
            seed = self.get_value('seed')
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) +
                         'all temp/csld ' +
                         ' {} {} {} {}'.format(T0, T1, Tdamp, seed)
                         )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) +
                         'all nve')
        elif self.Tcontrol_method == 'velocity rescaling':
            frequency = self.get_value('frequency')
            nevery = round(nsteps / (frequency / timestep))
            window = self.get_value('window')
            fraction = self.get_value('fraction')
            nfixes += 1
            lines.append(
                'fix                 {} '.format(nfixes) +
                'all temp/rescale ' +
                '{} {} {} {} {}'.format(nevery, T0, T1, window, fraction)
            )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) +
                         'all nve')
        elif self.Tcontrol_method == 'Langevin':
            seed = self.get_value('seed')
            nfixes += 1
            lines.append(
                'fix                 {} '.format(nfixes) +
                'all temp/langevin ' +
                '{} {} {} {} '.format(T0, T1, Tdamp, seed)
            )
            nfixes += 1
            lines.append('fix                 {} '.format(nfixes) +
                         'all nve')
        else:
            raise RuntimeError("Don't recognize temperature control " +
                               "'{}'".format(self.Tcontrol_method))

        # summary output written 100 times during run so we can see progress
        nevery = 10
        nfreq = int(nsteps / 100)
        nrepeat = int(nfreq / nevery)
        nfreq = nevery * nrepeat
        nfixes += 1
        lines.append(
            'fix                 {} '.format(nfixes) +
            'all ave/time ' +
            '{} {} {} {} file nve_summary.txt'.format(
                nevery, nrepeat, nfreq, properties)
        )
        # instantaneous output written every 10 steps for averaging
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
                '{} {} {} {} file nve.txt'.format(
                    nevery, nrepeat, nfreq, properties)
            )

        lines.append('run                 {}'.format(nsteps))
        lines.append('')
        for fix in range(1, nfixes+1):
            lines.append('unfix               {}'.format(fix))

        return lines

    def get_value(self, name):
        methodvar = name + '_method'
        if methodvar in self.__dict__:
            method = self.__dict__[methodvar]
        else:
            method = 'unknown'

        if method == 'is':
            result = self.__dict__[name]
        elif method == 'from variable':
            raise RuntimeError(
                'Variables not implemented yet! {}'.format(name)
            )
        else:
            raise RuntimeError(
                "Can't handle method '{}'!".format(method)
            )

        return self.parent.magnitude_in_lammps_units(result)
