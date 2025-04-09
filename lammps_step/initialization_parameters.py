# -*- coding: utf-8 -*-
"""Control parameters for the initialization in LAMMPS"""

import logging
import seamm

logger = logging.getLogger(__name__)

kspace_methods = {
    "automatic": "",
    "none": "none",
    "Ewald summation method": "ewald {kspace_accuracy}",
    "PPPM (Particle-particle particle-mesh) method": "pppm {kspace_accuracy}",
    "PPPM method for few charged atoms": "pppm/cg {kspace_accuracy} {smallq}",
    "PPPM method with a staggered mesh": "pppm/stagger {kspace_accuracy}",
    "PPPM method including dispersion terms": "pppm/disp {kspace_accuracy}",
    "MSM (Multilevel summation method)": "msm {kspace_accuracy}",
    "MSM method for few charged atoms": "msm/cg {kspace_accuracy} {smallq}",
}

charge_methods = {
    "default for forcefield": "default",
    "no charges": "none",
    "charge equilibration (QEq), older reaxff implementation": "qeq/reaxff",
    "charge equilibration (QEq) with point charges": "qeq/point",
    "charge equilibration (QEq) with shielded charges": "qeq/shielded",
    "charge equilibration with atom-condensed Kohn–Sham DFT (ACKS2)": "acks2/reaxff",
    "charge transfer with polarization current equilibration (QTPIE)": "qtpie/reaxff",
}


class InitializationParameters(seamm.Parameters):
    """The control parameters for the initialization in LAMMPS"""

    parameters = {
        "cutoff": {
            "default": 10.0,
            "kind": "float",
            "default_units": "Å",
            "enumeration": tuple(),
            "format_string": ".1f",
            "description": "Cutoff:",
            "help_text": (
                "The cutoff used for non-bonded and similar " "interactions."
            ),
        },
        "kspace_method": {
            "default": "automatic",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": tuple(kspace_methods.keys()),
            "format_string": "",
            "description": "K-space method:",
            "help_text": "The method for handling long-range interactions.",
        },
        "kspace_accuracy": {
            "default": 1.0e-05,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": ".2e",
            "description": "K-space accuracy:",
            "help_text": "The target accuracy for the k-space method.",
        },
        "kspace_smallq": {
            "default": 1.0e-05,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": ".1e",
            "description": "K-space negligable charge:",
            "help_text": (
                "The cutoff for the charge on an atom to be considered not zero."
            ),
        },
        "atomic charges": {
            "default": "default for forcefield",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": tuple(charge_methods.keys()),
            "format_string": "",
            "description": "Atomic charge handling:",
            "help_text": "The method for handling atomic charges.",
        },
        "qeq convergence": {
            "default": 1.0e-06,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": ".2e",
            "description": "QEq convergence:",
            "help_text": "The covergence goal for the charge equilibration.",
        },
        "qeq iterations": {
            "default": 100,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "d",
            "description": "QEq iterations:",
            "help_text": (
                "The maximun number of iterations for the charge equilibration."
            ),
        },
        "charged_atom_fraction_cutoff": {
            "default": 0.1,
            "kind": "float",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Cutoff for charged atoms:",
            "help_text": (
                "The upper limit of the fraction of charged atoms "
                "for using methods that exploit sparsity of charges."
            ),
        },
        "ewald_atom_cutoff": {
            "default": 1000,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "d",
            "description": "Atom limit for using Ewald summations:",
            "help_text": (
                "The upper limit of the number of atoms for using "
                "Ewald summations by default."
            ),
        },
        "msm_atom_cutoff": {
            "default": 5000,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "d",
            "description": "Lower atom limit for using the MSM method:",
            "help_text": (
                "The lower limit of the number of atoms for using "
                "the MSM method by default."
            ),
        },
        "tail_correction": {
            "default": "yes",
            "kind": "boolean",
            "format_string": "s",
            "enumeration": ("no", "yes"),
            "description": "Add tail corrections",
            "help_text": (
                "Add tail corrections for neglected long-range " "van der Waals terms."
            ),
        },
        "shift_nonbond": {
            "default": "no",
            "kind": "boolean",
            "format_string": "s",
            "enumeration": ("no", "yes"),
            "description": ("Shift nonbonded interactions to zero at the " "cutoff"),
            "help_text": (
                "Shift the nonbonded interactions to zero at the "
                "cutoff distance to avoid energy jumps."
            ),
        },
        "rigid_waters": {
            "default": "yes",
            "kind": "boolean",
            "format_string": "s",
            "enumeration": ("no", "yes"),
            "description": ("Keep the internal geometry of water molecules " "fixed"),
            "help_text": (
                "Use shake or rattle, as appropriate, to keep "
                "the water molecules rigid."
            ),
        },
        "fix_XH_bond_lengths": {
            "default": "none",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("none", "all", "CH"),
            "format_string": "",
            "description": "Fix X-H bond lengths",
            "help_text": "Fix the bond lengths for H bonded to X.",
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**InitializationParameters.parameters, **defaults}, data=data
        )
