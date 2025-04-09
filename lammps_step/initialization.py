# -*- coding: utf-8 -*-

"""A single-point initialization in LAMMPS"""

import logging
import pprint

import lammps_step
import seamm
import seamm_util
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger("lammps")
job = printing.getPrinter()
printer = printing.getPrinter("lammps")

msm_pair_styles = ["born", "buck", "", "lj/charmm", "lj/cut"]

thermo_variables = [
    "step",
    "elapsed",
    "elaplong",
    "dt",
    "time",
    "cpu",
    "tpcpu",
    "spcpu",
    "cpuremain",
    "part",
    "timeremain",
    "atoms",
    "temp",
    "press",
    "pe",
    "ke",
    "etotal",
    "enthalpy",
    "evdwl",
    "ecoul",
    "epair",
    "ebond",
    "eangle",
    "edihed",
    "eimp",
    "emol",
    "elong",
    "etail",
    "vol",
    "density",
    "lx",
    "ly",
    "lz",
    "xlo",
    "xhi",
    "ylo",
    "yhi",
    "zlo",
    "zhi",
    "xy",
    "xz",
    "yz",
    "xlat",
    "ylat",
    "zlat",
    "bonds",
    "angles",
    "dihedrals",
    "impropers",
    "pxx",
    "pyy",
    "pzz",
    "pxy",
    "pxz",
    "pyz",
    "fmax",
    "fnorm",
    "nbuild",
    "ndanger",
    "cella",
    "cellb",
    "cellc",
    "cellalpha",
    "cellbeta",
    "cellgamma",
]


class Initialization(seamm.Node):
    def __init__(self, flowchart=None, title="Initialization", extension=None):
        """Initialize the node"""

        logger.debug("Creating Initialization {}".format(self))

        super().__init__(flowchart=flowchart, title=title, extension=extension)

        self.description = []
        self.parameters = lammps_step.InitializationParameters()

    @property
    def header(self):
        """A printable header for this section of output"""
        return "Step {}: {}".format(".".join(str(e) for e in self._id), self.title)

    @property
    def version(self):
        """The semantic version of this module."""
        return lammps_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
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

        text = "Initialize the calculation with a cutoff of {cutoff}"
        if P["shift_nonbond"]:
            text += ", shifting the nonbond energies to 0 at the cutoff"
        text += ". If the system is periodic"
        if P["kspace_method"][0] == "$":
            text += " use the variable '{method}' to determine whether "
            text += "and how to accelerate the k-space summations."
        elif P["kspace_method"] == "none":
            text += " no k-space acceleration method will be used."
        elif P["kspace_method"] == "automatic":
            text += " the best k-space acceleration method for the "
            text += " molecular system will be chosen."
        else:
            text += " the {method} for k-space acceleration will be used."

        if P["kspace_method"] != "none":
            text += " The accuracy goal is {kspace_accuracy}."

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def get_input(self, extras=None):
        """Get the input for the initialization of LAMMPS"""

        self.description = []
        self.description.append("   " + self.header)

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Fix some things
        P["cutoff"] = P["cutoff"].to("angstrom").magnitude

        # Get the system
        system_db = self.get_variable("_system_db")
        configuration = system_db.system.configuration

        # See what type of forcefield we have and handle it
        ff = self.get_variable("_forcefield")
        if ff == "OpenKIM":
            return self.OpenKIM_input()

        ff_form = ff.ff_form

        # Valence forcefield...
        ffname = ff.current_forcefield
        n_atoms = configuration.n_atoms

        if ff_form != "reaxff":
            # And atom-type if necessary
            key = f"atom_types_{ffname}"
            if key not in configuration.atoms:
                logger.debug("Atom typing")
                ff.assign_forcefield(configuration)
            else:
                if any(typ is None for typ in configuration.atoms[key]):
                    ff.assign_forcefield(configuration)

        # Get the energy expression.
        style = (
            "LAMMPS-class2"
            if ff_form == "class2"
            else (
                "LAMMPS-dreiding"
                if ff_form == "dreiding"
                else "LAMMMPS-reaxff" if ff_form == "reaxff" else "LAMMPS"
            )
        )

        eex = ff.energy_expression(configuration, style=style)
        logger.debug("energy expression:\n" + pprint.pformat(eex))

        # Determine if we have any charges, and if so, if they are sparse
        key = f"charges_{ffname}"
        if key in configuration.atoms:
            charges = [*configuration.atoms[key]]
            n_charged_atoms = 0
            smallq = float(P["kspace_smallq"])
            for charge in charges:
                if abs(charge) > smallq:
                    n_charged_atoms += 1
            fraction_charged_atoms = n_charged_atoms / n_atoms
        else:
            n_charged_atoms = 0

        kspace_style = ""

        lines = []
        lines.append("")
        lines.append(f"# {self.header}")
        lines.append("")
        lines.append("units               real")

        periodicity = configuration.periodicity
        if periodicity == 0:
            lines.append("boundary            s s s")
            string = "Setup for a molecular (non-periodic) system."
        elif periodicity == 3:
            lines.append("boundary            p p p")
            tail_correction = (
                "yes" if P["tail_correction"] and not P["shift_nonbond"] else "no"
            )
            string = "Setup for a periodic (crystalline or fluid) system."
        else:
            raise RuntimeError(
                "The LAMMPS step can only handle 0-"
                " or 3-D periodicity at the moment!"
            )

        if ff_form == "reaxff":
            lines.append("atom_style          full")
        else:
            lines.append("atom_style          full")

        lines.append("newton              on")

        if n_atoms < 20:
            # LAMMPS has problems with bins for small systems
            lines.append("neighbor            2.0 nsq")
        lines.append("")
        if periodicity == 3:
            a, b, c, alpha, beta, gamma = configuration.cell.parameters
            if alpha != 90 or beta != 90 or gamma != 90:
                lines.append("box                 tilt large")
        lines.append(f"#    define the style of forcefield {ffname}")
        lines.append("")

        terms = eex["terms"]

        logging.debug("LAMMPS initialization, terms = \n" + pprint.pformat(terms))

        # Control of nonbonds

        if "nonbond parameters" in eex:
            # If using e.g. hybrid/overlay then must use explicit coeff lines here
            nonbond_forms = set([v[0] for v in eex["nonbond parameters"]])
            use_hybrid = len(nonbond_forms) > 1
            if use_hybrid and "hbond/dreiding/lj" in nonbond_forms:
                hybrid = (
                    "pair_style          hybrid/overlay hbond/dreiding/lj 2 6 6.5 90 "
                )
            else:
                hybrid = "pair_style          "

            nonbond_term = None
            if "pair" in terms:
                if len(terms["pair"]) != 1:
                    raise RuntimeError("Cannot handle multiple nonbond terms yet!")
                nonbond_term = terms["pair"][0]

            if nonbond_term is None:
                pprint.pprint(terms)
                raise RuntimeError("Cannot find nonbond term in forcefield!")

            if nonbond_term == "nonbond(9-6)":
                pair_style_base = "lj/class2"
                mixing = "sixthpower"
            elif nonbond_term == "nonbond(12-6)":
                pair_style_base = "lj/cut"
                # What type of mixing rule?
                modifiers = ff.ff["modifiers"]["nonbond(12-6)"]
                mixing = ""
                for section in modifiers:
                    for item in modifiers[section]:
                        if "combination" in item:
                            if mixing == "":
                                mixing = item.split()[1]
                            if mixing != item.split()[1]:
                                raise RuntimeError(
                                    "Conflicting combination rules in "
                                    "nonbond(12-6) section '" + section + "'"
                                )
                if mixing == "":
                    mixing = "geometric"
            elif nonbond_term == "buckingham":
                pair_style_base = "buck"
                mixing = None
            else:
                raise RuntimeError(
                    "Can't handle nonbond term {} yet!".format(nonbond_term)
                )

            shift = "yes" if P["shift_nonbond"] else "no"
            if P["kspace_method"] == "automatic":
                if periodicity == 3:
                    if n_charged_atoms == 0:
                        pair_style = pair_style_base
                        string += (
                            " The nonbonded interactions will be evaluated using "
                            "a cutoff of {cutoff} Å. Since there are no charges "
                            "on the atoms, no long-range coulomb method will be "
                            "used."
                        )
                    else:
                        string += (
                            " The nonbonded interactions will be evaluated using "
                            "a cutoff of {cutoff} Å, with the long-range terms "
                        )
                        pair_style = pair_style_base + "/coul/long"
                        if n_atoms < P["ewald_atom_cutoff"]:
                            kspace_style = "ewald {}".format(P["kspace_accuracy"])
                            string += (
                                "using the Ewald summation method with "
                                "an accuracy of {kspace_accuracy}."
                            )
                        elif fraction_charged_atoms < P["charged_atom_fraction_cutoff"]:
                            kspace_style = "pppm/cg {} {}".format(
                                P["kspace_accuracy"], P["charged_atom_fraction_cutoff"]
                            )
                            string += (
                                "using the PPPM method optimized for few "
                                "atoms with charges, with "
                                "an accuracy of {kspace_accuracy}."
                            )
                        else:
                            kspace_style = "pppm {}".format(P["kspace_accuracy"])
                            string += (
                                "using the PPPM method with "
                                "an accuracy of {kspace_accuracy}."
                            )
                    lines.append(f"{hybrid}          {pair_style} {P['cutoff']}")
                    if mixing is None:
                        lines.append(
                            f"pair_modify         tail {tail_correction} shift {shift}"
                        )
                    else:
                        lines.append(
                            "pair_modify         mix "
                            + mixing
                            + " tail {} shift {}".format(tail_correction, shift)
                        )
                    if shift:
                        string += (
                            " The van der Waals terms will be shifted "
                            "to zero energy at the cutoff distance."
                        )
                    if tail_correction:
                        string += (
                            " A long-range correction for the "
                            "van der Waals terms will be added."
                        )
                    if kspace_style != "":
                        lines.append("kspace_style        " + kspace_style)
                else:
                    kspace_style = ""
                    if n_charged_atoms == 0:
                        pair_style = pair_style_base
                        string += (
                            " The nonbonded interactions will be evaluated using "
                            "a simple cutoff of {cutoff} Å. Since there are no "
                            "charges on the atoms, no long-range coulomb method "
                            "will be used."
                        )
                    elif (
                        n_atoms > P["msm_atom_cutoff"]
                        and pair_style_base in msm_pair_styles
                    ):
                        pair_style = pair_style_base + "/coul/msm"
                        string += (
                            "The nonbonded interactions will be handled with "
                            " a cutoff of {cutoff} Å."
                        )
                        if fraction_charged_atoms < P["charged_atom_fraction_cutoff"]:
                            kspace_style = "msm/cg {kspace_accuracy} {kspace_smallq}"
                            string += (
                                " The MSM method will be used to handle "
                                "the longer range coulombic interactions, using "
                                "the approach tuned for systems with few charges."
                                "The accuracy goal is {kspace_accuracy}."
                            )
                        else:
                            kspace_style = "msm {kspace_accuracy}"
                            string += (
                                " The MSM method will be used to handle "
                                "the longer range coulombic interactions."
                                "The accuracy goal is {kspace_accuracy}."
                            )
                    else:
                        pair_style = pair_style_base + "/coul/cut"
                        string += (
                            "The nonbonded interactions will be handled with "
                            " a simple cutoff of {cutoff} Å."
                        )
                    if shift:
                        string += (
                            " The van der Waals terms will be shifted "
                            "to zero energy at the cutoff distance."
                        )
                    lines.append(f"{hybrid} {pair_style} {P['cutoff']}")
                    if mixing is None:
                        lines.append(f"pair_modify         shift {shift}")
                    else:
                        lines.append(f"pair_modify         mix {mixing} shift {shift}")
                    if kspace_style != "":
                        lines.append("kspace_style        " + kspace_style.format(**P))
                self.description.append(__(string, indent=7 * " ", **P))
            else:
                if periodicity == 3:
                    kspace_style = ""
                    if n_charged_atoms == 0 or P["kspace_style"] == "none":
                        pair_style = pair_style_base
                    elif fraction_charged_atoms < P["charged_atom_fraction_cutoff"]:
                        pair_style = pair_style_base + "/coul/long"
                        kspace_style = lammps_step.kspace_methods[
                            P["kspace_method"]
                        ].format(**P)
                    lines.append(f"{hybrid} {pair_style} {P['cutoff']}")
                    if mixing is None:
                        lines.append(
                            f"pair_modify         tail {tail_correction} shift {shift}"
                        )
                    else:
                        lines.append(
                            "pair_modify         mix "
                            + mixing
                            + " tail {} shift {}".format(tail_correction, shift)
                        )
                else:
                    if n_charged_atoms == 0:
                        pair_style = pair_style_base
                    else:
                        pair_style = pair_style_base + "/coul/cut"
                    lines.append(f"{hybrid} {pair_style} {P['cutoff']}")
                    if mixing is None:
                        lines.append(f"pair_modify         shift {shift}")
                    else:
                        lines.append(f"pair_modify         mix {mixing} shift {shift}")
                    if "msm" in lammps_step.kspace_methods[P["kspace_method"]]:
                        kspace_style = lammps_step.kspace_methods[
                            P["kspace_method"]
                        ].format(**P)

        if "n_bonds" in eex and eex["n_bonds"] > 0:
            forms = set([v[0] for v in eex["bond parameters"]])
            if len(forms) == 1:
                bond_style = lammps_step.bond_style[[*forms][0]]
                lines.append("bond_style          " + bond_style)
            else:
                line = "bond_style          hybrid"
                for term in forms:
                    line += " " + lammps_step.bond_style[term]
                lines.append(line)
        if "angle" in terms and eex["n_angles"] > 0:
            forms = set([v[0] for v in eex["angle parameters"]])
            if len(forms) == 1:
                angle_style = lammps_step.angle_style[[*forms][0]]
                if angle_style == "table":
                    lines.append(f"angle_style         {angle_style} linear 901")
                else:
                    lines.append("angle_style         " + angle_style)
            else:
                line = "angle_style         hybrid"
                for term in forms:
                    angle_style = lammps_step.angle_style[term]
                    if angle_style == "table":
                        line += f" {angle_style} linear 901"
                    else:
                        line += " " + angle_style
                lines.append(line)
        if "torsion" in terms and eex["n_torsions"] > 0:
            forms = set([v[0] for v in eex["torsion parameters"]])
            if len(forms) == 1:
                dihedral_style = lammps_step.dihedral_style[[*forms][0]]
                lines.append("dihedral_style      " + dihedral_style)
            else:
                line = "dihedral_style      hybrid"
                for term in forms:
                    line += " " + lammps_step.dihedral_style[term]
                lines.append(line)
        if "out-of-plane" in terms and eex["n_oops"] > 0:
            forms = set([v[0] for v in eex["oop parameters"]])
            if len(forms) == 1:
                improper_style = lammps_step.improper_style[[*forms][0]]
                lines.append("improper_style      " + improper_style)
            else:
                line = "improper_style      hybrid"
                for term in forms:
                    line += " " + lammps_step.improper_style[term]
                lines.append(line)

        lines.append("")

        if "nonbond parameters" in eex:
            # Handle 1-4s in nonbonds
            if "opls" in ffname:
                lines.append("special_bonds       lj/coul 0.0 0.0 0.5")
                lines.append("")
            elif "dreiding" in ffname:
                lines.append("special_bonds       dreiding")
                lines.append("")
            elif "cff" in ffname:
                lines.append("special_bonds       lj/coul 0.0 0.0 1.0")
                lines.append("")
            else:
                raise RuntimeError(
                    f"Special bonds not implemented for this force field: {ffname}"
                )

        if extras is not None and "read_data" in extras and extras["read_data"] is True:
            lines.append("read_data           structure.dat")
            lines.append("")

        if "nonbond parameters" in eex:
            if kspace_style != "":
                lines.append("kspace_style        " + kspace_style)

            # For hybrid nonbonds we need the explicit coeffs lines
            if use_hybrid:
                lines.append("")
                lines.append("#    nonbond parameters")
                for parameters in eex["nonbond parameters"]:
                    form, values, types, parameters_type, real_types = parameters
                    i, j = types
                    if form == "nonbond(9-6)":
                        lines.append(
                            f"pair_coeff    {i:6d} {j:6d} {pair_style} {values['eps']} "
                            f"{values['rmin']} # {types[0]} --> {real_types[0]}"
                        )
                    elif form == "nonbond(12-6)":
                        lines.append(
                            f"pair_coeff    {i:6d} {j:6d} {pair_style} {values['eps']} "
                            f"{values['sigma']} # {types[0]} --> {real_types[0]}"
                        )
                    elif form == "buckingham":
                        lines.append(
                            f"pair_coeff    {j:6d} {i:6d} {pair_style} {values['A']} "
                            f"{values['rho']} {values['C']} # {types[1]}-{types[0]} -->"
                            f" {real_types[1]}_{real_types[0]}"
                        )
                    elif form == "hbond/dreiding/lj":
                        lines.append(
                            f"pair_coeff    {j:6d} {i:6d} {form} {values['h_index']} "
                            f"{values['donor flag']} {values['eps']} {values['sigma']} "
                            f"{values['exponent']} "
                        )
                lines.append("")

        if ff_form == "reaxff":
            # Save the sum of the atomic energy terms for later
            self.parent._atomic_energy_sum = eex["Sum of atomic energies"]

            # Setup the Reax FF calculation
            lines.append("# Setting up ReaxFF")

            # How to handle charges...
            charge_method = P["atomic charges"]
            if charge_method not in lammps_step.charge_methods:
                tmp = "\n\t".join(lammps_step.charge_methods.keys())
                raise ValueError(
                    f"Do not recognize charge method '{charge_method}'. Should be "
                    f"one of\n\t{tmp}"
                )
            charge_method = lammps_step.charge_methods[charge_method]
            if charge_method == "default":
                charge_method = ff.charge_method

            if charge_method == "none":
                lines.append("pair_style          reaxff NULL checkqeq no enobonds yes")
            else:
                lines.append("pair_style          reaxff NULL enobonds yes")
            lines.append(
                f"pair_coeff          * * forcefield.dat {' '.join(eex['atom types'])}"
            )
            lines.append("")

            cutoff = P["cutoff"]
            convergence = P["qeq convergence"]
            niter = P["qeq iterations"]
            if charge_method in ("default", "qeq/reaxff"):
                lines.append(
                    f"fix                 charges all qeq/reaxff 1 0.0 {cutoff} "
                    f"{convergence} reaxff maxiter {niter}"
                )
            elif charge_method in ("qeq/point", "qeq/shielded"):
                lines.append(
                    f"fix                 charges all {charge_method} 1 {cutoff} "
                    f"{convergence} {niter} reaxff"
                )
            elif charge_method in ("qtpie/reaxff"):
                lines.append(
                    f"fix                 charges all {charge_method} 1 0.0 {cutoff} "
                    f"{convergence} reaxff Gaussians.dat maxiter {niter}"
                )
            elif charge_method in ("acks2/reaxff"):
                lines.append(
                    f"fix                 charges all {charge_method} 1 0.0 {cutoff} "
                    "{convergence} reaxff maxiter {niter}"
                )

        # Set up standard variables
        for variable in thermo_variables:
            lines.append("variable            {var} equal {var}".format(var=variable))

        # Special Dreiding hydrogen bond information
        self.parent.have_dreiding_hbonds = False
        if ff_form == "dreiding" and "hbond/dreiding/lj" in nonbond_forms:
            self.parent.have_dreiding_hbonds = True
            lines.append("")
            lines.append("compute             hb all pair hbond/dreiding/lj")
            lines.append("variable            N_hbond equal c_hb[1] #number hbonds")
            lines.append("variable            E_hbond equal c_hb[2] #hbond energy")

        return (lines, eex)

    def OpenKIM_input(self):
        """Create the initialization input for a calculation using OpenKIM."""
        # Get the configuration
        system_db = self.get_variable("_system_db")
        configuration = system_db.system.configuration

        # Get the (simple) energy expression for these systems
        eex = self.OpenKIM_energy_expression()

        lines = []
        lines.append("")
        lines.append(f"# {self.header}")
        lines.append("")
        lines.append("newton              on")
        lines.append("")
        potential = self.get_variable("_OpenKIM_Potential")
        if "reax" in potential.lower():
            lammps_step.set_lammps_unit_system("real")
            lines.append("atom_style          charge")
            lines.append(
                f"kim                 init {potential} real unit_conversion_mode"
            )
        else:
            lammps_step.set_lammps_unit_system("metal")
            lines.append(
                f"kim                 init {potential} metal unit_conversion_mode"
            )
        lines.append("")
        periodicity = configuration.periodicity
        if periodicity == 0:
            lines.append("boundary            s s s")
            string = "Setup for a molecular (non-periodic) system."
        elif periodicity == 3:
            lines.append("boundary            p p p")
            string = "Setup for a periodic (crystalline or fluid) system."
        else:
            raise RuntimeError(
                "The LAMMPS step can only handle 0-"
                " or 3-D periodicity at the moment!"
            )
        lines.append("")
        lines.append("read_data           structure.dat")
        lines.append("")
        lines.append(f'kim                 interactions {" ".join(eex["atom types"])}')

        # Set up standard variables
        for variable in thermo_variables:
            lines.append("variable            {var} equal {var}".format(var=variable))

        self.description.append(__(string, indent=self.indent + 4 * " "))

        return (lines, eex)

    def OpenKIM_energy_expression(self):
        """Create the (simple) energy expression for OpenKIM models."""
        eex = {}
        eex["terms"] = {"OpenKIM": []}

        # Get the configuration
        system_db = self.get_variable("_system_db")
        configuration = system_db.system.configuration
        atoms = configuration.atoms

        # The elements (1-based!) Probably not used...
        elements = atoms.symbols
        eex["elements"] = [""]
        eex["elements"].extend(elements)

        # The periodicity & cell parameters
        periodicity = eex["periodicity"] = configuration.periodicity
        if periodicity == 3:
            eex["cell"] = configuration.cell.parameters

        result = eex["atoms"] = []
        atom_types = eex["atom types"] = []
        masses = eex["masses"] = []

        coordinates = atoms.get_coordinates(fractionals=False)
        for element, xyz in zip(elements, coordinates):
            if element in atom_types:
                index = atom_types.index(element) + 1
            else:
                atom_types.append(element)
                index = len(atom_types)
                masses.append(
                    (seamm_util.element_data[element]["atomic weight"], element)
                )
            x, y, z = xyz
            result.append((x, y, z, index))

        eex["n_atoms"] = n_atoms = len(result)
        eex["n_atom_types"] = len(atom_types)

        # For ReaxFF we need charges
        potential = self.get_variable("_OpenKIM_Potential")
        if "reax" in potential.lower():
            if "charge" in atoms:
                eex["charges"] = atoms.get_column_data["charges"]
            else:
                eex["charges"] = [0.0] * n_atoms

        return eex
