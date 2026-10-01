# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialogs: create them and re-lay them out for every choice
that drives the layout, checking that the controls shown are those that the
parameters' rules say apply. Skipped when no display is available."""

import itertools

import pytest

SUBSTEPS = (
    "Initialization",
    "Energy",
    "Minimization",
    "Velocities",
    "NVE",
    "NVT",
    "NPT",
    "HeatFlux",
)

# Controls that a dialog creates but never lays out, so they are not checked for
# being shown when they apply. run_control is fixed at a fixed length of time in the
# dialog, and 'constrain X-H bonds' is not offered.
NEVER_SHOWN = {"run_control", "constrain X-H bonds"}


@pytest.fixture(scope="module")
def root():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw

    Pmw.initialise(root)
    yield root
    root.destroy()


def make(root, substep):
    import seamm

    flowchart = seamm.Flowchart(namespace="org.molssi.seamm.lammps", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.lammps.tk"
    )
    node = flowchart.create_node(substep)
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get(substep)
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    tk_node.create_dialog()
    return tk_node


def shown(widget):
    """Whether a widget is laid out, along with everything containing it in the
    dialog (a control can stay gridded inside a frame that is hidden)."""
    while widget is not None and widget != widget.winfo_toplevel():
        if widget.winfo_manager() == "":
            return False
        widget = widget.master
    return True


def values_of(tk_node):
    from lammps_step.tk_energy import widget_values

    return widget_values(tk_node)


def check(tk_node):
    """The shown controls are exactly those that apply, of the ones the dialog
    lays out."""
    tk_node.reset_dialog()
    P = tk_node.node.parameters
    values = values_of(tk_node)
    n_shown = 0
    for key in P:
        if key == "results" or key not in tk_node or key in NEVER_SHOWN:
            continue
        if shown(tk_node[key]):
            n_shown += 1
            assert P.applies(key, values), f"{key} is shown but does not apply"
        else:
            assert not P.applies(key, values), f"{key} applies but is not shown"
    assert n_shown > 0
    return values


def set_and_check(tk_node, settings):
    for key, value in settings.items():
        tk_node[key].set(value)
    return check(tk_node)


def combinations(P, keys):
    """Every combination of the choices of the given parameters."""
    for choice in itertools.product(*(P[key].enumeration for key in keys)):
        yield dict(zip(keys, choice))


SAMPLED = (
    "atomic positions",
    "com positions",
    "atomic velocities",
    "com velocities",
    "heat flux",
    "shear stress",
)


def check_trajectories(tk_node):
    """The trajectory tab of the dynamics steps."""
    P = tk_node.node.parameters
    for settings in combinations(
        P, ("trajectory", "trajectory save", "trajectory export", "trajectory extxyz")
    ):
        set_and_check(tk_node, settings)
    tk_node["trajectory"].set("never")
    for what in SAMPLED:
        if what not in tk_node or what in tk_node.trajectory_skip:
            continue
        for choice in P[what].enumeration:
            set_and_check(tk_node, {what: choice})
        tk_node[what].set("never")
        check(tk_node)


@pytest.mark.parametrize("substep", SUBSTEPS)
def test_layouts_follow_the_rules(root, substep):
    tk_node = make(root, substep)
    P = tk_node.node.parameters
    check(tk_node)

    if substep == "Initialization":
        for key in ("atomic charges", "kspace_method"):
            for choice in P[key].enumeration:
                set_and_check(tk_node, {key: choice})
        # The PPPM with dispersion handles the long-range dispersion itself
        values = set_and_check(
            tk_node, {"kspace_method": "PPPM method including dispersion terms"}
        )
        assert not shown(tk_node["tail_correction"])
        assert not P.applies("tail_correction", values)
        set_and_check(tk_node, {"kspace_method": "MSM method for few charged atoms"})
        assert shown(tk_node["kspace_smallq"])
    elif substep == "Velocities":
        for choice in P["method"].enumeration:
            set_and_check(tk_node, {"method": choice})
        set_and_check(tk_node, {"method": "scaling current velocities"})
        assert not shown(tk_node["seed"])
    elif substep == "Minimization":
        for settings in combinations(P, ("convergence", "minimizer")):
            set_and_check(tk_node, settings)
        for settings in combinations(
            P, ("optimize cell", "system type", "use_stress", "couple", "allow shear")
        ):
            set_and_check(tk_node, settings)
        # A solid under a general stress, with all components independent
        set_and_check(
            tk_node,
            {
                "optimize cell": "yes",
                "system type": "solid",
                "use_stress": "general stress",
                "couple": "none",
                "allow shear": "yes",
            },
        )
        for key in ("Sxx", "Syy", "Szz", "Syz", "Sxz", "Sxy"):
            assert shown(tk_node[key])
        assert not shown(tk_node["P"])
    elif substep in ("NVE", "NVT", "NPT", "HeatFlux"):
        check_trajectories(tk_node)
        if substep in ("NVT", "NPT"):
            for choice in P["thermostat"].enumeration:
                set_and_check(tk_node, {"thermostat": choice})
            set_and_check(tk_node, {"thermostat": "Nose-Hoover"})
        if substep == "NPT":
            for settings in combinations(
                P,
                (
                    "system type",
                    "barostat",
                    "Panneal",
                    "use_stress",
                    "couple",
                    "allow shear",
                ),
            ):
                set_and_check(tk_node, settings)
            values = set_and_check(
                tk_node,
                {
                    "system type": "solid",
                    "barostat": "Nose-Hoover",
                    "Panneal": "yes",
                    "use_stress": "general stress",
                    "couple": "x and y",
                    "allow shear": "no",
                },
            )
            for key in ("Sxx", "Szz"):
                for suffix in (",initial", ",final", " damp"):
                    assert shown(tk_node[key + suffix])
            assert not shown(tk_node["Syy,initial"])
            assert not shown(tk_node["Pinitial"])
            # A fluid always takes an isotropic pressure
            set_and_check(tk_node, {"system type": "fluid"})
            assert shown(tk_node["Pinitial"])
            assert not shown(tk_node["Sxx,initial"])


def test_variables_show_the_controls_they_might_need(root):
    """A variable counts as meeting a condition, since its value is known only when
    the flowchart runs."""
    tk_node = make(root, "NVT")
    set_and_check(tk_node, {"thermostat": "$thermostat"})
    for key in ("Tdamp", "Tchain", "seed", "frequency"):
        assert shown(tk_node[key])
