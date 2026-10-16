# -*- coding: utf-8 -*-

"""What the LAMMPS MD sub-steps will write to md.h5md, known from their parameters.

The NVE, NVT, NPT and Heat Flux sub-steps use these to tell analysis steps what a
run will provide (``md_provides``) before it runs, and to turn on sampling an
analysis needs (``request_sampling``). The names are those of the H5MD file (see
seamm_md_analysis, "The H5MD file").
"""

import logging

import seamm

logger = logging.getLogger(__name__)

#: The observables in the state files, and those of NPT runs
state_observables = (
    "temperature",
    "pressure",
    "pressure_tensor",
    "total_energy",
    "kinetic_energy",
    "potential_energy",
    "pair_energy",
)
npt_observables = ("volume", "density")

#: The sampling parameters for each name: (choice, rate, number of samples)
sampling_parameters = {
    "all/position": (
        "atomic positions",
        "atomic positions rate",
        "atomic positions number of samples",
    ),
    "all/velocity": (
        "atomic velocities",
        "atomic velocities rate",
        "atomic velocities number of samples",
    ),
    "molecules/position": (
        "com positions",
        "com positions rate",
        "com positions number of samples",
    ),
    "molecules/velocity": (
        "com velocities",
        "com velocities rate",
        "com velocities number of samples",
    ),
    "heat_flux": ("heat flux", "heat flux rate", "heat flux number of samples"),
    "pressure_tensor": (
        "pressure tensor",
        "pressure tensor rate",
        "pressure tensor number of samples",
    ),
}


def static_value(node, key):
    """A parameter's value before the run, or None if it cannot be known yet.

    A variable or expression is evaluated if the flowchart's variables already
    define it (at run time); otherwise it is unknown.
    """
    parameters = node.parameters
    if parameters is None or key not in parameters:
        return None
    parameter = parameters[key]
    try:
        if parameter.is_expr:
            return parameter.get(context=seamm.flowchart_variables._data)
        return parameter.get()
    except Exception:
        return None


def _interval(node, choice, rate, number, time):
    from seamm_md_analysis import Interval

    mode = static_value(node, choice)
    if mode is None:
        return Interval()
    if mode == "never":
        return None
    if "interval" in mode:
        return Interval(time=static_value(node, rate))
    n = static_value(node, number)
    if time is not None and n:
        return Interval(time=time / n)
    return Interval()


def md_provides(node):
    """What an NVE, NVT or NPT sub-step will write to md.h5md.

    Parameters
    ----------
    node : lammps_step.NVE
        The sub-step.

    Returns
    -------
    seamm_md_analysis.Provides
    """
    from seamm_md_analysis import Interval, Provides

    time = static_value(node, "time")
    timestep = static_value(node, "timestep")
    if isinstance(timestep, str):
        timestep = None  # "normal" etc. depend on the masses
    data = {}

    sampling = static_value(node, "sampling")
    if sampling is None or sampling != "none":
        interval = Interval(time=None if isinstance(sampling, str) else sampling)
        names = list(state_observables)
        if type(node).__name__ == "NPT":
            names += npt_observables
        for name in names:
            data[name] = interval
        # The user's extra properties are sampled with the state
        if hasattr(node, "extra_definitions"):
            for d in node.extra_definitions():
                data[f"user/{d['name']}"] = interval

    for name, (choice, rate, number) in sampling_parameters.items():
        interval = _interval(node, choice, rate, number, time)
        if interval is not None:
            data[name] = interval

    # The trajectory: in "all" unless atomic positions are sampled too
    trajectory = _interval(
        node, "trajectory", "trajectory rate", "trajectory number of samples", time
    )
    if trajectory is not None:
        group = "trajectory" if "all/position" in data else "all"
        data[f"{group}/position"] = trajectory
        if static_value(node, "trajectory velocities") and "all/velocity" not in data:
            data[f"{group}/velocity"] = trajectory
        if static_value(node, "trajectory forces"):
            data[f"{group}/force"] = trajectory

    parent = getattr(node, "parent", None)
    step = node.title if parent is None else f"{parent.title} / {node.title}"
    return Provides(data=data, timestep=timestep, step=step)


def request_sampling(node, need):
    """Turn on the sampling a need asks for, if possible.

    Parameters
    ----------
    node : lammps_step.NVE
        The sub-step to change.
    need : seamm_md_analysis.Need
        What is needed; its first name that this step can sample is used.

    Returns
    -------
    str or None
        A description of the change for the output, or None if the step could not
        provide it (e.g. the interval is in steps and the timestep is not known).
    """
    timestep = static_value(node, "timestep")
    if isinstance(timestep, str):
        timestep = None
    for name in need.names:
        if name not in sampling_parameters:
            continue
        choice, rate, _ = sampling_parameters[name]
        interval = need.max_time if need.sample_time is None else need.sample_time
        if need.max_steps is not None:
            if timestep is None:
                continue
            by_steps = need.max_steps * timestep
            if interval is None or by_steps < interval:
                interval = by_steps
        if interval is None:
            continue
        interval = interval.to("fs")
        node.parameters[choice].value = "by time interval"
        node.parameters[rate].value = f"{interval.magnitude:g}"
        node.parameters[rate].units = "fs"
        return f"turned on sampling of {name} every {interval:~P} in {node.title}"
    return None
