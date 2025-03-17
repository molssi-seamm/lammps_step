# -*- coding: utf-8 -*-

"""The graphical part of a Heat Flux step"""

import tkinter as tk

import lammps_step
import seamm_widgets as sw


class TkHeatFlux(lammps_step.TkNVE):
    """
    The graphical part of a Heat Flux step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in Heat Flux_parameters.py

    See Also
    --------
    HeatFlux, TkHeatFlux,
    HeatFluxParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )

    def create_dialog(self):
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the Heat Flux_parameters
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkHeatFlux.reset_dialog
        """

        frame = super().create_dialog(title="Heat Flux")

        # Shortcut for parameters
        P = self.node.parameters

        # Then create the widgets
        for key in ("time", "timestep", "heat flux", "use centroid stress", "sampling"):
            self[key] = P[key].widget(frame)

        # and lay them out
        self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Layout the widgets in the dialog.

        The widgets are chosen by default from the information in
        Heat Flux_parameter.

        This function simply lays them out row by row with
        aligned labels. You may wish a more complicated layout that
        is controlled by values of some of the control parameters.
        If so, edit or override this method

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None

        See Also
        --------
        TkHeatFlux.create_dialog
        """

        # Remove any widgets previously packed
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        # Shortcut for parameters
        P = self.node.parameters

        # keep track of the row in a variable, so that the layout is flexible
        # if e.g. rows are skipped to control such as "method" here
        row = 0
        widgets = []
        for key in ("time", "timestep", "heat flux", "use centroid stress", "sampling"):
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        # Align the labels
        sw.align_labels(widgets, sticky=tk.E)

        # And how to handle the structure
        if self.node.calculation == "Heat Flux":
            self["structure"].grid(row=row, column=0)
            row += 1

        # Setup the results if there are any
        have_results = (
            "results" in self.node.metadata and len(self.node.metadata["results"]) > 0
        )
        if have_results and "results" in P:
            self.setup_results()

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkHeatFlux.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
