# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Energy step"""

import lammps_step
import logging
import seamm
import tkinter.ttk as ttk

logger = logging.getLogger(__name__)


class TkEnergy(seamm.TkNode):

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
        my_logger=logger
    ):
        '''Initialize a node

        Keyword arguments:
        '''

        self.results_widgets = []

        # Call the constructor for the energy
        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
            my_logger=my_logger
        )

    def right_click(self, event):
        """Probably need to add our dialog...
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(
        self, title='Edit LAMMPS Energy Step', calculation='energy'
    ):
        """Create the dialog!"""

        super().create_dialog(title=title, widget='notebook', results_tab=True)

        self['message'] = ttk.Label(
            self['frame'],
            text='The LAMMPS energy step has no parameters\n'
            'All relevant parameters are set in the initialization step.'
        )
        self['message'].grid()

        self.setup_results(lammps_step.properties, calculation=calculation)
