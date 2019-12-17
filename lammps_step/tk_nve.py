# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Energy step"""

import lammps_step
import logging
import seamm_widgets as sw
import tkinter as tk
import tkinter.ttk as ttk

logger = logging.getLogger(__name__)


class TkNVE(lammps_step.TkEnergy):

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

    def create_dialog(
        self, title='Edit NVE dynamics parameters', calculation='nve'
    ):
        """Create the dialog!"""

        # Let parent classes do their thing.
        super().create_dialog(title=title, calculation=calculation)

        # Shortcut for parameters
        P = self.node.parameters

        # Frame to isolate widgets
        self['trj_frame'] = ttk.LabelFrame(
            self['frame'],
            borderwidth=4,
            relief='sunken',
            text='Trajectory',
            labelanchor='n',
            padding=10
        )

        self['time'] = P['time'].widget(self['trj_frame'])
        self['timestep'] = P['timestep'].widget(self['trj_frame'])
        self['sampling'] = P['sampling'].widget(self['trj_frame'])

        row = 0

        self['time'].grid(row=row, column=0, sticky=tk.W)
        row += 1
        self['timestep'].grid(row=row, column=0, sticky=tk.W)
        row += 1
        self['sampling'].grid(row=row, column=0, sticky=tk.W)
        row += 1

        sw.align_labels((self['time'], self['timestep'], self['sampling']))

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        frame = self['frame']
        for slave in frame.grid_slaves():
            slave.grid_forget()

        self['trj_frame'].grid(row=0, column=0)
        return 1

    def handle_dialog(self, result):
        if result == 'OK':
            # Shortcut for parameters
            P = self.node.parameters

            value, units = self['time'].get()
            P['time'].value = value
            P['time'].units = units

            tmp = self['timestep'].get()
            if tmp in P['timestep'].enumeration:
                P['timestep'].value = tmp
            else:
                P['timestep'].value = tmp[0]
                P['timestep'].units = tmp[1]

            tmp = self['sampling'].get()
            if tmp in P['sampling'].enumeration:
                P['sampling'].value = tmp
            else:
                P['sampling'].value = tmp[0]
                P['sampling'].units = tmp[1]

        # Let base classes reap their parameters
        super().handle_dialog(result)
