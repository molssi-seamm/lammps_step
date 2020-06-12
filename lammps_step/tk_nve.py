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

        # Metadata for the properties
        self.property_metadata = {}

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

        # Get the property metadata
        for item, data in lammps_step.properties.items():
            if ',' in item:
                continue
            if 'nve' in data["calculation"]:
                self.property_metadata[item] = data

    def create_dialog(
        self, title='Edit NVE dynamics parameters', calculation='nve'
    ):
        """Create the dialog!"""

        # Let parent classes do their thing.
        super().create_dialog(title=title, calculation=calculation)

        # Shortcut for parameters
        P = self.node.parameters

        # Frame to isolate widgets
        c_frame = self['control_frame'] = ttk.LabelFrame(
            self['frame'],
            borderwidth=4,
            relief='sunken',
            text='General Parameters',
            labelanchor='n',
            padding=10
        )

        for key in lammps_step.NVE_Parameters.parameters:
            if key == 'control_properties':
                self[key] = P[key].widget(
                    c_frame, metadata=self.property_metadata
                )
            else:
                self[key] = P[key].widget(c_frame)

        # make the control combobox wide enough
        self['run_control'].combobox.configure(width=40)

        # and binding to change as needed
        self['run_control'].combobox.bind(
            "<<ComboboxSelected>>", self.reset_control_frame
        )

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        frame = self['frame']
        # Clear the dialog
        for slave in frame.grid_slaves():
            slave.grid_forget()

        # Put in our control frame
        row = 0
        self['control_frame'].grid(row=row, column=0)
        row += 1

        # and the widgets in it
        self.reset_control_frame()

        return row

    def reset_control_frame(self, widget=None):
        """Layout the control widgets as needed for the current state"""

        run_control = self['run_control'].get()

        # Clear out the previous widgets
        c_frame = self['control_frame']
        for slave in c_frame.grid_slaves():
            slave.grid_forget()

        # And put them back in depending...
        row = 0
        self['run_control'].grid(row=row, column=0, sticky=tk.W)
        row += 1

        widgets = []
        if 'fixed length' in run_control:
            widgets.append(self['time'])
            self['time'].grid(row=row, column=0, sticky=tk.W)
            row += 1
        else:
            widgets.append(self['maximum_time'])
            self['maximum_time'].grid(row=row, column=0, sticky=tk.W)
            row += 1
        widgets.append(self['timestep'])
        self['timestep'].grid(row=row, column=0, sticky=tk.W)
        row += 1
        if 'fixed length' in run_control:
            widgets.append(self['sampling'])
            self['sampling'].grid(row=row, column=0, sticky=tk.W)
            row += 1
        else:
            self['control_properties'].grid(row=row, column=0, sticky=tk.NSEW)
            row += 1

        sw.align_labels(widgets)

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
