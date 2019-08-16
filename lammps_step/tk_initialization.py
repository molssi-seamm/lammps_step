# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Initialization step"""

import lammps_step
import seamm
import Pmw
import tkinter as tk
import tkinter.ttk as ttk


class TkInitialization(seamm.TkNode):

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50
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
            h=h
        )

    def right_click(self, event):
        """Probably need to add our dialog...
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(self):
        """Create the dialog!"""
        self.dialog = Pmw.Dialog(
            self.toplevel,
            buttons=('OK', 'Help', 'Cancel'),
            defaultbutton='OK',
            master=self.toplevel,
            title='Edit LAMMPS Initialization',
            command=self.handle_dialog
        )
        self.dialog.withdraw()

        self['hull'] = ttk.Frame(self.dialog.interior())
        self['hull'].pack(expand=tk.YES, fill=tk.BOTH)

        # Frame for options for all systems

        frame = self['frame'] = ttk.LabelFrame(
            self['hull'],
            text='For all systems',
            relief=tk.SUNKEN,
            borderwidth=5,
            labelanchor=tk.N
        )

        # Cutoff
        self['cutoff_label'] = ttk.Label(frame, text='Cutoff:')
        self['cutoff'] = ttk.Entry(frame, width=15)
        self['cutoff'].insert(0, self.node.cutoff)

        # Shift nonbond
        self.tk_var['shift_nonbond'] = tk.IntVar()
        self.tk_var['shift_nonbond'].set(self.node.shift_nonbond)
        self['shift_nonbond_label'] = ttk.Label(
            frame, text='Shift to 0 at cutoff:'
        )
        self['shift_nonbond'] = ttk.Checkbutton(
            frame, variable=self.tk_var['shift_nonbond']
        )

        # Frame for the periodic system options, i.e. kspace, etc.

        pframe = self['pframe'] = ttk.LabelFrame(
            self['hull'],
            text='For periodic systems',
            relief=tk.SUNKEN,
            borderwidth=5,
            labelanchor=tk.N
        )

        # Tail correction
        self.tk_var['tail_correction'] = tk.IntVar()
        self.tk_var['tail_correction'].set(self.node.use_tail_correction)
        self['tail_correction_label'] = ttk.Label(
            pframe, text='Correct for long range tail (periodic):'
        )
        self['tail_correction'] = ttk.Checkbutton(
            pframe, variable=self.tk_var['tail_correction']
        )

        # k-space methods
        self['kspace_method_label'] = ttk.Label(
            pframe, text="Long range method:"
        )
        self['kspace_method'] = ttk.Combobox(
            pframe,
            width=30,
            values=list(lammps_step.initialization.kspace_methods)
        )
        self['kspace_method'].bind(
            "<<ComboboxSelected>>", self.kspace_method_cb
        )
        self['kspace_method'].set(self.node.kspace_method)

        # accuracy
        self['accuracy_label'] = ttk.Label(pframe, text='Accuracy:')
        self['accuracy'] = ttk.Entry(pframe, width=15)
        self['accuracy'].insert(0, self.node.kspace_accuracy)

        # small charge cuttoff
        self['smallq_label'] = ttk.Label(pframe, text='Small charge cutoff:')
        self['smallq'] = ttk.Entry(pframe, width=15)
        self['smallq'].insert(0, self.node.kspace_smallq)

        # Grid in the static part of the dialog

        row = 0

        self['cutoff_label'].grid(row=row, column=0, sticky=tk.E)
        self['cutoff'].grid(row=row, column=1, sticky=tk.W)
        row += 1

        self['shift_nonbond_label'].grid(row=row, column=0, sticky=tk.E)
        self['shift_nonbond'].grid(row=row, column=1, sticky=tk.W)
        row += 1

        self['frame'].grid(row=0, column=0, sticky=tk.EW, pady=10)
        self['pframe'].grid(row=1, column=0, sticky=tk.EW, pady=10)

    def handle_dialog(self, result):
        if result is None or result == 'Cancel':
            self.dialog.deactivate(result)
            return

        if result == 'Help':
            # display help!!!
            return

        if result != "OK":
            self.dialog.deactivate(result)
            raise RuntimeError(
                "Don't recognize dialog result '{}'".format(result)
            )

        self.dialog.deactivate(result)

        self.node.cutoff = self['cutoff'].get()
        self.node.shift_nonbond = self.tk_var['shift_nonbond'].get()
        self.node.use_tail_correction = \
            self.tk_var['tail_correction'].get()
        self.node.kspace_method = self['kspace_method'].get()
        self.node.kspace_accuracy = self['accuracy'].get()
        self.node.kspace_smallq = self['smallq'].get()

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS initialization
        """

        if self.dialog is None:
            self.create_dialog()
            self.kspace_method_cb()

        self.dialog.activate(geometry='centerscreenfirst')

    def kspace_method_cb(self, event=None):
        """Grid the widgets into the dialog, depending on the current values
        of key variables. This provides a dyamic presentation to the user.
        """

        # Remove any widgets previously packed
        for slave in self['pframe'].grid_slaves():
            slave.grid_forget()

        row = 0
        self['kspace_method_label'].grid(row=row, column=0, sticky=tk.E)
        self['kspace_method'].grid(row=row, column=1, sticky=tk.W)
        row += 1

        method = self['kspace_method'].get()
        if method != 'none':
            self['accuracy_label'].grid(row=row, column=0, sticky=tk.E)
            self['accuracy'].grid(row=row, column=1, sticky=tk.W)
            row += 1

            if 'few charged' in method or method[0] == '$':
                self['smallq_label'].grid(row=row, column=0, sticky=tk.E)
                self['smallq'].grid(row=row, column=1, sticky=tk.W)
                row += 1

        if method[0] == '$' or 'dispersion' not in method:
            self['tail_correction_label'].grid(row=row, column=0, sticky=tk.E)
            self['tail_correction'].grid(row=row, column=1, sticky=tk.W)
            row += 1
