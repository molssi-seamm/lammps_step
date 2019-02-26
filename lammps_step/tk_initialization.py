# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS Initialization step"""

import molssi_workflow
import Pmw
import tkinter as tk
import tkinter.ttk as ttk


class TkInitialization(molssi_workflow.TkNode):
    def __init__(self, tk_workflow=None, node=None, canvas=None,
                 x=None, y=None, w=200, h=50):
        '''Initialize a node

        Keyword arguments:
        '''
        super().__init__(tk_workflow=tk_workflow, node=node,
                         canvas=canvas, x=x, y=y, w=w, h=h)

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
            command=self.handle_dialog)
        self.dialog.withdraw()

        frame = ttk.Frame(self.dialog.interior())
        frame.pack(expand=tk.YES, fill=tk.BOTH)
        self._widget['frame'] = frame

        row = 0
        if self.node.structure is None:
            if isinstance(self.node.previous(), molssi_workflow.StartNode):
                self.node.structure = 'initial'
            else:
                self.node.structure = 'current'

        # structure_label = ttk.Label(frame, text='Structure:')
        # structure = ttk.Combobox(
        #     frame,
        #     state='readonly',
        #     values=list(lammps_step.Initialization.structures))
        # structure.set(self.node.structure)
        # self._widget['structure'] = structure

        # structure_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
        # structure.grid(row=row, column=2, sticky=tk.W)
        # row += 1

        # Cutoff
        cutoff_label = ttk.Label(frame, text='Cutoff:')
        cutoff = ttk.Entry(frame, width=15)
        cutoff.insert(0, self.node.cutoff)
        self._widget['cutoff'] = cutoff

        cutoff_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
        cutoff.grid(row=row, column=2, sticky=tk.W)
        row += 1

        # Shift nonbond
        self._widget['shift_nonbond_var'] = tk.IntVar()
        shift_nonbond_label = ttk.Label(frame, text='Shift to 0 at cutoff:')
        shift_nonbond = ttk.Checkbutton(
            frame, variable=self._widget['shift_nonbond_var'])
        if self.node.shift_nonbond:
            self._widget['shift_nonbond_var'].set(1)
        else:
            self._widget['shift_nonbond_var'].set(0)
        self._widget['shift_nonbond'] = shift_nonbond

        shift_nonbond_label.grid(
            row=row, column=0, columnspan=2, sticky=tk.E)
        shift_nonbond.grid(row=row, column=2, sticky=tk.W)
        row += 1

        # Tail correction
        self._widget['tail_correction_var'] = tk.IntVar()
        tail_correction_label = ttk.Label(
            frame, text='Correct for long range tail (periodic):')
        tail_correction = ttk.Checkbutton(
            frame, variable=self._widget['tail_correction_var'])
        if self.node.use_tail_correction:
            self._widget['tail_correction_var'].set(1)
        else:
            self._widget['tail_correction_var'].set(0)
        self._widget['tail_correction'] = tail_correction

        tail_correction_label.grid(
            row=row, column=0, columnspan=2, sticky=tk.E)
        tail_correction.grid(row=row, column=2, sticky=tk.W)
        row += 1

        subframe = ttk.Frame(frame)
        self._widget['subframe'] = subframe
        subframe.grid(row=row, column=1, columnspan=5)

        frame.grid_columnconfigure(0, minsize=30)

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
                "Don't recognize dialog result '{}'".format(result))

        self.dialog.deactivate(result)

        w = self._widget
        self.node.structure = w['structure'].get()
        self.node.cutoff = w['cutoff'].get()
        self.node.shift_nonbond = w['shift_nonbond_var'].get()
        self.node.use_tail_correction = \
            w['tail_correction_var'].get()

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS energy
        calculation"""

        if self.dialog is None:
            self.create_dialog()
            self.reset_dialog()

        self.dialog.activate(geometry='centerscreenfirst')

    def reset_dialog(self, widget=None):
        pass
