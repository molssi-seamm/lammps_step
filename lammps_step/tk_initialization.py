# -*- coding: utf-8 -*-
"""The graphical part of a LAMMPS Initialization step"""

import molssi_workflow
import lammps_step
import tkinter as tk
import tkinter.ttk as ttk


class TkInitialization(molssi_workflow.TkNode):
    def __init__(self, node=None, canvas=None, x=None, y=None, w=None, h=None):
        '''Initialize a node

        Keyword arguments:
        '''
        super().__init__(node=node, canvas=canvas, x=x, y=y, w=w, h=h)

    def right_click(self, event):
        """Probably need to add our dialog...
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def edit(self):
        """Present a dialog for editing the input for the MOPAC initialization
        calculation"""

        # Create the dialog for editing this node
        self.dialog = tk.Toplevel(master=self.toplevel)
        self._tmp = {'dialog': self.dialog}
        self.dialog.transient(self.toplevel)
        self.dialog.title('LAMMPS Initialization')

        # Main frame holding the widgets
        frame = ttk.Frame(self.dialog)
        frame.pack(side='top', fill=tk.BOTH, expand=1)
        self._tmp['frame'] = frame

        row = 0
        if self.node.structure is None:
            if isinstance(self.node.previous(), molssi_workflow.StartNode):
                self.node.structure = 'initial'
            else:
                self.node.structure = 'current'

        structure_label = ttk.Label(frame, text='Structure:')
        structure = ttk.Combobox(
            frame,
            state='readonly',
            values=list(lammps_step.Initialization.structures))
        structure.set(self.node.structure)
        self._tmp['structure'] = structure

        structure_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
        structure.grid(row=row, column=2, sticky=tk.W)
        row += 1

        # Cutoff
        cutoff_label = ttk.Label(frame, text='Cutoff:')
        cutoff = ttk.Entry(frame, width=15)
        cutoff.insert(0, self.node.cutoff)
        self._tmp['cutoff'] = cutoff

        cutoff_label.grid(row=row, column=0, columnspan=2, sticky=tk.E)
        cutoff.grid(row=row, column=2, sticky=tk.W)
        row += 1

        # Shift nonbond
        self._tmp['shift_nonbond_var'] = tk.IntVar()
        shift_nonbond_label = ttk.Label(frame, text='Shift to 0 at cutoff:')
        shift_nonbond = ttk.Checkbutton(
            frame, variable=self._tmp['shift_nonbond_var'])
        if self.node.shift_nonbond:
            self._tmp['shift_nonbond_var'].set(1)
        else:
            self._tmp['shift_nonbond_var'].set(0)
        self._tmp['shift_nonbond'] = shift_nonbond

        shift_nonbond_label.grid(
            row=row, column=0, columnspan=2, sticky=tk.E)
        shift_nonbond.grid(row=row, column=2, sticky=tk.W)
        row += 1

        # Tail correction
        self._tmp['tail_correction_var'] = tk.IntVar()
        tail_correction_label = ttk.Label(
            frame, text='Correct for long range tail (periodic):')
        tail_correction = ttk.Checkbutton(
            frame, variable=self._tmp['tail_correction_var'])
        if self.node.use_tail_correction:
            self._tmp['tail_correction_var'].set(1)
        else:
            self._tmp['tail_correction_var'].set(0)
        self._tmp['tail_correction'] = tail_correction

        tail_correction_label.grid(
            row=row, column=0, columnspan=2, sticky=tk.E)
        tail_correction.grid(row=row, column=2, sticky=tk.W)
        row += 1

        subframe = ttk.Frame(frame)
        self._tmp['subframe'] = subframe
        subframe.grid(row=row, column=1, columnspan=5)

        frame.grid_columnconfigure(0, minsize=30)

        self.reset_dialog()

        # Button box with the OK, Help and Cancel buttons...
        button_box = ttk.Frame(self.dialog)
        button_box.pack(side='bottom', fill=tk.BOTH)

        ok_button = ttk.Button(button_box, text="OK", command=self.handle_ok)
        ok_button.pack(side='left')
        help_button = ttk.Button(
            button_box, text="Help", command=self.handle_help)
        help_button.pack(side='left')
        cancel_button = ttk.Button(
            button_box, text="Cancel", command=self.handle_cancel)
        cancel_button.pack(side='left')

    def reset_dialog(self, widget=None):
        frame = self._tmp['subframe']
        for slave in frame.grid_slaves():
            slave.destroy()

    def handle_ok(self):
        """Collect the changes from the dialog"""

        self.node.structure = self._tmp['structure'].get()
        self.node.cutoff = self._tmp['cutoff'].get()
        self.node.shift_nonbond = self._tmp['shift_nonbond_var'].get()
        self.node.use_tail_correction = \
            self._tmp['tail_correction_var'].get()

        self.dialog.destroy()
        self.dialog = None
        self._tmp = None

    def handle_help(self):
        print('Help')
        self.dialog.destroy()
        self.dialog = None
        self._tmp = None

    def handle_cancel(self):
        self.dialog.destroy()
        self.dialog = None
        self._tmp = None
