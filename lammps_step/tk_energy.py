# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Energy step"""

import lammps_step
import seamm_widgets as sw
import seamm
import Pmw
import tkinter as tk
import tkinter.ttk as ttk


class TkEnergy(seamm.TkNode):

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

        self.results_widgets = []

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
        P = self.node.parameters
        """Create the dialog!"""
        self.dialog = Pmw.Dialog(
            self.toplevel,
            buttons=('OK', 'Help', 'Cancel'),
            defaultbutton='OK',
            master=self.toplevel,
            title='Edit Energy parameters',
            command=self.handle_dialog
        )
        self.dialog.withdraw()

        # The tabbed notebook
        notebook = ttk.Notebook(self.dialog.interior())
        notebook.pack(side='top', fill=tk.BOTH, expand=tk.YES)
        self['notebook'] = notebook

        # Main frame holding the widgets
        frame = ttk.Frame(notebook)
        self['frame'] = frame
        notebook.add(frame, text='Parameters', sticky=tk.NW)

        self['message'] = ttk.Label(
            frame,
            text='The LAMMPS energy step has no parameters\n'
            'All relevant parameters are set in the initialization step.'
        )
        self['message'].grid()

        # Second tab for results
        rframe = self['results frame'] = ttk.Frame(notebook)
        notebook.add(rframe, text='Results', sticky=tk.NSEW)

        var = self.tk_var['create tables'] = tk.IntVar()
        if P['create tables'].value == 'yes':
            var.set(1)
        else:
            var.set(0)
        self['create tables'] = ttk.Checkbutton(
            rframe, text='Create tables if needed', variable=var
        )
        self['create tables'].grid(row=0, column=0, sticky=tk.W)

        self['results'] = sw.ScrolledColumns(
            rframe,
            columns=[
                'Result',
                'Save',
                'Variable name',
                'In table',
                'Column name',
            ]
        )
        self['results'].grid(row=1, column=0, sticky=tk.NSEW)
        rframe.columnconfigure(0, weight=1)
        rframe.rowconfigure(1, weight=1)

        self.setup_results('energy')

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS energy
        calculation"""

        if self.dialog is None:
            self.create_dialog()
            self.reset_dialog()

        self.dialog.activate(geometry='centerscreenfirst')

    def reset_dialog(self, widget=None):
        """Layout the dialog according to the current control
        parameters. For the energy, there are no parameters, so
        do nothing."""

        pass

    def setup_results(self, calculation='energy'):
        """Layout the results tab of the dialog"""
        results = self.node.parameters['results'].value

        self.results_widgets = []
        table = self['results']
        frame = table.interior()

        row = 0
        for key, entry in lammps_step.properties.items():
            if 'calculation' not in entry:
                continue
            if calculation not in entry['calculation']:
                continue
            if 'dimensionality' not in entry:
                continue
            if entry['dimensionality'] != 'scalar':
                continue

            widgets = []
            widgets.append(key)

            table.cell(row, 0, entry['description'])

            # variable
            var = self.tk_var[key] = tk.IntVar()
            var.set(0)
            w = ttk.Checkbutton(frame, variable=var)
            table.cell(row, 1, w)
            widgets.append(w)
            e = ttk.Entry(frame, width=15)
            e.insert(0, key.lower())
            table.cell(row, 2, e)
            widgets.append(e)

            if key in results:
                if 'variable' in results[key]:
                    var.set(1)
                    e.delete(0, tk.END)
                    e.insert(0, results[key]['variable'])

            # table
            w = ttk.Combobox(frame, width=10)
            table.cell(row, 3, w)
            widgets.append(w)
            e = ttk.Entry(frame, width=15)
            e.insert(0, key.lower())
            table.cell(row, 4, e)
            widgets.append(e)

            if key in results:
                if 'table' in results[key]:
                    w.set(results[key]['table'])
                    e.delete(0, tk.END)
                    e.insert(0, results[key]['column'])

            self.results_widgets.append(widgets)
            row += 1

        # And make the dialog wide enough
        frame.update_idletasks()
        width = frame.winfo_width() + 70  # extra space for frame, etc.
        height = frame.winfo_height()
        sw = frame.winfo_screenwidth()
        sh = frame.winfo_screenheight()

        mw = 0
        mh = 0
        for tab in self['notebook'].tabs():
            tab = frame.nametowidget(tab)
            w = tab.winfo_reqwidth()
            h = tab.winfo_reqheight()
            if w > mw:
                mw = w
            if h > mh:
                mh = h

        if width < mw:
            width = mw
        if width > sw:
            width = int(0.9 * sw)
        if height < mh:
            height = mh
        if height > sh:
            height = int(0.9 * sh)

        self.dialog.geometry('{}x{}'.format(width, height))

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

        # Shortcut for parameters
        P = self.node.parameters

        # and from the results tab...
        if self.tk_var['create tables'].get():
            P['create tables'].value = 'yes'
        else:
            P['create tables'].value = 'no'

        results = P['results'].value = {}
        for key, w_check, w_variable, w_table, w_column \
                in self.results_widgets:

            if self.tk_var[key].get():
                tmp = results[key] = dict()
                tmp['variable'] = w_variable.get()
            table = w_table.get()
            if table != '':
                if key not in results:
                    tmp = results[key] = dict()
                tmp['table'] = table
                tmp['column'] = w_column.get()
