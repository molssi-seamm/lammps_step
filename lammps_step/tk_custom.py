# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Custom step"""

import seamm
import Pmw
import tkinter as tk
import tkinter.ttk as ttk


class TkCustom(seamm.TkNode):

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
            master=self.toplevel,
            title='Edit Custom step',
            command=self.handle_dialog
        )
        self.dialog.withdraw()

        frame = ttk.Frame(self.dialog.interior())
        frame.pack(expand=tk.YES, fill=tk.BOTH)
        self._widget['frame'] = frame

        fixedFont = Pmw.logicalfont('Fixed')
        text = self._widget['text'] = Pmw.ScrolledText(
            frame,
            labelpos='n',
            label_text='Custom script for LAMMPS',
            text_wrap='none',
            text_font=fixedFont,
            text_padx=4,
            text_pady=4,
        )
        text.insert('1.0', self.node.text)
        text.pack(expand=tk.YES, fill=tk.BOTH)

    def edit(self):
        """Present a dialog for editing the input for the LAMMPS energy
        calculation"""

        if self.dialog is None:
            self.create_dialog()
            self.reset_dialog()

        self.dialog.activate(geometry='centerscreenfirst')

    def reset_dialog(self, widget=None):
        pass

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

        self.node.text = self._widget['text'].get('1.0', 'end')
