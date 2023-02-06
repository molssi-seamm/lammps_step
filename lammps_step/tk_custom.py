# -*- coding: utf-8 -*-

"""The graphical part of a LAMMPS Custom step"""

import seamm
import Pmw
import tkinter as tk


class TkCustom(seamm.TkNode):
    def __init__(
        self, tk_flowchart=None, node=None, canvas=None, x=None, y=None, w=200, h=50
    ):
        """Initialize a node

        Keyword arguments:
        """

        super().__init__(
            tk_flowchart=tk_flowchart, node=node, canvas=canvas, x=x, y=y, w=w, h=h
        )

    def right_click(self, event):
        """Probably need to add our dialog..."""

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(self):
        """Create the dialog!"""
        frame = super().create_dialog("Edit LAMMPS Custom Step")

        fixedFont = Pmw.logicalfont("Fixed")
        text = self._widget["text"] = Pmw.ScrolledText(
            frame,
            labelpos="n",
            label_text="Custom script for LAMMPS",
            text_wrap="none",
            text_font=fixedFont,
            text_padx=4,
            text_pady=4,
        )
        text.insert("1.0", self.node.parameters["script"].value)
        text.pack(expand=tk.YES, fill=tk.BOTH)

    def handle_dialog(self, result):
        """Do the right thing when the dialog is closed."""
        if result is None or result == "Cancel":
            self.dialog.deactivate(result)
            self["text"].delete(1.0, "end")
            self["text"].insert(1.0, self.node.parameters["script"].value)
        elif result == "Help":
            self.help()
        elif result == "OK":
            self.dialog.deactivate(result)
            # Capture the parameters from the widgets
            self.node.parameters["script"].value = (
                self["text"].get(1.0, tk.END).rstrip()
            )
        else:
            self.dialog.deactivate(result)
            raise RuntimeError("Don't recognize dialog result '{}'".format(result))
