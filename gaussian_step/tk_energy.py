# -*- coding: utf-8 -*-

"""The graphical part of a Gaussian Energy node"""

import logging
import tkinter as tk
import tkinter.ttk as ttk

import gaussian_step
import seamm
import seamm_widgets as sw

logger = logging.getLogger(__name__)


class TkEnergy(seamm.TkNode):
    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=120,
        y=20,
        w=200,
        h=50,
        my_logger=logger,
    ):
        """Initialize the graphical Tk Gaussian energy step

        Keyword arguments:
        """
        self.results_widgets = []

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
            my_logger=my_logger,
        )

    def right_click(self, event):
        """Probably need to add our dialog..."""

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)

    def create_dialog(self, title="Edit Gaussian Energy Step"):
        """Create the dialog!"""
        self.logger.debug("Creating the dialog")
        frame = super().create_dialog(title=title, widget="notebook", results_tab=True)

        # Create a frame for the calculation control
        self["calculation"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="Calculation",
            labelanchor="n",
            padding=10,
        )
        # Create a frame for the convergence control
        self["convergence frame"] = ttk.LabelFrame(
            frame,
            borderwidth=4,
            relief="sunken",
            text="SCF Convergence Control",
            labelanchor="n",
            padding=10,
        )

        # Create all the widgets
        P = self.node.parameters
        for key in (
            "level",
            "method",
            "basis",
            "advanced_method",
            "functional",
            "advanced_functional",
            "dispersion",
            "spin-restricted",
            "freeze-cores",
        ):
            self[key] = P[key].widget(self["calculation"])
        for key in (
            "maximum iterations",
            "convergence",
        ):
            self[key] = P[key].widget(self["convergence frame"])

        # bindings...
        for key in ("level", "method", "advanced_method"):
            self[key].bind("<<ComboboxSelected>>", self.reset_convergence)
            self[key].bind("<Return>", self.reset_convergence)
            self[key].bind("<FocusOut>", self.reset_convergence)

        # Top level needs to call reset_dialog
        if self.node.calculation == "energy":
            self.reset_dialog()

        self.logger.debug("Finished creating the dialog")

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        self["calculation"].grid(row=0, column=0)
        self.reset_calculation()
        self["convergence frame"].grid(row=1, column=0)
        self.reset_convergence()
        return 2

    def reset_calculation(self, widget=None):
        level = self["level"].get()

        if level == "recommended":
            long_method = self["method"].get()
            if self.is_expr(long_method):
                self.node.method = None
                meta = None
            else:
                self.node.method = gaussian_step.methods[long_method]["method"]
                meta = gaussian_step.methods[long_method]
            functional = self["functional"].get()
        else:
            long_method = self["advanced_method"].get()
            if self.is_expr(long_method):
                self.node.method = None
                meta = None
            else:
                self.node.method = gaussian_step.methods[long_method]["method"]
                meta = gaussian_step.methods[long_method]
            functional = self["advanced_functional"].get()

        # Set up the results table because it depends on the method
        self.results_widgets = []
        self.setup_results()

        frame = self["calculation"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        widgets2 = []
        row = 0
        self["level"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        row += 1
        if level == "recommended":
            self["method"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self["method"])
            row += 1
            self["basis"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self["basis"])
            row += 1
            if self.node.method is None or self.node.method == "DFT":
                self["functional"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["functional"])
                row += 1
            if meta is None or "freeze core?" in meta and meta["freeze core?"]:
                self["freeze-cores"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["freeze-cores"])
                row += 1
        else:
            self["advanced_method"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self["advanced_method"])
            row += 1
            self["basis"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self["basis"])
            row += 1
            if self.node.method is None or self.node.method == "DFT":
                self["advanced_functional"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["advanced_functional"])
                row += 1
            if meta is None or "freeze core?" in meta and meta["freeze core?"]:
                self["freeze-cores"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["freeze-cores"])
                row += 1
        if self.node.method is None or self.node.method == "DFT":
            dispersions = gaussian_step.dft_functionals[functional]["dispersion"]
            if len(dispersions) > 1:
                w = self["dispersion"]
                w.config(values=dispersions)
                if w.get() not in dispersions:
                    w.value(dispersions[1])
                w.grid(row=row, column=1, sticky=tk.W)
                widgets2.append(self["dispersion"])
                row += 1
                sw.align_labels(widgets2, sticky=tk.E)
            frame.columnconfigure(0, minsize=30)
        self["spin-restricted"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        widgets.append(self["spin-restricted"])
        row += 1
        sw.align_labels(widgets, sticky=tk.E)

        return row

    def reset_convergence(self, widget=None):
        """Layout the convergence widgets as needed for the current state"""

        frame = self["convergence frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        row = 0

        for key in (
            "maximum iterations",
            "convergence",
        ):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        frame.columnconfigure(0, minsize=150)
        sw.align_labels(widgets, sticky=tk.E)
