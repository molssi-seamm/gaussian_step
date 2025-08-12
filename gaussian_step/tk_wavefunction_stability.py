# -*- coding: utf-8 -*-

"""The graphical part of a Gaussian wavefunction stability analysis node"""

import logging
import tkinter as tk

import gaussian_step
import seamm_widgets as sw

logger = logging.getLogger("Gaussian")


class TkWavefunctionStability(gaussian_step.TkEnergy):
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
        """Initialize the graphical Tk Gaussian wavefucntion stability step

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

    def create_dialog(self, title="Edit Gaussian Wavefunction Stability Step"):
        """Create the dialog!"""
        self.logger.debug("Creating the dialog")

        # Let parent classes do their thing.
        super().create_dialog(title=title)

        P = self.node.parameters

        # Create the rest of the widgets
        for key in ("stability analysis", "test spin multiplicity"):
            self[key] = P[key].widget(self["calculation"])

        # Top level needs to call reset_dialog
        if self.node.calculation == "wavefunction stability":
            self.reset_dialog()

        self.logger.debug("Finished creating the dialog")

    def reset_calculation(self, widget=None):
        level = self["level"].get()

        long_method = self["method"].get()
        self["advanced_method"].set(long_method)
        if self.is_expr(long_method):
            self.node.method = None
            meta = None
        else:
            if long_method in gaussian_step.methods:
                self.node.method = gaussian_step.methods[long_method]["method"]
                meta = gaussian_step.methods[long_method]
            else:
                # See if it matches the keyword part
                for key, mdata in gaussian_step.methods.items():
                    if long_method == mdata["method"]:
                        long_method = key
                        self["method"].set(long_method)
                        meta = mdata
                        self.node.method = meta["method"]
                        break
                else:
                    self.node.method = long_method
                    meta = None

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

        for key in (
            "initial checkpoint",
            "checkpoint",
            "geometry",
            "initial wavefunction",
            "stability analysis",
            "test spin multiplicity",
        ):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        self["method"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        widgets.append(self["method"])
        row += 1

        if level == "recommended":
            if self.node.method is None or self.node.method == "DFT":
                self["functional"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["functional"])
                row += 1
                self["integral grid"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["integral grid"])
                row += 1
        else:
            if self.node.method is None or self.node.method == "DFT":
                self["advanced_functional"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["advanced_functional"])
                row += 1
                self["integral grid"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["integral grid"])
                row += 1

        if meta is None or "nobasis" not in meta or not meta["nobasis"]:
            self["basis"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self["basis"])
            row += 1

        for key in (
            "spin-restricted",
            "use symmetry",
        ):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        width0 = sw.align_labels(widgets, sticky=tk.E)
        width1 = sw.align_labels(widgets2, sticky=tk.E)
        frame.columnconfigure(0, minsize=width0 - width1 + 50)

        return row
