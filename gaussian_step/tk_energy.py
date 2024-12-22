# -*- coding: utf-8 -*-

"""The graphical part of a Gaussian Energy node"""

import logging
import tkinter as tk
import tkinter.ttk as ttk

import gaussian_step
import seamm
import seamm_widgets as sw

logger = logging.getLogger("Gaussian")


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

        P = self.node.parameters

        # The option to just write input
        self["input only"] = P["input only"].widget(frame)
        self["file handling"] = P["file handling"].widget(frame)

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

        # Create the rest of the widgets
        for key in (
            "level",
            "method",
            "basis",
            "advanced_method",
            "functional",
            "advanced_functional",
            "integral grid",
            "dispersion",
            "spin-restricted",
            "freeze-cores",
            "use symmetry",
            "bond orders",
            "apply bond orders",
            "calculate gradient",
            "print basis set",
            "save basis set",
            "basis set file",
        ):
            self[key] = P[key].widget(self["calculation"])

        # bindings...
        for key in (
            "level",
            "method",
            "advanced_method",
            "bond orders",
            "save basis set",
        ):
            self[key].bind("<<ComboboxSelected>>", self.reset_calculation)
            self[key].bind("<Return>", self.reset_calculation)
            self[key].bind("<FocusOut>", self.reset_calculation)

        for key in (
            "maximum iterations",
            "convergence",
        ):
            self[key] = P[key].widget(self["convergence frame"])

        # Create the structure-handling widgets
        sframe = self["structure frame"] = ttk.LabelFrame(
            frame, text="Configuration Handling", labelanchor=tk.N
        )
        row = 0
        widgets = []
        for key in ("structure handling", "system name", "configuration name"):
            self[key] = P[key].widget(sframe)
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1
        sw.align_labels(widgets, sticky=tk.E)

        # A tab for output of orbitals, etc.
        notebook = self["notebook"]
        self["output frame"] = oframe = ttk.Frame(notebook)
        notebook.insert(self["results frame"], oframe, text="Output", sticky="new")

        # Frame to isolate widgets
        p_frame = self["plot frame"] = ttk.LabelFrame(
            self["output frame"],
            borderwidth=4,
            relief="sunken",
            text="Plots",
            labelanchor="n",
            padding=10,
        )

        for key in gaussian_step.EnergyParameters.output_parameters:
            self[key] = P[key].widget(p_frame)

        # Set the callbacks for changes
        for widget in ("orbitals", "region"):
            w = self[widget]
            w.combobox.bind("<<ComboboxSelected>>", self.reset_plotting)
            w.combobox.bind("<Return>", self.reset_plotting)
            w.combobox.bind("<FocusOut>", self.reset_plotting)
        p_frame.grid(row=0, column=0, sticky="new")
        oframe.columnconfigure(0, weight=1)

        # and lay them out
        self.reset_plotting()

        # Top level needs to call reset_dialog
        if self.node.calculation == "energy":
            self.reset_dialog()

        self.logger.debug("Finished creating the dialog")

    def reset_dialog(self, widget=None):
        """Layout the widgets as needed for the current state"""

        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        input_only = self["input only"].get().lower() == "yes"
        row = 0
        # Whether to just write input
        self["input only"].grid(row=row, column=0, sticky=tk.W)
        row += 1
        # And how to handle files
        if not input_only:
            self["file handling"].grid(row=row, column=0, sticky=tk.W)
            row += 1

        self["calculation"].grid(row=row, column=0)
        row += 1
        self.reset_calculation()
        self["convergence frame"].grid(row=row, column=0)
        row += 1
        self.reset_convergence()
        self["structure frame"].grid(
            row=row, column=0, columnspan=2, sticky=tk.N, pady=5
        )
        row += 1
        return row

    def reset_calculation(self, widget=None):
        level = self["level"].get()

        if level == "recommended":
            long_method = self["method"].get()
            functional = self["functional"].get()
            widget = self["method"]
        else:
            long_method = self["advanced_method"].get()
            functional = self["advanced_functional"].get()
            widget = self["advanced_method"]
        bond_orders = self["bond orders"].get()
        save_basis_set = self["save basis set"].get().lower() != "no"
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
                        widget.set(long_method)
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
        if level == "recommended":
            self["method"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self["method"])
            row += 1
            if meta is None or "nobasis" not in meta or not meta["nobasis"]:
                self["basis"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                widgets.append(self["basis"])
                row += 1
            if self.node.method is None or self.node.method == "DFT":
                self["functional"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["functional"])
                row += 1
                self["integral grid"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["integral grid"])
                row += 1
            if meta is None or "freeze core?" in meta and meta["freeze core?"]:
                self["freeze-cores"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["freeze-cores"])
                row += 1
        else:
            self["advanced_method"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self["advanced_method"])
            row += 1
            if meta is None or "nobasis" not in meta or not meta["nobasis"]:
                self["basis"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
                widgets.append(self["basis"])
                row += 1
            if self.node.method is None or self.node.method == "DFT":
                self["advanced_functional"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["advanced_functional"])
                row += 1
                self["integral grid"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["integral grid"])
                row += 1
            if meta is None or "freeze core?" in meta and meta["freeze core?"]:
                self["freeze-cores"].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self["freeze-cores"])
                row += 1
        if self.node.method is None or self.node.method == "DFT":
            if functional in gaussian_step.dft_functionals:
                dispersions = gaussian_step.dft_functionals[functional]["dispersion"]
                if len(dispersions) > 1:
                    w = self["dispersion"]
                    w.config(values=dispersions)
                    if w.get() not in dispersions:
                        w.value(dispersions[1])
                    w.grid(row=row, column=1, sticky=tk.W)
                    widgets2.append(self["dispersion"])
                    row += 1
        for key in ("spin-restricted", "use symmetry", "bond orders"):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        if bond_orders != "none":
            for key in ("apply bond orders",):
                self[key].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self[key])
                row += 1

        self["use symmetry"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        widgets.append(self["use symmetry"])
        row += 1

        if self.__class__.__name__ == "TkEnergy":
            self["calculate gradient"].grid(
                row=row, column=0, columnspan=2, sticky=tk.W
            )
            widgets.append(self["calculate gradient"])
            row += 1

        for key in ("print basis set", "save basis set"):
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        if save_basis_set:
            for key in ("basis set file",):
                self[key].grid(row=row, column=1, sticky=tk.EW)
                widgets2.append(self[key])
                row += 1

        sw.align_labels(widgets, sticky=tk.E)
        sw.align_labels(widgets2, sticky=tk.E)
        frame.columnconfigure(0, minsize=50)

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

    def reset_plotting(self, widget=None):
        frame = self["plot frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        plot_orbitals = self["orbitals"].get() == "yes"
        region = self["region"].get()

        widgets = []

        row = 0
        for key in (
            "total density",
            "total spin density",
            "orbitals",
        ):
            # "difference density",
            self[key].grid(row=row, column=0, columnspan=4, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        if plot_orbitals:
            key = "selected orbitals"
            self[key].grid(row=row, column=1, columnspan=4, sticky=tk.EW)
            row += 1

        key = "region"
        self[key].grid(row=row, column=0, columnspan=4, sticky=tk.EW)
        widgets.append(self[key])
        row += 1

        if region == "explicit":
            key = "nx"
            self[key].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
            widgets.append(self[key])
            self["ny"].grid(row=row, column=2, sticky=tk.EW)
            self["nz"].grid(row=row, column=3, sticky=tk.EW)

        sw.align_labels(widgets, sticky=tk.E)
        frame.columnconfigure(0, minsize=10)
        frame.columnconfigure(4, weight=1)

        return row
