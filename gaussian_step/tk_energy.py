# -*- coding: utf-8 -*-

"""The graphical part of a Gaussian Energy node"""

import logging
import tkinter.ttk as ttk

import gaussian_step
from gaussian_step.energy_parameters import find_method
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
            "geometry",
            "initial checkpoint",
            "checkpoint",
            "initial wavefunction",
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
            "functional",
            "advanced_functional",
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
            frame, text="Configuration Handling", labelanchor="n"
        )
        row = 0
        widgets = []
        for key in (
            "save standard orientation",
            "structure handling",
            "system name",
            "configuration name",
        ):
            self[key] = P[key].widget(sframe)
            self[key].grid(row=row, column=0, sticky="ew")
            widgets.append(self[key])
            row += 1
        sw.align_labels(widgets, sticky="e")

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

        # Whether to write input only decides whether file handling applies
        for event in ("<<ComboboxSelected>>", "<Return>", "<FocusOut>"):
            self["input only"].bind(event, self.reset_dialog)

        # Top level needs to call reset_dialog
        if self.node.calculation == "energy":
            self.reset_dialog()

        self.logger.debug("Finished creating the dialog")

    def reset_dialog(self, widget=None, row=0):
        """Layout the widgets as needed for the current state"""

        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        # Whether to just write input
        self["input only"].grid(row=row, column=0, columnspan=2, sticky="w")
        row += 1
        # And how to handle files
        if self.node.parameters.applies("file handling", self._widget_values()):
            self["file handling"].grid(row=row, column=0, sticky="w")
            row += 1

        self["calculation"].grid(row=row, column=0, columnspan=2)
        row += 1
        self.reset_calculation()

        self["convergence frame"].grid(row=row, column=0, sticky="new")
        self.reset_convergence()
        self["structure frame"].grid(
            row=row, column=1, columnspan=1, sticky="new", pady=5
        )
        row += 1

        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)

        return row

    def reset_calculation(self, widget=None):
        """Lay out the calculation frame for the current choices.

        Which controls are shown, and what they offer, come from the parameters'
        rules (gaussian_step.EnergyParameters and its subclasses), which the
        flowchart builder uses too. Indentation and order are kept here.
        """
        P = self.node.parameters
        values = self._widget_values()

        # The method, written out in full if given by its Gaussian keyword. The
        # node's method sets up the results, which depend on it.
        key = P.method_parameter(values) or "advanced_method"
        method_string = self[key].get()
        if self.is_expr(method_string):
            self.node.method = None
        else:
            method, meta, name = find_method(method_string)
            if meta is not None and name != method_string:
                self[key].set(name)
            self.node.method = method

        # Values that the others imply
        values = self._widget_values()
        for key, value in P.implied(values).items():
            if key in self and values.get(key) != value:
                self[key].set(value)
        values = self._widget_values()

        def applies(key):
            return P.applies(key, values)

        # Set up the results table because it depends on the method
        self.results_widgets = []
        self.setup_results()

        frame = self["calculation"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        widgets2 = []
        row = 0

        def add_full(key):
            nonlocal row
            if applies(key):
                self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
                widgets.append(self[key])
                row += 1

        def add_indented(key, sticky="ew"):
            nonlocal row
            if applies(key):
                self[key].grid(row=row, column=1, sticky=sticky)
                widgets2.append(self[key])
                row += 1

        self["level"].grid(row=row, column=0, columnspan=2, sticky="ew")
        row += 1

        for key in self._header_keys():
            add_full(key)

        # The method (one or the other, depending on the level), then indented
        # the functional and grid for DFT, the dispersion correction if the
        # functional offers a choice, and freezing the core for correlated methods.
        add_full("method")
        add_full("advanced_method")
        add_indented("functional")
        add_indented("advanced_functional")
        add_indented("integral grid")
        if applies("dispersion"):
            self._filter_dispersions(values)
            add_indented("dispersion", sticky="w")
        add_indented("freeze-cores")

        # The basis, unless the method does not use one
        add_full("basis")

        for key in ("spin-restricted", "use symmetry", "bond orders"):
            add_full(key)
        add_indented("apply bond orders")

        for key in ("calculate gradient", "print basis set", "save basis set"):
            add_full(key)
        add_indented("basis set file")

        width0 = sw.align_labels(widgets, sticky="e")
        width1 = sw.align_labels(widgets2, sticky="e")
        frame.columnconfigure(0, minsize=width0 - width1 + 50)

        return row

    def _header_keys(self):
        """The full-width controls at the top of the calculation frame."""
        return ("initial checkpoint", "checkpoint", "geometry", "initial wavefunction")

    def _filter_dispersions(self, values=None):
        """Offer only the dispersion corrections the functional has, keeping the
        selection one of them (preferring a correction to none)."""
        P = self.node.parameters
        values = self._widget_values() if values is None else values
        allowed = P.choices("dispersion", values)
        if allowed is None:
            allowed = P["dispersion"].enumeration
        allowed = list(allowed)
        w = self["dispersion"]
        w.config(values=allowed)
        current = w.get()
        if current not in allowed and not self.is_expr(current):
            w.set(allowed[1] if len(allowed) > 1 else allowed[0])

    def _widget_values(self):
        """The dialog's current values, {name: value}, for the parameters' rules."""
        values = {}
        for key in self.node.parameters:
            if key == "results" or key not in self:
                continue
            try:
                value = self[key].get()
            except Exception:
                continue
            values[key] = value[0] if isinstance(value, tuple) else value
        return values

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
            self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
            widgets.append(self[key])
            row += 1

        frame.columnconfigure(0, minsize=150)
        sw.align_labels(widgets, sticky="e")

    def reset_plotting(self, widget=None):
        frame = self["plot frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        P = self.node.parameters
        values = self._widget_values()

        widgets = []

        row = 0
        for key in (
            "total density",
            "total spin density",
            "orbitals",
            "save wfx",
        ):
            # "difference density",
            self[key].grid(row=row, column=0, columnspan=4, sticky="ew")
            widgets.append(self[key])
            row += 1

        if P.applies("selected orbitals", values):
            key = "selected orbitals"
            self[key].grid(row=row, column=1, columnspan=4, sticky="ew")
            row += 1

        key = "region"
        self[key].grid(row=row, column=0, columnspan=4, sticky="ew")
        widgets.append(self[key])
        row += 1

        if P.applies("nx", values):
            key = "nx"
            self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
            widgets.append(self[key])
        if P.applies("ny", values):
            self["ny"].grid(row=row, column=2, sticky="ew")
        if P.applies("nz", values):
            self["nz"].grid(row=row, column=3, sticky="ew")

        sw.align_labels(widgets, sticky="e")
        frame.columnconfigure(0, minsize=10)
        frame.columnconfigure(4, weight=1)

        return row
