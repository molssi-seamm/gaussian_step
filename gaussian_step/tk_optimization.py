# -*- coding: utf-8 -*-

"""The graphical part of a Gaussian Optimization node"""

import logging
import pprint
import tkinter.messagebox as msg
import tkinter.ttk as ttk

import gaussian_step

# import seamm
import seamm_widgets as sw

logger = logging.getLogger("Gaussian")


class TkOptimization(gaussian_step.TkEnergy):
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
        """Initialize the graphical Tk Gaussian optimization step

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

    def create_dialog(self, title="Edit Gaussian Optimization Step"):
        """Create the edit dialog!

        This is reasonably complicated, so a bit of description
        is in order. The superclass Energy creates the dialog
        along with the calculation parameters in a 'calculation'
        frame..

        This method adds a second frame for controlling the optimizer.

        The layout is handled in part by the Energy superclass, which
        handles the calculation frame. Our part is handled by two
        methods:

        * reset_dialog does the general layout of the main frames.
        * reset_optimization handles the layout of the optimization
          section.
        """

        logger.debug("TkOptimization.create_dialog")

        # Let parent classes do their thing.
        super().create_dialog(title=title)

        # Shortcut for parameters
        P = self.node.parameters

        logger.debug("Parameters:\n{}".format(pprint.pformat(P.to_dict())))

        # Frame to isolate widgets
        opt_frame = self["optimization"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Geometry Optimization",
            labelanchor="n",
            padding=10,
        )

        for key in gaussian_step.OptimizationParameters.parameters:
            self[key] = P[key].widget(opt_frame)

        # bindings...
        for key in ("target", "hessian"):
            self[key].bind("<<ComboboxSelected>>", self.reset_optimization)
            self[key].bind("<Return>", self.reset_optimization)
            self[key].bind("<FocusOut>", self.reset_optimization)

        for key in ("geometry convergence",):
            self[key].bind("<<ComboboxSelected>>", self.check_grid)
            self[key].bind("<Return>", self.check_grid)
            self[key].bind("<FocusOut>", self.check_grid)

        # Top level needs to call reset_dialog
        if self.node.calculation == "optimization":
            self.reset_dialog()

    def check_grid(self, widget=None):
        """Check the grid for the optimization parameters"""
        convergence = self["geometry convergence"].get()
        grid = self["integral grid"].get()

        if "tight" in convergence.lower():
            if grid != "SuperFine":
                if msg.askyesno(
                    "Tighten integral grid?",
                    f"The integral grid is currently '{grid}'.\n"
                    "Change to 'SuperFine'? (recommmended)",
                    parent=self["frame"],
                ):
                    self["integral grid"].set("SuperFine")

    def reset_dialog(self, widget=None, row=0):
        """Layout the widgets, letting our parents go first."""
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        input_only = self["input only"].get().lower() == "yes"

        # Whether to just write input
        self["input only"].grid(row=row, column=0, sticky="w")
        row += 1

        # And how to handle files
        if not input_only:
            self["file handling"].grid(row=row, column=0, columnspan=2, sticky="w")
            row += 1

        self["calculation"].grid(row=row, column=0, columnspan=2)
        row += 1
        self.reset_calculation()

        self["convergence frame"].grid(row=row, column=0, columnspan=2)
        row += 1
        self.reset_convergence()

        self["optimization"].grid(row=row, column=0, sticky="new")
        self.reset_optimization()
        self["structure frame"].grid(row=row, column=1, sticky="new")
        row += 1

        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)

        return row

    def reset_optimization(self, widget=None):
        frame = self["optimization"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        target = self["target"].get()
        hessian = self["hessian"].get()

        widgets = []
        widgets2 = []
        row = 0

        self["target"].grid(row=row, column=0, columnspan=2, sticky="ew")
        widgets.append(self["target"])
        row += 1

        if target not in ("minimum", "transition state"):
            self["saddle order"].grid(row=row, column=1, sticky="ew")
            widgets2.append(self["saddle order"])
            row += 1

        if target not in ("minimum"):
            self["ignore curvature error"].grid(row=row, column=1, sticky="ew")
            widgets2.append(self["ignore curvature error"])
            row += 1

        for key in (
            "geometry convergence",
            "coordinates",
            "max geometry steps",
            "hessian",
        ):
            self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
            widgets.append(self[key])
            row += 1

        if hessian == "calculate":
            for key in ("recalc hessian",):
                self[key].grid(row=row, column=1, sticky="ew")
                widgets2.append(self[key])
                row += 1

        for key in ("ignore unconverged optimization",):
            self[key].grid(row=row, column=0, columnspan=2, sticky="ew")
            widgets.append(self[key])
            row += 1

        w1 = sw.align_labels(widgets, sticky="e")
        w2 = sw.align_labels(widgets2, sticky="e")
        frame.columnconfigure(0, minsize=w1 - w2 + 30)
