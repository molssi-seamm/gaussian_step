# -*- coding: utf-8 -*-

"""The graphical part of a Gaussian Optimization node"""

import logging
import pprint
import tkinter as tk
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

        # Top level needs to call reset_dialog
        if self.node.calculation == "optimization":
            self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Layout the widgets, letting our parents go first."""
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        row = 0
        self["calculation"].grid(row=row, column=0)
        row += 1
        self.reset_calculation()
        self["convergence frame"].grid(row=row, column=0)
        row += 1
        self.reset_convergence()
        self["optimization"].grid(row=row, column=0)
        row += 1
        self.reset_optimization()
        self["structure frame"].grid(row=row, column=0)

        return row

    def reset_optimization(self, widget=None):
        frame = self["optimization"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        widgets = []
        # widgets2 = []
        row = 0

        self["geometry convergence"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        widgets.append(self["geometry convergence"])
        row += 1

        self["coordinates"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        widgets.append(self["coordinates"])
        row += 1

        self["max geometry steps"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        widgets.append(self["max geometry steps"])
        row += 1

        self["recalc hessian"].grid(row=row, column=0, columnspan=2, sticky=tk.EW)
        widgets.append(self["recalc hessian"])
        row += 1

        self["ignore unconverged optimization"].grid(
            row=row, column=0, columnspan=2, sticky=tk.EW
        )
        widgets.append(self["ignore unconverged optimization"])
        row += 1

        sw.align_labels(widgets, sticky=tk.E)
