# -*- coding: utf-8 -*-

"""The graphical part of a Thermodynamics step"""

import pprint  # noqa: F401
import tkinter as tk
import tkinter.ttk as ttk

import gaussian_step  # noqa: F401, E999
from seamm_util import ureg, Q_, units_class  # noqa: F401, E999
import seamm_widgets as sw


class TkThermodynamics(gaussian_step.TkOptimization):
    """
    The graphical part of a Thermodynamics step in a flowchart.

    Attributes
    ----------
    tk_flowchart : TkFlowchart = None
        The flowchart that we belong to.
    node : Node = None
        The corresponding node of the non-graphical flowchart
    canvas: tkCanvas = None
        The Tk Canvas to draw on
    dialog : Dialog
        The Pmw dialog object
    x : int = None
        The x-coordinate of the center of the picture of the node
    y : int = None
        The y-coordinate of the center of the picture of the node
    w : int = 200
        The width in pixels of the picture of the node
    h : int = 50
        The height in pixels of the picture of the node
    self[widget] : dict
        A dictionary of tk widgets built using the information
        contained in Thermodynamics_parameters.py

    See Also
    --------
    Thermodynamics, TkThermodynamics,
    ThermodynamicsParameters,
    """

    def __init__(
        self,
        tk_flowchart=None,
        node=None,
        canvas=None,
        x=None,
        y=None,
        w=200,
        h=50,
    ):
        """
        Initialize a graphical node.

        Parameters
        ----------
        tk_flowchart: Tk_Flowchart
            The graphical flowchart that we are in.
        node: Node
            The non-graphical node for this step.
        namespace: str
            The stevedore namespace for finding sub-nodes.
        canvas: Canvas
           The Tk canvas to draw on.
        x: float
            The x position of the nodes center on the canvas.
        y: float
            The y position of the nodes cetner on the canvas.
        w: float
            The nodes graphical width, in pixels.
        h: float
            The nodes graphical height, in pixels.

        Returns
        -------
        None
        """
        self.dialog = None

        super().__init__(
            tk_flowchart=tk_flowchart,
            node=node,
            canvas=canvas,
            x=x,
            y=y,
            w=w,
            h=h,
        )

    def create_dialog(self, title="Edit Gaussian Thermodynamics Step"):
        """
        Create the dialog. A set of widgets will be chosen by default
        based on what is specified in the Thermodynamics_parameters
        module.

        Parameters
        ----------
        None

        Returns
        -------
        None

        See Also
        --------
        TkThermodynamics.reset_dialog
        """

        # Let parent classes do their thing.
        super().create_dialog(title=title)

        # Shortcut for parameters
        P = self.node.parameters

        # Frame to isolate widgets
        thermo_frame = self["thermodynamics"] = ttk.LabelFrame(
            self["frame"],
            borderwidth=4,
            relief="sunken",
            text="Thermodynamics",
            labelanchor="n",
            padding=10,
        )

        for key in gaussian_step.ThermodynamicsParameters.parameters:
            self[key] = P[key].widget(thermo_frame)

        # bindings...
        for key in ("optimize first",):
            self[key].bind("<<ComboboxSelected>>", self.reset_dialog)
            self[key].bind("<Return>", self.reset_dialog)
            self[key].bind("<FocusOut>", self.reset_dialog)

        # Top level needs to call reset_dialog
        if self.node.calculation == "thermodynamics":
            self.reset_dialog()

    def reset_dialog(self, widget=None):
        """Layout the widgets, letting our parents go first."""
        frame = self["frame"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        optimize_first = self["optimize first"].get() != "no"

        row = 0
        # Whether to just write input
        self["input only"].grid(row=row, column=0, sticky=tk.W)
        row += 1

        self["thermodynamics"].grid(row=row, column=0)
        self.reset_thermodynamics()

        self["calculation"].grid(row=row, column=1)
        self.reset_calculation()
        row += 1

        if optimize_first:
            self["optimization"].grid(row=row, column=0)
            self.reset_optimization()

        self["convergence frame"].grid(row=row, column=1)
        self.reset_convergence()
        row += 1

        self["structure frame"].grid(row=row, column=0, columnspan=2)

        return row

    def reset_thermodynamics(self, widget=None):
        """Layout the widgets in the dialog.

        The widgets are chosen by default from the information in
        Thermodynamics_parameter.

        This function simply lays them out row by row with
        aligned labels. You may wish a more complicated layout that
        is controlled by values of some of the control parameters.
        If so, edit or override this method

        Parameters
        ----------
        widget : Tk Widget = None

        Returns
        -------
        None

        See Also
        --------
        TkThermodynamics.create_dialog
        """

        # Remove any widgets previously packed
        frame = self["thermodynamics"]
        for slave in frame.grid_slaves():
            slave.grid_forget()

        # keep track of the row in a variable, so that the layout is flexible
        # if e.g. rows are skipped to control such as "method" here
        row = 0
        widgets = []
        for key in ("optimize first",):
            self[key].grid(row=row, column=0, sticky=tk.EW)
            widgets.append(self[key])
            row += 1

        # Align the labels
        sw.align_labels(widgets, sticky=tk.E)

    def right_click(self, event):
        """
        Handles the right click event on the node.

        Parameters
        ----------
        event : Tk Event

        Returns
        -------
        None

        See Also
        --------
        TkThermodynamics.edit
        """

        super().right_click(event)
        self.popup_menu.add_command(label="Edit..", command=self.edit)

        self.popup_menu.tk_popup(event.x_root, event.y_root, 0)
