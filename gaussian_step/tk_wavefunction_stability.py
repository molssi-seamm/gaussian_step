# -*- coding: utf-8 -*-

"""The graphical part of a Gaussian wavefunction stability analysis node"""

import logging

import gaussian_step

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

    def _header_keys(self):
        """The full-width controls at the top of the calculation frame. The rest
        of the layout, and which controls apply, is as for the Energy (see
        WavefunctionStabilityParameters for the rules)."""
        return (
            *super()._header_keys(),
            "stability analysis",
            "test spin multiplicity",
        )
