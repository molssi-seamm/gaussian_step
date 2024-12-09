# -*- coding: utf-8 -*-
"""Global control parameters for Gaussian
"""

import logging

import gaussian_step

# import seamm

logger = logging.getLogger("Gaussian")


class OptimizationParameters(gaussian_step.EnergyParameters):
    """The control parameters for the energy."""

    parameters = {
        "target": {
            "default": "minimum",
            "kind": "string",
            "default_units": "",
            "enumeration": ("minimum", "transition state", "saddle point"),
            "format_string": "",
            "description": "Optimization target:",
            "help_text": (
                "The type of structure that is the target for the optimization."
            ),
        },
        "saddle order": {
            "default": 2,
            "kind": "integer",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Order of saddle point:",
            "help_text": "Number of directions down from the saddle point.",
        },
        "max geometry steps": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": ("default", "6*nAtoms", "9*nAtoms"),
            "format_string": "",
            "description": "Maximum steps:",
            "help_text": (
                "The maximum number of steps to take in the optimization. "
                "'default' is based on the system size, giving a reasonable "
                "limit in most cases."
            ),
        },
        "geometry convergence": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": [x for x in gaussian_step.optimization_convergence],
            "format_string": "",
            "description": "Convergence criteria:",
            "help_text": "The criteria to use for convergence.",
        },
        "hessian": {
            "default": "default guess",
            "kind": "string",
            "default_units": "",
            "enumeration": (
                "default guess",
                "read from checkpoint",
                "calculate",
                "from property on configuration",
            ),
            # "from previous step",
            "format_string": "",
            "description": "Hessian:",
            "help_text": "How to get the (initial) Hessian.",
        },
        "hessian step": {
            "default": "",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "Step:",
            "help_text": "Which step to get the Hessian from.",
        },
        "recalc hessian": {
            "default": "at beginning",
            "kind": "integer",
            "default_units": "",
            "enumeration": (
                "every step",
                "at beginning",
                "HF at beginning",
                "5",
                "10",
            ),
            "format_string": "",
            "description": "Recalculate Hessian:",
            "help_text": (
                "How often to recalculate the Hessian (in steps). Smaller "
                "values help convergence but are expensive."
            ),
        },
        "coordinates": {
            "default": "redundant",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "redundant",
                "cartesian",
                "generalized internal (GIC)",
            ),
            "format_string": "s",
            "description": "Type of coordinates:",
            "help_text": "The typ of coordinates to use in the minimization.",
        },
        "ignore curvature error": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Ignore incorrect curvature for TS:",
            "help_text": (
                "Whether to ignore incorrect curvvatures during TS optimizations"
            ),
        },
        "ignore unconverged optimization": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Ignore lack of convergence:",
            "help_text": (
                "Whether to ignore lack of convergence in the optimization. Otherwise, "
                "an error is thrown."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**OptimizationParameters.parameters, **defaults}, data=data
        )

        # Do any local editing of defaults
        tmp = self["configuration name"]
        tmp._data["enumeration"] = ["optimized with {model}", *tmp.enumeration[1:]]
        tmp.default = "keep current name"

        tmp = self["configuration name"]
        tmp._data["enumeration"] = ["optimized with {model}", *tmp.enumeration]
        tmp.default = "optimized with {model}"
