# -*- coding: utf-8 -*-
"""
Control parameters for the Thermodynamics step in a SEAMM flowchart
"""

import logging

import gaussian_step

logger = logging.getLogger(__name__)


class ThermodynamicsParameters(gaussian_step.OptimizationParameters):
    """
    The control parameters for Thermodynamics.

    You need to replace the "time" entry in dictionary below these comments with the
    definitions of parameters to control this step. The keys are parameters for the
    current plugin,the values are dictionaries as outlined below.

    Examples
    --------
    ::

        parameters = {
            "time": {
                "default": 100.0,
                "kind": "float",
                "default_units": "ps",
                "enumeration": tuple(),
                "format_string": ".1f",
                "description": "Simulation time:",
                "help_text": ("The time to simulate in the dynamics run.")
            },
        }

    parameters : {str: {str: str}}
        A dictionary containing the parameters for the current step.
        Each key of the dictionary is a dictionary that contains the
        the following keys:

    parameters["default"] :
        The default value of the parameter, used to reset it.

    parameters["kind"] : enum()
        Specifies the kind of a variable. One of  "integer", "float", "string",
        "boolean", or "enum"

        While the "kind" of a variable might be a numeric value, it may still have
        enumerated custom values meaningful to the user. For instance, if the parameter
        is a convergence criterion for an optimizer, custom values like "normal",
        "precise", etc, might be adequate. In addition, any parameter can be set to a
        variable of expression, indicated by having "$" as the first character in the
        field. For example, $OPTIMIZER_CONV.

    parameters["default_units"] : str
        The default units, used for resetting the value.

    parameters["enumeration"]: tuple
        A tuple of enumerated values.

    parameters["format_string"]: str
        A format string for "pretty" output.

    parameters["description"]: str
        A short string used as a prompt in the GUI.

    parameters["help_text"]: str
        A longer string to display as help for the user.

    See Also
    --------
    Thermodynamics, TkThermodynamics, ThermodynamicsParameters, ThermodynamicsStep
    """

    parameters = {
        "optimize first": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Optimize structure first:",
            "help_text": (
                "Whether to optimize the structure before calculating the frequencies."
            ),
        },
        "T": {
            "default": 298.15,
            "kind": "float",
            "default_units": "K",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Temperature:",
            "help_text": "The temperature for the calculation.",
        },
        "P": {
            "default": 1.0,
            "kind": "float",
            "default_units": "atm",
            "enumeration": tuple(),
            "format_string": ".2f",
            "description": "Pressure:",
            "help_text": "The pressure for the calculation.",
        },
    }

    # Rules shared by the dialog and the flowchart builder (see seamm.Parameters):
    # the optimization's parameters apply only if optimizing first, and the
    # temperature and pressure are not used yet.
    unused = {
        **gaussian_step.OptimizationParameters.unused,
        "T": (
            "the Gaussian step does not use it yet; the thermodynamic functions are "
            "at Gaussian's default of 298.15 K"
        ),
        "P": (
            "the Gaussian step does not use it yet; the thermodynamic functions are "
            "at Gaussian's default of 1 atm"
        ),
    }

    def __init__(self, defaults={}, data=None):
        """
        Initialize the parameters, by default with the parameters defined above

        Parameters
        ----------
        defaults: dict
            A dictionary of parameters to initialize. The parameters
            above are used first and any given will override/add to them.
        data: dict
            A dictionary of keys and a subdictionary with value and units
            for updating the current, default values.

        Returns
        -------
        None
        """

        logger.debug("ThermodynamicsParameters.__init__")

        # The optimization's parameters apply only when optimizing first. This is
        # in the definitions, not edited afterwards, so that from_dict() keeps it.
        optimization = {}
        for key, definition in gaussian_step.OptimizationParameters.parameters.items():
            if key in self.unused:
                continue
            conditions = definition.get("applies_when") or {}
            optimization[key] = {
                **definition,
                "applies_when": {"optimize first": "yes", **conditions},
            }

        super().__init__(
            defaults={
                **ThermodynamicsParameters.parameters,
                **optimization,
                **defaults,
            },
            data=data,
        )
