# -*- coding: utf-8 -*-
"""Global control parameters for Gaussian
"""

import logging

from gaussian_step import methods, dft_functionals
import seamm

logger = logging.getLogger(__name__)


class EnergyParameters(seamm.Parameters):
    """The control parameters for the energy."""

    parameters = {
        "level": {
            "default": "recommended",
            "kind": "string",
            "format_string": "s",
            "enumeration": ("recommended", "advanced"),
            "description": "The level of disclosure in the interface",
            "help_text": (
                "How much detail to show in the GUI. Currently 'recommended' "
                "or 'advanced', which shows everything."
            ),
        },
        "basis": {
            "default": "6-31G**",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "6-31G",
                "6-31G*",
                "6-31G**",
                "cc-pVDZ",
                "cc-pVTZ",
                "cc-pVQZ",
                "Def2SV",
                "Def2SVP",
                "Def2SVPP",
                "Def2TZP",
            ),
            "format_string": "s",
            "description": "Basis:",
            "help_text": ("The basis set to use."),
        },
        "method": {
            "default": "DFT: Kohn-Sham density functional theory",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in methods if methods[x]["level"] == "normal"],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use."),
        },
        "advanced_method": {
            "default": "DFT: Kohn-Sham density functional theory",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in methods],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use."),
        },
        "functional": {
            "default": "B3LYP hybrid GGA exchange-correlation functional",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [
                x for x in dft_functionals if dft_functionals[x]["level"] == "normal"
            ],
            "format_string": "s",
            "description": "DFT Functional:",
            "help_text": ("The exchange-correlation functional to use."),
        },
        "advanced_functional": {
            "default": "B3LYP hybrid GGA exchange-correlation functional",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in dft_functionals],
            "format_string": "s",
            "description": "DFT Functional:",
            "help_text": ("The exchange-correlation functional to use."),
        },
        "dispersion": {
            "default": "GD3BJ",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ["none", "GD3BJ", "GD3", "DG2"],
            "format_string": "s",
            "description": "Dispersion correction:",
            "help_text": ("The dispersion correction to use."),
        },
        "spin-restricted": {
            "default": "default",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("default", "yes", "no"),
            "format_string": "s",
            "description": "Spin-restricted:",
            "help_text": (
                "Whether to restrict the spin (RHF, ROHF, RKS) or not "
                "(UHF, UKS)."
                " Default is restricted for singlets, unrestricted otherwise."
            ),
        },
        "convergence": {
            "default": "default",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "s",
            "description": "Energy convergence criterion:",
            "help_text": (
                "Criterion for convergence of the RMS of the density (10^-N) and "
                "maximum change in the density matrix (10^-(N+2))."
            ),
        },
        "maximum iterations": {
            "default": "default",
            "kind": "integer",
            "default_units": "",
            "enumeration": ("default",),
            "format_string": "s",
            "description": "Maximum iterations:",
            "help_text": "Maximum number of SCF iterations.",
        },
        "ignore convergence": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Ignore lack of convergence:",
            "help_text": (
                "Whether to ignore lack of convergence in the SCF. Otherwise, "
                "an error is thrown."
            ),
        },
        "freeze-cores": {
            "default": "yes",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Freeze core orbitals:",
            "help_text": (
                "Whether to freeze the core orbitals in correlated " "methods"
            ),
        },
        "results": {
            "default": {},
            "kind": "dictionary",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "",
            "description": "results",
            "help_text": ("The results to save to variables or in " "tables. "),
        },
        "create tables": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Create tables as needed:",
            "help_text": (
                "Whether to create tables as needed for "
                "results being saved into tables."
            ),
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={**EnergyParameters.parameters, **defaults}, data=data
        )
