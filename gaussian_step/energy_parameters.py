# -*- coding: utf-8 -*-
"""Global control parameters for Gaussian
"""

import logging

from gaussian_step import methods, dft_functionals
import seamm

logger = logging.getLogger("Gaussian")


class EnergyParameters(seamm.Parameters):
    """The control parameters for the energy."""

    parameters = {
        "input only": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": (
                "yes",
                "no",
            ),
            "format_string": "s",
            "description": "Write the input files and stop:",
            "help_text": "Don't run MOPAC. Just write the input files.",
        },
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
            "default": "B3LYP : hybrid functional of Becke and Lee, Yang, and Parr",
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
            "default": "B3LYP : hybrid functional of Becke and Lee, Yang, and Parr",
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
        "integral grid": {
            "default": "UltraFine",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": (
                "96,32,64",
                "SuperFine",
                "UltraFine",
                "Fine",
                "SG1",
                "Coarse",
            ),
            "format_string": "s",
            "description": "Numerical grid:",
            "help_text": (
                "The grid to use for numerical integrations in e.g. DFT. "
                "The 'UltraFine' is the normal default for Gaussian."
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
        "use symmetry": {
            "default": "yes",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("yes", "loose", "identify only", "no"),
            "format_string": "s",
            "description": "Use symmetry:",
            "help_text": "Whether to use symmetry, and if so how much.",
        },
        "calculate gradient": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Calculate gradient:",
            "help_text": "Whether to calculate the gradient:",
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
        "bond orders": {
            "default": "Wiberg",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("Wiberg", "none"),
            "format_string": "s",
            "description": "Calculate bond orders:",
            "help_text": "Whether to calculate the bond orders and if so how.",
        },
        "apply bond orders": {
            "default": "yes",
            "kind": "bool",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Apply bond orders to structure:",
            "help_text": (
                "Whether to use the calculated bond orders to update the structure"
            ),
        },
        "print basis set": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Print basis set:",
            "help_text": "Whether to print the basis set to the output.",
        },
        "save basis set": {
            "default": "no",
            "kind": "str",
            "default_units": "",
            "enumeration": ("yes", "no", "append to"),
            "format_string": "s",
            "description": "Save basis set:",
            "help_text": "Whether to save the basis set to a file.",
        },
        "basis set file": {
            "default": "basis.gbs",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "File for basis set:",
            "help_text": "The file for the output basis set.",
        },
        "file handling": {
            "default": "remove checkpoint files",
            "kind": "string",
            "format_string": "s",
            "enumeration": ("keep all", "remove all", "remove checkpoint files"),
            "description": "File handling",
            "help_text": (
                "How to handle files after a successful calculation."
                "or 'advanced', which shows everything."
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

    output_parameters = {
        "total density": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Plot total density:",
            "help_text": "Whether to plot the total charge density.",
        },
        "total spin density": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Plot total spin density:",
            "help_text": "Whether to plot the total spin density.",
        },
        "difference density": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Plot difference density:",
            "help_text": "Whether to plot the difference density.",
        },
        "orbitals": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Plot orbitals:",
            "help_text": "Whether to plot orbitals.",
        },
        "selected orbitals": {
            "default": "HOMO, LUMO",
            "kind": "string",
            "default_units": "",
            "enumeration": ("HOMO, LUMO", "-1, HOMO, LUMO, +1", "all"),
            "format_string": "",
            "description": "Selected orbitals:",
            "help_text": "Which orbitals to plot.",
        },
        "region": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": ("default", "explicit"),
            "format_string": "",
            "description": "Region:",
            "help_text": "The region for the plots",
        },
        "nx": {
            "default": 50,
            "kind": "integer",
            "default_units": "",
            "enumeration": None,
            "format_string": "",
            "description": "Grid:",
            "help_text": "Number of grid points in first direction",
        },
        "ny": {
            "default": 50,
            "kind": "integer",
            "default_units": "",
            "enumeration": None,
            "format_string": "",
            "description": "x",
            "help_text": "Number of grid points in second direction",
        },
        "nz": {
            "default": 50,
            "kind": "integer",
            "default_units": "",
            "enumeration": None,
            "format_string": "",
            "description": "x",
            "help_text": "Number of grid points in first direction",
        },
    }

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={
                **EnergyParameters.parameters,
                **EnergyParameters.output_parameters,
                **seamm.standard_parameters.structure_handling_parameters,
                **defaults,
            },
            data=data,
        )

        # Do any local editing of defaults
        tmp = self["configuration name"]
        tmp._data["enumeration"] = ["single-point with {model}", *tmp.enumeration[1:]]
        tmp.default = "keep current name"

        tmp = self["configuration name"]
        tmp._data["enumeration"] = ["single-point with {model}", *tmp.enumeration]
        tmp.default = "single-point with {model}"
