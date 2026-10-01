# -*- coding: utf-8 -*-
"""Global control parameters for Gaussian"""

import logging

from gaussian_step import methods, dft_functionals
import seamm

logger = logging.getLogger("Gaussian")

# All the dispersion corrections that any functional offers
dispersion_corrections = []
for _data in dft_functionals.values():
    for _dispersion in _data["dispersion"]:
        if _dispersion not in dispersion_corrections:
            dispersion_corrections.append(_dispersion)


def find_method(method_string):
    """The method given by its full name or its Gaussian keyword.

    Parameters
    ----------
    method_string : str
        The full name, e.g. "HF: Hartree-Fock self consistent field (SCF)", or the
        keyword, e.g. "HF".

    Returns
    -------
    (str, dict or None, str)
        The keyword, the method's metadata (None if the method is not known) and
        the full name.
    """
    if method_string in methods:
        return methods[method_string]["method"], methods[method_string], method_string
    for name, data in methods.items():
        if method_string == data["method"]:
            return data["method"], data, name
    return method_string, None, method_string


def find_functional(functional):
    """The full name of a functional given by its full name, the first part of it,
    or its Gaussian name; None if not known."""
    if functional in dft_functionals:
        return functional
    for name in dft_functionals:
        if functional == name.split(":")[0].strip():
            return name
    for name, data in dft_functionals.items():
        if functional == data["name"]:
            return name
    return None


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
                "6-311G",
                "6-311G*",
                "6-311G**",
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
        "checkpoint": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": ("default", "job:gaussian.chk"),
            "format_string": "s",
            "description": "Checkpoint file:",
            "help_text": "The checkpoint file to use.",
        },
        "initial checkpoint": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": ("default", "job:gaussian.chk"),
            "format_string": "s",
            "description": "Initial checkpoint file:",
            "help_text": "The initial checkpoint file to start with.",
        },
        "geometry": {
            "default": "current",
            "kind": "string",
            "default_units": "",
            "enumeration": ("current", "from configuration", "from checkpoint file"),
            "format_string": "s",
            "description": "Geometry:",
            "help_text": "Where to get the geometry.",
        },
        "initial wavefunction": {
            "default": "default",
            "kind": "string",
            "default_units": "",
            "enumeration": ("default", "Read", "Harris", "Core"),
            "format_string": "s",
            "description": "Initial wavefunction:",
            "help_text": "The starting wavefunction for the calculation.",
        },
        "method": {
            "applies_when": {"level": "recommended"},
            "default": "DFT: Kohn-Sham density functional theory",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in methods if methods[x]["level"] == "normal"],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use."),
        },
        "advanced_method": {
            "applies_when": {"level": {"not": "recommended"}},
            "default": "DFT: Kohn-Sham density functional theory",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in methods],
            "format_string": "s",
            "description": "Method:",
            "help_text": ("The computational method to use."),
        },
        "functional": {
            "applies_when": {"level": "recommended"},
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
            "applies_when": {"level": {"not": "recommended"}},
            "default": "B3LYP : hybrid functional of Becke and Lee, Yang, and Parr",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": [x for x in dft_functionals],
            "format_string": "s",
            "description": "DFT Functional:",
            "help_text": ("The exchange-correlation functional to use."),
        },
        "dispersion": {
            "default": "none",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": dispersion_corrections,
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
            "applies_when": {"bond orders": {"not": "none"}},
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
        "save standard orientation": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Save standard orientation:",
            "help_text": "Keep the standard orientation rather than input orientation.",
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
            "applies_when": {"save basis set": {"not": "no"}},
            "default": "basis.gbs",
            "kind": "string",
            "default_units": "",
            "enumeration": tuple(),
            "format_string": "s",
            "description": "File for basis set:",
            "help_text": "The file for the output basis set.",
        },
        "file handling": {
            "applies_when": {"input only": {"not": "yes"}},
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
        "save wfx": {
            "default": "no",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "",
            "description": "Write a wavefunction (wfx) file:",
            "help_text": (
                "Whether to write an AIMPAC wavefunction (.wfx) file for the "
                "calculation. This is required by a following Atomic Charges step "
                "using DDEC6 (Chargemol). The file is written to the job directory "
                "as '<step number>.wfx'."
            ),
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
            "applies_when": {"orbitals": "yes"},
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
            "applies_when": {"region": "explicit"},
            "default": 50,
            "kind": "integer",
            "default_units": "",
            "enumeration": None,
            "format_string": "",
            "description": "Grid:",
            "help_text": "Number of grid points in first direction",
        },
        "ny": {
            "applies_when": {"region": "explicit"},
            "default": 50,
            "kind": "integer",
            "default_units": "",
            "enumeration": None,
            "format_string": "",
            "description": "x",
            "help_text": "Number of grid points in second direction",
        },
        "nz": {
            "applies_when": {"region": "explicit"},
            "default": 50,
            "kind": "integer",
            "default_units": "",
            "enumeration": None,
            "format_string": "",
            "description": "x",
            "help_text": "Number of grid points in first direction",
        },
    }

    # Rules shared by the dialog and the flowchart builder (see seamm.Parameters).
    # The simple conditions are the "applies_when" entries above.

    methods_by_level = True
    """Whether the method is 'method' or 'advanced_method' depending on the level.
    False for steps that always use 'method'."""

    unused = {
        "ignore convergence": (
            "the Gaussian step does not use it; an SCF that does not converge is "
            "an error"
        ),
        "difference density": "the Gaussian step cannot plot it yet",
    }
    """Parameters that this kind of step never uses, with the reason."""

    def _value(self, key, values):
        return values.get(key, self[key].value)

    def method_parameter(self, values=None):
        """Which parameter gives the method: 'method' or 'advanced_method', or None
        if that depends on a variable."""
        if values is None:
            values = self.current_values()
        if not self.methods_by_level:
            return "method"
        level = self._value("level", values)
        if self._is_expr(level):
            return None
        return "method" if level == "recommended" else "advanced_method"

    def functional_parameter(self, values=None):
        """Which parameter gives the functional: 'functional' or
        'advanced_functional', or None if that depends on a variable."""
        if values is None:
            values = self.current_values()
        level = self._value("level", values)
        if self._is_expr(level):
            return None
        return "functional" if level == "recommended" else "advanced_functional"

    def _method(self, values):
        """The method's keyword, metadata and name; the keyword is None if the
        method is given by a variable, and the metadata None if it is not known."""
        key = self.method_parameter(values)
        if key is None:
            return None, None, None
        value = self._value(key, values)
        if self._is_expr(value):
            return None, None, value
        return find_method(value)

    def _is_dft(self, values):
        keyword, _, _ = self._method(values)
        return keyword is None or keyword == "DFT"

    def _functional(self, values):
        """The functional as given, and its full name (None if it is not known)."""
        key = self.functional_parameter(values)
        if key is None:
            return None, None
        value = self._value(key, values)
        if self._is_expr(value):
            return value, None
        return value, find_functional(value)

    def _dispersions(self, values):
        """The dispersion corrections the functional offers, or None if the
        functional is given by a variable."""
        value, name = self._functional(values)
        if value is None or self._is_expr(value):
            return None
        if name is None:
            return ()
        return tuple(dft_functionals[name]["dispersion"])

    def applies(self, key, values=None, _seen=None):
        """As seamm.Parameters.applies, plus: the functionals and integral grid
        only for DFT; the dispersion correction only for a DFT functional that
        offers one; freezing the core only for methods that can; the basis only
        for methods that use one; and nothing in this step's 'unused'."""
        if values is None:
            values = self.current_values()
        if key in self.unused:
            return False
        if not super().applies(key, values, _seen):
            return False
        if key in ("functional", "advanced_functional", "integral grid"):
            return self._is_dft(values)
        if key == "dispersion":
            if not self._is_dft(values):
                return False
            dispersions = self._dispersions(values)
            return dispersions is None or len(dispersions) > 1
        if key == "freeze-cores":
            _, data, _ = self._method(values)
            return data is None or bool(data.get("freeze core?", False))
        if key == "basis":
            _, data, _ = self._method(values)
            return data is None or not data.get("nobasis", False)
        return True

    def not_applicable_reason(self, key, values=None):
        """Why a parameter does not apply, for the builder's messages."""
        if values is None:
            values = self.current_values()
        if key in self.unused:
            return self.unused[key]
        reason = super().not_applicable_reason(key, values)
        if reason or self.applies(key, values):
            return reason
        _, _, method = self._method(values)
        if key in ("functional", "advanced_functional", "integral grid", "dispersion"):
            if not self._is_dft(values):
                return f"it applies only to DFT, and the method is {method!r}"
        if key == "dispersion":
            functional, name = self._functional(values)
            if name is None:
                return f"the functional {functional!r} is not one the step knows"
            return f"the functional {name!r} has no dispersion correction to choose"
        if key == "freeze-cores":
            return f"the method {method!r} has no frozen-core option"
        if key == "basis":
            return f"the method {method!r} does not use a basis set"
        return ""

    def choices(self, key, values=None):
        """The dispersion corrections are those the functional offers."""
        if values is None:
            values = self.current_values()
        if key == "dispersion":
            dispersions = self._dispersions(values)
            if dispersions:
                return dispersions
        return super().choices(key, values)

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        # The dispersion choices used to include 'DG2', a misspelling of 'GD2', which
        # Gaussian rejects; flowcharts saved then may hold it.
        if data is not None and isinstance(data.get("dispersion"), dict):
            if data["dispersion"].get("value") == "DG2":
                data["dispersion"]["value"] = "GD2"

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
        tmp = self["system name"]
        tmp._data["enumeration"] = ["single-point with {model}", *tmp.enumeration]
        tmp.default = "keep current name"

        tmp = self["configuration name"]
        tmp._data["enumeration"] = ["single-point with {model}", *tmp.enumeration]
        tmp.default = "single-point with {model}"
