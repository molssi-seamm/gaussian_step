# -*- coding: utf-8 -*-
"""Global control parameters for Gaussian"""

import logging

from .energy_parameters import EnergyParameters

# import gaussian_step
# import seamm

logger = logging.getLogger("Gaussian")


class WavefunctionStabilityParameters(EnergyParameters):
    """The control parameters for the energy."""

    parameters = {
        "stability analysis": {
            "default": "reoptimize if unstable",
            "kind": "enumeration",
            "default_units": "",
            "enumeration": ("check for stability", "reoptimize if unstable"),
            "format_string": "s",
            "description": "Wavefunction analysis:",
            "help_text": (
                "Whether to check the stability of the wavefunction, and optionally "
                "to reoptimize it with lowered symmetry."
            ),
        },
        "test spin multiplicity": {
            "default": "yes",
            "kind": "boolean",
            "default_units": "",
            "enumeration": ("yes", "no"),
            "format_string": "s",
            "description": "Test spin multiplicity:",
            "help_text": (
                "Whether to check for lower energy solutions with different spin "
                "multiplicity, e.g. triplets for a singlets"
            ),
        },
    }

    # Rules shared by the dialog and the flowchart builder (see seamm.Parameters).
    # The method is always 'method', whatever the level (the dialog copies it to
    # 'advanced_method', which the calculation reads at the advanced level), and
    # the analysis does not use these Energy settings.
    methods_by_level = False

    unused = {
        **EnergyParameters.unused,
        "advanced_method": "this step always takes the method from 'method'",
        "calculate gradient": "only the Energy sub-step uses it",
        "dispersion": "the stability analysis does not use a dispersion correction",
        "freeze-cores": "the stability analysis is only for HF and DFT",
        "bond orders": "the stability analysis does not calculate them",
        "apply bond orders": "the stability analysis does not calculate bond orders",
        "print basis set": "the stability analysis does not print the basis set",
        "save basis set": "the stability analysis does not save the basis set",
        "basis set file": "the stability analysis does not save the basis set",
    }

    def implied(self, values=None):
        """'advanced_method' follows 'method', which is what the calculation reads
        at the advanced level."""
        if values is None:
            values = self.current_values()
        result = super().implied(values)
        method = self._value("method", values)
        if self._value("advanced_method", values) != method:
            result["advanced_method"] = method
        return result

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        # Only HF and DFT, at either level. This is in the definition, not edited
        # afterwards, so that from_dict() keeps it.
        method = {
            **EnergyParameters.parameters["method"],
            "enumeration": (
                "DFT: Kohn-Sham density functional theory",
                "HF: Hartree-Fock self consistent field (SCF)",
            ),
            "applies_when": None,
        }

        super().__init__(
            defaults={
                **WavefunctionStabilityParameters.parameters,
                **EnergyParameters.parameters,
                "method": method,
                **defaults,
            },
            data=data,
        )
