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

    def __init__(self, defaults={}, data=None):
        """Initialize the instance, by default from the default
        parameters given in the class"""

        super().__init__(
            defaults={
                **WavefunctionStabilityParameters.parameters,
                **EnergyParameters.parameters,
                **defaults,
            },
            data=data,
        )

        # Do any local editing of defaults
        tmp = self["method"]
        tmp._data["enumeration"] = (
            "DFT: Kohn-Sham density functional theory",
            "HF: Hartree-Fock self consistent field (SCF)",
        )
