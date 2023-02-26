# -*- coding: utf-8 -*-

"""Setup and run Gaussian"""

import logging
import textwrap

from tabulate import tabulate

import gaussian_step
import seamm
import seamm.data
from seamm_util import units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("gaussian")


class Energy(seamm.Node):
    def __init__(self, flowchart=None, title="Energy", extension=None):
        """Initialize the node"""

        logger.debug("Creating Energy {}".format(self))

        super().__init__(
            flowchart=flowchart, title=title, extension=extension, logger=logger
        )

        self._method = None

        self._calculation = "energy"
        self._model = None
        self._metadata = gaussian_step.metadata
        self.parameters = gaussian_step.EnergyParameters()

        self.description = "A single point energy calculation"

    @property
    def header(self):
        """A printable header for this section of output"""
        return "Step {}: {}".format(".".join(str(e) for e in self._id), self.title)

    @property
    def method(self):
        """The method ... HF, DFT, ... used."""
        return self._method

    @method.setter
    def method(self, value):
        self._method = value

    @property
    def version(self):
        """The semantic version of this module."""
        return gaussian_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return gaussian_step.__git_revision__

    def description_text(self, P=None, calculation_type="Single-point energy"):
        """Prepare information about what this node will do"""

        if not P:
            P = self.parameters.values_to_dict()

        if P["level"] == "recommended":
            method = P["method"]
        else:
            method = P["advanced_method"]

        if self.is_expr(method):
            text = f"{calculation_type} using method given by {method}."
        elif (
            method in gaussian_step.methods
            and gaussian_step.methods[method]["method"] == "DFT"
        ):
            if P["level"] == "recommended":
                functional = P["functional"]
            else:
                functional = P["advanced_functional"]
            basis = P["basis"]
            text = f"{calculation_type} using {method} using {functional}"
            if (
                len(gaussian_step.dft_functionals[functional]["dispersion"]) > 1
                and P["dispersion"] != "none"
            ):
                text += f" with the {P['dispersion']} dispersion correction"
            text += f", using the {basis} basis set."
        else:
            text = f"{calculation_type} using {method}."

        # Spin
        if P["spin-restricted"] == "default":
            text += (
                " The spin will be restricted to a pure eigenstate for singlets and "
                "unrestricted for other states in which case the result may not be "
                "a pure eigenstate."
            )
        elif P["spin-restricted"] == "yes":
            text += " The spin will be restricted to a pure eigenstate."
        elif self.is_expr(P["spin-restricted"]):
            text += " Whether the spin will be restricted to a pure "
            text += "eigenstate will be determined by {P['spin-restricted']}"
        else:
            text += (
                " The spin will not be restricted and the result  may not be a "
                "proper eigenstate."
            )

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def get_input(self, calculation_type="energy", restart=None):
        """Get the input for an energy calculation for Gaussian"""

        _, configuration = self.get_system_configuration(None)

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        self.description = []
        self.description.append(__(self.description_text(PP), **PP, indent=self.indent))

        keywords = []

        # Figure out what we are doing!
        if P["level"] == "recommended":
            method_string = P["method"]
        else:
            method_string = P["advanced_method"]

        method_data = gaussian_step.methods[method_string]
        method = method_data["method"]

        # How to handle spin restricted.
        multiplicity = configuration.spin_multiplicity
        spin_restricted = P["spin-restricted"]
        if spin_restricted == "default":
            if multiplicity == 1:
                restricted = True
            else:
                restricted = False
        elif spin_restricted == "yes":
            restricted = True
        else:
            restricted = False

        basis = P["basis"]
        if method == "DFT":
            if P["level"] == "recommended":
                functional = P["functional"]
            else:
                functional = P["advanced_functional"]
            functional_data = gaussian_step.dft_functionals[functional]
            if restricted:
                if multiplicity == 1:
                    keywords.append(f"R{functional_data['name']}/{basis}")
                else:
                    keywords.append(f"RO{functional_data['name']}/{basis}")
            else:
                keywords.append(f"U{functional_data['name']}/{basis}")
            if len(functional_data["dispersion"]) > 1 and P["dispersion"] != "none":
                keywords.append(f"EmpiricalDispersion={P['dispersion']}")
        else:
            if restricted:
                if multiplicity == 1:
                    keywords.append(f"RHF/{basis}")
                else:
                    keywords.append(f"ROHF/{basis}")
            else:
                keywords.append(f"UHF/{basis}")

        if "freeze core" in method_data:
            if method_data["freeze core?"] and P["freeze-cores"] == "no":
                keywords.append("FULL")

        if P["maximum iterations"] != "default":
            keywords.append(f"MaxCycle={P['maximum iterations']}")
        if P["convergence"] != "default":
            keywords.append("Conver={P['convergence']}")

        return keywords

    def analyze(self, indent="", data={}, out=[], table=None):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """

        # P = self.parameters.current_values_to_dict(
        #     context=seamm.flowchart_variables._data
        # )

        text = ""
        if table is None:
            table = {
                "Property": [],
                "Value": [],
                "Units": [],
            }

        metadata = gaussian_step.metadata["results"]
        if "Total Energy" not in data:
            text += "Gaussian did not produce the energy. Something is wrong!"

        for key in ("Total Energy", "Virial Ratio", "RMS Density"):
            tmp = data[key]
            mdata = metadata[key]
            table["Property"].append(key)
            table["Value"].append(f"{tmp:{mdata['format']}}")
            if "units" in mdata:
                table["Units"].append(mdata["units"])
            else:
                table["Units"].append("")

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "decimal", "left"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        text_lines.append("Results".center(length))
        text_lines.append(tmp)

        if text != "":
            text = str(__(text, **data, indent=self.indent + 4 * " "))
            text += "\n\n"
        text += textwrap.indent("\n".join(text_lines), self.indent + 7 * " ")

        printer.normal(text)
