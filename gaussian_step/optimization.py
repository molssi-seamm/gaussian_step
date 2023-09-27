# -*- coding: utf-8 -*-

"""Setup and run Gaussian"""

import logging
import textwrap

from tabulate import tabulate

import gaussian_step
import seamm
import seamm.data
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __
from seamm_util import Q_

logger = logging.getLogger("Gaussian")
job = printing.getPrinter()
printer = printing.getPrinter("gaussian")


class Optimization(gaussian_step.Energy):
    def __init__(
        self,
        flowchart=None,
        title="Optimization",
        extension=None,
        module=__name__,
        logger=logger,
    ):
        """Initialize the node"""

        logger.debug("Creating Optimization {}".format(self))

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._calculation = "optimization"
        self._model = None
        self._metadata = gaussian_step.metadata
        self.parameters = gaussian_step.OptimizationParameters()

        self.description = "A geometry optimization"

    def description_text(self, P=None, calculation="Geometry optimization"):
        """Prepare information about what this node will do"""

        if not P:
            P = self.parameters.values_to_dict()

        text = super().description_text(P=P, calculation=calculation)

        coordinates = P["coordinates"]
        added = f"\nThe geometry optimization will use {coordinates} coordinates,"
        added += " a {geometry convergence} convergence criterion, "
        if P["max geometry steps"] == "default":
            added += (
                "and the default maximum number of steps, which is based on the "
                "system size."
            )
        else:
            added += "and no more than {max geometry steps} steps."

        if P["recalc hessian"] != "never":
            added += " The Hessian will be recalculated every {recalc hessian}"
            added += " steps. Note that calculating the second derivatives is "
            added += "quite expensive!"

        return text + "\n" + __(added, **P, indent=4 * " ").__str__()

    def run(self, keywords=set()):
        """Run an optimization calculation with Gaussian"""
        _, configuration = self.get_system_configuration()

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        subkeywords = []
        convergence = gaussian_step.optimization_convergence[P["geometry convergence"]]
        if convergence != "":
            subkeywords.append(convergence)
        max_steps = P["max geometry steps"]
        if max_steps != "default":
            if "nAtoms" in max_steps:
                n_atoms = configuration.n_atoms
                max_steps = max_steps.replace("nAtoms", str(n_atoms))
                max_steps = eval(max_steps)
            else:
                max_steps = int(max_steps)
            subkeywords.append(f"MaxCycles={max_steps}")

        if P["recalc hessian"] == "every step":
            subkeywords.append("CalcAll")
        elif P["recalc hessian"] == "at beginning":
            subkeywords.append("CalcFC")
        elif P["recalc hessian"] == "HF at beginning":
            subkeywords.append("CalcHFFC")
        elif P["recalc hessian"] == "never":
            pass
        else:
            subkeywords.append(f"RecalcFC={P['recalc hessian']}")

        coordinates = P["coordinates"]
        if "GIC" in coordinates:
            subkeywords.append("GIC")
        elif coordinates in ("redundant", "cartesian"):
            subkeywords.append(coordinates.capitalize())
        else:
            raise RuntimeError(
                f"Don't recognize optimization coordinates '{coordinates}'"
            )

        if len(subkeywords) == 1:
            keywords.add(f"Opt={subkeywords[0]}")
        elif len(subkeywords) > 1:
            keywords.add(f"Opt=({','.join(subkeywords)})")

        super().run(keywords=keywords)

    def analyze(self, indent="", data={}, out=[], table=None):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        text = ""

        if table is None:
            table = {
                "Property": [],
                "Value": [],
                "Units": [],
            }

        # metadata = gaussian_step.metadata["results"]
        if "Total Energy" not in data:
            text += "Gaussian did not produce the energy. Something is wrong!"

        # Get the system & configuration
        _, configuration = self.get_system_configuration(None)

        if configuration.n_atoms == 1:
            text += "System is an atom, so nothing to optimize."
        else:
            # Information about the optimization
            n_steps = data["Optimization Number of geometries"][0]
            if data["Geometry Optimization Converged"]:
                text += f"The geometry optimization converged in {n_steps} steps."
            else:
                text += (
                    f"Warning: The geometry optimization did not converge in {n_steps} "
                    "steps."
                )
                table2 = {}
                for key in (
                    "Maximum Force",
                    "RMS Force",
                    "Maximum Displacement",
                    "RMS Displacement",
                ):
                    table2[key] = [f"{v:.6f}" for v in data[key + " Trajectory"]]
                    table2[key].append("-")
                    table2[key].append(f"{data[key + ' Threshold']:.6f}")
                tmp = tabulate(
                    table2,
                    headers="keys",
                    tablefmt="rounded_outline",
                    colalign=("decimal", "decimal", "decimal", "decimal"),
                    disable_numparse=True,
                )
                length = len(tmp.splitlines()[0])
                text_lines = []
                text_lines.append("Convergence".center(length))
                text_lines.append(tmp)

                printer.normal(__(text, indent=self.indent + 4 * " "))
                printer.normal("")
                text = ""
                printer.normal(
                    textwrap.indent("\n".join(text_lines), self.indent + 7 * " ")
                )

        # Update the structure
        if "Current cartesian coordinates" in data:
            factor = Q_(1, "a0").to("Å").magnitude
            system_db = self.get_variable("_system_db")
            configuration = system_db.system.configuration
            xs = []
            ys = []
            zs = []
            it = iter(data["Current cartesian coordinates"])
            for x in it:
                xs.append(factor * x)
                ys.append(factor * next(it))
                zs.append(factor * next(it))
            configuration.atoms["x"][0:] = xs
            configuration.atoms["y"][0:] = ys
            configuration.atoms["z"][0:] = zs
            text += " Updated the system with the structure from Gaussian."

        if text != "":
            text = str(__(text, **data, indent=self.indent + 4 * " "))
            text += "\n\n"
        printer.normal(text)

        super().analyze(data=data)

        if configuration.n_atoms > 1:
            if (
                not data["Geometry Optimization Converged"]
                and not P["ignore unconverged optimization"]
            ):
                raise RuntimeError("Gaussian geometry optimization failed to converge.")
