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
        added = "\nThe geometry optimization is targeting "

        target = P["target"].lower()
        if self.is_expr(target):
            added += "a minimum or saddle point, depending on {target},"
        elif "min" in target:
            added += "the minimum"
        elif "trans" in target or target == "ts":
            added += "a transition state"
        elif "saddle" in target:
            added += "a saddle point with {saddle order} downhill directions"
        added += f" using {coordinates} coordinates,"
        added += " a {geometry convergence} convergence criterion, "
        if P["max geometry steps"] == "default":
            added += (
                "and the default maximum number of steps, which is based on the "
                "system size."
            )
        else:
            added += "and no more than {max geometry steps} steps."

        if P["recalc hessian"] == "never":
            pass
        elif self.is_expr(P["recalc hessian"]):
            added += " Whether and how to calculate the Hessian will be determined by "
            added += "{recalc hessian}."
        elif "every" in P["recalc hessian"]:
            added += " The Hessian will be recalculated every step."
        elif P["recalc hessian"] == "at beginning":
            added += " The Hessian will be calculated once at the beginning."
        elif P["recalc hessian"] == "HF at beginning":
            added += (
                " The Hartree-Fock Hessian will be calculated once at the beginning."
            )
        else:
            added += " The Hessian will be recalculated every {recalc hessian} steps."
        if P["recalc hessian"] != "never":
            added += " Note that calculating the second derivatives is quite expensive!"

        if (
            isinstance(P["input only"], bool)
            and P["input only"]
            or P["input only"] == "yes"
        ):
            if type(self) is Optimization:
                added += (
                    "\n\nThe input file will be written. No calculation will be run."
                )

        return text + "\n" + __(added, **P, indent=4 * " ").__str__()

    def run(self, keywords=None):
        """Run an optimization calculation with Gaussian"""
        if keywords is None:
            keywords = set()

        _, configuration = self.get_system_configuration()

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Set the attribute for writing just the input
        self.input_only = P["input only"]

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
            # Also need to use an IOP to set the max. Odd.
            # https://mattermodeling.stackexchange.com/questions/5087/ (continued)
            #       why-does-gaussian-ignore-the-opt-maxcycles-keyword-for-optimizations
            keywords.add(f"iop(1/152={max_steps})")

        # Handle the target for the optimization
        target = P["target"].lower()
        if "min" in target:
            pass
        elif "trans" in target or target == "ts":
            subkeywords.append("TS")
        elif "saddle" in target:
            subkeywords.append(f"Saddle={P['saddle order']}")

        if "min" not in target:
            if P["ignore curvature error"]:
                subkeywords.append("NoEigenTest")

        # Handle options for the calculation of the Hessian
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

    def analyze(self, indent="", data={}, out=[], table=None, P=None):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """
        if P is None:
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
        if "energy" not in data:
            text += "Gaussian did not produce the energy. Something is wrong!"

        # Get the system & configuration
        _, configuration = self.get_system_configuration(None)

        if configuration.n_atoms == 1:
            text += "System is an atom, so nothing to optimize."
        else:
            # Information about the optimization
            if "Optimization Number of geometries" in data:
                n_steps = data["Optimization Number of geometries"][0]
            else:
                n_steps = -1
            data["nsteps"] = n_steps
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

            # If calculating 2nd derivatives each step have the vibrations
            if "vibfreqs" in data:
                imaginary = [-v for v in data["vibfreqs"] if v < 0]
                data["number of saddle modes"] = len(imaginary)
                target = P["target"].lower()
                if "trans" in target or target == "ts":
                    if len(imaginary) == 1:
                        text += (
                            " The structure is a transition state, as requested, "
                            "with one mode with negative curvature of "
                            f"{imaginary[0]:.2f} cm^-1."
                        )
                        data["ts frequency"] = round(imaginary[0], 2)
                    elif len(imaginary) == 0:
                        text += (
                            " Optimization to a transition state was requested, "
                            "however, the structure has no modes with negative "
                            "curvature."
                        )
                    else:
                        freqs = ", ".join([f"{v:.2}" for v in imaginary])
                        text += (
                            " A transition state was requested, but the structure is "
                            f"a saddle point with {len(imaginary)} modes with negative "
                            f"curvature: {freqs}"
                        )
                else:
                    if len(imaginary) == 1:
                        text += (
                            " The structure is a transition state "
                            "with one mode with negative curvature of "
                            f"{imaginary[0]:.2f} cm^-1."
                        )
                    elif len(imaginary) == 0:
                        pass
                    else:
                        freqs = ", ".join([f"{v:.2}" for v in imaginary])
                        text += (
                            " The structure is "
                            f"a saddle point with {len(imaginary)} modes with negative "
                            f"curvature: {freqs}"
                        )
        if text != "":
            text = str(__(text, **data, indent=self.indent + 4 * " "))
            text += "\n\n"
            printer.normal(text)

        super().analyze(data=data, P=P)

        if configuration.n_atoms > 1:
            if (
                not data["Geometry Optimization Converged"]
                and not P["ignore unconverged optimization"]
            ):
                raise RuntimeError("Gaussian geometry optimization failed to converge.")
