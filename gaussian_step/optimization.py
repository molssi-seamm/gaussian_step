# -*- coding: utf-8 -*-

"""Setup and run Gaussian"""

import logging
import textwrap

from tabulate import tabulate

import gaussian_step
import seamm
import seamm.data
from seamm_util import Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __


try:
    from itertools import batched
except ImportError:
    from itertools import islice

    def batched(iterable, n):
        "Batch data into tuples of length n. The last batch may be shorter."
        # batched('ABCDEFG', 3) --> ABC DEF G
        if n < 1:
            raise ValueError("n must be at least one")
        it = iter(iterable)
        while batch := tuple(islice(it, n)):
            yield batch


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

        if P is None:
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

        if P["hessian"] == "default guess":
            added += " The Hessian will be initialized with a default guess."
        elif self.is_expr(P["hessian"]):
            added += " How the Hessian is initialized will be determined by {hessian}."
        elif P["hessian"] == "read from checkpoint":
            added += " The initial Hessian will be read from the checkpoint file."
        elif P["hessian"] == "from property on configuration":
            added += " The initial Hessian will be read from the configuration."
        elif P["hessian"] == "calculate":
            if self.is_expr(P["recalc hessian"]):
                added += " When to calculate the Hessian will be determined by "
                added += "{recalc hessian}."
            elif "every" in P["recalc hessian"]:
                added += " The Hessian will be recalculated every step."
            elif P["recalc hessian"] == "at beginning":
                added += " The Hessian will be calculated once at the beginning."
            elif P["recalc hessian"] == "HF at beginning":
                added += (
                    " The Hartree-Fock Hessian will be calculated once at the "
                    "beginning."
                )
            else:
                added += (
                    " The Hessian will be recalculated every {recalc hessian} steps."
                )

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

    def format_force_constants(
        self,
        n_atoms,
        force_constants,
        fc_units,
        energy=0.0,
        energy_units="E_h",
        gradients=None,
        gradient_units="E_h/Å",
    ):
        """Format the force constants for the Gaussian input.

        Parameters
        ----------
        n_atoms : int
            The number of atoms in the system
        force_constants : [float]
            The force constants in lower triangular order or as a square matrix
        fc_units: str
            The units of the force constants
        energy : float
            The energy, optional
        energy_units : str
            The units of the energy, defaults to Hartrees
        gradients : [float]
            The gradients, optional
        gradient_units : str
            The units of the gradients, defaults to Hartrees/Å

        Returns
        -------
        []
            The section lines

        Note
        ----
        From the Gaussian documentation:

        FCCards
        Requests that read the energy (although value is not used), cartesian forces and
        force constants from the input stream, as written out by Punch=Derivatives. The
        format for this input is:

        Energy  	Format (D24.16)
        Cartesian forces  	Lines of format (6F12.8)
        Force constants	Lines of format (6F12.8)

        The force constants are in lower triangular form: ((F(J,I),J=1,I),I=1,3Natoms),
        where 3Natoms is the number of Cartesian coordinates.
        """

        lines = []

        # The energy is not used, but is required
        tmp = Q_(energy, energy_units).m_as("E_h")
        lines.append(f"{tmp:24.16e}")

        # The gradients
        if gradients is None:
            tmp = [0.0] * 3 * n_atoms
        elif (
            len(gradients) == n_atoms
            and isinstance(gradients[0], (list, tuple))
            and len(gradients[0]) == 3
        ):
            tmp = [v for row in gradients for v in row]
        elif isinstance(gradients, (list, tuple)) and len(gradients) == 3 * n_atoms:
            tmp = gradients
        else:
            raise ValueError(
                f"Can't handle gradients of type {type(gradients)}. Wrong dimensions?"
            )

        fac = Q_(1.0, gradient_units).m_as("E_h/Å")
        for row in batched(tmp, 6):
            lines.append("".join([f"{fac * v:12.8f}" for v in row]))

        # The force constants
        ndof = 3 * n_atoms
        ntri = ndof * (ndof + 1) // 2
        if len(force_constants) == ntri:
            tmp = force_constants
        elif (
            len(force_constants) == ndof
            and isinstance(force_constants[0], (list, tuple))
            and len(force_constants[0]) == ndof
        ):
            # Square matrix
            tmp = [v for i, row in enumerate(force_constants) for v in row[: i + 1]]
        else:
            raise ValueError(
                f"Can't handle force constants of type {type(force_constants)}. "
                "Wrong dimensions?"
            )

        fac = Q_(1.0, fc_units).m_as("E_h/a_0^2")
        for row in batched(tmp, 6):
            lines.append("".join([f"{fac * v:12.8f}" for v in row]))

        # Section ends in blank line
        lines.append("")

        return lines

    def run(self, keywords=None, extra_sections={}):
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

        # Handle options for the Hessian
        if P["hessian"] == "default guess":
            pass
        elif P["hessian"] == "read from checkpoint":
            subkeywords.append("ReadFC")
        elif P["hessian"] == "from property on configuration":
            subkeywords.append("FCCards")
            properties = configuration.properties.get(
                "force constants*", include_system_properties=True
            )

            print("properties:")
            print("\t" + "\n\t".join(properties.keys()))

            possibilities = [*properties]
            if len(possibilities) == 0:
                raise ValueError(
                    "No force constants found in the configuration properties."
                )
            elif len(possibilities) == 1:
                d2E = properties[possibilities[0]]["value"]
                units = configuration.properties.units(possibilities[0])
            else:
                printer.normal(
                    "   Warning: found more than one set of force constants:"
                )
                printer.normal("\t" + "\n\t".join(possibilities))
                print.normal(f"    using {possibilities[-1]}.")
                d2E = properties[possibilities[-1]]["value"]
                units = configuration.properties.units(possibilities[-1])

            extra_sections["Initial force constants"] = self.format_force_constants(
                configuration.n_atoms,
                d2E,
                units,
            )

        elif P["hessian"] == "calculate":
            if P["recalc hessian"] == "every step":
                subkeywords.append("CalcAll")
            elif P["recalc hessian"] == "at beginning":
                subkeywords.append("CalcFC")
            elif P["recalc hessian"] == "HF at beginning":
                subkeywords.append("CalcHFFC")
            else:
                try:
                    tmp = int(P["recalc hessian"])
                    if tmp < 1:
                        raise ValueError(
                            "The Hessian recalculation interval must be > 0"
                        )
                except Exception:
                    raise ValueError(
                        "Calculating the Hessian must either be a keyword or integer."
                    )
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

        super().run(keywords=keywords, extra_sections=extra_sections)

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
            if "N steps optimization" in data:
                n_steps = data["N steps optimization"]
            else:
                n_steps = -1
            data["N steps optimization"] = n_steps
            if data["optimization is converged"]:
                text += f"The geometry optimization converged in {n_steps} steps."
            else:
                text += (
                    f"Warning: The geometry optimization did not converge in {n_steps} "
                    "steps."
                )

                printer.normal(__(text, indent=self.indent + 4 * " "))
                printer.normal("")
                text = ""

                table2 = {}
                for key in (
                    "maximum atom force",
                    "RMS atom force",
                    "maximum atom displacement",
                    "RMS atom displacement",
                ):
                    if key + " trajectory" not in data:
                        continue
                    table2[key] = [f"{v:.6f}" for v in data[key + " trajectory"]]
                    table2[key].append("-")
                    table2[key].append(f"{data[key + ' Threshold']:.6f}")
                if table2 != {}:
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
                    printer.normal(
                        textwrap.indent("\n".join(text_lines), self.indent + 7 * " ")
                    )

            # If calculating 2nd derivatives each step has the vibrations
            if "vibrational frequencies" in data:
                imaginary = [-v for v in data["vibrational frequencies"] if v < 0]
                data["N saddle modes"] = len(imaginary)
                target = P["target"].lower()
                if "trans" in target or target == "ts":
                    if len(imaginary) == 1:
                        text += (
                            " The structure is a transition state, as requested, "
                            "with one mode with negative curvature of "
                            f"{imaginary[0]:.2f} cm^-1."
                        )
                        data["TS frequency"] = round(imaginary[0], 2)
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
                not data["optimization is converged"]
                and not P["ignore unconverged optimization"]
            ):
                raise RuntimeError("Gaussian geometry optimization failed to converge.")
