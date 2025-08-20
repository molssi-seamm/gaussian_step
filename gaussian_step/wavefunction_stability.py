# -*- coding: utf-8 -*-

"""Setup and run Gaussian"""

import logging
import pprint  # noqa: F401
import textwrap

from tabulate import tabulate

import gaussian_step
import seamm
from seamm_util import units_class, Q_
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger("Gaussian")
job = printing.getPrinter()
printer = printing.getPrinter("gaussian")


class WavefunctionStability(gaussian_step.Energy):
    def __init__(
        self,
        flowchart=None,
        title="Wavefunction Stability",
        extension=None,
        module=__name__,
        logger=logger,
    ):
        """Initialize the node"""

        logger.debug("Creating WavefunctionStability {}".format(self))

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._method = None

        self._calculation = "wavefunction stability"
        self._model = None
        self._metadata = gaussian_step.metadata
        self.parameters = gaussian_step.WavefunctionStabilityParameters()

        self.description = "A wavefunction stability calculation"

    def description_text(self, P=None, calculation="Wavefunction stability analysis"):
        """Prepare information about what this node will do"""

        if not P:
            P = self.parameters.values_to_dict()

        method, method_data, method_string = self.get_method(P)

        tmp = method_string.split(":", 1)
        if len(tmp) > 1:
            method_string = f"{tmp[0].strip()} ({tmp[1].strip()})"

        if self.is_expr(method):
            text = f"{calculation} using method given by {method}"
        elif method == "DFT":
            functional = self.get_functional(P)

            tmp = functional.split(":", 1)
            if len(tmp) > 1:
                functional_string = f"{tmp[0].strip()} ({tmp[1].strip()})"
            else:
                functional_string = functional

            text = (
                f"{calculation} using {method_string} using {functional_string} with "
                "the {integral grid} grid for the numerical integration"
            )
        else:
            text = f"{calculation} using {method_string}"

        text += f" and the {P['basis']} basis set."

        # Spin
        if P["spin-restricted"] == "default":
            text += (
                " Initially, the spin will be restricted to a pure eigenstate for "
                "singlets and unrestricted for other states."
            )
        elif P["spin-restricted"] == "yes":
            text += " Initially the spin will be restricted to a pure eigenstate."
        elif self.is_expr(P["spin-restricted"]):
            text += " Whether the spin will be restricted initially to a pure "
            text += f"eigenstate will be determined by {P['spin-restricted']}"
        else:
            text += (
                " Initially, the spin will not be restricted and the result may not be"
                " a proper eigenstate."
            )
        text += (
            " The stability analysis may result in a lower energy wavefunction, which "
            "may not be proper eigenstate of spin, i.e. will be a UHF wavefunction."
        )

        if (
            isinstance(P["input only"], bool)
            and P["input only"]
            or P["input only"] == "yes"
        ):
            if type(self) is WavefunctionStability:
                text += (
                    "\n\nThe input file will be written. No calculation will be run."
                )
        else:
            text += "\n\n"
            text += seamm.standard_parameters.structure_handling_description(
                P, model="<model>"
            )

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def run(self, keywords=None, extra_sections={}):
        """Run a stability analysis calculation."""
        if keywords is None:
            keywords = set()

        _, starting_configuration = self.get_system_configuration(None)

        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )
        # Have to fix formatting for printing...
        PP = dict(P)
        for key in PP:
            if isinstance(PP[key], units_class):
                PP[key] = "{:~P}".format(PP[key])

        # Set the attribute for writing just the input
        self.input_only = P["input only"]

        # Print what we are doing
        printer.important(__(self.description_text(PP), indent=self.indent))

        # The stability keyword
        stability = P["stability analysis"].lower()
        if "check" in stability:
            keywords.add("Stable=RExt")
        elif "opt" in stability:
            keywords.add("Stable=Opt")

        # Sort out the checkpoint files
        initial_chkpt = P["initial checkpoint"]
        if initial_chkpt == "default":
            step_no = int(self._id[-1])
            if step_no == 1:
                initial_chkpt = None
            else:
                initial_chkpt = self.file_path(
                    f"{step_no - 1}.chk", relative_to=self.wd.parent
                )
                if not initial_chkpt.exists():
                    initial_chkpt = None
        else:
            initial_chkpt = self.file_path(initial_chkpt, relative_to=self.wd.parent)
            if not initial_chkpt.exists():
                tmp = P["initial checkpoint"]
                raise ValueError(
                    f"The requested initial checkpoint file '{tmp}' ({initial_chkpt}) "
                    "does not exist, so stopping."
                )

        chkpt = P["checkpoint"]
        if chkpt == "default":
            step_no = self._id[-1]
            chkpt = self.file_path(f"{step_no}.chk", relative_to=self.wd.parent)
        else:
            chkpt = self.file_path(chkpt, relative_to=self.wd.parent)

        wavefunction = P["initial wavefunction"].strip().title()
        if wavefunction == "Default":
            if initial_chkpt is not None:
                keywords.add("Guess=Read")
        else:
            keywords.add(f"Guess={wavefunction}")

        geometry = P["geometry"].strip().lower()
        if "check" in geometry or ("current" in geometry and initial_chkpt is not None):
            keywords.add("Geom=AllCheck")

        # Figure out what we are doing!
        method, method_data, method_string = self.get_method(P)

        # Citations
        if "citations" in method_data:
            for level, citation in method_data["citations"]:
                self.references.cite(
                    raw=self._bibliography[citation],
                    alias=citation,
                    module="gaussian_step",
                    level=level,
                    note=method_string,
                )

        if P["use symmetry"] == "loose":
            keywords.add("Symmetry=(Loose)")
        elif P["use symmetry"] == "identify only":
            keywords.add("NoSymmetry")
        elif P["use symmetry"] == "no":
            keywords.add("Symmetry=None")

        if P["maximum iterations"] != "default":
            keywords.add(f"MaxCycle={P['maximum iterations']}")
        if P["convergence"] != "default":
            keywords.add(f"Conver={P['convergence']}")

        keywords_save = set(keywords)

        # How to handle spin restricted.
        multiplicity = starting_configuration.spin_multiplicity
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
            functional = self.get_functional(P)
            functional_data = gaussian_step.dft_functionals[functional]
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"{functional_data['name']}/{basis}")
                else:
                    keywords.add(f"RO{functional_data['name']}/{basis}")
            else:
                keywords.add(f"U{functional_data['name']}/{basis}")
            # Set the numerical integration grid
            grid = P["integral grid"]
            if grid == "96,32,64":
                grid = "-96032"
            keywords.add(f"Integral(Grid={grid})")

            # Citations
            if "citations" in functional_data:
                for level, citation in functional_data["citations"]:
                    self.references.cite(
                        raw=self._bibliography[citation],
                        alias=citation,
                        module="gaussian_step",
                        level=level,
                        note=functional,
                    )
        elif method == "HF":
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"RHF/{basis}")
                else:
                    keywords.add(f"ROHF/{basis}")
            else:
                keywords.add(f"UHF/{basis}")

        if self._timing_data is not None:
            try:
                self._timing_data[6] = starting_configuration.to_smiles(
                    canonical=True, flavor="openbabel"
                )
            except Exception:
                self._timing_data[6] = ""
            try:
                self._timing_data[7] = starting_configuration.to_smiles(
                    canonical=True, hydrogens=True, flavor="openbabel"
                )
            except Exception:
                self._timing_data[7] = ""
            try:
                self._timing_data[8] = starting_configuration.formula[0]
            except Exception:
                self._timing_data[7] = ""
            try:
                self._timing_data[9] = str(starting_configuration.charge)
            except Exception:
                self._timing_data[9] = ""
            try:
                self._timing_data[10] = str(starting_configuration.spin_multiplicity)
            except Exception:
                self._timing_data[10] = ""

        data = self.run_gaussian(
            keywords,
            extra_sections=extra_sections,
            old_chkpt=initial_chkpt,
            chkpt=chkpt,
        )

        if not self.input_only:
            # Follow instructions for where to put the coordinates,
            system, configuration = self.get_system_configuration(
                P=P, same_as=starting_configuration, model=self.model
            )

            self.analyze(data=data, P=P)

            S2 = f"{data['S**2']:.4f}" if "S**2" in data else ""
            S2_annihilated = (
                f"{data['S**2 after annihilation']:.4f}"
                if "S**2 after annihilation" in data
                else ""
            )

            # Get the properties to save
            spin_data = [
                [
                    multiplicity,
                    data["energy"],
                    self.chkpt,
                    S2,
                    data["ideal S**2"],
                    S2_annihilated,
                    configuration.properties.get(),
                ]
            ]

        if P["test spin multiplicity"]:
            text = "Testing whether the "
            if multiplicity > 2:
                multiplicities = (multiplicity + 2, multiplicity - 2)
                text += (
                    f"{self.spin_state(multiplicity + 2)} or "
                    f"{self.spin_state(multiplicity - 2)} have a lower energy."
                )
            else:
                multiplicities = (multiplicity + 2,)
                text += f"{self.spin_state(multiplicity + 2)} has a lower energy."
            printer.normal(__(text, indent=self.indent + 4 * " "))

            chkpt_stem = self.chkpt.stem
            chkpt_save = self.chkpt
            id_save = self._id
            for multiplicity in multiplicities:
                configuration.spin_multiplicity = multiplicity
                printer.normal("")
                spin_state = self.spin_state(multiplicity).title()
                text = f"{spin_state} spin state:"
                printer.normal(__(text, indent=self.indent + 4 * " "))

                # Work in a subdirectory
                self._id = (*id_save, f"spin_{multiplicity}")

                # New chkpt file
                chkpt = self.chkpt.with_stem(f"{chkpt_stem}_{multiplicity}")
                self.chkpt = chkpt

                keywords = set(keywords_save)

                # How to handle spin restricted.
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
                    functional = self.get_functional(P)
                    functional_data = gaussian_step.dft_functionals[functional]
                    if restricted:
                        if multiplicity == 1:
                            keywords.add(f"{functional_data['name']}/{basis}")
                        else:
                            keywords.add(f"RO{functional_data['name']}/{basis}")
                    else:
                        keywords.add(f"U{functional_data['name']}/{basis}")
                    # Set the numerical integration grid
                    grid = P["integral grid"]
                    if grid == "96,32,64":
                        grid = "-96032"
                    keywords.add(f"Integral(Grid={grid})")

                    # Citations
                    if "citations" in functional_data:
                        for level, citation in functional_data["citations"]:
                            self.references.cite(
                                raw=self._bibliography[citation],
                                alias=citation,
                                module="gaussian_step",
                                level=level,
                                note=functional,
                            )
                elif method == "HF":
                    if restricted:
                        if multiplicity == 1:
                            keywords.add(f"RHF/{basis}")
                        else:
                            keywords.add(f"ROHF/{basis}")
                    else:
                        keywords.add(f"UHF/{basis}")

                if self._timing_data is not None:
                    try:
                        self._timing_data[6] = starting_configuration.to_smiles(
                            canonical=True, flavor="openbabel"
                        )
                    except Exception:
                        self._timing_data[6] = ""
                    try:
                        self._timing_data[7] = starting_configuration.to_smiles(
                            canonical=True, hydrogens=True, flavor="openbabel"
                        )
                    except Exception:
                        self._timing_data[7] = ""
                    try:
                        self._timing_data[8] = starting_configuration.formula[0]
                    except Exception:
                        self._timing_data[7] = ""
                    try:
                        self._timing_data[9] = str(starting_configuration.charge)
                    except Exception:
                        self._timing_data[9] = ""
                    try:
                        self._timing_data[10] = str(multiplicity)
                    except Exception:
                        self._timing_data[10] = ""

                try:
                    data = self.run_gaussian(
                        keywords,
                        extra_sections=extra_sections,
                        spin_multiplicity=multiplicity,
                        old_chkpt=chkpt_save,
                        chkpt=chkpt,
                    )

                    self.analyze(data=data, P=P)
                except RuntimeError:
                    text = "Gaussian failed! Skipping this state and hoping..."
                    printer.normal(__(text, indent=self.indent + 8 * " "))
                else:
                    S2 = f"{data['S**2']:.4f}" if "S**2" in data else ""
                    S2_annihilated = (
                        f"{data['S**2 after annihilation']:.4f}"
                        if "S**2 after annihilation" in data
                        else ""
                    )
                    spin_data.append(
                        [
                            multiplicity,
                            data["energy"],
                            self.chkpt,
                            S2,
                            data["ideal S**2"],
                            S2_annihilated,
                            configuration.properties.get(),
                        ]
                    )

                # Set the id back to its original value
                self._id = id_save

            # Find the lowest energy structure, and make that the default chkpoint
            spin_data.sort(key=lambda E: E[1])

            multiplicity = spin_data[0][0]
            best_chk = spin_data[0][2]
            if best_chk != chkpt_save:
                chkpt_save.rename(
                    f"{chkpt_stem}_{starting_configuration.spin_multiplicity}"
                )
                best_chk.rename(chkpt_save)

            # Update the configuration with the correct data
            for name, tmp in spin_data[0][6].items():
                configuration.properties.put(name, tmp["value"])
            configuration.spin_multiplicity = multiplicity

            # And print a table of results
            state = self.spin_state(multiplicity)
            Emin = spin_data[0][1]
            text = f"The {state} state has the lowest energy: {Emin:.6f} E_h.\n"
            printer.normal(__(text, indent=self.indent + 4 * " "))
            factor = Q_("E_h").m_as("kJ/mol")

            table = {
                "Spin State": [self.spin_state(v[0]).title() for v in spin_data],
                "Energy": [f"{factor*(v[1]-Emin):.1f}" for v in spin_data],
                "Units": ["kJ/mol"] * len(spin_data),
                "S\N{SUPERSCRIPT TWO}": [v[3] for v in spin_data],
                "annihilated S\N{SUPERSCRIPT TWO}": [v[5] for v in spin_data],
                "Ideal S\N{SUPERSCRIPT TWO}": [f"{v[4]:.4f}" for v in spin_data],
            }

            tmp = tabulate(
                table,
                headers="keys",
                tablefmt="rounded_outline",
                colalign=("center", "decimal", "center", "decimal", "decimal"),
                disable_numparse=True,
            )
            length = len(tmp.splitlines()[0])
            text_lines = []
            text_lines.append("Spin states examined".center(length))
            text_lines.append(tmp)
            text = textwrap.indent("\n".join(text_lines), self.indent + 4 * " ")
            printer.normal("")
            printer.normal(text)
            printer.normal("")

    def analyze(self, indent="", data={}, table=None, P=None):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """

        stability = P["stability analysis"].lower()

        if "wavefunction is stable" not in data:
            text = (
                "The wavefunction stability analysis did not work for some reason. "
                "Please check the Gaussian output in output.txt for more information."
            )
            printer.normal(__(text, indent=self.indent + 4 * " "))
        else:
            spin_state = data["requested spin state"]
            chg = data["charge"]
            if chg == 0:
                text = (
                    f"The stability analysis found that the initial {spin_state} "
                    "wavefunction "
                )
                text = (
                    f"The stability analysis found that the initial {spin_state} "
                    f"wavefunction for a system with charge of {chg} "
                )

            status = data["wavefunction stability"][0]
            if status == "stable":
                text += "was stable."
            else:
                text += f"had an {status} instability."
            printer.normal(__(text, indent=self.indent + 4 * " "))
            printer.normal("")
            printer.normal(
                __(
                    data["wavefunction stability text"][0],
                    indent=self.indent + 8 * " ",
                )
            )

            if "opt" in stability:
                # Requested to optimize the wavefunction
                nsteps = len(data["wavefunction stability"]) - 1
                if nsteps > 0:
                    spin_state = data["requested spin state"]
                    if data["wavefunction is stable"]:
                        if nsteps == 1:
                            text = (
                                "The analysis took one step to find a "
                                f"stable {spin_state} wavefunction."
                            )
                        else:
                            text = (
                                f"The analysis took another {nsteps} steps to find "
                                f"a stable {spin_state} wavefunction."
                            )
                    else:
                        if nsteps == 1:
                            text = (
                                "The analysis took one step but could not find a "
                                f"stable {spin_state} wavefunction."
                            )
                        else:
                            text = (
                                f"The analysis took another {nsteps} steps but could "
                                f"not find a stable {spin_state} wavefunction."
                            )
                    printer.normal(__(text, indent=self.indent + 4 * " "))
                    step = 0
                    for lines, status in zip(
                        data["wavefunction stability text"][1:],
                        data["wavefunction stability"][1:],
                    ):
                        step += 1
                        printer.normal(
                            __(f"\nStep {step}:", indent=self.indent + 8 * " ")
                        )

                        text = "The stability analysis found that the wavefunction "
                        if status == "stable":
                            text += "was stable."
                        else:
                            text += f"had an {status} instability."

                        printer.normal(__(text, indent=self.indent + 12 * " "))
                        printer.normal(__(lines, indent=self.indent + 12 * " "))

        printer.normal("")
        if data["wavefunction is stable"]:
            text = "The final wavefunction is stable, i.e. a minimum."
        else:
            text = (
                "The final wavefunction is NOT stable, i.e. there is a lower energy "
                "solution than the one found."
            )
        printer.normal(__(text, indent=self.indent + 4 * " "))
        printer.normal("")

        super().analyze(indent, data, table, P)
