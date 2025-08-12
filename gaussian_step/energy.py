# -*- coding: utf-8 -*-

"""Setup and run Gaussian"""

import csv
import logging
import textwrap

import numpy as np
from openbabel import openbabel
from tabulate import tabulate

import gaussian_step
from .substep import Substep
import seamm
import seamm.data
from seamm_util import Q_, units_class
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

logger = logging.getLogger("Gaussian")
job = printing.getPrinter()
printer = printing.getPrinter("gaussian")

superscript = {
    "1": "\N{SUPERSCRIPT ONE}",
    "2": "\N{SUPERSCRIPT TWO}",
    "3": "\N{SUPERSCRIPT THREE}",
    "4": "\N{SUPERSCRIPT FOUR}",
    "5": "\N{SUPERSCRIPT FIVE}",
    "6": "\N{SUPERSCRIPT SIX}",
    "7": "\N{SUPERSCRIPT SEVEN}",
    "8": "\N{SUPERSCRIPT EIGHT}",
    "9": "\N{SUPERSCRIPT NINE}",
}


class Energy(Substep):
    def __init__(
        self,
        flowchart=None,
        title="Energy",
        extension=None,
        module=__name__,
        logger=logger,
    ):
        """Initialize the node"""

        logger.debug("Creating Energy {}".format(self))

        super().__init__(
            flowchart=flowchart,
            title=title,
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._method = None

        self._calculation = "energy"
        self._model = None
        self._metadata = gaussian_step.metadata
        self.parameters = gaussian_step.EnergyParameters()

        self.description = "A single point energy calculation"

    def description_text(self, P=None, calculation="Single-point energy"):
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
            if (
                functional in gaussian_step.dft_functionals
                and len(gaussian_step.dft_functionals[functional]["dispersion"]) > 1
                and P["dispersion"] != "none"
            ):
                text += f" with the {P['dispersion']} dispersion correction"
        else:
            text = f"{calculation} using {method_string}"

        if method_data != {}:
            if "nobasis" not in method_data or not method_data["nobasis"]:
                text += f" and the {P['basis']} basis set"
            if "freeze core?" in method_data and method_data["freeze core?"]:
                if P["freeze-cores"] == "no":
                    text += " with no core orbitals frozen."
                elif P["freeze-cores"] == "yes":
                    text += " with the core orbitals frozen."
                else:
                    text += (
                        f" with core orbitals frozen depending on {P['freeze-cores']}."
                    )
            else:
                text += "."
        else:
            text += f". If the method uses a basis set, it will be the {P['basis']}"
            text += " basis set."
            text += " If the method supports freezing core orbitals, it will be run "
            if P["freeze-cores"] == "no":
                text += " with no core orbitals frozen."
            elif P["freeze-cores"] == "yes":
                text += " with the core orbitals frozen."
            else:
                text += f" with core orbitals frozen depending on {P['freeze-cores']}."
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
            text += f"eigenstate will be determined by {P['spin-restricted']}"
        else:
            text += (
                " The spin will not be restricted and the result may not be a "
                "proper eigenstate."
            )

        if (
            isinstance(P["input only"], bool)
            and P["input only"]
            or P["input only"] == "yes"
        ):
            if type(self) is Energy:
                text += (
                    "\n\nThe input file will be written. No calculation will be run."
                )
        else:
            # Plotting
            plots = []
            for key in ("total density", "total spin density"):
                plt = P[key] if isinstance(P[key], bool) else P[key] == "yes"
                if plt:
                    plots.append(key)

            # if P["difference density"]:
            #     plots.append("difference density")

            key = "orbitals"
            plt = P[key] if isinstance(P[key], bool) else P[key] == "yes"
            if plt:
                if len(plots) > 0:
                    text += f"\nThe {', '.join(plots)} and orbitals "
                    text += f"{P['selected orbitals']} will be plotted."
                else:
                    text += f"\nThe orbitals {P['selected orbitals']} will be plotted."

            text += (
                " The final structure and any charges, etc. will "
                f"{P['structure handling'].lower()} "
            )

            confname = P["configuration name"]
            if confname == "use SMILES string":
                text += "using SMILES as its name."
            elif confname == "use Canonical SMILES string":
                text += "using canonical SMILES as its name."
            elif confname == "keep current name":
                text += "keeping the current name."
            elif confname == "optimized with {model}":
                text += "with 'optimized with <model>' as its name."
            elif confname == "use configuration number":
                text += "using the index of the configuration (1, 2, ...) as its name."
            else:
                confname = confname.replace("{model}", "<model>")
                text += f"with '{confname}' as its name."

        return self.header + "\n" + __(text, **P, indent=4 * " ").__str__()

    def run(self, keywords=None, extra_sections={}):
        """Run a single-point Gaussian calculation."""
        if keywords is None:
            keywords = set()

        _, starting_configuration = self.get_system_configuration(None)

        atnos = set(starting_configuration.atoms.atomic_numbers)
        max_atno = max(atnos)

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

        # If doing a single point, add the correct keyword for the job type
        if self.__class__.__name__ == "Energy":
            if P["calculate gradient"]:
                keywords.add("FORCE")
            else:
                keywords.add("SP")

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

        # And possible frozen core
        if "freeze core?" in method_data:
            if method_data["freeze core?"] and P["freeze-cores"] == "no":
                method_options = ["FULL"]
            else:
                method_options = ["FC"]
        else:
            method_options = []

        if len(method_options) == 0:
            mopts = ""
        else:
            mopts = "=(" + ", ".join(method_options) + ")"

        basis = P["basis"]
        if method == "DFT":
            functional = self.get_functional(P)
            functional_data = gaussian_step.dft_functionals[functional]
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"{functional_data['name']}{mopts}/{basis}")
                else:
                    keywords.add(f"RO{functional_data['name']}{mopts}/{basis}")
            else:
                keywords.add(f"U{functional_data['name']}{mopts}/{basis}")
            if len(functional_data["dispersion"]) > 1 and P["dispersion"] != "none":
                keywords.add(f"EmpiricalDispersion={P['dispersion']}")
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
                    keywords.add(f"RHF{mopts}/{basis}")
                else:
                    keywords.add(f"ROHF{mopts}/{basis}")
            else:
                keywords.add(f"UHF{mopts}/{basis}")
        elif method[0:2] == "MP":
            if "(" in method:
                method, option = method.split("(")
                method_options.append(option.rstrip(")"))
            if len(method_options) == 0:
                mopts = ""
            else:
                mopts = "=(" + ", ".join(method_options) + ")"
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"R{method}{mopts}/{basis}")
                else:
                    keywords.add(f"RO{method}{mopts}/{basis}")
            else:
                keywords.add(f"U{method}{mopts}/{basis}")
        elif method in ("QCISD", "QCISD(T)"):
            if "(" in method:
                method, option = method.split("(")
                method_options.append(option.rstrip(")"))
            if len(method_options) == 0:
                mopts = ""
            else:
                mopts = "=(" + ", ".join(method_options) + ")"
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"R{method}{mopts}/{basis}")
                else:
                    keywords.add(f"RO{method}{mopts}/{basis}")
            else:
                keywords.add(f"U{method}{mopts}/{basis}")
        elif method in ("CCSD", "CCSD(T)"):
            if "(" in method:
                method, option = method.split("(")
                method_options.append(option.rstrip(")"))
            if len(method_options) == 0:
                mopts = ""
            else:
                mopts = "=(" + ", ".join(method_options) + ")"
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"R{method}{mopts}/{basis}")
                else:
                    keywords.add(f"RO{method}{mopts}/{basis}")
            else:
                keywords.add(f"U{method}{mopts}/{basis}")
        elif method in ("CBS-4M", "CBS-QB3"):
            if self.gversion == "g09":
                if self.__class__ == Energy:
                    raise RuntimeError(
                        "G09 does not appear to be able to run the CBS methods without "
                        "optimizing the structure during the calculation."
                    )
                else:
                    if restricted and multiplicity != 1:
                        keywords.add(f"RO{method}{mopts}")
                    else:
                        keywords.add(f"{method}{mopts}")
            else:
                if self.__class__ == Energy:
                    mopts.append("NoOpt")
                    mopts = "=(" + ", ".join(method_options) + ")"
                    if restricted and multiplicity != 1:
                        keywords.add(f"RO{method}{mopts}")
                    else:
                        keywords.add(f"{method}{mopts}")
                else:
                    if restricted and multiplicity != 1:
                        keywords.add(f"RO{method}{mopts}")
                    else:
                        keywords.add(f"{method}{mopts}")
            if max_atno > 36:
                raise RuntimeError(
                    f"{method} cannot handle systems with atoms heavier than Kr (36)"
                )
        elif method == "CBS-APNO":
            if self.gversion == "g09":
                if self.__class__ == Energy:
                    raise RuntimeError(
                        "G09 does not appear to be able to run the CBS methods without "
                        "optimizing the structure during the calculation."
                    )
                else:
                    keywords.add(f"{method}{mopts}")
            else:
                if self.__class__ == Energy:
                    mopts.append("NoOpt")
                    mopts = "=(" + ", ".join(method_options) + ")"
                    keywords.add(f"{method}{mopts}")
                else:
                    keywords.add(f"{method}{mopts}")
            if max_atno > 18:
                raise RuntimeError(
                    f"{method} cannot handle systems with atoms heavier than Ar (18)"
                )
        elif method in (
            "G1",
            "G2",
            "G3",
            "G4",
            "G2MP2",
            "G3B3",
            "G3MP2",
            "G3MP2B3",
            "G4MP2",
        ):
            if self.gversion == "g09":
                if self.__class__ == Energy:
                    raise RuntimeError(
                        "G09 does not appear to be able to run the Gn methods without "
                        "optimizing the structure during the calculation."
                    )
                else:
                    keywords.add(f"{method}{mopts}")
            else:
                if self.__class__ == Energy:
                    mopts.append("NoOpt")
                    mopts = "=(" + ", ".join(method_options) + ")"
                    keywords.add(f"{method}{mopts}")
                else:
                    keywords.add(f"{method}{mopts}")
            if max_atno > 36:
                raise RuntimeError(
                    f"{method} cannot handle systems with atoms heavier than Kr (36)"
                )
            if method[:2] in ("G3", "G4"):
                for atno in range(21, 31):
                    if atno in atnos:
                        raise RuntimeError(
                            f"{method} can only handle elements 1-20 (H-Ca) and "
                            "31-36 (Ga-Kr)"
                        )
        elif method in ("AM1", "PM3", "PM3MM", "PM6", "PDDG", "PM7", "PM7MOPAC"):
            if restricted and multiplicity != 1:
                keywords.add(f"RO{method}{mopts}")
            else:
                keywords.add(f"{method}{mopts}")
        else:
            keywords.add(f"{method}{mopts}/{basis}")

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

        if self.__class__.__name__ == "Energy":
            if P["calculate gradient"]:
                keywords.add("FORCE")
            else:
                keywords.add("SP")

        if P["bond orders"] == "Wiberg":
            keywords.add("Pop=NBORead")
            extra_sections["NBO input"] = "$nbo bndidx $end\n"

        if P["print basis set"] or P["save basis set"].lower() != "no":
            keywords.add("GFInput")

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

    def analyze(self, indent="", data={}, table=None, P=None):
        """Parse the output and generating the text output and store the
        data in variables for other stages to access
        """
        if P is None:
            P = self.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )

        if "energy" not in data:
            text = "Gaussian did not produce the energy. Something is wrong!"
            printer.normal(__(text, indent=self.indent + 4 * " "))

        text = ""

        # Calculate the enthalpy of formation, if possible
        tmp_text = self.calculate_enthalpy_of_formation(data)
        if tmp_text != "":
            path = self.wd / "Thermochemistry.txt"
            path.write_text(tmp_text)

        if table is None:
            table = {
                "Property": [],
                "Value": [],
                "Units": [],
            }

        metadata = gaussian_step.metadata["results"]

        # The semiempirical energy is the enthalpy, not energy!
        key = "energy"
        if key in data:
            # # Figure out the method.
            method, method_data, method_string = self.get_method(P)

            if method in ("AM1", "PM3", "PM3MM", "PM6", "PDDG", "PM7", "PM7MOPAC"):
                tmp = data[key]
                mdata = metadata[key]
                table["Property"].append(
                    "\N{GREEK CAPITAL LETTER DELTA}fH\N{SUPERSCRIPT ZERO}"
                )
                table["Value"].append(f"{tmp:{mdata['format']}}")
                if "units" in mdata:
                    table["Units"].append(mdata["units"])
                else:
                    table["Units"].append("")
                table["Property"].append("")
                tmp = Q_(float(data[key]), "hartree").m_as("kcal/mol")
                table["Value"].append(f"{tmp:.2f}")
                table["Units"].append("kcal/mol")
                table["Property"].append("")
                tmp = Q_(float(data[key]), "hartree").m_as("kJ/mol")
                table["Value"].append(f"{tmp:.2f}")
                table["Units"].append("kJ/mol")

                keys = []
            else:
                if "DfH0" in data:
                    tmp = data["DfH0"]
                    table["Property"].append(
                        "\N{GREEK CAPITAL LETTER DELTA}fH\N{SUPERSCRIPT ZERO}"
                    )
                    table["Value"].append(f"{Q_(tmp, 'kJ/mol').m_as('kcal/mol'):.2f}")
                    table["Units"].append("kcal/mol")
                    table["Property"].append("")
                    table["Value"].append(f"{tmp:.2f}")
                    table["Units"].append("kJ/mol")
                keys = [
                    "H atomization",
                    "DfE0",
                    "E atomization",
                    "energy",
                ]

        keys.extend(
            [
                "E_0",
                "H",
                "F",
                "S",
                "T",
                "P",
                "virial ratio",
                "RMS density difference",
                "E qcisd_t",
                "E qcisd",
                "E ccsd_t",
                "E cc",
                "E mp5",
                "E mp4",
                "E mp4sdq",
                "E mp4dq",
                "E mp3",
                "E mp2",
                "E scf",
            ]
        )
        for key in keys:
            if key in data:
                tmp = data[key]
                mdata = metadata[key]
                table["Property"].append(key)
                table["Value"].append(f"{tmp:{mdata['format']}}")
                if "units" in mdata:
                    table["Units"].append(mdata["units"])
                else:
                    table["Units"].append("")

        keys = [
            ("S**2", "S\N{SUPERSCRIPT TWO}"),
            ("S**2 after annihilation", "S\N{SUPERSCRIPT TWO} after annihilation"),
            ("ideal S**2", "ideal S\N{SUPERSCRIPT TWO}"),
            ("symmetry group", "Symmetry"),
            ("symmetry group used", "Symmetry used"),
            ("state", "Electronic state"),
        ]
        if data["method"][0] == "U":
            for letter, symbol in (("alpha", "α"), ("beta", "β")):
                keys.extend(
                    [
                        (f"E {letter} homo", f"{symbol}-HOMO Energy"),
                        (f"E {letter} lumo", f"{symbol}-LUMO Energy"),
                        (f"E {letter} gap", f"{symbol}-Gap"),
                        (f"{letter} HOMO symmetry", f"{symbol}-HOMO Symmetry"),
                        (f"{letter} LUMO symmetry", f"{symbol}-LUMO Symmetry"),
                    ]
                )
        else:
            keys.extend(
                [
                    ("E alpha homo", "HOMO Energy"),
                    ("E alpha lumo", "LUMO Energy"),
                    ("E alpha gap", "Gap"),
                    ("alpha HOMO symmetry", "HOMO Symmetry"),
                    ("alpha LUMO symmetry", "LUMO Symmetry"),
                ]
            )
        keys.extend(
            [
                ("dipole moment magnitude", "Dipole moment"),
            ]
        )
        for key, name in keys:
            if key in data:
                tmp = data[key]
                if key == "state":
                    tmp = superscript[tmp[0]] + tmp[1:]
                mdata = metadata[key]
                table["Property"].append(name)
                table["Value"].append(f"{tmp:{mdata['format']}}")
                if "units" in mdata:
                    table["Units"].append(mdata["units"])
                else:
                    table["Units"].append("")

        for key, name in (
            ("cpu time", "CPU time"),
            ("elapsed time", "Elapsed time"),
        ):
            if key in data:
                tmp = data[key]
                table["Property"].append(name)
                if ":" in tmp:
                    units = ""
                else:
                    units = "s"
                    tmp = f"{float(tmp):.2f}"
                table["Value"].append(tmp)
                table["Units"].append(units)

        tmp = tabulate(
            table,
            headers="keys",
            tablefmt="rounded_outline",
            colalign=("center", "decimal", "left"),
            disable_numparse=True,
        )
        length = len(tmp.splitlines()[0])
        text_lines = []
        if "method" in data:
            if data["method"].startswith("U"):
                spin_text = "U-"
            elif data["method"].startswith("RO"):
                spin_text = "RO-"
            else:
                spin_text = "R-"
        else:
            spin_text = ""
        spin_state = data["requested spin state"]
        chg = data["charge"]
        if chg == 0:
            header = f"Results for {spin_text}{self.model} for the {spin_state} state"
        else:
            header = (
                f"Results for {spin_text}{self.model} for the {spin_state} state, "
                f"charge {chg}"
            )
        text_lines.append(header.center(length))
        text_lines.append(method_string.center(length))
        if method == "DFT":
            functional = self.get_functional(P)
            text_lines.append(functional.center(length))
        text_lines.append(tmp)

        if text != "":
            text = str(__(text, **data, indent=self.indent + 4 * " "))
            text += "\n\n"
        text += textwrap.indent("\n".join(text_lines), self.indent + 7 * " ")

        if "Composite/summary" in data:
            text += "\n\n\n"
            text += textwrap.indent(data["Composite/summary"], self.indent + 4 * " ")

        # Handle the basis set as requested
        if P["print basis set"] and "basis set" in data:
            name = data["basis set name"]
            text += "\n\n"
            text += textwrap.indent(
                f"{name} basis set:\n" + data["basis set"], self.indent + 4 * " "
            )

        save_basis_set = P["save basis set"].lower()
        if (
            "basis set" in data
            and save_basis_set == "yes"
            or "append" in save_basis_set
        ):
            filename = P["basis set file"].strip()
            path = self.file_path(filename)

            if "append" in save_basis_set:
                with path.open("a") as fd:
                    fd.write(data["basis set"])
            else:
                path.write_text(data["basis set"])

        # Update the structure. Gaussian may have reoriented.
        system, configuration = self.get_system_configuration()

        if "atomcoords" in data and P["save standard orientation"]:
            coords = data["atomcoords"][-1]
            xs = [xyz[0] for xyz in coords]
            ys = [xyz[1] for xyz in coords]
            zs = [xyz[2] for xyz in coords]
            configuration.atoms["x"][0:] = xs
            configuration.atoms["y"][0:] = ys
            configuration.atoms["z"][0:] = zs
        elif "Current cartesian coordinates" in data:
            factor = Q_(1, "a0").to("Å").magnitude
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

        if "atomcharges/mulliken" in data:
            text_lines = ["\n"]
            symbols = configuration.atoms.asymmetric_symbols
            atoms = configuration.atoms
            symmetry = configuration.symmetry

            # Add to atoms (in coordinate table)
            if "charge" not in atoms:
                atoms.add_attribute(
                    "charge", coltype="float", configuration_dependent=True
                )
            if symmetry.n_symops == 1:
                chgs = data["atomcharges/mulliken"]
            else:
                chgs, delta = symmetry.symmetrize_atomic_scalar(data["ATOM_CHARGES"])
                delta = np.array(delta)
                max_delta = np.max(abs(delta))
                text_lines.append(
                    "The maximum difference of the charges of symmetry related atoms "
                    f"was {max_delta:.4f}\n"
                )
            atoms["charge"][0:] = chgs

            # Print the charges and dump to a csv file
            chg_tbl = {
                "Atom": [*range(1, len(symbols) + 1)],
                "Element": symbols,
                "Charge": [],
            }
            with open(self.wd / "atom_properties.csv", "w", newline="") as fd:
                writer = csv.writer(fd)
                if "atomspins/mulliken" in data:
                    # Sum to atom spins...
                    spins = data["atomspins/mulliken"]

                    # Add to atoms (in coordinate table)
                    if "spin" not in atoms:
                        atoms.add_attribute(
                            "spin", coltype="float", configuration_dependent=True
                        )
                        if symmetry.n_symops == 1:
                            atoms["spin"][0:] = spins
                        else:
                            spins, delta = symmetry.symmetrize_atomic_scalar(spins)
                            atoms["spins"][0:] = spins
                            delta = np.array(delta)
                            max_delta = np.max(abs(delta))
                            text_lines.append(
                                " The maximum difference of the spins of symmetry "
                                f"related atoms was {max_delta:.4f}.\n"
                            )

                    header = "        Atomic charges and spins"
                    chg_tbl["Spin"] = []
                    writer.writerow(["Atom", "Element", "Charge", "Spin"])
                    for atom, symbol, q, s in zip(
                        range(1, len(symbols) + 1),
                        symbols,
                        chgs,
                        spins,
                    ):
                        q = f"{q:.3f}"
                        s = f"{s:.3f}"

                        writer.writerow([atom, symbol, q, s])

                        chg_tbl["Charge"].append(q)
                        chg_tbl["Spin"].append(s)
                else:
                    header = "        Atomic charges"
                    writer.writerow(["Atom", "Element", "Charge"])
                    for atom, symbol, q in zip(
                        range(1, len(symbols) + 1),
                        symbols,
                        chgs,
                    ):
                        q = f"{q:.2f}"
                        writer.writerow([atom, symbol, q])

                        chg_tbl["Charge"].append(q)
            if len(symbols) <= int(self.parent.options["max_atoms_to_print"]):
                text_lines.append(header)
                if "Spin" in chg_tbl:
                    text_lines.append(
                        tabulate(
                            chg_tbl,
                            headers="keys",
                            tablefmt="rounded_outline",
                            colalign=("center", "center", "decimal", "decimal"),
                            disable_numparse=True,
                        )
                    )
                else:
                    text_lines.append(
                        tabulate(
                            chg_tbl,
                            headers="keys",
                            tablefmt="rounded_outline",
                            colalign=("center", "center", "decimal"),
                            disable_numparse=True,
                        )
                    )
                text += "\n\n"
                text += textwrap.indent("\n".join(text_lines), self.indent + 7 * " ")

        # Bond orders, if calculated
        if "Wiberg bond order matrix" in data:
            text += self._bond_orders(
                P, data["Wiberg bond order matrix"], configuration
            )

        printer.normal(text)

        # Write the structure locally for use in density and orbital plots
        obConversion = openbabel.OBConversion()
        obConversion.SetOutFormat("sdf")
        obMol = configuration.to_OBMol(properties="*")
        title = f"SEAMM={system.name}/{configuration.name}"
        obMol.SetTitle(title)
        # Try to get the rotated coordinates
        if "atomcoords" in data:
            coords = data["atomcoords"][-1]
            for i, atom in enumerate(openbabel.OBMolAtomIter(obMol)):
                atom.SetVector(coords[i][0], coords[i][1], coords[i][2])
        sdf = obConversion.WriteString(obMol)
        path = self.wd / "structure.sdf"
        path.write_text(sdf)

        text = self.make_plots(data)
        if text != "":
            printer.normal(__(text, indent=self.indent + 4 * " "))

        text = seamm.standard_parameters.set_names(system, configuration, P, **data)

        printer.normal("")
        printer.normal(__(text, indent=self.indent + 4 * " "))
        printer.normal("")

        # Put any requested results into variables or tables
        self.store_results(data=data, create_tables=True, configuration=configuration)

    def _bond_orders(self, P, bond_order_matrix, configuration):
        """Analyze and print the bond orders, and optionally use for the bonding
        in the structure.

        Parameters
        ----------
        P : {}
            The control options for the step
        bond_order_matrix : [natoms * [natoms * float]]
            The square bond order matrix.
        configuration : molsystem.Configuration
            The configuration to put the bonds on, if requested.
        """
        text = ""
        n_atoms = configuration.n_atoms

        if n_atoms == 1:
            return "\n\n        No bonds, since there is only one atom.\n"

        bond_i = []
        bond_j = []
        bond_order = []
        bond_order_str = []
        orders = []
        for j, row in enumerate(bond_order_matrix):
            for i, order in enumerate(row):
                if i == j:
                    break
                if order > 0.4:
                    bond_i.append(i)
                    bond_j.append(j)
                    if order > 1.3 and order < 1.7:
                        bond_order.append(5)
                        bond_order_str.append("aromatic")
                    else:
                        bond_order.append(round(order))
                        bond_order_str.append(str(round(order)))
                    orders.append(order)

        if len(bond_order) > 0:
            symbols = configuration.atoms.symbols
            text_lines = []
            if len(symbols) <= int(self.parent.options["max_atoms_to_print"]):
                name = configuration.atoms.names
                table = {
                    "i": [name[i] for i in bond_i],
                    "j": [name[j] for j in bond_j],
                    "bond order": [f"{o:6.3f}" for o in orders],
                    "bond multiplicity": bond_order_str,
                }
                tmp = tabulate(
                    table,
                    headers="keys",
                    tablefmt="rounded_outline",
                    disable_numparse=True,
                    colalign=("center", "center", "center", "center"),
                )
                length = len(tmp.splitlines()[0])
                text_lines.append("\n")
                text_lines.append("Bond Orders".center(length))
                text_lines.append(tmp)
                text += "\n\n"
                text += textwrap.indent("\n".join(text_lines), 12 * " ")

            if P["apply bond orders"]:
                ids = configuration.atoms.ids
                iatoms = [ids[i] for i in bond_i]
                jatoms = [ids[j] for j in bond_j]
                configuration.new_bondset()
                configuration.bonds.append(i=iatoms, j=jatoms, bondorder=bond_order)
                text2 = (
                    "\nReplaced the bonds in the configuration with those from the "
                    "calculated bond orders.\n"
                )

                text += str(__(text2, indent=8 * " "))

        return text
