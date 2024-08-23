# -*- coding: utf-8 -*-

"""Setup and run Gaussian"""

import csv
import logging
from pathlib import Path
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

        if P["level"] == "recommended":
            method = P["method"]
        else:
            method = P["advanced_method"]

        if not self.is_expr(method) and method not in gaussian_step.methods:
            # See if it matches the keyword part
            for key, mdata in gaussian_step.methods.items():
                if method == mdata["method"]:
                    method = key
                    if P["level"] == "recommended":
                        self.parameters["method"].value = method
                    else:
                        self.parameters["advanced_method"].value = method

        if self.is_expr(method):
            text = f"{calculation} using method given by {method}."
        elif (
            method in gaussian_step.methods
            and gaussian_step.methods[method]["method"] == "DFT"
        ):
            if P["level"] == "recommended":
                functional = P["functional"]
            else:
                functional = P["advanced_functional"]
            basis = P["basis"]
            text = f"{calculation} using {method} using {functional}"
            if (
                functional in gaussian_step.dft_functionals
                and len(gaussian_step.dft_functionals[functional]["dispersion"]) > 1
                and P["dispersion"] != "none"
            ):
                text += f" with the {P['dispersion']} dispersion correction"
            text += f", using the {basis} basis set."
        else:
            text = f"{calculation} using {method}."

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
            if P["total density"]:
                plots.append("total density")
            # if P["difference density"]:
            #     plots.append("difference density")
            if P["total spin density"]:
                plots.append("spin density")
            if P["orbitals"]:
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

    def run(self, keywords=None, extra_lines=None):
        """Run a single-point Gaussian calculation."""
        if keywords is None:
            keywords = set()
        if extra_lines is None:
            extra_lines = []

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

        # If doing a single point, add the correct keyword for the job type
        if self.__class__.__name__ == "Energy":
            if P["calculate gradient"]:
                keywords.add("FORCE")
            else:
                keywords.add("SP")

        # Figure out what we are doing!
        if P["level"] == "recommended":
            method_string = P["method"]
        else:
            method_string = P["advanced_method"]

        # If we don't recognize the string presume (hope?) it is a Gaussian method
        if method_string in gaussian_step.methods:
            method_data = gaussian_step.methods[method_string]
            method = method_data["method"]
        else:
            # See if it matches the keyword part
            for key, mdata in gaussian_step.methods.items():
                if method_string == mdata["method"]:
                    method_string = key
                    method_data = mdata
                    method = method_data["method"]
                    break
            else:
                method_data = {}
                method = method_string

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
            if P["level"] == "recommended":
                functional = P["functional"]
            else:
                functional = P["advanced_functional"]
            functional_data = gaussian_step.dft_functionals[functional]
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"R{functional_data['name']}/{basis}")
                else:
                    keywords.add(f"RO{functional_data['name']}/{basis}")
            else:
                keywords.add(f"U{functional_data['name']}/{basis}")
            if len(functional_data["dispersion"]) > 1 and P["dispersion"] != "none":
                keywords.add(f"EmpiricalDispersion={P['dispersion']}")
        elif method == "HF":
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"RHF/{basis}")
                else:
                    keywords.add(f"ROHF/{basis}")
            else:
                keywords.add(f"UHF/{basis}")
            # keywords.add("Guess=Mix")
        elif method[0:2] == "MP":
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"R{method}/{basis}")
                else:
                    keywords.add(f"RO{method}/{basis}")
            else:
                keywords.add(f"U{method}/{basis}")
        elif method in ("CCSD", "CCSD(T)"):
            if restricted:
                if multiplicity == 1:
                    keywords.add(f"R{method}/{basis}")
                else:
                    keywords.add(f"RO{method}/{basis}")
            else:
                keywords.add(f"U{method}/{basis}")
        elif method in ("CBS-4M", "CBS-QB3"):
            if self.gversion == "g09":
                if self.__class__ == Energy:
                    raise RuntimeError(
                        "G09 does not appear to be able to run the CBS methods without "
                        "optimizing the structure during the calculation."
                    )
                else:
                    if restricted and multiplicity != 1:
                        keywords.add(f"RO{method}")
                    else:
                        keywords.add(f"{method}")
            else:
                if self.__class__ == Energy:
                    if restricted and multiplicity != 1:
                        keywords.add(f"RO{method}=NoOpt")
                    else:
                        keywords.add(f"{method}=(NoOpt)")
                else:
                    if restricted and multiplicity != 1:
                        keywords.add(f"RO{method}")
                    else:
                        keywords.add(f"{method}")
        elif method == "CBS-APNO":
            if self.gversion == "g09":
                if self.__class__ == Energy:
                    raise RuntimeError(
                        "G09 does not appear to be able to run the CBS methods without "
                        "optimizing the structure during the calculation."
                    )
                else:
                    keywords.add(f"{method}")
            else:
                if self.__class__ == Energy:
                    keywords.add(f"{method}=NoOpt")
                else:
                    keywords.add(f"{method}")
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
                        "G09 does not appear to be able to run the CBS methods without "
                        "optimizing the structure during the calculation."
                    )
                else:
                    keywords.add(f"{method}")
            else:
                if self.__class__ == Energy:
                    keywords.add(f"{method}=(NoOpt)")
                else:
                    keywords.add(f"{method}")
        else:
            keywords.add(f"{method}/{basis}")

        if P["use symmetry"] == "loose":
            keywords.add("Symmetry=Loose")
        elif P["use symmetry"] == "identify only":
            keywords.add("NoSymmetry")
        elif P["use symmetry"] == "no":
            keywords.add("Symmetry=None")

        if "freeze core" in method_data:
            if method_data["freeze core?"] and P["freeze-cores"] == "no":
                keywords.add("FULL")

        if P["maximum iterations"] != "default":
            keywords.add(f"MaxCycle={P['maximum iterations']}")
        if P["convergence"] != "default":
            keywords.add("Conver={P['convergence']}")

        if P["bond orders"] == "Wiberg":
            keywords.add("Pop=NBORead")
            extra_lines.append("$nbo bndidx $end")

        data = self.run_gaussian(keywords, extra_lines=extra_lines)

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

        text = ""
        if table is None:
            table = {
                "Property": [],
                "Value": [],
                "Units": [],
            }

        metadata = gaussian_step.metadata["results"]
        if "energy" not in data:
            text += "Gaussian did not produce the energy. Something is wrong!"
            printer.normal(__(text, indent=self.indent + 4 * " "))

        for key in (
            "energy",
            "Virial Ratio",
            "RMS Density",
            "Cluster Energy with triples",
            "Cluster Energy",
            "MP5 Energy",
            "MP4 Energy",
            "MP4SDQ Energy",
            "MP4DQ Energy",
            "MP3 Energy",
            "MP2 Energy",
        ):
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
            ("metadata/symmetry_detected", "Symmetry"),
            ("metadata/symmetry_used", "Symmetry used"),
        ]
        if "E(β-homo)" in data:
            for letter in ("α", "β"):
                keys.extend(
                    [
                        (f"E({letter}-homo)", f"{letter}-HOMO Energy"),
                        (f"E({letter}-lumo)", f"{letter}-LUMO Energy"),
                        (f"E({letter}-gap)", f"{letter}-Gap"),
                        (f"Sym({letter}-homo)", f"{letter}-HOMO Symmetry"),
                        (f"Sym({letter}-lumo)", f"{letter}-LUMO Symmetry"),
                    ]
                )
        else:
            keys.extend(
                [
                    ("E(homo)", "HOMO Energy"),
                    ("E(lumo)", "LUMO Energy"),
                    ("E(gap)", "Gap"),
                    ("Sym(homo)", "HOMO Symmetry"),
                    ("Sym(lumo)", "LUMO Symmetry"),
                ]
            )
        keys.extend(
            [
                ("dipole_moment_magnitude", "Dipole moment"),
            ]
        )
        for key, name in keys:
            if key in data:
                tmp = data[key]
                mdata = metadata[key]
                table["Property"].append(name)
                table["Value"].append(f"{tmp:{mdata['format']}}")
                if "units" in mdata:
                    table["Units"].append(mdata["units"])
                else:
                    table["Units"].append("")

        for key, name in (
            ("metadata/cpu_time", "CPU time"),
            ("metadata/wall_time", "Wall-clock time"),
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
        text_lines.append("Results".center(length))
        text_lines.append(tmp)

        if text != "":
            text = str(__(text, **data, indent=self.indent + 4 * " "))
            text += "\n\n"
        text += textwrap.indent("\n".join(text_lines), self.indent + 7 * " ")

        if "Composite/summary" in data:
            text += "\n\n\n"
            text += textwrap.indent(data["Composite/summary"], self.indent + 4 * " ")

        # Update the structure. Gaussian may have reoriented.
        system, configuration = self.get_system_configuration(None)
        directory = Path(self.directory)
        if "Current cartesian coordinates" in data:
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
            with open(directory / "atom_properties.csv", "w", newline="") as fd:
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
        obMol = configuration.to_OBMol(properties="all")
        title = f"SEAMM={system.name}/{configuration.name}"
        obMol.SetTitle(title)
        sdf = obConversion.WriteString(obMol)
        path = directory / "structure.sdf"
        path.write_text(sdf)

        text = self.make_plots(data)
        if text != "":
            printer.normal(__(text, indent=self.indent + 4 * " "))

        text = (
            "The structure and charges, etc. were placed in configuration "
            f"'{system.name}/{configuration.name}'."
        )
        printer.normal("")
        printer.normal(__(text, indent=self.indent + 4 * " "))
        printer.normal("")

        # Put any requested results into variables or tables
        self.store_results(data=data, create_tables=True)

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
