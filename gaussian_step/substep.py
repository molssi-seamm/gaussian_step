# -*- coding: utf-8 -*-

"""Setup and run Gaussian"""

from collections import Counter
import configparser
import csv
from datetime import datetime, timezone
import gzip
import importlib
import logging
from math import isnan
import os
from pathlib import Path
import pkg_resources
import platform
import pprint
import re
import shutil
import string
import time

import cclib
from cpuinfo import get_cpu_info
import numpy as np
import pandas

import gaussian_step
import seamm
import seamm_exec
import seamm.data
from seamm_util import Configuration, Q_
import seamm_util.printing as printing

logger = logging.getLogger("Gaussian")
job = printing.getPrinter()
printer = printing.getPrinter("gaussian")


def humanize(memory, suffix="B", kilo=1024):
    """
    Scale memory to its proper format e.g:

        1253656 => '1.20 MiB'
        1253656678 => '1.17 GiB'
    """
    if kilo == 1000:
        units = ["", "k", "M", "G", "T", "P"]
    elif kilo == 1024:
        units = ["", "Ki", "Mi", "Gi", "Ti", "Pi"]
    else:
        raise ValueError("kilo must be 1000 or 1024!")

    for unit in units:
        if memory < 10 * kilo:
            return f"{int(memory)}{unit}{suffix}"
        memory /= kilo


def dehumanize(memory, suffix="B"):
    """
    Unscale memory from its human readable form e.g:

        '1.20 MB' => 1200000
        '1.17 GB' => 1170000000
    """
    units = {
        "": 1,
        "k": 1000,
        "M": 1000**2,
        "G": 1000**3,
        "P": 1000**4,
        "Ki": 1024,
        "Mi": 1024**2,
        "Gi": 1024**3,
        "Pi": 1024**4,
    }

    tmp = memory.split()
    if len(tmp) == 1:
        return memory
    elif len(tmp) > 2:
        raise ValueError("Memory must be <number> <units>, e.g. 1.23 GB")

    amount, unit = tmp
    amount = float(amount)

    for prefix in units:
        if prefix + suffix == unit:
            return int(amount * units[prefix])

    raise ValueError(f"Don't recognize the units on '{memory}'")


class Substep(seamm.Node):
    def __init__(
        self,
        flowchart=None,
        title="no title",
        extension=None,
        logger=logger,
        module=__name__,
    ):
        """Initialize the node"""

        logger.debug("Creating Energy {}".format(self))

        super().__init__(
            flowchart=flowchart, title=title, extension=extension, logger=logger
        )

        self._input_only = False
        self._timing_data = []
        self._timing_path = Path("~/.seamm.d/timing/gaussian.csv").expanduser()

        # Set up the timing information
        self._timing_header = [
            "node",  # 0
            "cpu",  # 1
            "cpu_version",  # 2
            "cpu_count",  # 3
            "cpu_speed",  # 4
            "date",  # 5
            "SMILES",  # 6
            "H_SMILES",  # 7
            "formula",  # 8
            "net_charge",  # 9
            "spin_multiplicity",  # 10
            "model",  # 11
            "keywords",  # 12
            "symmetry",  # 13
            "symmetry_used",  # 14
            "nbf",  # 15
            "nproc",  # 16
            "time",  # 17
        ]
        try:
            self._timing_path.parent.mkdir(parents=True, exist_ok=True)

            self._timing_data = 18 * [""]
            self._timing_data[0] = platform.node()
            tmp = get_cpu_info()
            if "arch" in tmp:
                self._timing_data[1] = tmp["arch"]
            if "cpuinfo_version_string" in tmp:
                self._timing_data[2] = tmp["cpuinfo_version_string"]
            if "count" in tmp:
                self._timing_data[3] = str(tmp["count"])
            if "hz_advertized_friendly" in tmp:
                self._timing_data[3] = tmp["hz_advertized_friendly"]

            if not self._timing_path.exists():
                with self._timing_path.open("w", newline="") as fd:
                    writer = csv.writer(fd)
                    writer.writerow(self._timing_header)
        except Exception:
            self._timing_data = None

    @property
    def version(self):
        """The semantic version of this module."""
        return gaussian_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return gaussian_step.__git_revision__

    @property
    def global_options(self):
        return self.parent.global_options

    @property
    def gversion(self):
        return self.parent.gversion

    @property
    def input_only(self):
        """Whether to write the input only, not run MOPAC."""
        return self._input_only

    @input_only.setter
    def input_only(self, value):
        self._input_only = value

    @property
    def is_runable(self):
        """Indicate whether this not runs or just adds input."""
        return True

    @property
    def options(self):
        return self.parent.options

    def calculate_enthalpy_of_formation(self, data):
        """Calculate the enthalpy of formation from the results of a calculation.

        This uses tabulated values of the enthalpy of formation of the atoms for
        the elements and tabulated energies calculated for atoms with the current
        method.

        Parameters
        ----------
        data : dict
            The results of the calculation.
        """
        # Read the tabulated values from either user or data directory
        personal_file = Path("~/.seamm.d/data/atom_energies.csv").expanduser()
        if personal_file.exists():
            personal_table = pandas.read_csv(personal_file, index_col=False)
        else:
            personal_table = None

        path = Path(pkg_resources.resource_filename(__name__, "data/"))
        csv_file = path / "atom_energies.csv"
        table = pandas.read_csv(csv_file, index_col=False)

        self.logger.debug(f"self.model = {self.model}")
        self.logger.debug("\n\t".join(table.columns))

        # Check if have the data
        atom_energies = None
        correction_energy = None
        if self.model.startswith("U") or self.model.startswith("R"):
            column = self.model[1:]
        else:
            column = self.model

        if self.logger.isEnabledFor(logging.DEBUG):
            self.logger.debug(f"Looking for '{column}' in")
            for i, col in enumerate(table.columns, start=1):
                self.logger.debug(f"\t{i:5d}: {col}")

        column2 = column + " correction"
        if personal_table is not None and column in personal_table.columns:
            atom_energies = personal_table[column].to_list()
            if column2 in personal_table.columns:
                correction_energy = personal_table[column2].to_list()
        elif column in table.columns:
            atom_energies = table[column].to_list()
            if column2 in table.columns:
                correction_energy = table[column2].to_list()

        DfH0gas = None
        if personal_table is not None and "ΔfH°gas" in personal_table.columns:
            DfH0gas = personal_table["ΔfH°gas"].to_list()
        elif "ΔfH°gas" in table.columns:
            DfH0gas = table["ΔfH°gas"].to_list()

        if DfH0gas is None or atom_energies is None:
            return

        # Get the atomic numbers and counts
        _, configuration = self.get_system_configuration(None)
        counts = Counter(configuration.atoms.atomic_numbers)
        Ecorr = 0.0
        for atno, count in counts.items():
            if isnan(atom_energies[atno - 1]):
                # Don't have the data for this element
                return
            Ecorr += count * (DfH0gas[atno - 1] - atom_energies[atno - 1])
            if correction_energy is not None and not isnan(correction_energy[atno - 1]):
                Ecorr += count * correction_energy[atno - 1]

        data["DfH0"] = Q_(data["H"], "E_h").m_as("kJ/mol") + Ecorr

        return

    def get_functional(self, P=None):
        """Work out the DFT functional"""
        if P is None:
            P = self.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )
        if P["level"] == "recommended":
            functional = P["functional"]
        else:
            functional = P["advanced_functional"]
        found = functional in gaussian_step.dft_functionals
        if not found:
            # Might be first part of name, or Gaussian-encoded name
            for _key, _data in gaussian_step.dft_functionals.items():
                if functional == _key.split(":")[0].strip():
                    functional = _key
                    found = True
                    break
        if not found:
            # Might be the internal Gaussian name
            for _key, _data in gaussian_step.dft_functionals.items():
                if functional == _data["name"]:
                    functional = _key
                    found = True
                    break
        if not found:
            raise ValueError(f"Don't recognize functional '{functional}'")

        return functional

    def get_method(self, P=None):
        """The method ... HF, DFT, ... used."""
        # Figure out the method.
        if P is None:
            P = self.parameters.current_values_to_dict(
                context=seamm.flowchart_variables._data
            )
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
            for _key, _mdata in gaussian_step.methods.items():
                if method_string == _mdata["method"]:
                    method_string = _key
                    method_data = _mdata
                    method = method_data["method"]
                    break
            else:
                method_data = {}
                method = method_string

        return method, method_data

    def make_plots(self, data):
        """Create the density and orbital plots if requested.

        Parameters
        ----------
        data : dict()
             Dictionary of results from the calculation (results.tag file)
        """
        text = "\n\n"

        directory = Path(self.directory)
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Get the configuration and basic information
        system, configuration = self.get_system_configuration(None)

        periodicity = configuration.periodicity
        if periodicity != 0:
            raise NotImplementedError("Periodic cube files not implemented yet!")

        # Have the needed data?
        if "homos" not in data:
            return ""

        spin_polarized = len(data["homos"]) == 2

        # Prepare to run
        executor = self.parent.flowchart.executor

        # Read configuration file for Gaussian
        seamm_options = self.global_options
        ini_dir = Path(seamm_options["root"]).expanduser()
        full_config = configparser.ConfigParser()
        full_config.read(ini_dir / "gaussian.ini")
        executor_type = executor.name
        if executor_type not in full_config:
            raise RuntimeError(
                f"No section for '{executor_type}' in MOPAC ini file "
                f"({ini_dir / 'mopac.ini'})"
            )
        config = dict(full_config.items(executor_type))

        # Set up the environment
        if config["root-directory"] != "":
            env = {"g09root": config["root-directory"]}
        else:
            env = {}

        if config["setup-environment"] != "":
            cmd = [". {setup-environment} && cubegen"]
        else:
            cmd = ["cubegen"]

        npts = "-2"

        keys = []
        if P["total density"]:
            keys.append("total density")
        if spin_polarized and P["total spin density"]:
            keys.append("spin density")

        n_errors = 0
        for key in keys:
            if key == "total density":
                args = f"1 Density=SCF gaussian.fchk Total_Density.cube {npts} h"
            elif key == "spin density":
                args = f"1 Spin=SCF gaussian.fchk Spin_Density.cube {npts} h"

            # And run CUBEGEN
            result = executor.run(
                cmd=[*cmd, args],
                config=config,
                directory=self.directory,
                files={},
                return_files=["*"],
                in_situ=True,
                shell=True,
                env=env,
            )
            if not result:
                self.logger.error("There was an error running CubeGen")
                n_errors += 1
                printer.important(f"There was an error calling CUBEGEN, {cmd} {args}")

        # Any requested orbitals
        if P["orbitals"]:
            n_orbitals = data["nmo"]
            # and work out the orbitals
            txt = P["selected orbitals"]
            for spin, homo in enumerate(data["homos"]):
                if txt == "all":
                    orbitals = [*range(n_orbitals)]
                else:
                    orbitals = []
                    for chunk in txt.split(","):
                        chunk = chunk.strip()
                        if ":" in chunk or ".." in chunk:
                            if ":" in chunk:
                                first, last = chunk.split(":")
                            elif ".." in chunk:
                                first, last = chunk.split("..")
                            first = first.strip().upper()
                            last = last.strip().upper()

                            if first == "HOMO":
                                first = homo
                            elif first == "LUMO":
                                first = homo + 1
                            else:
                                first = int(
                                    first.removeprefix("HOMO").removeprefix("LUMO")
                                )
                                if first < 0:
                                    first = homo + first
                                else:
                                    first = homo + 1 + first

                            if last == "HOMO":
                                last = homo
                            elif last == "LUMO":
                                last = homo + 1
                            else:
                                last = int(
                                    last.removeprefix("HOMO").removeprefix("LUMO")
                                )
                                if last < 0:
                                    last = homo + last
                                else:
                                    last = homo + 1 + last

                            orbitals.extend(range(first, last + 1))
                        else:
                            first = chunk.strip().upper()

                            if first == "HOMO":
                                first = homo
                            elif first == "LUMO":
                                first = homo + 1
                            else:
                                first = int(
                                    first.removeprefix("HOMO").removeprefix("LUMO")
                                )
                                if first < 0:
                                    first = homo + first
                                else:
                                    first = homo + 1 + first
                            orbitals.append(first)

                # Remove orbitals out of limits
                tmp = orbitals
                orbitals = []
                for x in tmp:
                    if x >= 0 and x < n_orbitals:
                        orbitals.append(x)

                if spin_polarized:
                    l1 = ("A", "B")[spin]
                    l2 = ("α-", "β-")[spin]
                else:
                    l1 = ""
                    l2 = ""
                for mo in orbitals:
                    if mo == homo:
                        filename = f"{l2}HOMO.cube"
                    elif mo < homo:
                        filename = f"{l2}HOMO-{homo - mo}.cube"
                    elif mo == homo + 1:
                        filename = f"{l2}LUMO.cube"
                    else:
                        filename = f"{l2}LUMO+{mo - homo - 1}.cube"
                    args = f"1 {l1}MO={mo + 1} gaussian.fchk {filename} {npts} h"

                    # And run CUBEGEN
                    result = executor.run(
                        cmd=[*cmd, args],
                        config=config,
                        directory=self.directory,
                        files={},
                        return_files=["*"],
                        in_situ=True,
                        shell=True,
                        env=env,
                    )
                    if not result:
                        self.logger.error("There was an error running CubeGen")
                        n_errors += 1
                        printer.important(
                            f"There was an error calling CUBEGEN, {cmd} {args}"
                        )

        # Finally rename and gzip the cube files
        n_processed = 0
        paths = directory.glob("*.cube")
        for path in paths:
            out = path.with_suffix(".cube.gz")
            with path.open("rb") as f_in:
                with gzip.open(out, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)
            n_processed += 1
            path.unlink()
        if n_errors > 0:
            text += (
                f"Created {n_processed} density and orbital cube files, but there were "
                f"{n_errors} errors trying to create cube files."
            )
        else:
            text += f"Created {n_processed} density and orbital cube files."

        return text

    def parse_fchk(self, path, data={}):
        """Process the data of a formatted Chk file given as lines of data.

        Parameters
        ----------
        path : pathlib.Path
            The path to the checkpoint file
        """
        lines = path.read_text().splitlines()

        it = iter(lines)
        # Ignore first potentially truncated title line
        next(it)

        # Type line (A10,A30,A30)
        line = next(it)
        data["calculation"] = line[0:10].strip()
        data["method"] = line[10:40].strip()

        # The rest of the file consists of a line defining the data.
        # If the data is a scalar, it is on the control line, otherwise it follows
        translation = self.metadata["translation"]
        while True:
            try:
                line = next(it)
            except StopIteration:
                break
            try:
                key = line[0:40].strip()
                if key in translation:
                    key = translation[key]
                code = line[43]
                is_array = line[47:49] == "N="
                if is_array:
                    count = int(line[49:61].strip())
                    value = []
                    if code == "I":
                        i = 0
                        while i < count:
                            line = next(it)
                            for pos in range(0, 6 * 12, 12):
                                value.append(int(line[pos : pos + 12].strip()))
                                i += 1
                                if i == count:
                                    break
                    elif code == "R":
                        i = 0
                        while i < count:
                            line = next(it)
                            for pos in range(0, 5 * 16, 16):
                                text = line[pos : pos + 16].strip()
                                # Fortran drops E in format for large exponents...
                                text = re.sub(r"([0-9])-", r"\1E-", text)
                                value.append(float(text))
                                i += 1
                                if i == count:
                                    break
                    elif code == "C":
                        value = ""
                        i = 0
                        while i < count:
                            line = next(it)
                            for pos in range(0, 5 * 12, 12):
                                value += line[pos : pos + 12]
                                i += 1
                                if i == count:
                                    break
                                value = value.rstrip()
                    elif code == "H":
                        value = ""
                        i = 0
                        while i < count:
                            line = next(it)
                            for pos in range(0, 9 * 8, 8):
                                value += line[pos : pos + 8]
                                i += 1
                                if i == count:
                                    break
                                value = value.rstrip()
                    elif code == "L":
                        i = 0
                        while i < count:
                            line = next(it)
                            for pos in range(72):
                                value.append(line[pos] == "T")
                                i += 1
                                if i == count:
                                    break
                else:
                    if code == "I":
                        value = int(line[49:].strip())
                    elif code == "R":
                        value = float(line[49:].strip())
                    elif code == "C":
                        value = line[49:].strip()
                    elif code == "L":
                        value = line[49] == "T"
                data[key] = value
            except Exception:
                pass
        return data

    def parse_output(self, path, data={}):
        """Process the output.

        Parameters
        ----------
        path : pathlib.Path
            The Gaussian log file.
        """
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        lines = path.read_text().splitlines()

        # Did it end properly?
        data["success"] = "Normal termination" in lines[-1]

        # Find the date and version of Gaussian
        # Gaussian 09:  EM64M-G09RevE.01 30-Nov-2015
        it = iter(lines)
        for line in it:
            if "Cite this work" in line:
                for line in it:
                    if "**********************" in line:
                        line = next(it)
                        if "Gaussian" in line:
                            try:
                                _, version, revision, date = line.split()
                                _, month, year = date.split("-")
                                revision = revision.split("Rev")[1]
                                data["G revision"] = revision
                                data["G version"] = f"G{version.strip(':')}"
                                data["G month"] = month
                                data["G year"] = year
                            except Exception as e:
                                self.logger.warning(
                                    f"Could not find the Gaussian citation: {e}"
                                )
                            break
                break

        # And the optimization steps, if any.
        #
        # Need to be careful about end of the first (and presumably only?) optimization.
        # The FORCE calculation prints out the same information about convergence, but
        # may indicate no convergence. This can confuse this code unless we look for the
        # end of the optimization step and quit then
        it = iter(lines)
        n_steps = 0
        max_force = []
        rms_force = []
        max_displacement = []
        rms_displacement = []
        converged = None
        for line in it:
            if line == "         Item               Value     Threshold  Converged?":
                n_steps += 1
                converged = True

                tmp1, tmp2, value, threshold, criterion = next(it).split()
                if tmp1 == "Maximum" and tmp2 == "Force":
                    max_force.append(float(value))
                    data["maximum atom force threshold"] = float(threshold)
                    if criterion != "YES":
                        converged = False

                tmp1, tmp2, value, threshold, criterion = next(it).split()
                if tmp1 == "RMS" and tmp2 == "Force":
                    rms_force.append(float(value))
                    data["RMS atom force threshold"] = float(threshold)
                    if criterion != "YES":
                        converged = False

                tmp1, tmp2, value, threshold, criterion = next(it).split()
                if tmp1 == "Maximum" and tmp2 == "Displacement":
                    max_displacement.append(float(value))
                    data["maximum atom displacement threshold"] = float(threshold)
                    if criterion != "YES":
                        converged = False

                tmp1, tmp2, value, threshold, criterion = next(it).split()
                if tmp1 == "RMS" and tmp2 == "Displacement":
                    rms_displacement.append(float(value))
                    data["RMS atom displacement threshold"] = float(threshold)
                    if criterion != "YES":
                        converged = False
            elif line == " Optimization completed.":
                line = next(it)
                if line == "    -- Stationary point found.":
                    converged = True
                else:
                    self.logger.warning(f"Optimization completed: {line}")
                break
            elif line == "    -- Stationary point found.":
                converged = True
                break

        data["N steps optimization"] = n_steps

        if converged is not None:
            data["optimization is converged"] = converged
            data["maximum atom force"] = max_force[-1]
            data["RMS atom force"] = rms_force[-1]
            data["maximum atom displacement"] = max_displacement[-1]
            data["RMS atom displacement"] = rms_displacement[-1]
        data["maximum atom force Trajectory"] = max_force
        data["RMS atom force Trajectory"] = rms_force
        data["maximum atom displacement Trajectory"] = max_displacement
        data["RMS atom displacement Trajectory"] = rms_displacement

        # CBS calculations

        # Complete Basis Set (CBS) Extrapolation:
        # M. R. Nyden and G. A. Petersson, JCP 75, 1843 (1981)
        # G. A. Petersson and M. A. Al-Laham, JCP 94, 6081 (1991)
        # G. A. Petersson, T. Tensfeldt, and J. A. Montgomery, JCP 94, 6091 (1991)
        # J. A. Montgomery, J. W. Ochterski, and G. A. Petersson, JCP 101, 5900 (1994)
        #
        # Temperature=               298.150000 Pressure=                       1.000000
        # E(ZPE)=                      0.050496 E(Thermal)=                     0.053508
        # E(SCF)=                    -78.059017 DE(MP2)=                       -0.281841
        # DE(CBS)=                    -0.071189 DE(MP34)=                      -0.024136
        # DE(Int)=                     0.021229 DE(Empirical)=                 -0.075463
        # CBS-4 (0 K)=               -78.439921 CBS-4 Energy=                 -78.436908
        # CBS-4 Enthalpy=            -78.435964 CBS-4 Free Energy=            -78.460753

        method, method_data = self.get_method(P)
        self.logger.debug(f"Checking for CBS extrapolation for {method}")
        if method[0:4] == "CBS-":
            # Need last section
            if method == "CBS-4M":
                match = "CBS-4 Enthalpy="
            else:
                match = f"{method} Enthalpy="
            self.logger.debug(f"Looking for '{match}'")
            text = []
            found = False
            for line in reversed(lines):
                if found:
                    text.append(line)
                    if "Complete Basis Set" in line:
                        break
                elif match in line:
                    found = True
                    text.append(line)

            self.logger.debug(f"Found CBS extrapolation: {found}")
            if found:
                text = text[::-1]
                it = iter(text)
                next(it)
                citations = []
                for line in it:
                    tmp = line.strip()
                    if tmp == "":
                        break
                    citations.append(tmp)
                data["citations"] = citations

                for line in it:
                    line = line.strip()
                    if len(line) > 40:
                        part = [line[0:37], line[38:]]
                    else:
                        part = [line]
                    for p in part:
                        if "=" not in p:
                            continue
                        key, value = p.split("=", 1)
                        key = key.strip()
                        value = float(value.strip())
                        if "(0 K)" in key:
                            key = "E 0"
                        elif "Free Energy" in key:
                            key = "F"
                        elif "Energy" in key:
                            key = "energy"
                        elif "Enthalpy" in key:
                            key = "H"
                        data[key] = value
                data["model"] = method
                data["Composite/summary"] = "\n".join(text)

        # Gn calculations. No header!!!!!

        # Temperature=              298.150000 Pressure=                      1.000000
        # E(ZPE)=                     0.050251 E(Thermal)=                    0.053306
        # E(CCSD(T))=               -78.321715 E(Empiric)=                   -0.041682
        # DE(Plus)=                  -0.005930 DE(2DF)=                      -0.076980
        # E(Delta-G3XP)=             -0.117567 DE(HF)=                       -0.008255
        # G4(0 K)=                  -78.521880 G4 Energy=                   -78.518825
        # G4 Enthalpy=              -78.517880 G4 Free Energy=              -78.542752

        if method in (
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
            # Need last section
            if method == "G2":
                match = "G2MP2 Enthalpy="
            elif method == "G3B3":
                match = "G3 Enthalpy="
                method = "G3"
            elif method == "G3MP2B3":
                match = "G3MP2 Enthalpy="
                method = "G3MP2"
            else:
                match = f"{method} Enthalpy="
            text = []
            found = False
            for line in reversed(lines):
                if found:
                    if line.strip() == "":
                        break
                    text.append(line)
                elif match in line:
                    found = True
                    text.append(line)

            if found:
                translation = self.metadata["translation"]
                text = text[::-1]
                for line in text:
                    line = line.strip()
                    if len(line) > 36:
                        part = [line[0:36], line[37:]]
                    else:
                        part = [line]
                    for p in part:
                        if "=" not in p:
                            continue
                        key, value = p.split("=", 1)
                        key = key.strip()
                        value = float(value.strip())
                        if key.startswith(method):
                            key = key.split(" ", 1)[1]
                        elif key == "E(Empiric)":
                            key = "E(empirical)"
                        key = "Composite/" + key
                        if key in translation:
                            key = translation[key]
                        data[key] = value

                data["model"] = method
                tmp = " " * 20 + f"{method[0:2]} composite method extrapolation\n\n"
                data["Composite/summary"] = tmp + "\n".join(text)
                data["Total Energy"] = data["Composite/Free Energy"]

        # The Wiberg bond orders ... which look like this:

        # Wiberg bond index matrix in the NAO basis:
        #
        #     Atom    1       2       3       4       5       6       7       8       9
        #     ---- ------  ------  ------  ------  ------  ------  ------  ------  -----
        #   1.  C  0.0000  1.8962  0.0134  0.1327  0.9261  0.9249  0.0071  0.0005  0.00x
        #   2.  C  1.8962  0.0000  1.1131  0.0127  0.0044  0.0049  0.9112  0.0029  0.01x
        #  ...
        #  10.  H  0.0002  0.0022  0.0049  0.9269  0.0000  0.0002  0.0015  0.0171  0.00x
        #
        #     Atom   10
        #     ---- ------
        #   1.  C  0.0002
        #   2.  C  0.0022
        #  ...

        it = iter(lines)
        for line in it:
            n_atoms = None
            if line.startswith(" Wiberg bond index matrix in the NAO basis:"):
                bond_orders = []
                next(it)
                # Read each chunk of output
                while True:
                    # Skip the two header lines
                    next(it)
                    next(it)
                    count = 0
                    # And add the data to the bond_order matrix
                    for line in it:
                        line = line.strip()
                        if line == "":
                            if n_atoms is None:
                                n_atoms = count
                            break
                        count += 1
                        vals = [float(val) for val in line.split()[2:]]
                        if len(bond_orders) < count:
                            bond_orders.append(vals)
                        else:
                            bond_orders[count - 1].extend(vals)
                    if len(bond_orders[0]) >= n_atoms:
                        break

                data["Wiberg bond order matrix"] = bond_orders

        return data

    def process_data(self, data):
        """Massage the cclib data to a more easily used form."""
        self.logger.debug(pprint.pformat(data))
        # Convert numpy arrays to Python lists
        new = {}
        for key, value in data.items():
            if isinstance(value, np.ndarray):
                new[key] = value.tolist()
            elif isinstance(value, list):
                if len(value) > 0 and isinstance(value[0], np.ndarray):
                    new[key] = [i.tolist() for i in value]
                else:
                    new[key] = value
            elif isinstance(value, dict):
                for k, v in value.items():
                    newkey = f"{key}/{k}"
                    if isinstance(v, np.ndarray):
                        new[newkey] = v.tolist()
                    else:
                        new[newkey] = v
            else:
                new[key] = value

        for key in ("metadata/cpu_time", "metadata/wall_time"):
            if key in new:
                time = new[key][0]
                for tmp in new[key][1:]:
                    time += tmp
                new[key] = str(time).lstrip("0:")
                if "." in new[key]:
                    new[key] = new[key].rstrip("0")

        # Pull out the HOMO and LUMO energies as scalars
        if "homos" in new and "moenergies" in new:
            homos = new["homos"]
            if len(homos) == 2:
                for i, letter in enumerate(["α", "β"]):
                    Es = new["moenergies"][i]
                    homo = homos[i]
                    new[f"{letter} HOMO orbital number"] = homo + 1
                    new[f"E {letter} homo"] = Es[homo]
                    if homo > 0:
                        new[f"E {letter} nhomo"] = Es[homo - 1]
                    if homo + 1 < len(Es):
                        new[f"E {letter} lumo"] = Es[homo + 1]
                        new[f"E {letter} gap"] = Es[homo + 1] - Es[homo]
                    if homo + 2 < len(Es):
                        new[f"E {letter} slumo"] = Es[homo + 2]
                    if "mosyms" in new:
                        syms = new["mosyms"][i]
                        new[f"{letter} HOMO symmetry"] = syms[homo]
                        if homo > 0:
                            new[f"{letter} NHOMO symmetry"] = syms[homo - 1]
                        if homo + 1 < len(syms):
                            new[f"{letter} LUMO symmetry)"] = syms[homo + 1]
                        if homo + 2 < len(syms):
                            new[f"{letter} SLUMO symmetry"] = syms[homo + 2]
            else:
                Es = new["moenergies"][0]
                homo = homos[0]
                new["HOMO orbital number"] = homo + 1
                new["E homo"] = Es[homo]
                if homo > 0:
                    new["E nhomo"] = Es[homo - 1]
                if homo + 1 < len(Es):
                    new["E lumo"] = Es[homo + 1]
                    new["E gap"] = Es[homo + 1] - Es[homo]
                if homo + 2 < len(Es):
                    new["E slumo"] = Es[homo + 2]
                if "mosyms" in new:
                    syms = new["mosyms"][0]
                    new["HOMO symmetry"] = syms[homo]
                    if homo > 0:
                        new["NHOMO symmetry"] = syms[homo - 1]
                    if homo + 1 < len(syms):
                        new["LUMO symmetry"] = syms[homo + 1]
                    if homo + 2 < len(syms):
                        new["SLUMO symmetry"] = syms[homo + 2]

        # moments
        if "moments" in new:
            moments = new["moments"]
            new["multipole reference point"] = moments[0]
            new["dipole moment"] = moments[1]
            new["dipole moment magnitude"] = np.linalg.norm(moments[1])
            if len(moments) > 2:
                new["quadrupole moment"] = moments[2]
            if len(moments) > 3:
                new["octapole moment"] = moments[3]
            if len(moments) > 4:
                new["hexadecapole moment"] = moments[4]
            del new["moments"]

        for key in ("symmetry group", "symmetry group used"):
            if key in new:
                new[key] = new[key].capitalize()

        return new

    def run_gaussian(self, keywords, extra_lines=None):
        """Run Gaussian.

        Parameters
        ----------
        None

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        # Check for successful run, don't rerun
        success_file = directory / "success.dat"
        if not success_file.exists():
            # Get the system & configuration
            system, configuration = self.get_system_configuration(None)

            # Access the options
            options = self.options
            seamm_options = self.global_options

            # Get the computational environment and set limits
            ce = seamm_exec.computational_environment()

            # How many threads to use
            n_cores = ce["NTASKS"]
            self.logger.debug("The number of cores available is {}".format(n_cores))

            if seamm_options["parallelism"] not in ("openmp", "any"):
                n_threads = 1
            else:
                if options["ncores"] == "available":
                    n_threads = n_cores
                else:
                    n_threads = int(options["ncores"])
                if n_threads > n_cores:
                    n_threads = n_cores
                if n_threads < 1:
                    n_threads = 1
                if seamm_options["ncores"] != "available":
                    n_threads = min(n_threads, int(seamm_options["ncores"]))
            ce["NTASKS"] = n_threads
            self.logger.debug(f"Gaussian will use {n_threads} threads.")

            # How much memory to use
            if seamm_options["memory"] == "all":
                mem_limit = ce["MEM_PER_NODE"]
            elif seamm_options["memory"] == "available":
                # For the default, 'available', use in proportion to number of
                # cores used
                mem_limit = ce["MEM_PER_CPU"] * n_threads
            else:
                mem_limit = dehumanize(seamm_options["memory"])

            if options["memory"] == "all":
                memory = ce["MEM_PER_NODE"]
            elif options["memory"] == "available":
                # For the default, 'available', use in proportion to number of
                # cores used
                memory = ce["MEM_PER_CPU"] * n_threads
            else:
                memory = dehumanize(options["memory"])

            memory = min(memory, mem_limit)
            ce["MEM_PER_NODE"] = memory

            # Apply a minimum of 800 MB
            min_memory = dehumanize("800 MB")
            if min_memory > memory:
                memory = min_memory

            # Gaussian allows no decimal points.
            memory = humanize(memory, kilo=1000)

            lines = []
            lines.append("%Chk=gaussian")
            lines.append(f"%Mem={memory}")
            lines.append(f"%NProcShared={n_threads}")

            lines.append("# " + " ".join(keywords))

            lines.append(" ")
            lines.append(f"{system.name}/{configuration.name}")
            lines.append(" ")
            lines.append(f"{configuration.charge}    {configuration.spin_multiplicity}")

            # Atoms with coordinates
            symbols = configuration.atoms.symbols
            XYZs = configuration.atoms.coordinates
            for symbol, xyz in zip(symbols, XYZs):
                x, y, z = xyz
                lines.append(f"{symbol:2}   {x:10.6f} {y:10.6f} {z:10.6f}")
            lines.append(" ")

            if extra_lines is not None:
                lines.extend(extra_lines)

            files = {"input.dat": "\n".join(lines)}
            self.logger.debug("input.dat:\n" + files["input.dat"])

            printer.important(
                self.indent + f"    Gaussian will use {n_threads} OpenMP threads and "
                f"up to {memory} of memory.\n"
            )
            if self.input_only:
                # Just write the input files and stop
                for filename in files:
                    path = directory / filename
                    path.write_text(files[filename])
                data = {}
            else:
                executor = self.parent.flowchart.executor

                # Read configuration file for Gaussian if it exists
                executor_type = executor.name
                full_config = configparser.ConfigParser()
                ini_dir = Path(seamm_options["root"]).expanduser()
                path = ini_dir / "gaussian.ini"

                # If the config file doesn't exist, get the default
                if not path.exists():
                    resources = importlib.resources.files("gaussian_step") / "data"
                    ini_text = (resources / "gaussian.ini").read_text()
                    txt_config = Configuration(path)
                    txt_config.from_string(ini_text)
                    txt_config.save()

                full_config.read(ini_dir / "gaussian.ini")

                # Getting desperate! Look for an executable in the path
                if (
                    executor_type not in full_config
                    or "root-directory" not in full_config[executor_type]
                    or "setup-environment" not in full_config[executor_type]
                ):
                    # See if we can find the Gaussian environment variables
                    if "g16root" in os.environ:
                        g_ver = "g16"
                        root_directory = os.environ["g16root"]
                        if "GAUSS_BSDDIR" in os.environ:
                            setup_directory = Path(os.environ["GAUSS_BSDDIR"])
                        else:
                            setup_directory = Path(root_directory) / g_ver / "bsd"
                    elif "g09root" in os.environ:
                        g_ver = "g09"
                        root_directory = os.environ["g09root"]
                        if "GAUSS_BSDDIR" in os.environ:
                            setup_directory = Path(os.environ["GAUSS_BSDDIR"])
                        else:
                            setup_directory = Path(root_directory) / g_ver / "bsd"
                    else:
                        root_directory = None
                        exe_path = shutil.which("g16")
                        if exe_path is None:
                            exe_path = shutil.which("g09")
                        if exe_path is None:
                            raise RuntimeError(
                                f"No section for '{executor_type}' in Gaussian ini file"
                                f" ({ini_dir / 'gaussian.ini'}), nor in the defaults, "
                                "nor in the path!"
                            )
                        g_ver = exe_path.name
                        root_directory = str(exe_path.parent.parent)
                        setup_directory = Path(root_directory) / g_ver / "bsd"
                    setup_environment = str(setup_directory / f"{g_ver}.profile")

                    txt_config = Configuration(path)

                    if not txt_config.section_exists(executor_type):
                        txt_config.add_section(executor_type)

                    txt_config.set_value(executor_type, "installation", "local")
                    txt_config.set_value(executor_type, "code", g_ver)
                    txt_config.set_value(
                        executor_type, "root-directory", root_directory
                    )
                    txt_config.set_value(
                        executor_type, "setup-environment", setup_environment
                    )
                    txt_config.save()
                    full_config.read(ini_dir / "gaussian.ini")

                config = dict(full_config.items(executor_type))
                # Use the matching version of the seamm-gaussian image by default.
                config["version"] = self.version

                g_ver = config["code"]

                # Setup the calculation environment definition
                if config["root-directory"] != "":
                    env = {f"{g_ver}root": config["root-directory"]}
                else:
                    env = {}

                if config["setup-environment"] != "":
                    cmd = f". {config['setup-environment']} ; {g_ver}"
                else:
                    cmd = g_ver

                cmd += " < input.dat > output.txt ; formchk gaussian.chk"

                return_files = [
                    "output.txt",
                    "gaussian.chk",
                    "gaussian.fchk",
                ]

                self.logger.debug(f"{cmd=}")
                self.logger.debug(f"{env=}")

                self._timing_data[12] = " ".join(keywords)
                self._timing_data[5] = datetime.now(timezone.utc).isoformat()
                t0 = time.time_ns()

                result = executor.run(
                    cmd=[cmd],
                    config=config,
                    directory=self.directory,
                    files=files,
                    return_files=return_files,
                    in_situ=True,
                    shell=True,
                    env=env,
                )

                t = (time.time_ns() - t0) / 1.0e9
                if self._timing_data is not None:
                    self._timing_data[17] = f"{t:.3f}"
                    self._timing_data[16] = str(n_threads)

                if not result:
                    self.logger.error("There was an error running Gaussian")
                    return None

                # self.logger.debug("\n" + pprint.pformat(result))

        if not self.input_only:
            # And output
            path = directory / "output.txt"
            if path.exists():
                try:
                    data = vars(cclib.io.ccread(path))
                    data = self.process_data(data)
                except Exception as e:
                    self.logger.warning(f"\ncclib raised exception {e}\nPlease report!")
                    data = {}
            else:
                data = {}

            # Switch to standard names
            translation = self.metadata["translation"]
            keys = [*data.keys()]
            for key in keys:
                if key in translation:
                    to = translation[key]
                    if to in data and data[to] != data[key]:
                        self.logger.warning(
                            f"Overwriting {to} with {key}\n"
                            f"\t{data[key]} --> {data[to]}"
                        )
                    data[to] = data[key]
                    del data[key]

            # Adding timing information from SEAMM
            data["SEAMM elapsed time"] = round(t, 1)
            data["SEAMM np"] = n_threads

            self.logger.debug("after cclib")
            self.logger.debug(f"{pprint.pformat(data)}")
            self.logger.debug(80 * "*")

            success = "success" if "success" in data else False

            if not success:
                raise RuntimeError("Gaussian did not complete successfully")

            # Get the data from the formatted checkpoint file
            data = self.parse_fchk(directory / "gaussian.fchk", data)

            # Debug output
            if self.logger.isEnabledFor(logging.DEBUG):
                keys = "\n".join(data.keys())
                self.logger.debug("After parse_fchk")
                self.logger.debug(f"Data keys:\n{keys}")
                self.logger.debug(f"Data:\n{pprint.pformat(data)}")

            # And parse a bit more out of the output
            if path.exists():
                data = self.parse_output(path, data)

            # The model chemistry
            if "model" in data:
                self.model = data["model"]
            elif "method" in data:
                if data["method"].endswith("DFT"):
                    model = data["method"] + "/" + data["density functional"]
                else:
                    model = data["method"]
                if "basis set name" in data:
                    self.model = model + "/" + data["basis set name"]
                else:
                    self.model = model
            else:
                self.model = "unknown"
            self.logger.debug(f"model = {self.model}")
            data["model"] = "Gaussian/" + self.model

            # Capitalize symmetry names
            for key in ("symmetry group", "symmetry group used"):
                if key in data:
                    data[key] = data[key].capitalize()

            # Check for QCI and change the energy names
            if "QCISD" in data["model"]:
                if "E cc" in data:
                    data["E qcisd"] = data["E cc"]
                    del data["E cc"]
                if "E ccsd_t" in data:
                    data["E qcisd_t"] = data["E ccsd_t"]
                    del data["E ccsd_t"]

            # printer.normal(f"\n\n\n\n\nData:\n{pprint.pformat(data)}\n\n\n\n\n")

            # Debug output
            if self.logger.isEnabledFor(logging.DEBUG):
                keys = "\n".join(data.keys())
                self.logger.debug(f"Data keys:\n{keys}")
                self.logger.debug(f"Data:\n{pprint.pformat(data)}")

            # If ran successfully, put out the success file
            if data["success"]:
                success_file.write_text("success")

                if self._timing_data is not None:
                    self._timing_data[11] = data["model"]
                    if "metadata/symmetry_detected" in data:
                        self._timing_data[13] = str(data["metadata/symmetry_detected"])
                    else:
                        self._timing_data[13] = ""
                    if "metadata/symmetry_used" in data:
                        self._timing_data[14] = str(data["metadata/symmetry_used"])
                    else:
                        self._timing_data[14] = ""
                    if "Number of basis functions" in data:
                        self._timing_data[15] = str(data["Number of basis functions"])
                    else:
                        self._timing_data[15] = ""
                    try:
                        with self._timing_path.open("a", newline="") as fd:
                            writer = csv.writer(fd)
                            writer.writerow(self._timing_data)
                    except Exception:
                        pass

        # Add other citations here or in the appropriate place in the code.
        # Add the bibtex to data/references.bib, and add a self.reference.cite
        # similar to the above to actually add the citation to the references.
        if "G version" in data:
            try:
                template = string.Template(self._bibliography[data["G version"]])
                citation = template.substitute(
                    month=data["G month"],
                    version=data["G revision"],
                    year=data["G year"],
                )
                self.references.cite(
                    raw=citation,
                    alias="Gaussian",
                    module="gaussian_step",
                    level=1,
                    note="The principle Gaussian citation.",
                )
            except Exception:
                pass

        return data
