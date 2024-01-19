# -*- coding: utf-8 -*-

"""Setup and run Gaussian"""

import configparser
import gzip
import logging
from pathlib import Path
import pprint
import re
import shutil
import string

import cclib
import numpy as np

import gaussian_step
import seamm
import seamm_exec
import seamm.data
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
    def method(self):
        """The method ... HF, DFT, ... used."""
        return self._method

    @method.setter
    def method(self, value):
        self._method = value

    @property
    def options(self):
        return self.parent.options

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
        data["basis"] = line[40:70].strip()

        # The rest of the file consists of a line defining the data.
        # If the data is a scalar, it is on the control line, otherwise it follows
        while True:
            try:
                line = next(it)
            except StopIteration:
                break
            key = line[0:40].strip()
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
                                logger.warning(
                                    f"Could not find the Gaussian citation: {e}"
                                )
                            break
                break

        # And the optimization steps, if any.
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
                    data["Maximum Force Threshold"] = float(threshold)
                    if criterion != "YES":
                        converged = False

                tmp1, tmp2, value, threshold, criterion = next(it).split()
                if tmp1 == "RMS" and tmp2 == "Force":
                    rms_force.append(float(value))
                    data["RMS Force Threshold"] = float(threshold)
                    if criterion != "YES":
                        converged = False

                tmp1, tmp2, value, threshold, criterion = next(it).split()
                if tmp1 == "Maximum" and tmp2 == "Displacement":
                    max_displacement.append(float(value))
                    data["Maximum Displacement Threshold"] = float(threshold)
                    if criterion != "YES":
                        converged = False

                tmp1, tmp2, value, threshold, criterion = next(it).split()
                if tmp1 == "RMS" and tmp2 == "Displacement":
                    rms_displacement.append(float(value))
                    data["RMS Displacement Threshold"] = float(threshold)
                    if criterion != "YES":
                        converged = False

        if converged is not None:
            data["Geometry Optimization Converged"] = converged
            data["Maximum Force"] = max_force[-1]
            data["RMS Force"] = rms_force[-1]
            data["Maximum Displacement"] = max_displacement[-1]
            data["RMS Displacement"] = rms_displacement[-1]
            data["Maximum Force Trajectory"] = max_force
            data["RMS Force Trajectory"] = rms_force
            data["Maximum Displacement Trajectory"] = max_displacement
            data["RMS Displacement Trajectory"] = rms_displacement

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

        if P["method"][0:4] == "CBS-":
            # Need last section
            if P["method"] in gaussian_step.methods:
                method = gaussian_step.methods[P["method"]]["method"]
            else:
                method = P["method"]

            match = f"{method} Enthalpy="
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
                        if method in key:
                            key = key.split(" ", 1)[1]
                        data[f"Composite/{key}"] = value
                data["Composite/model"] = method
                data["Composite/summary"] = "\n".join(text)
                data["Total Energy"] = data["Composite/Free Energy"]

        # Gn calculations. No header!!!!!

        # Temperature=              298.150000 Pressure=                      1.000000
        # E(ZPE)=                     0.050251 E(Thermal)=                    0.053306
        # E(CCSD(T))=               -78.321715 E(Empiric)=                   -0.041682
        # DE(Plus)=                  -0.005930 DE(2DF)=                      -0.076980
        # E(Delta-G3XP)=             -0.117567 DE(HF)=                       -0.008255
        # G4(0 K)=                  -78.521880 G4 Energy=                   -78.518825
        # G4 Enthalpy=              -78.517880 G4 Free Energy=              -78.542752

        if P["method"][0:2] in ("G1", "G2", "G3", "G4"):
            # Need last section
            method = P["method"][0:2]
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
                        if method in key:
                            key = key.split(" ", 1)[1]
                        elif key == "E(Empiric)":
                            key = "E(empirical)"
                        data[f"Composite/{key}"] = value

                data["Composite/model"] = method
                tmp = " " * 20 + f"{method[0:2]} composite method extrapolation\n\n"
                data["Composite/summary"] = tmp + "\n".join(text)
                data["Total Energy"] = data["Composite/Free Energy"]

        return data

    def process_data(self, data):
        """Massage the cclib data to a more easily used form."""
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
                    new[f"N({letter}-homo)"] = homo + 1
                    new[f"E({letter}-homo)"] = Es[homo]
                    if homo > 0:
                        new[f"E({letter}-homo-1)"] = Es[homo - 1]
                    if homo + 1 < len(Es):
                        new[f"E({letter}-lumo)"] = Es[homo + 1]
                        new[f"E({letter}-gap)"] = Es[homo + 1] - Es[homo]
                    if homo + 2 < len(Es):
                        new[f"E({letter}-lumo+1)"] = Es[homo + 2]
                    if "mosyms" in new:
                        syms = new["mosyms"][i]
                        new[f"Sym({letter}-homo)"] = syms[homo]
                        if homo > 0:
                            new[f"Sym({letter}-homo-1)"] = syms[homo - 1]
                        if homo + 1 < len(syms):
                            new[f"Sym({letter}-lumo)"] = syms[homo + 1]
                        if homo + 2 < len(syms):
                            new[f"Sym({letter}-lumo+1)"] = syms[homo + 2]
            else:
                Es = new["moenergies"][0]
                homo = homos[0]
                new["N(homo)"] = homo + 1
                new["E(homo)"] = Es[homo]
                if homo > 0:
                    new["E(homo-1)"] = Es[homo - 1]
                if homo + 1 < len(Es):
                    new["E(lumo)"] = Es[homo + 1]
                    new["E(gap)"] = Es[homo + 1] - Es[homo]
                if homo + 2 < len(Es):
                    new["E(lumo+1)"] = Es[homo + 2]
                if "mosyms" in new:
                    syms = new["mosyms"][0]
                    new["Sym(homo)"] = syms[homo]
                    if homo > 0:
                        new["Sym(homo-1)"] = syms[homo - 1]
                    if homo + 1 < len(syms):
                        new["Sym(lumo)"] = syms[homo + 1]
                    if homo + 2 < len(syms):
                        new["Sym(lumo+1)"] = syms[homo + 2]

        # moments
        if "moments" in new:
            moments = new["moments"]
            new["multipole_reference"] = moments[0]
            new["dipole_moment"] = moments[1]
            new["dipole_moment_magnitude"] = np.linalg.norm(moments[1])
            if len(moments) > 2:
                new["quadrupole_moment"] = moments[2]
            if len(moments) > 3:
                new["octapole_moment"] = moments[3]
            if len(moments) > 4:
                new["hexadecapole_moment"] = moments[4]
            del new["moments"]

        for key in ("metadata/symmetry_detected", "metadata/symmetry_used"):
            if key in new:
                new[key] = new[key].capitalize()

        return new

    def run_gaussian(self, keywords):
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
        success = directory / "success.dat"
        if not success.exists():
            # Get the system & configuration
            system, configuration = self.get_system_configuration(None)

            # Access the options
            options = self.options
            seamm_options = self.global_options

            # Get the computational environment and set limits
            ce = seamm_exec.computational_environment()

            # How many threads to use
            n_cores = ce["NTASKS"]
            self.logger.info("The number of cores available is {}".format(n_cores))

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
            self.logger.info(f"Gaussian will use {n_threads} threads.")

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

            # keywords.add("FormCheck=ForceCart")
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

            files = {"input.dat": "\n".join(lines)}
            logger.info("input.dat:\n" + files["input.dat"])

            if self.input_only:
                # Just write the input files and stop
                for filename in files:
                    path = directory / filename
                    path.write_text(files[filename])
                data = {}
            else:
                # Run Gaussian
                executor = self.parent.flowchart.executor

                # Read configuration file for Gaussian
                ini_dir = Path(seamm_options["root"]).expanduser()
                full_config = configparser.ConfigParser()
                full_config.read(ini_dir / "gaussian.ini")
                executor_type = executor.name
                if executor_type not in full_config:
                    raise RuntimeError(
                        f"No section for '{executor_type}' in the Gaussian ini file "
                        f"({ini_dir / 'mopac.ini'})"
                    )
                config = dict(full_config.items(executor_type))

                # Set up the environment
                if config["root-directory"] != "":
                    env = {"g09root": config["root-directory"]}
                else:
                    env = {}

                if config["setup-environment"] != "":
                    cmd = [". {setup-environment} ; {code}"]
                else:
                    cmd = ["{code}"]

                cmd.append(" < input.dat > output.txt ; formchk gaussian.chk")

                printer.important(
                    self.indent + f"    Gaussian will use {n_threads} OpenMP threads "
                    f"and up to {memory} of memory.\n"
                )

                return_files = [
                    "output.txt",
                    "gaussian.chk",
                    "gaussian.fchk",
                ]
                result = executor.run(
                    cmd=cmd,
                    ce=ce,
                    config=config,
                    directory=self.directory,
                    files=files,
                    return_files=return_files,
                    in_situ=True,
                    shell=True,
                    env=env,
                )

                if not result:
                    self.logger.error("There was an error running Gaussian")
                    return None

                self.logger.debug("\n" + pprint.pformat(result))

        if not self.input_only:
            # And output
            path = directory / "output.txt"
            if path.exists():
                data = vars(cclib.io.ccread(path))
                data = self.process_data(data)
            else:
                data = {}

            # Get the data from the formatted checkpoint file
            data = self.parse_fchk(directory / "gaussian.fchk", data)

            # And parse a bit more out of the output
            if path.exists():
                data = self.parse_output(path, data)

            # Debug output
            if self.logger.isEnabledFor(logging.INFO):
                keys = "\n".join(data.keys())
                logger.info(f"Data keys:\n{keys}")
            if self.logger.isEnabledFor(logging.DEBUG):
                keys = "\n".join(data.keys())
                logger.info(f"Data keys:\n{keys}")
                logger.debug(f"Data:\n{pprint.pformat(data)}")

            # The model chemistry
            # self.model = f"{data['metadata/functional']}/{data['metadata/basis_set']}"
            if "Composite/model" in data:
                self.model = data["Composite/model"]
            else:
                self.model = (
                    f"{data['metadata/methods'][-1]}/{data['method']}/"
                    f"{data['metadata/basis_set']}"
                )
            logger.info(f"model = {self.model}")

            # If ran successfully, put out the success file
            if data["success"]:
                success.write_text("success")

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
