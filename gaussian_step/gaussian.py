# -*- coding: utf-8 -*-

"""Non-graphical part of the Gaussian step in a SEAMM flowchart
"""

import logging
from pathlib import Path
import pprint  # noqa: F401
import string
import sys

import psutil

import gaussian_step
import seamm
from seamm_util import ureg, Q_  # noqa: F401
import seamm_util.printing as printing
from seamm_util.printing import FormattedText as __

# In addition to the normal logger, two logger-like printing facilities are
# defined: "job" and "printer". "job" send output to the main job.out file for
# the job, and should be used very sparingly, typically to echo what this step
# will do in the initial summary of the job.
#
# "printer" sends output to the file "step.out" in this steps working
# directory, and is used for all normal output from this step.

logger = logging.getLogger(__name__)
job = printing.getPrinter()
printer = printing.getPrinter("Gaussian")


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


class Gaussian(seamm.Node):
    """
    The non-graphical part of a Gaussian step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : GaussianParameters
        The control parameters for Gaussian.

    See Also
    --------
    TkGaussian,
    Gaussian, GaussianParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Gaussian",
        namespace="org.molssi.seamm.gaussian",
        extension=None,
        logger=logger,
    ):
        """A step for Gaussian in a SEAMM flowchart.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        namespace : str
            The namespace for the plug-ins of the subflowchart
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Gaussian {self}")
        self.subflowchart = seamm.Flowchart(
            parent=self, name="Gaussian", namespace=namespace
        )  # yapf: disable

        super().__init__(
            flowchart=flowchart,
            title="Gaussian",
            extension=extension,
            module=__name__,
            logger=logger,
        )  # yapf: disable

        self.parameters = gaussian_step.GaussianParameters()
        self._data = {}

    @property
    def version(self):
        """The semantic version of this module."""
        return gaussian_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return gaussian_step.__git_revision__

    def set_id(self, node_id):
        """Set the id for node to a given tuple"""
        self._id = node_id

        # and set our subnodes
        self.subflowchart.set_ids(self._id)

        return self.next()

    def create_parser(self):
        """Setup the command-line / config file parser"""
        # parser_name = 'gaussian-step'
        parser_name = self.step_type
        parser = self.flowchart.parser

        # Remember if the parser exists ... this type of step may have been
        # found before
        parser_exists = parser.exists(parser_name)

        # Create the standard options, e.g. log-level
        result = super().create_parser(name=parser_name)

        if parser_exists:
            return result

        # Options for Gaussian
        parser.add_argument(
            parser_name,
            "--gaussian-path",
            default="",
            help="the path to the Gaussian executable",
        )

        parser.add_argument(
            parser_name,
            "--gaussian-exe",
            default="g16",
            help="the Gaussian executable",
        )

        parser.add_argument(
            parser_name,
            "--gaussian-root",
            default="",
            help="The location of the root direction for the Gaussian installation",
        )

        parser.add_argument(
            parser_name,
            "--gaussian-environment",
            default="",
            help="A file to source to setup the Gaussian environment",
        )

        parser.add_argument(
            parser_name,
            "--ncores",
            default="4",
            help="How many threads to use in Gaussian",
        )

        parser.add_argument(
            parser_name,
            "--memory",
            default="available",
            help=(
                "The maximum amount of memory to use for Gaussian, which can be "
                "'all' or 'available', or a number, which may use k, Ki, "
                "M, Mi, etc. suffixes. Default: available."
            ),
        )

        return result

    def description_text(self, P=None):
        """Create the text description of what this step will do.
        The dictionary of control values is passed in as P so that
        the code can test values, etc.

        Parameters
        ----------
        P: dict
            An optional dictionary of the current values of the control
            parameters.
        Returns
        -------
        str
            A description of the current step.
        """
        self.subflowchart.root_directory = self.flowchart.root_directory

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        text = self.header + "\n\n"
        while node is not None:
            try:
                text += __(node.description_text(), indent=3 * " ").__str__()
            except Exception as e:
                print(f"Error describing gaussian flowchart: {e} in {node}")
                logger.critical(f"Error describing gaussian flowchart: {e} in {node}")
                raise
            except:  # noqa: E722
                print(
                    "Unexpected error describing gaussian flowchart: {} in {}".format(
                        sys.exc_info()[0], str(node)
                    )
                )
                logger.critical(
                    "Unexpected error describing gaussian flowchart: {} in {}".format(
                        sys.exc_info()[0], str(node)
                    )
                )
                raise
            text += "\n"
            node = node.next()

        return text

    def run(self):
        """Run a Gaussian step.

        Parameters
        ----------
        None

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        printer.important(self.header)
        printer.important("")

        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        # Get the system & configuration
        system, configuration = self.get_system_configuration(None)

        # Access the options
        options = self.options
        seamm_options = self.global_options

        # Work out how many cores and how much memory to use
        n_cores = psutil.cpu_count(logical=False)
        self.logger.info("The number of cores is {}".format(n_cores))

        # How many threads to use
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
        self.logger.info(f"Psi4 will use {n_threads} threads.")

        # How much memory to use
        svmem = psutil.virtual_memory()

        if seamm_options["memory"] == "all":
            mem_limit = svmem.total
        elif seamm_options["memory"] == "available":
            # For the default, 'available', use in proportion to number of
            # cores used
            mem_limit = svmem.total * (n_threads / n_cores)
        else:
            mem_limit = dehumanize(seamm_options["memory"])

        if options["memory"] == "all":
            memory = svmem.total
        elif options["memory"] == "available":
            # For the default, 'available', use in proportion to number of
            # cores used
            memory = svmem.total * (n_threads / n_cores)
        else:
            memory = dehumanize(options["memory"])

        memory = min(memory, mem_limit)

        # Apply a minimum of 800 MB
        min_memory = dehumanize("800 MB")
        if min_memory > memory:
            memory = min_memory

        # Gaussian allows no decimal points.
        memory = humanize(memory, kilo=1000)

        # The node after this one, to return at end
        next_node = super().run(printer)

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        lines = []
        lines.append("%Chk=CheckPoint")
        lines.append(f"%Mem={memory}")
        lines.append(f"%NProcShared={n_threads}")
        keywords = set()
        while node is not None:
            keywords = keywords.union(node.get_input())
            node = node.next()
        keywords.add("FormCheck=ForceCart")

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

        exe = options["gaussian_exe"]
        exe_path = options["gaussian_path"]
        if exe_path != "":
            exe = f"{exe_path}/{exe}"

        printer.important(
            self.indent + f"    Gaussian will use {n_threads} OpenMP threads and "
            f"up to {memory} of memory.\n"
        )

        if options["gaussian_root"] != "":
            env = {"g09root": options["gaussian_root"]}
        else:
            env = {}

        if options["gaussian_environment"] != "":
            cmd = f". {options['gaussian_environment']} ; {exe}"
        else:
            cmd = exe

        cmd += " < input.dat"

        local = seamm.ExecLocal()
        result = local.run(
            shell=True,
            cmd=cmd,
            files=files,
            env=env,
            return_files=[
                "CheckPoint.chk",
                "Test.FChk",
            ],
            in_situ=True,
            directory=directory,
        )
        # result = local.run(
        #     cmd=[exe],
        #     files=files,
        #     input_data=files["input.dat"],
        #     return_files=[
        #         "CheckPoint.chk",
        #         "Test.FChk",
        #     ],
        #     in_situ=True,
        #     directory=directory,
        # )

        if result is None:
            logger.error("There was an error running Gaussian")
            return None

        logger.debug("\n" + pprint.pformat(result))

        logger.info("stdout:\n" + result["stdout"])
        if result["stderr"] != "":
            logger.warning("stderr:\n" + result["stderr"])

        # Analyze the results
        self.analyze(
            fchk=result["Test.FChk"]["data"].splitlines(),
            output=result["stdout"].splitlines(),
        )

        # Add other citations here or in the appropriate place in the code.
        # Add the bibtex to data/references.bib, and add a self.reference.cite
        # similar to the above to actually add the citation to the references.
        if "G version" in self._data:
            data = self._data
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

        return next_node

    def analyze(self, indent="", fchk=[], output=[], **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        self.parse_fchk(fchk)
        self.parse_output(output)

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        # Loop over the subnodes, asking them to do their analysis
        while node is not None:
            for value in node.description:
                printer.important(value)
            node.analyze(data=self._data)
            printer.normal("")
            node = node.next()

        # Put any requested results into variables or tables
        self.store_results(data=self._data, create_tables=True)

        # Update the structure
        if "Current cartesian coordinates" in self._data:
            factor = Q_(1, "a0").to("Å").magnitude
            system_db = self.get_variable("_system_db")
            configuration = system_db.system.configuration
            xs = []
            ys = []
            zs = []
            it = iter(self._data["Current cartesian coordinates"])
            for x in it:
                xs.append(factor * x)
                ys.append(factor * next(it))
                zs.append(factor * next(it))
            configuration.atoms["x"][0:] = xs
            configuration.atoms["y"][0:] = ys
            configuration.atoms["z"][0:] = zs
            printer.important(
                self.indent + "    Updated the system with the structure from Gaussian",
            )
            printer.important("")

    def parse_fchk(self, lines):
        """Process the data of a formatted Chk file given as lines of data.

        Parameters
        ----------
        lines : iterable
            The data of the file as an iterable.
        """
        data = self._data

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
                            value.append(float(line[pos : pos + 16].strip()))
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

    def parse_output(self, lines):
        """Process the output given as lines of data.

        Parameters
        ----------
        lines : iterable
            The data of the file as an iterable.
        """
        data = self._data

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

        data["Geometry Optimization Converged"] = converged
        data["Maximum Force"] = max_force[-1]
        data["RMS Force"] = rms_force[-1]
        data["Maximum Displacement"] = max_displacement[-1]
        data["RMS Displacement"] = rms_displacement[-1]
        data["Maximum Force Trajectory"] = max_force
        data["RMS Force Trajectory"] = rms_force
        data["Maximum Displacement Trajectory"] = max_displacement
        data["RMS Displacement Trajectory"] = rms_displacement
