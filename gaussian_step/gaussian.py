# -*- coding: utf-8 -*-

"""Non-graphical part of the Gaussian step in a SEAMM flowchart"""

import logging
import pkg_resources
from pathlib import Path
import pprint  # noqa: F401
import sys

import gaussian_step
import molsystem
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

logger = logging.getLogger("Gaussian")
job = printing.getPrinter()
printer = printing.getPrinter("Gaussian")

# Add this module's properties to the standard properties
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


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
        )

        super().__init__(
            flowchart=flowchart,
            title="Gaussian",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self.parameters = gaussian_step.GaussianParameters()
        self._gversion = "g09"
        self._data = {}

    @property
    def version(self):
        """The semantic version of this module."""
        return gaussian_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return gaussian_step.__git_revision__

    @property
    def gversion(self):
        """The Gaussian version to target."""
        return self._gversion

    @gversion.setter
    def gversion(self, version):
        self._gversion = version

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
            "--max-atoms-to-print",
            default=25,
            help="Maximum number of atoms to print charges, etc.",
        )

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
        # Create the directory
        directory = Path(self.directory)
        directory.mkdir(parents=True, exist_ok=True)

        # The node after this one, to return at end
        next_node = super().run(printer)

        printer.important(self.header)
        printer.important("")

        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        while node is not None:
            if node.is_runable:
                node.run()
            node = node.next()

        # Handle any cleanup of files requested
        node = self.subflowchart.get_node("1").next()
        while node is not None:
            if node.is_runable:
                node.cleanup()
            node = node.next()

        return next_node

    def analyze(self, indent="", fchk=[], output=[], configuration=None, **kwargs):
        """Do any analysis of the output from this step.

        Also print important results to the local step.out file using
        "printer".

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        """
        # Get the first real node
        node = self.subflowchart.get_node("1").next()

        # Loop over the subnodes, asking them to do their analysis
        while node is not None:
            for value in node.description:
                printer.important(value)
            node.analyze(data=self._data, configuration=configuration)
            printer.normal("")
            node = node.next()

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
