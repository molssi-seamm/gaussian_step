# -*- coding: utf-8 -*-

"""Non-graphical part of the Thermodynamics step in a Gaussian flowchart
"""

import logging
from pathlib import Path
import pkg_resources
import pprint  # noqa: F401

import gaussian_step  # noqa: E999
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
printer = printing.getPrinter("gaussian")

# Add this module's properties to the standard properties
path = Path(pkg_resources.resource_filename(__name__, "data/"))
csv_file = path / "properties.csv"
if path.exists():
    molsystem.add_properties_from_file(csv_file)


class Thermodynamics(gaussian_step.Optimization):
    """
    The non-graphical part of a Thermodynamics step in a flowchart.

    Attributes
    ----------
    parser : configargparse.ArgParser
        The parser object.

    options : tuple
        It contains a two item tuple containing the populated namespace and the
        list of remaining argument strings.

    subflowchart : seamm.Flowchart
        A SEAMM Flowchart object that represents a subflowchart, if needed.

    parameters : ThermodynamicsParameters
        The control parameters for Thermodynamics.

    See Also
    --------
    TkThermodynamics,
    Thermodynamics, ThermodynamicsParameters
    """

    def __init__(
        self,
        flowchart=None,
        title="Thermodynamics",
        extension=None,
        module=__name__,
        logger=logger,
    ):
        """A substep for Thermodynamics in a subflowchart for Gaussian.

        You may wish to change the title above, which is the string displayed
        in the box representing the step in the flowchart.

        Parameters
        ----------
        flowchart: seamm.Flowchart
            The non-graphical flowchart that contains this step.

        title: str
            The name displayed in the flowchart.
        extension: None
            Not yet implemented
        logger : Logger = logger
            The logger to use and pass to parent classes

        Returns
        -------
        None
        """
        logger.debug(f"Creating Thermodynamics {self}")

        super().__init__(
            flowchart=flowchart,
            title="Thermodynamics",
            extension=extension,
            module=__name__,
            logger=logger,
        )

        self._calculation = "thermodynamics"
        self._model = None
        self._metadata = gaussian_step.metadata
        self.parameters = gaussian_step.ThermodynamicsParameters()

    @property
    def header(self):
        """A printable header for this section of output"""
        return "Step {}: {}".format(".".join(str(e) for e in self._id), self.title)

    @property
    def version(self):
        """The semantic version of this module."""
        return gaussian_step.__version__

    @property
    def git_revision(self):
        """The git version of this module."""
        return gaussian_step.__git_revision__

    def description_text(self, P=None, calculation="Thermodynamics"):
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
        if P is None:
            P = self.parameters.values_to_dict()

        if (
            isinstance(P["optimize first"], bool)
            and P["optimize first"]
            or P["optimize first"] == "yes"
        ):
            calculation = " The structure will be optimized and then the "
        else:
            calculation = "The "

        calculation += "thermodynamic functions will be calculated at {T} and {P}"

        if P["optimize first"]:
            text = gaussian_step.Optimization.description_text(
                self, P=P, calculation=calculation
            )
        else:
            text = gaussian_step.Energy.description_text(
                self, P=P, calculation=calculation
            )

        if (
            isinstance(P["input only"], bool)
            and P["input only"]
            or P["input only"] == "yes"
        ):
            if type(self) is Thermodynamics:
                added = (
                    "\n\nThe input file will be written. No calculation will be run."
                )
                text += "\n" + __(added, **P, indent=4 * " ").__str__()

        return text

    def run(self, keywords=None):
        """Run a Thermodynamics step.

        Parameters
        ----------
        keywords: set()
            Any existing keywords.

        Returns
        -------
        seamm.Node
            The next node object in the flowchart.
        """
        if keywords is None:
            keywords = set()

        # Get the values of the parameters, dereferencing any variables
        P = self.parameters.current_values_to_dict(
            context=seamm.flowchart_variables._data
        )

        # Set the attribute for writing just the input
        self.input_only = P["input only"]

        # Get the current system and configuration (ignoring the system...)
        # _, configuration = self.get_system_configuration(None)

        # Figure out what we are doing!
        method, method_data = self.get_method(P)

        if method not in gaussian_step.composite_methods:
            keywords.add("Freq")

        if P["optimize first"]:
            gaussian_step.Optimization.run(self, keywords=keywords)
        else:
            gaussian_step.Energy.run(self, keywords=keywords)

    def analyze(self, indent="", data={}, out=[], table=None, P=None):
        """Analyze and print the thermodynamic and vibrational data.

        Parameters
        ----------
        indent: str
            An extra indentation for the output
        data: {}
            The data from the Gaussian calculation.
        out: []
            ?
        table: {} = None
            The tabular data for the main output.
        P: {}
            The control parameters.
        """
        if P["optimize first"]:
            gaussian_step.Optimization.analyze(self, data=data, P=P)
        else:
            gaussian_step.Energy.analyze(self, data=data, P=P)
