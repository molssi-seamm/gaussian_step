# -*- coding: utf-8 -*-

"""
gaussian_step
A SEAMM plug-in for Gaussian
"""

# Bring up the classes so that they appear to be directly in
# the gaussian_step package.

from .metadata import methods  # noqa: F401
from .metadata import Gn_composite_methods, CBS_composite_methods  # noqa: F401
from .metadata import composite_methods  # noqa: F401
from .metadata import dft_functionals  # noqa: F401
from .metadata import optimization_convergence  # noqa: F401
from .metadata import metadata  # noqa: F401

from .gaussian import Gaussian  # noqa: F401, E501
from .gaussian_parameters import GaussianParameters  # noqa: F401, E501
from .gaussian_step import GaussianStep  # noqa: F401, E501
from .tk_gaussian import TkGaussian  # noqa: F401, E501

from .energy import Energy  # noqa: F401
from .energy_parameters import EnergyParameters  # noqa: F401
from .energy_step import EnergyStep  # noqa: F401
from .tk_energy import TkEnergy  # noqa: F401

from .optimization import Optimization  # noqa: F401
from .optimization_parameters import OptimizationParameters  # noqa: F401
from .optimization_step import OptimizationStep  # noqa: F401
from .tk_optimization import TkOptimization  # noqa: F401

from .thermodynamics_step import ThermodynamicsStep  # noqa: F401
from .thermodynamics import Thermodynamics  # noqa: F401
from .thermodynamics_parameters import ThermodynamicsParameters  # noqa: F401
from .tk_thermodynamics import TkThermodynamics  # noqa: F401

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@vt.edu"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
