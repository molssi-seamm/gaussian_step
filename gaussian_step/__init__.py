# -*- coding: utf-8 -*-

"""
gaussian_step
A SEAMM plug-in for Gaussian
"""

# Bring up the classes so that they appear to be directly in
# the gaussian_step package.

from gaussian_step.gaussian import Gaussian  # noqa: F401, E501
from gaussian_step.gaussian_parameters import GaussianParameters  # noqa: F401, E501
from gaussian_step.gaussian_step import GaussianStep  # noqa: F401, E501
from gaussian_step.tk_gaussian import TkGaussian  # noqa: F401, E501

# Handle versioneer
from ._version import get_versions

__author__ = "Paul Saxe"
__email__ = "psaxe@vt.edu"
versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions
