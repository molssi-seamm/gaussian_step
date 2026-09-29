# -*- coding: utf-8 -*-

"""Installer for the Gaussian plug-in.

Gaussian is licensed software that users install themselves, so this installer
does not install the code. It gives the SEAMM installation a ``gaussian.ini`` to
edit, from the template in ``data/gaussian.ini``, which says how to run Gaussian.
"""

import importlib
import logging
import shutil

import seamm_installer

logger = logging.getLogger(__name__)


class Installer(seamm_installer.InstallerBase):
    """Give the installation a gaussian.ini describing how to run Gaussian."""

    def __init__(self, logger=logger):
        super().__init__(logger=logger)

        logger.debug("Initializing the Gaussian installer object.")

        self.section = "gaussian-step"
        self.executables = ["g16"]
        self.resource_path = importlib.resources.files("gaussian_step") / "data"

    def exe_version(self, config):
        """Return the name and version of Gaussian.

        Gaussian has no option that prints its version, so this reports which
        Gaussian ('g16' or 'g09') the configuration uses and whether it is found.
        """
        code = config.get("code", "g16")
        found = shutil.which(code) is not None or config.get("root-directory", "")
        return "Gaussian", code if found else f"{code} (not found)"
