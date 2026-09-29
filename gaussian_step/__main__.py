# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""Handle the installation of the Gaussian step."""

from .installer import Installer


def run():
    """Give the SEAMM installation a gaussian.ini describing how to run Gaussian."""
    installer = Installer()
    installer.run()


if __name__ == "__main__":
    run()
