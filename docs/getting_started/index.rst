***************
Getting Started
***************

Installation
============
The Gaussian step is installed with the `SEAMM Manager`_, and is probably already part of
your SEAMM installation. To add it, or bring it up to date::

  seamm-manager install gaussian-step
  seamm-manager update gaussian-step

or use the Manager's window. Gaussian itself is licensed software that you install
yourself, so installing the step does not install Gaussian. Instead it creates
``~/SEAMM/gaussian.ini`` -- or ``gaussian.ini`` in whichever SEAMM installation you are
working on -- which tells SEAMM where Gaussian is and how to run it. An existing file
is never changed.

.. _SEAMM Manager: https://molssi-seamm.github.io/getting_started/installation/seamm-manager.html

Configuring Gaussian
====================
SEAMM needs to know which Gaussian to run (``g16`` or ``g09``), where it is installed,
and the profile that sets up its environment. The first time it runs Gaussian, SEAMM
looks for them itself: from ``$g16root`` or ``$g09root`` if either is set, otherwise
from ``g16`` or ``g09`` on your ``PATH``, and saves what it finds in the ``[local]``
section of ``~/SEAMM/gaussian.ini``. To give them yourself, edit that section:

.. code-block:: ini

    [local]
    installation = local
    code = g16
    root-directory = /opt/gaussian
    setup-environment = /opt/gaussian/g16/bsd/g16.profile

``root-directory`` is the directory that contains the ``g16`` (or ``g09``) directory.
The file itself has comments explaining each option.

That should be enough to get started. For more detail about the functionality in this plug-in, see the :ref:`User Guide <user-guide>`.
