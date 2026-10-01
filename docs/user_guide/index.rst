.. _user-guide:

**********
User Guide
**********
The Gaussian plug-in provides an interface to the Gaussian quantum chemistry code. 

..
   The following sections cover accessing and controlling this functionality.

   .. toctree::
      :maxdepth: 2
      :titlesonly:

Settings that depend on each other
==================================

Which settings apply depends on others: the functional and integration grid only for
DFT, a dispersion correction only from the functional's own list, freezing the core only
for methods that can, and the basis set only for methods that use one. The step's
dialogs show only the settings that apply with the current choices, and the same rules
are used when a flowchart is built or edited without the editor (``seamm-flowchart`` or
SEAMM's MCP server): a setting that would have no effect is refused, with the reason,
and a value that contradicts another is refused too. See "Flowcharts without the editor"
in SEAMM's user guide.

Each sub-step offers its own naming choice for the system and configuration (for
example "optimized with {model}") alongside the usual ones. Flowcharts saved with the
dispersion correction ``DG2``, a misspelling of ``GD2`` that Gaussian rejects, are read
as ``GD2``.


Index
=====

* :ref:`genindex`
