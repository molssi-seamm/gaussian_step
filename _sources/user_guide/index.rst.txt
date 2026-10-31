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

The double hybrid B2PLYP is offered with Grimme's D3(BJ) dispersion only, which runs
Gaussian's built-in ``B2PLYPD3``, and its results are named B2PLYP-D3(BJ), as ORCA's.
With ``EmpiricalDispersion=GD3`` or ``GD3BJ``, Gaussian 09 E.01 uses the parameter s6 =
1.0 for B2PLYP instead of its 0.64, giving a dispersion energy 36 % of the r⁻⁶ term too
large (1.4 kJ/mol for octane, against Grimme's own code), and ignores the IOps that set
the parameters; a flowchart asking for GD3 with B2PLYP stops with an error.

Results stored for other steps
==============================

Whatever results are selected, the sub-steps always store the energy as ``total
energy#Gaussian#<model>`` and, for Thermodynamics, the zero-point energy (``ZPE``), the
absolute enthalpy (``H``) at the temperature (``T``) and, at 298.15 K, ``H(298.15 K)``:
the names that SEAMM's quantum-chemistry steps and the Thermochemistry step share.

They also store a fingerprint of the settings of the calculation, so that one can tell
whether results on different structures -- the energies of the species of a reaction,
say -- come from the same kind of calculation: ``settings fingerprint#Gaussian#<model>``,
a hash of the route less the task (``Opt``, ``Freq``, ...), output, starting guess and
Link 0 lines and less what belongs to the structure (its title, charge, multiplicity and
geometry, and the R or U the step chose from the multiplicity); ``element
fingerprints#Gaussian#<model>``; and ``calculation#Gaussian#<model>``, the whole record
-- the full Gaussian version, the fingerprints, a fingerprint of the geometry and the
text that was hashed -- for finding the cause of a mismatch. Thermodynamics stores them
as ``frequency fingerprint`` etc. too. The Gaussian release and revision letter count
(G09 E), not the number of the revision.


Index
=====

* :ref:`genindex`
