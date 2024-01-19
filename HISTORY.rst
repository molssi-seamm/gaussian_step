=======
History
=======

2024.1.19: Switched to new way to run Gaussian, added option to just write input file
    * Switched to using the new way to run executables, which supports containers.
    * Added an option to just write the input file, without running
      Gaussian. This is useful for debugging, and for running Gaussian
      on a remote server.

2023.10.25 Bug fixes: variable for functional, and parsing FChk file
    * Fixed a problem with handling the functional if it was a variable rather than a
      specific functional.
    * Fixed a problem parsing the FChk file. For exponents > 99 the FORTRAN format used
      in Gaussian grops the "E", resulting in numbers like 0.947-104 that caused a
      problem when trying to read them.
      
2023.10.22 Bug fixes: orbital plots and output
    * The plots of the HOMO and LUMO were shifted by one orbital due to some code
      counting from 1 and other, from 0. Sigh.
    * The output to Job.out was inadvertently truncated.

2023.10.7 Added structure file for plots of density and orbitals.
    * Always write the current structure as 'structure.sdf' in the directory where the
      cube files for orbitals and densities are written. The Dashboard picks up this
      file to render the structure along with the surfaces.
      
2023.9.27 Added composite and other methods, DFT functionals
    * Now support HF, DFT, MP4, CCD & CCSD, CBS-x, and Gn methods
    * Added PBE, PBE-98, PBE0, and HSE06 functionals
    * Added analysis of HOMO/LUMO gap energy
    * Added plotting of orbitals and densities
    * Added otuput of atomic charges and spins, and placing them on the configuration.
    * Added ability to control the system/configuration update

2023.2.26.1 Moved Gaussian output to output.txt
    * Capturing stdout prevent users from seeing the output during a calculation.
      This fixes that.
      
2023.2.26 Initial version with energy and optimization 
    * Support running the energy or optimization with HF, DFT, MP2 and MP3 though
      testing has not yet been thorough.
    * The DFT functional supported are at the moment limited.
      
2023.2.24
    * Plug-in created using the SEAMM plug-in cookiecutter.
