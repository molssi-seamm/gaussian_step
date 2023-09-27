=======
History
=======

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
