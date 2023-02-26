"""This file contains metadata describing the results from Gaussian
"""

metadata = {}

methods = {
    "DFT: Kohn-Sham density functional theory": {
        "method": "DFT",
        "calculation": ["energy", "gradients"],
        "level": "normal",
        "gradients": "analytic",
    },
    "HF: Hartree-Fock self consistent field (SCF)": {
        "method": "HF",
        "calculation": ["energy", "gradients"],
        "level": "normal",
        "gradients": "analytic",
    },
    "MP2: 2nd-order Møller–Plesset perturbation theory": {
        "method": "MP2",
        "calculation": ["energy", "gradients"],
        "level": "normal",
        "gradients": "analytic",
        "freeze core?": True,
    },
    "MP3: 3rd-order Møller–Plesset perturbation theory": {
        "method": "MP3",
        "calculation": ["energy", "gradients"],
        "level": "normal",
        "gradients": "analytic",
        "freeze core?": True,
    },
}

dft_functionals = {
    "B3LYP hybrid GGA exchange-correlation functional": {
        "name": "B3LYP",
        "dispersion": ["none", "GD3BJ", "GD3", "GD2"],
        "level": "normal",
    },
}

optimization_convergence = {
    "default": "",
    "tight": "Tight",
    "very tight": "VeryTight",
    "loose": "Loose",
}

"""Properties that Gaussian produces.
`metadata["results"]` describes the results that this step can produce. It is a
dictionary where the keys are the internal names of the results within this step, and
the values are a dictionary describing the result. For example::

    metadata["results"] = {
        "total_energy": {
            "calculation": [
                "energy",
                "optimization",
            ],
            "description": "The total energy",
            "dimensionality": "scalar",
            "methods": [
                "ccsd",
                "ccsd(t)",
                "dft",
                "hf",
            ],
            "property": "total energy#QuickMin#{model}",
            "type": "float",
            "units": "E_h",
        },
    }

Fields
______

calculation : [str]
    Optional metadata describing what subtype of the step produces this result.
    The subtypes are completely arbitrary, but often they are types of calculations
    which is why this is name `calculation`. To use this, the step or a substep
    define `self._calculation` as a value. That value is used to select only the
    results with that value in this field.

description : str
    A human-readable description of the result.

dimensionality : str
    The dimensions of the data. The value can be "scalar" or an array definition
    of the form "[dim1, dim2,...]". Symmetric tringular matrices are denoted
    "triangular[n,n]". The dimensions can be integers, other scalar
    results, or standard parameters such as `n_atoms`. For example, '[3]',
    [3, n_atoms], or "triangular[n_aos, n_aos]".

methods : str
    Optional metadata like the `calculation` data. `methods` provides a second
    level of filtering, often used for the Hamiltionian for *ab initio* calculations
    where some properties may or may not be calculated depending on the type of
    theory.

property : str
    An optional definition of the property for storing this result. Must be one of
    the standard properties defined either in SEAMM or in this steps property
    metadata in `data/properties.csv`.

type : str
    The type of the data: string, integer, or float.

units : str
    Optional units for the result. If present, the value should be in these units.
"""
metadata["results"] = {
    "Virial Ratio": {
        "calculation": ["energy", "optimization", "thermodynamics", "vibrations"],
        "description": "the virial ratio",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "format": ".4f",
    },
    "Total Energy": {
        "calculation": ["energy", "optimization", "thermodynamics", "vibrations"],
        "description": "total energy including all terms",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "units": "E_h",
        "format": ".6f",
    },
    "RMS Density": {
        "calculation": ["energy", "optimization", "thermodynamics", "vibrations"],
        "description": "the RMS density difference in the SCF",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "format": ".2e",
    },
    "Geometry Optimization Converged": {
        "calculation": ["optimization"],
        "description": "whether the geometry optimzation converged",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "boolean",
        "format": "s",
    },
    "RMS Force": {
        "calculation": ["optimization"],
        "description": "the RMS force on the atoms",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "units": "E_h/Å",
        "format": ".6f",
    },
    "Maximum Force": {
        "calculation": ["optimization"],
        "description": "the maximum force on an atom",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "units": "E_h/Å",
        "format": ".6f",
    },
    "RMS Displacement": {
        "calculation": ["optimization"],
        "description": "the RMS displacement of the atoms",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "units": "Å",
        "format": ".6f",
    },
    "Maximum Displacement": {
        "calculation": ["optimization"],
        "description": "the maximum displacement of an atom",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "units": "Å",
        "format": ".6f",
    },
    "RMS Force Threshold": {
        "calculation": ["optimization"],
        "description": "the RMS force threshold for the atoms",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "units": "E_h/Å",
        "format": ".6f",
    },
    "Maximum Force Threshold": {
        "calculation": ["optimization"],
        "description": "the maximum force threshold for an atom",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "units": "E_h/Å",
        "format": ".6f",
    },
    "RMS Displacement Threshold": {
        "calculation": ["optimization"],
        "description": "the RMS displacement threshold for the atoms",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "units": "Å",
        "format": ".6f",
    },
    "Maximum Displacement Threshold": {
        "calculation": ["optimization"],
        "description": "the maximum displacement threshold for an atom",
        "dimensionality": "scalar",
        "methods": ["DFT", "HF", "MP2", "MP3", "MP4"],
        "type": "float",
        "units": "Å",
        "format": ".6f",
    },
    "Gaussian Version": {
        "description": "the version of Gaussian used",
        "dimensionality": "scalar",
        "type": "string",
        "format": "s",
    },
}
