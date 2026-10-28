# -*- coding: utf-8 -*-

"""Fingerprints of the settings of Gaussian calculations.

What Gaussian's input holds that does not belong to the settings is removed here,
and :func:`seamm_util.fingerprint.calculation_record` hashes the rest: the route,
with its keywords in any order and in upper case, since Gaussian's input is
case-insensitive.

Removed:

* the Link 0 lines (``%Chk``, ``%OldChk``, ``%Mem``, ``%NProcShared``, ...);
* the task (``SP``, ``Force``, ``Opt``, ``Freq``), output (``Pop``, ``GFInput``,
  ``Output=wfx``, ``Punch``, ...) and the starting guess and geometry from a
  checkpoint (``Guess``, ``Geom``), which do not change the energy at a geometry;
* the title, charge, multiplicity and geometry, which belong to the species, and
  the sections after them (initial force constants, NBO input, the wfx file);
* the R or U the step puts before the method when it chooses from the
  multiplicity (restricted for singlets, unrestricted otherwise). When the user
  chose, it stays: ROHF is not the same calculation as UHF.

The step only uses basis sets by name, so there are no per-element blocks.
"""

import re

from seamm_util.fingerprint import calculation_record

#: Route keywords that do not change the energy at a given geometry
_DROPPED_KEYWORDS = {
    "SP",
    "FORCE",
    "OPT",
    "FREQ",
    "IRC",
    "POP",
    "GFINPUT",
    "GFPRINT",
    "OUTPUT",
    "PUNCH",
    "GEOM",
    "GUESS",
    "FORMCHECK",
    "TEST",
}


def gaussian_version(data):
    """The Gaussian version from the step's parsed data, e.g. "G09 E.01", and the
    part that is hashed, e.g. "G09 E": the release and the revision letter, not
    the revision's number."""
    release = str(data.get("G version", "") or "").strip()
    revision = str(data.get("G revision", "") or "").strip()
    version = f"{release} {revision}".strip()
    hashed = f"{release} {revision.split('.')[0]}".strip()
    return version, hashed


def _route_words(text):
    """The words of a route, keeping parenthesized options together, e.g.
    'Integral(Grid=UltraFine)' or 'Opt=(Tight, MaxCycles=50)'."""
    words = []
    word = ""
    depth = 0
    for c in text:
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        if c.isspace() and depth == 0:
            if word:
                words.append(word)
            word = ""
        else:
            word += c
    if word:
        words.append(word)
    return words


def settings_text(input_text, spin_prefix=None):
    """The route keywords of a Gaussian input that make up its settings.

    Parameters
    ----------
    input_text : str
        The contents of the input file.
    spin_prefix : str or None
        The R, RO or U that the step put before the method because of the
        multiplicity, to remove; None if the user chose it.

    Returns
    -------
    [str]
        The route, as one normalized line.
    """
    route = []
    in_route = False
    for line in input_text.splitlines():
        stripped = line.strip()
        if not in_route:
            if stripped.startswith("%") or stripped == "":
                continue
            if stripped.startswith("#"):
                in_route = True
            else:
                continue
        if stripped == "":
            break  # the end of the route; the rest belongs to the species
        if stripped.startswith("#"):
            stripped = re.sub(r"^#[NPT]?", "", stripped, flags=re.IGNORECASE)
        route.append(stripped)

    words = set()
    for word in _route_words(" ".join(route)):
        word = word.upper()
        name = re.split(r"[=(]", word, maxsplit=1)[0]
        if name in _DROPPED_KEYWORDS:
            continue
        if spin_prefix and "/" in word and word.startswith(spin_prefix.upper()):
            word = word[len(spin_prefix) :]
        words.add(word)
    return ["# " + " ".join(sorted(words))] if words else []


def gaussian_calculation_record(
    input_text, data, spin_prefix=None, atomic_numbers=None, coordinates=None
):
    """The record of a Gaussian calculation's settings, with its fingerprints.

    Parameters
    ----------
    input_text : str
        The input file as run.
    data : dict
        The step's parsed results, for the Gaussian version.
    spin_prefix : str or None
        See :func:`settings_text`.
    atomic_numbers, coordinates : optional
        The geometry, for its fingerprint.

    Returns
    -------
    dict
        See :func:`seamm_util.fingerprint.calculation_record`.
    """
    version, hashed = gaussian_version(data)
    return calculation_record(
        "Gaussian",
        version,
        {"input.dat": settings_text(input_text, spin_prefix=spin_prefix)},
        atomic_numbers=atomic_numbers,
        coordinates=coordinates,
        hashed_version=hashed,
    )
