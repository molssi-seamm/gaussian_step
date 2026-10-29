# -*- coding: utf-8 -*-

"""The fingerprints of the settings of Gaussian calculations
(gaussian_step.fingerprint), and that the Energy step stores them."""

import os
from pathlib import Path

import pytest

from gaussian_step.fingerprint import (
    gaussian_calculation_record,
    gaussian_version,
    settings_text,
)

SOURCE = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "run_path"

INPUT = """\
%OldChk=/somewhere/1.chk
%Chk=gaussian.chk
%Mem=2GB
%NProcShared=4
# UB3LYP/6-31G* FORCE Integral(Grid=UltraFine) Pop=NBORead
# Punch=(Coord,Derivatives)

Energy of ethyl/ethyl

0    2
C    0.000000   0.000000   0.000000

$nbo bndidx $end
"""
VERSION = {"G version": "G09", "G revision": "E.01"}


def fingerprint(text, data=VERSION, spin_prefix="U"):
    return gaussian_calculation_record(text, data, spin_prefix=spin_prefix)[
        "fingerprint"
    ]


def test_what_is_hashed():
    assert settings_text(INPUT, spin_prefix="U") == [
        "# B3LYP/6-31G* INTEGRAL(GRID=ULTRAFINE)"
    ]
    # A spin the user chose is kept
    assert settings_text(INPUT) == ["# INTEGRAL(GRID=ULTRAFINE) UB3LYP/6-31G*"]


def test_same_settings_written_differently():
    """Order and case, the task, Link 0, the guess, the species and the spin the
    step chose from the multiplicity do not count."""
    other = """\
%Mem=1GB
#p Opt=(Tight, MaxCycles=50) integral(grid=ultrafine) RB3LYP/6-31G* Guess=Read

Closed shell

0 1
C 0 0 0

"""
    assert fingerprint(INPUT) == fingerprint(other, spin_prefix="R")
    # The revision's number does not count, its letter does
    assert fingerprint(INPUT) == fingerprint(
        other, data={"G version": "G09", "G revision": "E.02"}, spin_prefix="R"
    )
    assert fingerprint(INPUT) != fingerprint(
        other, data={"G version": "G09", "G revision": "D.01"}, spin_prefix="R"
    )
    assert fingerprint(INPUT) != fingerprint(
        other, data={"G version": "G16", "G revision": "E.01"}, spin_prefix="R"
    )


@pytest.mark.parametrize(
    "change",
    [
        ("UltraFine", "SuperFine"),
        ("6-31G*", "6-31G**"),
        ("UB3LYP", "UPBE1PBE"),
        ("FORCE", "FORCE SCF=Tight"),
    ],
)
def test_settings_that_count(change):
    old, new = change
    assert fingerprint(INPUT) != fingerprint(INPUT.replace(old, new, 1))


def test_forced_spin_counts():
    rohf = INPUT.replace("UB3LYP", "ROB3LYP")
    assert fingerprint(rohf, spin_prefix=None) != fingerprint(INPUT, spin_prefix=None)


def test_gaussian_version():
    assert gaussian_version(VERSION) == ("G09 E.01", "G09 E")
    assert gaussian_version({}) == ("", "")


SPEC = """\
title: Gaussian fingerprint test
steps:
- Water: {}
- Gaussian:
    steps:
    - Energy:
        basis: 6-31G*
"""


def test_energy_step_stores_fingerprints(tmp_path, monkeypatch):
    import molsystem
    from seamm_exec.testing import fake_program, run_spec

    bin = tmp_path / "bin"
    code = fake_program(
        bin / "g09",
        files=[DATA / "fort.7", DATA / "gaussian.chk"],
        stdout=DATA / "output.txt",
        need_input=False,
    )
    fake_program(bin / "formchk", files=[DATA / "gaussian.fchk"], need_input=False)
    monkeypatch.setenv("PATH", f"{bin}{os.pathsep}{os.environ['PATH']}")
    job = run_spec(
        tmp_path,
        SPEC,
        inis={
            "gaussian": {"code": code, "root-directory": "", "setup-environment": ""}
        },
        args=["--ncores", "1"],
        source=SOURCE,
    )

    db = molsystem.SystemDB(filename=str(job / "seamm.db"))
    try:
        configuration = db.system.configuration
        properties = configuration.properties
        stored = properties.get("settings fingerprint#Gaussian#*")
        assert len(stored) == 1
        ((name, value),) = stored.items()
        model = name.split("#", 2)[2]
        record = properties.get(f"calculation#Gaussian#{model}")
        record = record[f"calculation#Gaussian#{model}"]["value"]
        energy = properties.get(f"total energy#Gaussian#{model}")
        energy = energy[f"total energy#Gaussian#{model}"]["value"]
    finally:
        db.close()
    assert record["fingerprint"] == value["value"]
    assert record["version"] == "G09 E.01"
    assert record["files"] == {"input.dat": "# B3LYP/6-31G* INTEGRAL(GRID=ULTRAFINE)"}
    assert record["geometry fingerprint"].startswith("sha256:")
    assert energy == pytest.approx(-76.4087114098, abs=1.0e-6)


def test_b2plyp_is_offered_with_gd3():
    """B2PLYP, the double hybrid of CBH-ANL, with Grimme's D3 (GD3, zero damping)."""
    from pathlib import Path

    import gaussian_step

    name = "B2PLYP : double hybrid of Grimme, with MP2 correlation"
    data = gaussian_step.dft_functionals[name]
    assert data["name"] == "B2PLYP"
    assert data["dispersion"] == ["none", "GD3BJ", "GD3"]
    bib = (Path(gaussian_step.__file__).parent / "data" / "references.bib").read_text()
    for _, key in data["citations"]:
        assert "{" + key + "," in bib
    # the dispersion is a setting
    route = (
        "# B2PLYP/cc-pVTZ Opt Freq Integral(Grid=UltraFine)\n\nt\n\n0 1\nC 0 0 0\n\n"
    )
    with_gd3 = route.replace("Freq", "Freq EmpiricalDispersion=GD3")
    assert fingerprint(route, spin_prefix="R") != fingerprint(with_gd3, spin_prefix="R")


def test_dispersion_in_the_model():
    """The dispersion correction is part of the model, as ORCA's D3ZERO/D3BJ, and
    the atoms' references are those without it."""
    from gaussian_step.substep import model_with_dispersion, model_without_dispersion

    gd3 = {"Opt", "EmpiricalDispersion=GD3"}
    assert model_with_dispersion("B2PLYP-FC/cc-pVTZ", gd3) == "B2PLYP-D3(0)-FC/cc-pVTZ"
    assert (
        model_with_dispersion("B3LYP/6-31G(d)", ["EmpiricalDispersion=GD3BJ"])
        == "B3LYP-D3(BJ)/6-31G(d)"
    )
    assert model_with_dispersion("B3LYP/6-31G(d)", ["Opt"]) == "B3LYP/6-31G(d)"
    assert model_with_dispersion("unknown", gd3) == "unknown"
    assert model_without_dispersion("B2PLYP-D3(0)-FC/cc-pVTZ") == "B2PLYP-FC/cc-pVTZ"
    assert model_without_dispersion("B3LYP-D3(BJ)/6-31G(d)") == "B3LYP/6-31G(d)"
