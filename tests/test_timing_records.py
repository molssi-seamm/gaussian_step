# -*- coding: utf-8 -*-
"""The timing records Gaussian runs write (seamm_exec campaign 2026-10-05)."""

from types import SimpleNamespace

from gaussian_step.substep import task_kind, timing_descriptors

LOG = """\
    25 basis functions,    42 primitive gaussians,    25 cartesian basis functions
     5 alpha electrons        5 beta electrons
 SCF Done:  E(RB3LYP) =  -76.4089533203     A.U. after   11 cycles
 SCF Done:  E(RB3LYP) =  -76.4089612345     A.U. after    6 cycles
 Job cpu time:       0 days  0 hours  0 minutes 32.5 seconds.
 Elapsed time:       0 days  0 hours  0 minutes  8.2 seconds.
 Normal termination of Gaussian 09 at Mon Oct  5 10:07:06 2026.
"""


def test_task_kind():
    assert task_kind("B3LYP/6-31G** SCF=Tight") == "energy"
    assert task_kind("Force B3LYP/6-31G**") == "gradient"
    assert task_kind("Opt=Redundant B3LYP/6-31G**") == "opt"
    assert task_kind("Opt Freq B3LYP/6-31G**") == "opt+freq"
    assert task_kind("Freq=Raman B3LYP/6-31G**") == "freq"


def test_descriptors():
    conf = SimpleNamespace(
        atoms=SimpleNamespace(atomic_numbers=[8, 1, 1]),
        charge=0,
        spin_multiplicity=1,
        periodicity=0,
    )
    data = {
        "model": "Gaussian/B3LYP/6-31G(d,p)",
        "metadata/symmetry_detected": "C2v",
        "metadata/symmetry_used": "C2v",
        "success": True,
    }
    d = timing_descriptors("Opt=Redundant B3LYP/6-31G**", data, LOG, conf)
    assert d["method"] == "B3LYP" and d["basis"] == "6-31G(d,p)"
    assert d["task"] == "opt"
    assert d["symmetry_used"] == "C2v"
    assert d["nbf"] == 25  # from the log, since data had none
    assert d["n_electrons"] == 10
    assert d["scf_runs"] == 2 and d["scf_cycles"] == 17
    assert abs(d["cpu_seconds"] - 32.5) < 1e-9
    assert abs(d["code_seconds"] - 8.2) < 1e-9
    assert d["terminated_normally"] is True
    assert d["n_atoms"] == 3 and d["n_heavy"] == 1
