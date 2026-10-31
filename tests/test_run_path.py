# -*- coding: utf-8 -*-

"""Gaussian's run path, end to end: a flowchart built from a spec and run by
``run_flowchart``, through the Energy sub-step -- the input, ``gaussian.ini``, the
task, the output and its analysis (seamm_exec.testing).

* With a fake ``g09`` and ``formchk`` replaying a real run's output, so the path
  runs in CI.
* With the real Gaussian when it is installed (skipped otherwise, as in CI).
"""

import configparser
import os
from pathlib import Path

import pytest

from seamm_exec.testing import fake_program, run_spec, table_value

SOURCE = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "run_path"

SPEC = """\
title: Gaussian run-path test
steps:
- Water: {}
- Gaussian:
    steps:
    - Energy:
        basis: 6-31G*
"""

ENERGY = -76.408711  # E_h, B3LYP/6-31G*, Gaussian 09 E.01


def installed_gaussian():
    """The [local] section of the installation's gaussian.ini, if it has a code."""
    root = os.environ.get("SEAMM_ROOT", "").strip() or "~/SEAMM"
    config = configparser.ConfigParser()
    try:
        config.read(Path(root).expanduser() / "gaussian.ini")
        local = dict(config["local"])
    except (configparser.Error, KeyError):
        return None
    return local if local.get("code") else None


def check_job(job):
    run = job / "2" / "1"
    text = (run / "input.dat").read_text()
    assert "B3LYP/6-31G*" in text
    assert "Normal termination" in (run / "output.txt").read_text()
    assert (run / "gaussian.fchk").exists()
    assert table_value(job / "job.out", "energy") == pytest.approx(ENERGY, abs=1.0e-6)


def test_run_path_with_a_fake_gaussian(tmp_path, monkeypatch):
    # Gaussian reads input.dat on its standard input and writes the punch and
    # checkpoint files (a placeholder, since formchk is faked too);
    # the step then runs formchk, found on the PATH.
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
    check_job(job)
    # The fake ran, not a Gaussian on the PATH: the log is the recorded one.
    log = (job / "2" / "1" / "output.txt").read_text()
    assert "Normal termination of Gaussian 09 at Wed Oct 28 14:19:51 2026" in log


@pytest.mark.skipif(installed_gaussian() is None, reason="Gaussian is not installed")
def test_run_path_with_the_real_gaussian(tmp_path):
    job = run_spec(
        tmp_path,
        SPEC,
        inis={"gaussian": installed_gaussian()},
        args=["--ncores", "1"],
        source=SOURCE,
    )
    check_job(job)


B2PLYP_SPEC = """\
title: Gaussian B2PLYP-D3(BJ)
steps:
- Water: {}
- Gaussian:
    steps:
    - Energy:
        functional: "B2PLYP : double hybrid of Grimme, with MP2 correlation"
        dispersion: GD3BJ
        basis: 6-31G
"""


def test_b2plyp_d3bj_route_uses_the_builtin_keyword(tmp_path, monkeypatch):
    """B2PLYP with D3(BJ) runs Gaussian's B2PLYPD3: with EmpiricalDispersion=GD3BJ,
    Gaussian 09 E.01 uses s6 = 1.0 instead of B2PLYP's 0.64."""
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
        B2PLYP_SPEC,
        inis={
            "gaussian": {"code": code, "root-directory": "", "setup-environment": ""}
        },
        args=["--ncores", "1"],
        source=SOURCE,
    )
    route = (job / "2" / "1" / "input.dat").read_text().upper()
    assert "B2PLYPD3" in route
    assert "EMPIRICALDISPERSION" not in route


@pytest.mark.skipif(installed_gaussian() is None, reason="Gaussian is not installed")
def test_b2plyp_d3bj_with_the_real_gaussian(tmp_path):
    """Gaussian's B2PLYP-D3(BJ) dispersion is Grimme's, with B2PLYP's s6 = 0.64:
    simple-dftd3 1.7.0 gives -0.00023429 E_h for this water (and -0.0003197 with the
    s6 = 1.0 that EmpiricalDispersion=GD3BJ used)."""
    job = run_spec(
        tmp_path,
        B2PLYP_SPEC,
        inis={"gaussian": installed_gaussian()},
        args=["--ncores", "1"],
        source=SOURCE,
    )
    log = (job / "2" / "1" / "output.txt").read_text()
    assert "Normal termination" in log
    # The dispersion energy is the change in the nuclear repulsion
    before = after = None
    for x in log.splitlines():
        if "nuclear repulsion energy" in x and before is None:
            before = float(x.split()[3])
        elif "Nuclear repulsion after empirical dispersion term" in x:
            after = float(x.split("=")[1].split()[0])
    energy = after - before
    assert energy == pytest.approx(-0.00023429, abs=1.0e-8)
    assert "B2PLYP-D3(BJ)" in (job / "job.out").read_text()
