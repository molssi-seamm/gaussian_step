#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""gaussian.ini: the shipped template and how the step reads it."""

import configparser
import importlib.resources

import pytest

from gaussian_step import substep

TEMPLATE = importlib.resources.files("gaussian_step") / "data" / "gaussian.ini"


def _step():
    """A substep without a flowchart; _gaussian_config needs no state."""
    return substep.Substep.__new__(substep.Substep)


@pytest.fixture
def no_gaussian(monkeypatch):
    monkeypatch.delenv("g16root", raising=False)
    monkeypatch.delenv("g09root", raising=False)
    monkeypatch.delenv("GAUSS_BSDDIR", raising=False)
    monkeypatch.setattr(substep.shutil, "which", lambda name: None)


def test_template_is_shipped():
    config = configparser.ConfigParser()
    config.read_string(TEMPLATE.read_text())
    assert config["local"]["installation"] == "local"
    assert "root-directory" not in config["local"]  # found at run time
    assert (TEMPLATE.parent / "configuration.txt").is_file()


def test_missing_file_gets_the_template_and_a_clear_error(tmp_path, no_gaussian):
    with pytest.raises(RuntimeError, match="SEAMM cannot find Gaussian"):
        _step()._gaussian_config("local", tmp_path)
    assert "setup-environment" in (tmp_path / "gaussian.ini").read_text()


def test_gaussian_found_from_its_environment_variable(
    tmp_path, no_gaussian, monkeypatch
):
    monkeypatch.setenv("g16root", "/opt/gaussian")
    config = _step()._gaussian_config("local", tmp_path)
    assert config["code"] == "g16"
    assert config["root-directory"] == "/opt/gaussian"
    assert config["setup-environment"] == "/opt/gaussian/g16/bsd/g16.profile"
    saved = configparser.ConfigParser()
    saved.read(tmp_path / "gaussian.ini")
    assert saved["local"]["root-directory"] == "/opt/gaussian"


def test_given_configuration_is_used_as_it_is(tmp_path, no_gaussian):
    mine = (
        "[local]\ninstallation = local\ncode = g09\n"
        "root-directory = /Users/me/Software\n"
        "setup-environment = /Users/me/Software/g09/bsd/g09.profile\n"
    )
    (tmp_path / "gaussian.ini").write_text(mine)
    config = _step()._gaussian_config("local", tmp_path)
    assert config["root-directory"] == "/Users/me/Software"
    assert (tmp_path / "gaussian.ini").read_text() == mine
