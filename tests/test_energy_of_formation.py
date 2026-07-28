#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for Substep.calculate_energy_of_formation, the production
formation-energy method (0 K, ZPE-free, via seamm_thermochemistry) that
`energy.py`'s `analyze()` now calls in place of the old
`calculate_enthalpy_of_formation` (298 K, required a frequency calc).

A missing/broken formation-energy dependency must never break a Gaussian
job, so most of this suite is about the graceful-degradation paths, not
just the happy path.
"""

import logging
import sys

import pytest

from gaussian_step.substep import Substep

seamm_thermochemistry = pytest.importorskip("seamm_thermochemistry")

pytestmark = pytest.mark.skipif(
    not seamm_thermochemistry.DEFAULT_DB_PATH.exists(),
    reason=f"{seamm_thermochemistry.DEFAULT_DB_PATH} not built -- run "
    "seamm_thermochemistry/scripts/build_prototype_db.py",
)


class _FakeAtoms:
    def __init__(self, atomic_numbers):
        self.atomic_numbers = atomic_numbers


class _FakeConfiguration:
    def __init__(self, atomic_numbers):
        self.atoms = _FakeAtoms(atomic_numbers)


class _FakeSelf:
    def __init__(self, model, atomic_numbers):
        self.model = model
        self.logger = logging.getLogger("test_energy_of_formation")
        self._configuration = _FakeConfiguration(atomic_numbers)

    def get_system_configuration(self, _arg):
        return None, self._configuration


ETHANE_ATOMIC_NUMBERS = [6, 6, 1, 1, 1, 1, 1, 1]
ENERGY_E_H = -79.5


def test_energy_of_formation_from_energy_alone():
    # No "H", no "ZPE" -- exactly the case the old method couldn't handle.
    fake_self = _FakeSelf("CBS-QB3", ETHANE_ATOMIC_NUMBERS)
    data = {"energy": ENERGY_E_H}

    text = Substep.calculate_energy_of_formation(fake_self, data)

    assert "E atomization" in data
    assert "DfE" in data
    assert "DfH0" not in data  # no ZPE available
    assert "Energy of Formation" in text
    assert f"{data['DfE']:.2f}" in text


def test_energy_of_formation_adds_zpe_when_available():
    fake_self = _FakeSelf("CBS-QB3", ETHANE_ATOMIC_NUMBERS)
    zpe_hartree = 0.075
    data = {"energy": ENERGY_E_H, "ZPE": zpe_hartree}

    Substep.calculate_energy_of_formation(fake_self, data)

    assert "DfH0" in data
    from seamm_util import Q_

    zpe_kJ = Q_(zpe_hartree, "E_h").m_as("kJ/mol")
    assert data["DfH0"] == pytest.approx(data["DfE"] + zpe_kJ)


def test_unknown_method_reports_cleanly_no_exception():
    fake_self = _FakeSelf("this-method-does-not-exist", ETHANE_ATOMIC_NUMBERS)
    data = {"energy": ENERGY_E_H}

    text = Substep.calculate_energy_of_formation(fake_self, data)

    assert "E atomization" not in data
    assert "DfE" not in data
    assert "cannot calculate" in text.lower()


def test_element_without_0K_anchor_reports_cleanly():
    # Uranium: in the VASP table but has no experimental DfH0 at all (0K
    # or 298K) in the master reference sheet -- and isn't a gaussian_step
    # composition anyone would run, but exercises the same failure path a
    # heavier organic element would hit before its 0K anchor is populated.
    fake_self = _FakeSelf("CBS-QB3", [92])
    data = {"energy": ENERGY_E_H}

    text = Substep.calculate_energy_of_formation(fake_self, data)

    assert "DfE" not in data
    assert "cannot calculate" in text.lower()


def test_graceful_when_seamm_thermochemistry_not_installed(monkeypatch):
    # The well-known trick: sys.modules[name] = None makes `import name`
    # raise ImportError without actually uninstalling anything.
    monkeypatch.setitem(sys.modules, "seamm_thermochemistry", None)

    fake_self = _FakeSelf("CBS-QB3", ETHANE_ATOMIC_NUMBERS)
    data = {"energy": ENERGY_E_H}

    text = Substep.calculate_energy_of_formation(fake_self, data)

    assert "DfE" not in data
    assert "not installed" in text


def test_graceful_when_database_not_built(monkeypatch, tmp_path):
    monkeypatch.setattr(
        seamm_thermochemistry.db, "DEFAULT_DB_PATH", tmp_path / "does_not_exist.db"
    )

    fake_self = _FakeSelf("CBS-QB3", ETHANE_ATOMIC_NUMBERS)
    data = {"energy": ENERGY_E_H}

    text = Substep.calculate_energy_of_formation(fake_self, data)

    assert "DfE" not in data
    assert "not built" in text


def test_result_differs_from_legacy_298K_method_as_expected():
    """Confirms the acknowledged, intentional behavior change: switching
    the production path from the 298 K enthalpy convention to the 0 K
    energy convention changes the reported number for the same molecule.
    """
    fake_self = _FakeSelf("CBS-QB3", ETHANE_ATOMIC_NUMBERS)

    old_data = {"energy": ENERGY_E_H, "H": -79.4}
    Substep.calculate_enthalpy_of_formation(fake_self, old_data)

    new_data = {"energy": ENERGY_E_H}
    Substep.calculate_energy_of_formation(fake_self, new_data)

    assert old_data["DfH0"] != pytest.approx(new_data["DfE"])
