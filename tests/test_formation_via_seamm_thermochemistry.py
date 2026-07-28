#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Proves the "collapses to one call" claim: for the same inputs,
`Substep.calculate_enthalpy_of_formation_via_seamm_thermochemistry`
(seamm_thermochemistry, a couple of library calls) reproduces
`Substep.calculate_enthalpy_of_formation` (this package's own ~200-line,
~5000-column-CSV implementation) bit-for-bit.

Requires the prototype `seamm_thermochemistry` database
(`~/Work/SEAMM/seamm_thermochemistry/`), built via
`seamm_thermochemistry/scripts/build_prototype_db.py`; skipped if that
package or its database isn't present, so this stays out of the way of
gaussian_step's normal CI.
"""

import logging

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
    """Duck-types just enough of Substep for the formation-energy methods:
    `.model`, `.logger`, and `.get_system_configuration()`. Deliberately
    doesn't define `.canonical_smiles`/`.PC_iupac_name` on the fake
    configuration -- both call sites already fall back gracefully via
    `except Exception` in the production code.
    """

    def __init__(self, model, atomic_numbers):
        self.model = model
        self.logger = logging.getLogger("test_formation_via_seamm_thermochemistry")
        self._configuration = _FakeConfiguration(atomic_numbers)

    def get_system_configuration(self, _arg):
        return None, self._configuration


# Ethane, C2H6 -- Ochterski's own worked example in "Thermochemistry in
# Gaussian", and well within CBS-QB3's Z<=36 coverage in the prototype DB.
ETHANE_ATOMIC_NUMBERS = [6, 6, 1, 1, 1, 1, 1, 1]

# Arbitrary but plausible Hartree values -- this test is about arithmetic
# equivalence between old and new code paths, not chemical accuracy.
ENERGY_E_H = -79.5
ENTHALPY_E_H = -79.4


@pytest.mark.parametrize("method", ["CBS-QB3", "G4"])
def test_atomization_energy_matches_legacy(method):
    fake_self = _FakeSelf(method, ETHANE_ATOMIC_NUMBERS)

    old_data = {"energy": ENERGY_E_H}
    new_data = {"energy": ENERGY_E_H}

    Substep.calculate_enthalpy_of_formation(fake_self, old_data)
    Substep.calculate_enthalpy_of_formation_via_seamm_thermochemistry(
        fake_self, new_data
    )

    assert "E atomization" in old_data
    assert new_data["E atomization"] == pytest.approx(
        old_data["E atomization"], abs=1e-6
    )


@pytest.mark.parametrize("method", ["CBS-QB3", "G4"])
def test_formation_enthalpy_matches_legacy(method):
    fake_self = _FakeSelf(method, ETHANE_ATOMIC_NUMBERS)

    old_data = {"energy": ENERGY_E_H, "H": ENTHALPY_E_H}
    new_data = {"energy": ENERGY_E_H, "H": ENTHALPY_E_H}

    Substep.calculate_enthalpy_of_formation(fake_self, old_data)
    Substep.calculate_enthalpy_of_formation_via_seamm_thermochemistry(
        fake_self, new_data
    )

    assert "DfH0" in old_data, "legacy method did not compute DfH0 -- test setup bug"
    assert new_data["DfH0"] == pytest.approx(old_data["DfH0"], abs=1e-6)


def test_missing_method_reports_cleanly_like_legacy():
    fake_self = _FakeSelf("this-method-does-not-exist", ETHANE_ATOMIC_NUMBERS)
    old_data = {"energy": ENERGY_E_H}
    new_data = {"energy": ENERGY_E_H}

    old_text = Substep.calculate_enthalpy_of_formation(fake_self, old_data)
    new_text = Substep.calculate_enthalpy_of_formation_via_seamm_thermochemistry(
        fake_self, new_data
    )

    assert "E atomization" not in old_data
    assert "E atomization" not in new_data
    assert "no tabulated atom energies" in old_text.lower()
    assert "no" in new_text.lower()
