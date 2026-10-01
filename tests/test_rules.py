# -*- coding: utf-8 -*-

"""The rules the dialogs and the flowchart builder share (see seamm.Parameters):
settings that have no effect are refused with the reason, implied values are
filled in, and choices are narrowed."""

import pytest

import seamm
from seamm.builder import FlowchartBuilder, FlowchartBuildError, set_parameters

HF = "HF: Hartree-Fock self consistent field (SCF)"
PM6 = "PM6: PM6 semiempirical HF"
MP2 = "MP2: 2nd-order Møller–Plesset perturbation theory"
HSE06 = "HSE06 : hybrid functional of Heyd, Scuseria and Ernzerhof"
WB97XD = "wB97XD : hybrid functional of Chai and Head-Gordon, with dispersion"


def node(substep):
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm.gaussian", directory=".")
    result = flowchart.create_node(substep)
    flowchart.add_node(result)
    return result


def test_functional_refused_for_hf():
    energy = node("Energy")
    with pytest.raises(FlowchartBuildError, match="it applies only to DFT"):
        set_parameters(energy, method=HF, functional=HSE06)
    # The step is left as it was
    assert energy.parameters["method"].value != HF
    set_parameters(energy, method=HF)
    assert not energy.parameters.applies("integral grid")


def test_level_chooses_the_method_parameter():
    energy = node("Energy")
    with pytest.raises(FlowchartBuildError, match="'level' is not 'recommended'"):
        set_parameters(energy, advanced_method=MP2)
    set_parameters(energy, level="advanced", advanced_method=MP2, freeze_cores="no")
    with pytest.raises(FlowchartBuildError, match="'level' is 'recommended'"):
        set_parameters(energy, method=HF)


def test_basis_and_frozen_core_follow_the_method():
    energy = node("Energy")
    with pytest.raises(FlowchartBuildError, match="does not use a basis set"):
        set_parameters(energy, method=PM6, basis="cc-pVDZ")
    with pytest.raises(FlowchartBuildError, match="has no frozen-core option"):
        set_parameters(energy, method=HF, freeze_cores="no")
    set_parameters(energy, method=MP2, freeze_cores="no", basis="cc-pVDZ")


def test_dispersion_narrowed_by_the_functional():
    energy = node("Energy")
    P = energy.parameters
    set_parameters(energy, functional=HSE06, dispersion="PFD")
    assert P.choices("dispersion") == ("none", "PFD")
    with pytest.raises(FlowchartBuildError, match="choose one of: 'none', 'PFD'"):
        set_parameters(energy, dispersion="GD3BJ")
    with pytest.raises(FlowchartBuildError, match="has no dispersion correction"):
        set_parameters(energy, functional=WB97XD, dispersion="none")
    # A functional given by a variable might offer any of them
    set_parameters(energy, functional="$functional", dispersion="GD3BJ")


def test_sub_controls():
    energy = node("Energy")
    with pytest.raises(FlowchartBuildError, match="'save basis set' is not 'no'"):
        set_parameters(energy, basis_set_file="my.gbs")
    set_parameters(energy, save_basis_set="yes", basis_set_file="my.gbs")
    with pytest.raises(FlowchartBuildError, match="'bond orders' is not 'none'"):
        set_parameters(energy, bond_orders="none", apply_bond_orders="no")
    with pytest.raises(FlowchartBuildError, match="'region' is 'explicit'"):
        set_parameters(energy, nx=20)
    set_parameters(energy, region="explicit", nx=20)
    with pytest.raises(FlowchartBuildError, match="'input only' is not 'yes'"):
        set_parameters(energy, input_only="yes", file_handling="keep all")


def test_optimization():
    opt = node("Optimization")
    with pytest.raises(FlowchartBuildError, match="always calculates the gradient"):
        set_parameters(opt, calculate_gradient="no")
    with pytest.raises(FlowchartBuildError, match="'target' is not"):
        set_parameters(opt, saddle_order=3)
    set_parameters(opt, target="saddle point", saddle_order=3)
    with pytest.raises(FlowchartBuildError, match="'hessian' is 'calculate'"):
        set_parameters(opt, recalc_hessian="every step")
    set_parameters(opt, hessian="calculate", recalc_hessian="every step")


def test_thermodynamics_optimize_first():
    thermo = node("Thermodynamics")
    with pytest.raises(FlowchartBuildError, match="'optimize first' is 'yes'"):
        set_parameters(thermo, optimize_first="no", target="transition state")
    set_parameters(thermo, target="transition state")
    with pytest.raises(FlowchartBuildError, match="does not use it yet"):
        set_parameters(thermo, T=400)


def test_stability_method_is_implied_for_the_advanced_level():
    stability = node("Wavefunction Stability")
    P = stability.parameters
    set_parameters(stability, level="advanced", method=HF)
    # The calculation reads 'advanced_method' at the advanced level
    assert P["advanced_method"].value == HF
    with pytest.raises(FlowchartBuildError, match="always takes the method"):
        set_parameters(stability, advanced_method=MP2)
    with pytest.raises(FlowchartBuildError, match="dispersion"):
        set_parameters(stability, dispersion="GD3BJ")


def test_flowchart_builder():
    fb = FlowchartBuilder(title="Gaussian rules")
    gaussian = fb.add("Gaussian")
    gaussian.add("Energy", method=HF, basis="cc-pVDZ")
    with pytest.raises(FlowchartBuildError, match="it applies only to DFT"):
        gaussian.add("Optimization", method=PM6, integral_grid="SuperFine")


@pytest.mark.parametrize("substep", ("Thermodynamics", "Wavefunction Stability"))
def test_rules_survive_from_dict(substep):
    """The builder restores a step with from_dict() after refusing a setting, and
    flowcharts are read with it, so the rules must be in the definitions."""
    P = node(substep).parameters
    expected = P.applicable({**P.current_values(), "optimize first": "no"})
    P.from_dict(P.to_dict())
    assert P.applicable({**P.current_values(), "optimize first": "no"}) == expected
    if substep == "Wavefunction Stability":
        assert len(P["method"].enumeration) == 2
        assert P.applies("method", {**P.current_values(), "level": "advanced"})
