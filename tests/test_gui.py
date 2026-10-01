# -*- coding: utf-8 -*-

"""Smoke test of the Tk dialogs: create them and re-lay them out for every choice
that drives the layout, checking that the controls shown are exactly those that
the parameters' rules say apply. Skipped when no display is available."""

import pytest

SUBSTEPS = ("Energy", "Optimization", "Thermodynamics", "Wavefunction Stability")

# The controls in the results tab, which is laid out by seamm itself
NOT_LAID_OUT = ("results", "create tables")


@pytest.fixture()
def root():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw

    Pmw.initialise(root)
    yield root
    root.destroy()


def make(root, substep):
    import seamm

    flowchart = seamm.Flowchart(namespace="org.molssi.seamm.gaussian", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.gaussian.tk"
    )
    node = flowchart.create_node(substep)
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get(substep)
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    tk_node.create_dialog()
    # The results tab is rebuilt on every layout and is slow; it is not under test.
    tk_node.setup_results = lambda *args, **kwargs: None
    return tk_node


def shown(tk_node, key):
    """Whether a control is shown: gridded, in frames that are themselves shown."""
    top = (tk_node["frame"], tk_node["output frame"])
    widget = tk_node[key]
    while widget not in top:
        if widget.winfo_manager() == "":
            return False
        widget = widget.master
    return True


def check(tk_node):
    """The shown controls are exactly those that apply."""
    tk_node.reset_dialog()
    tk_node.reset_plotting()
    P = tk_node.node.parameters
    values = tk_node._widget_values()
    for key in P:
        if key in NOT_LAID_OUT or key not in tk_node:
            continue
        if shown(tk_node, key):
            assert P.applies(key, values), f"{key} is shown but does not apply"
        else:
            assert not P.applies(key, values), f"{key} applies but is not shown"
    return values


def set_and_check(tk_node, key, value):
    tk_node[key].set(value)
    return check(tk_node)


@pytest.mark.parametrize("substep", SUBSTEPS)
def test_layouts_follow_the_rules(root, substep):
    tk_node = make(root, substep)
    P = tk_node.node.parameters
    check(tk_node)

    for value in ("yes", "no"):
        set_and_check(tk_node, "input only", value)
        assert shown(tk_node, "file handling") == (value == "no")

    for level in P["level"].enumeration:
        set_and_check(tk_node, "level", level)
        method_key = P.method_parameter(tk_node._widget_values())
        functional_key = P.functional_parameter(tk_node._widget_values())
        for method in P[method_key].enumeration:
            values = set_and_check(tk_node, method_key, method)
            if substep == "Wavefunction Stability":
                # The calculation reads 'advanced_method' at the advanced level
                assert values["advanced_method"] == method
            if "DFT" in method:
                for functional in P[functional_key].enumeration:
                    values = set_and_check(tk_node, functional_key, functional)
                    # Only the dispersion corrections the functional offers
                    if shown(tk_node, "dispersion"):
                        offered = tk_node["dispersion"].combobox.cget("values")
                        assert tuple(offered) == P.choices("dispersion", values)
                        assert values["dispersion"] in offered
        # A method given by a variable: anything might apply
        values = set_and_check(tk_node, method_key, "$method")
        assert shown(tk_node, "basis")
        set_and_check(tk_node, method_key, P[method_key].default)
    set_and_check(tk_node, "level", "recommended")

    for key in ("bond orders", "save basis set", "orbitals", "region"):
        if not P.applies(key, tk_node._widget_values()):
            continue
        for value in P[key].enumeration:
            set_and_check(tk_node, key, value)
        set_and_check(tk_node, key, P[key].default)

    if "optimize first" in P:
        for value in ("no", "yes"):
            set_and_check(tk_node, "optimize first", value)
            assert shown(tk_node, "target") == (value == "yes")

    if "target" in P:
        for key in ("target", "hessian"):
            for value in P[key].enumeration:
                set_and_check(tk_node, key, value)
            set_and_check(tk_node, key, P[key].default)


def test_dispersion_follows_the_functional(root):
    """Changing to a functional without the current dispersion correction picks
    one it has; one with no choice hides the control."""
    tk_node = make(root, "Energy")
    tk_node["dispersion"].set("GD3BJ")
    values = set_and_check(
        tk_node,
        "functional",
        "HSE06 : hybrid functional of Heyd, Scuseria and Ernzerhof",
    )
    assert values["dispersion"] == "PFD"
    set_and_check(
        tk_node,
        "functional",
        "wB97XD : hybrid functional of Chai and Head-Gordon, with dispersion",
    )
    assert not shown(tk_node, "dispersion")


def test_method_keyword_is_written_out(root):
    tk_node = make(root, "Energy")
    set_and_check(tk_node, "method", "HF")
    assert tk_node["method"].get() == "HF: Hartree-Fock self consistent field (SCF)"
    assert tk_node.node.method == "HF"
    assert not shown(tk_node, "functional")
