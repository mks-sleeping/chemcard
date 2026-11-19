"""Level 2: chemistry lint for single species."""

import pytest

from chemcard.chem.lint import lint_species, medium_from_conditions


def codes(smiles, **kw):
    return {(i.level, i.code) for i in lint_species(smiles, **kw)}


def test_missing_h_on_charged_atom_is_an_error_with_a_fix():
    issues = lint_species("C[N+](C)C")
    assert issues[0].level == "error" and issues[0].code == "unexpected_radical"
    assert "Did you mean C[NH+](C)C?" in issues[0].message


def test_radicals_allowed_when_declared():
    assert codes("C[CH2]", radicals_expected=True) == set()
    assert ("error", "unexpected_radical") in codes("C[CH2]")


@pytest.mark.parametrize("smiles", ["C[NH3+]", "C[NH+](C)C", "CC(=O)O", "C[N+](=O)[O-]", "c1ccccc1"])
def test_ordinary_species_are_clean(smiles):
    assert codes(smiles) == set()


def test_unstabilised_cations_warn():
    assert ("warning", "methyl_cation") in codes("[CH3+]")
    assert ("warning", "primary_cation") in codes("CC[CH2+]")


@pytest.mark.parametrize("smiles", ["C=C[CH2+]", "c1ccccc1[CH2+]", "CO[CH2+]", "CC(C)(C)[C+](C)C"])
def test_stabilised_cations_do_not_warn(smiles):
    assert not {c for _, c in codes(smiles)} & {"methyl_cation", "primary_cation"}


def test_doubly_charged_warns():
    assert ("warning", "multiply_charged") in codes("[O-]C(=O)CC(=O)[O-]")


@pytest.mark.parametrize("smiles", ["[OH-]", "C[O-]", "CC([O-])=C", "CC(=O)[CH2-]"])
def test_strong_base_in_acid(smiles):
    assert ("warning", "base_in_acid") in codes(smiles, medium="acidic")
    assert ("warning", "base_in_acid") not in codes(smiles, medium="basic")


@pytest.mark.parametrize("smiles", ["C[OH2+]", "CC(=[OH+])C"])
def test_strong_acid_in_base(smiles):
    assert ("warning", "acid_in_base") in codes(smiles, medium="basic")


def test_cyanide_is_not_a_carbanion_problem():
    assert codes("[C-]#N", medium="acidic") == set()


@pytest.mark.parametrize(
    "conditions, medium",
    [
        ("H⁺", "acidic"),
        ("H⁺, −H₂O", "acidic"),
        ("TsOH, 甲苯", "acidic"),
        ("NaOH, H₂O", "basic"),
        ("LDA, THF, −78 °C", "basic"),
        ("pH 2.5–4.5", "acidic"),
        ("pH 10", "basic"),
        ("NaOH, Na₂HPO₄ (pH 7)", None),
        ("H₂O", None),
        ("MeOH", None),
        (None, None),
    ],
)
def test_medium_from_conditions(conditions, medium):
    assert medium_from_conditions(conditions) == medium
