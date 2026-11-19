"""Acceptance tests for the chemistry core (v0.1)."""

import pytest

from chemcard.chem.render import MAX_HEAVY_ATOMS, InvalidSmiles, render_smiles

CAFFEINE = "Cn1cnc2c1c(=O)n(C)c(=O)n2C"


@pytest.mark.parametrize(
    "smiles",
    [
        CAFFEINE,
        "c1ccccc1",                # benzene
        "CC(=O)Oc1ccccc1C(=O)O",   # aspirin
        "C[C@H](N)C(=O)O",         # L-alanine (stereo)
        "F/C=C/F",                 # E-alkene
        "C[NH3+]",                 # charge
        "*O[Si](C)(C)C(C)(C)C",    # R-O-TBS (dummy atom)
    ],
)
def test_renders_valid_svg(smiles):
    s = render_smiles(smiles)
    assert s.svg.lstrip().startswith(("<svg", "<?xml"))
    assert "<script" not in s.svg.lower()
    assert s.width > 0 and s.height > 0


def test_formula():
    assert render_smiles(CAFFEINE).formula == "C8H10N4O2"
    assert render_smiles("*O[Si](C)(C)C(C)(C)C").formula == "C6H15OSiR"


def test_parse_error_explains_reason():
    with pytest.raises(InvalidSmiles, match="unclosed ring"):
        render_smiles("C1CC")


def test_dummy_atom_drawn_as_R():
    svg = render_smiles("*O[Si](C)(C)C(C)(C)C").svg
    assert "*" not in svg.split("<!-- END OF HEADER -->")[-1]


@pytest.mark.parametrize("bad", ["", "C1CC", "xyz", "C" * 600])
def test_invalid_input(bad):
    with pytest.raises(InvalidSmiles):
        render_smiles(bad)


def test_too_many_atoms():
    with pytest.raises(InvalidSmiles):
        render_smiles("C" * (MAX_HEAVY_ATOMS + 1))


def test_deterministic():
    assert render_smiles(CAFFEINE).svg == render_smiles(CAFFEINE).svg


def test_fixed_bond_length():
    # with a fixed bond length, a long chain is wider than ethanol
    assert render_smiles("CCCCCCCCCC").width > render_smiles("CCO").width


@pytest.mark.parametrize("smiles", ["[H+]", "[H-]", "[Br-]", "[Na+]"])
def test_lone_ions_are_visible(smiles):
    # CoordGen returns NaN coordinates for a lone H+, which drew a 0x0 image
    s = render_smiles(smiles)
    assert s.width > 10 and s.height > 10


@pytest.mark.parametrize("smiles", ["CN1C2CCC1CC(=O)C2", "C1CC2CCC1C2"])  # tropinone, norbornane
def test_bridged_bicyclics_render(smiles):
    s = render_smiles(smiles)
    assert s.width > 0 and s.height > 0
    assert render_smiles(smiles).svg == s.svg
