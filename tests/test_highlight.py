"""Highlights (model → figure) and selections (user → model)."""

import pytest
from rdkit import Chem

from chemcard.chem.highlight import HighlightSpec as H
from chemcard.chem.highlight import InvalidHighlight, describe_selection, resolve
from chemcard.chem.render import render_smiles
from chemcard.server import Highlight, Species, StepHighlight, mcp, render_structure_svg, show_reaction_step, show_structure
from chemcard.server import describe_selection_tool

ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"


@pytest.fixture
def anyio_backend():
    return "asyncio"


def elements(smiles, idx):
    return sorted(Chem.MolFromSmiles(smiles).GetAtomWithIdx(i).GetSymbol() for i in idx)


def test_smarts_highlights_atoms_and_bonds_in_canonical_indices():
    r = resolve(ASPIRIN, [H("red", smarts="C(=O)[OX2H1]")])
    assert r.canonical_smiles == ASPIRIN
    assert elements(r.canonical_smiles, r.atom_colors) == ["C", "O", "O"]
    assert set(r.atom_colors.values()) == {"red"} and len(r.bond_colors) == 2


def test_atom_maps_survive_canonicalisation():
    # the map sits on N; the canonical SMILES puts N first, so indices must be remapped
    r = resolve("OCC[NH2:7]", [H("blue", atoms=(7,))])
    assert r.canonical_smiles == "NCCO"
    assert r.atom_colors == {0: "blue"}


def test_later_highlight_wins_on_overlap():
    r = resolve(ASPIRIN, [H("red", smarts="C=O"), H("blue", smarts="C(=O)[OX2H1]")])
    assert "red" in r.atom_colors.values() and "blue" in r.atom_colors.values()


@pytest.mark.parametrize(
    "spec, message",
    [
        (H("red", smarts="C(=O"), "invalid SMARTS"),
        (H("red", smarts="[N+]"), "matches nothing"),
        (H("red", atoms=(3,)), "not found"),
        (H("red"), "exactly one"),
        (H("pink", smarts="C=O"), "unknown colour"),
    ],
)
def test_bad_highlights_explain_themselves(spec, message):
    with pytest.raises(InvalidHighlight, match=message):
        resolve(ASPIRIN, [spec])


def test_drawing_reports_one_click_target_per_atom():
    s = render_smiles(ASPIRIN)
    assert len(s.atoms) == Chem.MolFromSmiles(ASPIRIN).GetNumAtoms()
    assert [a["el"] for a in s.atoms[:3]] == ["C", "C", "O"]
    assert all(0 <= a["x"] <= s.width and 0 <= a["y"] <= s.height for a in s.atoms)


def test_highlight_changes_the_drawing():
    r = resolve(ASPIRIN, [H("red", smarts="C=O")])
    plain = render_smiles(ASPIRIN).svg
    marked = render_smiles(r.canonical_smiles, r.atom_colors, r.bond_colors).svg
    assert marked != plain


def test_describe_selection():
    d = describe_selection(ASPIRIN, [2, 1])  # the ester C=O, clicked O first
    assert d["fragment_smiles"] == "C=O"
    assert d["mapped_smiles"] == "C[C:2](Oc1ccccc1C(=O)O)=[O:1]"
    assert [(a["element"], a["map"]) for a in d["atoms"]] == [("O", 1), ("C", 2)]


def test_selection_numbers_can_be_fed_back_as_highlight():
    d = describe_selection(ASPIRIN, [2, 1])
    r = resolve(d["mapped_smiles"], [H("green", atoms=(1, 2))])
    assert r.canonical_smiles == ASPIRIN
    assert set(r.atom_colors) == {1, 2}


# ---- through the tools ----------------------------------------------------

def test_show_structure_reports_highlight_and_click_targets():
    r = show_structure(ASPIRIN, highlights=[Highlight(smarts="C=O", color="orange")])
    assert not r.is_error
    marks = r.structured_content["highlight"]
    assert set(marks["atoms"].values()) == {"orange"}
    assert len(r.meta["atoms"]) == 13


def test_bad_highlight_is_a_tool_error():
    r = show_structure(ASPIRIN, highlights=[Highlight(smarts="[N+]")])
    assert r.is_error and "matches nothing" in r.content[0].text


def test_fallback_redraws_identically_with_highlight():
    r = show_structure(ASPIRIN, highlights=[Highlight(smarts="C=O", color="red")])
    again = render_structure_svg(r.structured_content["canonical_smiles"], r.structured_content["highlight"])
    assert again.structured_content["svg"] == r.meta["svg"]


def test_reaction_step_highlight_targets_one_species():
    r = show_reaction_step(
        reactants=[Species(smiles="C[NH2:1]"), Species(smiles="CC=O")],
        products=[Species(smiles="CC(O)NC")],
        highlights=[
            StepHighlight(side="reactants", index=1, atoms=[1], color="blue"),
            StepHighlight(side="reactants", index=2, smarts="C=O", color="red"),
        ],
    )
    assert not r.is_error
    assert r.structured_content["reactants"][0]["highlight"]["atoms"] == {"1": "blue"}
    assert set(r.structured_content["reactants"][1]["highlight"]["atoms"].values()) == {"red"}
    assert "highlight" not in r.structured_content["products"][0]


def test_reaction_step_highlight_out_of_range():
    r = show_reaction_step(
        reactants=[Species(smiles="CC")], products=[Species(smiles="CC")],
        highlights=[StepHighlight(side="products", index=2, smarts="C")],
    )
    assert r.is_error and "only 1 products" in r.content[0].text


def test_describe_selection_tool():
    r = describe_selection_tool(ASPIRIN, [2, 1])
    assert r.structured_content["fragment_smiles"] == "C=O"
    assert describe_selection_tool(ASPIRIN, [99]).is_error


@pytest.mark.anyio
async def test_card_only_tools_are_hidden_from_model():
    tools = {t.name: t for t in await mcp.list_tools()}
    assert tools["describe_selection"].meta["ui"]["visibility"] == ["app"]
    assert tools["render_structure_svg"].meta["ui"]["visibility"] == ["app"]
