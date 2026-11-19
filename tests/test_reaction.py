"""Reaction / mechanism step: A + B → C."""

import pytest

from chemcard.chem.reaction import balance_note
from chemcard.server import CARD_URI, Species, mcp, show_reaction_step


@pytest.fixture
def anyio_backend():
    return "asyncio"


def test_balanced_step_has_no_note():
    assert balance_note(["CC(C)=O", "[OH-]"], ["[CH2-]C(C)=O", "O"]) is None


def test_omitted_water_is_reported():
    # hemiaminal → iminium without H+ or H2O: net loss of OH-, so one H, one O and a charge
    note = balance_note(["CN1C(O)CCC1O"], ["C[N+]1=CCCC1O"])
    assert "reactant side has extra H O" in note
    assert "charge +0 → +1" in note


def test_r_groups_are_ignored():
    assert balance_note(["*C(=O)O", "CO"], ["*C(=O)OC", "O"]) is None


def test_step_result_layout():
    r = show_reaction_step(
        reactants=[Species(smiles="O=CCCC=O", label="丁二醛"), Species(smiles="CN")],
        products=[Species(smiles="CN1C(O)CCC1O")],
        conditions="H₂O",
        caption="第 1 步",
    )
    assert not r.is_error
    d = r.structured_content
    assert d["kind"] == "reaction"
    assert [s["label"] for s in d["reactants"]] == ["丁二醛", None]
    assert d["conditions"] == "H₂O" and d["caption"] == "第 1 步" and d["reversible"] is False
    assert d["warnings"] == []
    # SVGs only in _meta, one per species, never in what the model reads
    assert len(r.meta["reactants"]) == 2 and len(r.meta["products"]) == 1
    assert "<svg" not in str(d) and "<svg" not in r.content[0].text
    assert r.content[0].text.startswith("Displayed step: 丁二醛 + CH5N → ")


def test_imbalance_is_told_to_the_model():
    r = show_reaction_step(reactants=[Species(smiles="CN1C(O)CCC1O")], products=[Species(smiles="C[N+]1=CCCC1O")])
    assert "extra H O" in r.content[0].text and r.structured_content["warnings"]


def test_invalid_species_names_its_position():
    r = show_reaction_step(reactants=[Species(smiles="CC")], products=[Species(smiles="C1CC")])
    assert r.is_error
    assert "products[1]" in r.content[0].text


@pytest.mark.anyio
async def test_reaction_tool_bound_to_card():
    tools = {t.name: t for t in await mcp.list_tools()}
    t = tools["draw_reaction_step"]
    assert t.meta["ui"]["resourceUri"] == CARD_URI
    assert t.annotations.read_only_hint is True
    assert set(t.input_schema["required"]) == {"reactants", "products"}


def test_result_tells_the_model_to_explain_before_next_step():
    r = show_reaction_step(reactants=[Species(smiles="CC")], products=[Species(smiles="CC")])
    assert "Only after that, draw the next step." in r.content[0].text


def test_step_without_highlights_gets_a_tip():
    from chemcard.server import StepHighlight

    plain = show_reaction_step(reactants=[Species(smiles="CC=O")], products=[Species(smiles="CC=O")])
    assert "no highlights" in plain.content[0].text
    marked = show_reaction_step(
        reactants=[Species(smiles="CC=O")], products=[Species(smiles="CC=O")],
        highlights=[StepHighlight(side="reactants", index=1, smarts="C=O")],
    )
    assert "no highlights" not in marked.content[0].text
