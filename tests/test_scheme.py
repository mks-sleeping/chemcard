"""Curved arrows: resolving the model's specs and drawing the reactants together."""

import math
import re

import pytest

from chemcard.chem.mechanism import Step, lint_step
from chemcard.chem.scheme import ArrowSpec as A
from chemcard.chem.scheme import InvalidArrow, compose, resolve_arrows
from chemcard.server import Arrow, Species, StepHighlight, render_reactant_scheme, show_reaction_step

SN2 = (["[OH-:1]", "[CH3:2][Br:3]"], [A("lp:1", "bond:1-2"), A("bond:2-3", "atom:3")])
E2 = (["CC[O-:1]", "C[CH:2]([H:5])[CH2:3][Br:4]"],
      [A("lp:1", "atom:5"), A("bond:2-5", "bond:2-3"), A("bond:3-4", "atom:4")])


# ---- resolving ---------------------------------------------------------------

def test_sn2_resolves_to_species_and_canonical_atoms():
    canonicals, arrows = resolve_arrows(*SN2)
    assert canonicals == ["[OH-]", "CBr"]
    assert arrows[0] == {"from": ["lp", 0, 0], "to": ["bond", 0, 0, 1, 0], "style": "full"}
    assert arrows[1] == {"from": ["bond", 1, 0, 1, 1], "to": ["atom", 1, 1], "style": "full"}


def test_explicit_h_is_kept_for_proton_transfer():
    canonicals, arrows = resolve_arrows(*E2)
    assert "[H]" in canonicals[1]
    assert arrows[0]["to"][0] == "atom"


@pytest.mark.parametrize(
    "smiles, specs, message",
    [
        (["[OH-:1]", "[CH3:1]Br"], [A("lp:1", "atom:1")], "more than one reactant"),
        (["[OH-:1]", "CBr"], [A("lp:1", "atom:9")], "not found"),
        (["[OH-:1]", "[CH3:2]Br"], [A("lp:2", "atom:1")], "no lone pair"),
        (["[OH-:1]", "[CH3:2][Br:3]"], [A("bond:1-2", "atom:3")], "not bonded"),
        (["[OH-:1]", "[CH3:2]Br"], [A("pair:1", "atom:2")], "cannot read"),
        (["[OH-:1]", "[CH3:2]Br"], [A("lp:1", "atom:2", "double")], "style"),
        (["[OH-:1]", "[CH3:2]Br"], [A("lp:1", "atom:2")] * 9, "at most"),
    ],
)
def test_bad_arrows_explain_themselves(smiles, specs, message):
    with pytest.raises(InvalidArrow, match=message):
        resolve_arrows(smiles, specs)


def test_carbanion_has_a_lone_pair():
    resolve_arrows(["N#[C-:1]", "C[C:2](C)=O"], [A("lp:1", "atom:2")])


# ---- drawing -----------------------------------------------------------------

@pytest.mark.parametrize("case", [SN2, E2, (["[Br:1][Br:2]"], [A("bond:1-2", "atom:1", "half"), A("bond:1-2", "atom:2", "half")])])
def test_compose_draws_every_arrow(case):
    canonicals, arrows = resolve_arrows(*case)
    sc = compose(canonicals, arrows)
    assert sc.svg.count("class='mechanism-arrow'") == len(arrows)
    assert sc.svg.count(">+</text>") == len(canonicals) - 1
    assert not re.search(r"\bnan\b", sc.svg, re.IGNORECASE)
    assert sc.width > 0 and sc.height > 0
    # the canvas grows to hold the arrows: every path point is inside the viewBox
    vx, vy, vw, vh = (float(v) for v in re.search(r"viewBox='([^']*)'", sc.svg).group(1).split())
    for x, y in re.findall(r"(-?\d+\.\d) (-?\d+\.\d)", " ".join(re.findall(r"d='([^']*)'", sc.svg.split("mechanism-arrow", 1)[1]))):
        assert vx - 1 <= float(x) <= vx + vw + 1 and vy - 1 <= float(y) <= vy + vh + 1


def test_click_targets_know_their_molecule():
    canonicals, arrows = resolve_arrows(*SN2)
    sc = compose(canonicals, arrows)
    assert [(a["species"], a["index"], a["el"]) for a in sc.atoms] == [(0, 0, "O"), (1, 0, "C"), (1, 1, "Br")]
    assert len(sc.species) == 2 and sc.species[0]["cx"] < sc.species[1]["cx"]


def test_reacting_atoms_face_each_other():
    # OH- on the left attacks C: C must sit on the left of CH3Br (backside, away from Br)
    canonicals, arrows = resolve_arrows(*SN2)
    sc = compose(canonicals, arrows)
    c, br = sc.atoms[1], sc.atoms[2]
    assert c["x"] < br["x"]


def test_compose_is_deterministic():
    canonicals, arrows = resolve_arrows(*E2)
    assert compose(canonicals, arrows).svg == compose(canonicals, arrows).svg


def test_explicit_h_keeps_step_balanced():
    step = Step(["CC[O-]", "C[CH]([H])CBr"], ["CC=C", "CCO", "[Br-]"], "NaOEt")
    assert not [i for i in lint_step(step) if i.code == "unbalanced"]


# ---- through the tool ----------------------------------------------------------

def sn2_step(**kw):
    return show_reaction_step(
        reactants=[Species(smiles="[OH-:1]"), Species(smiles="[CH3:2][Br:3]")],
        products=[Species(smiles="CO"), Species(smiles="[Br-]")],
        arrows=[Arrow(start="lp:1", end="bond:1-2"), Arrow(start="bond:2-3", end="atom:3")],
        **kw,
    )


def test_tool_returns_arrows_and_combined_drawing():
    r = sn2_step()
    assert not r.is_error
    assert r.structured_content["arrows"][0]["from"] == ["lp", 0, 0]
    assert "mechanism-arrow" in r.meta["reactant_scheme"]["svg"]
    assert "mechanism-arrow" not in str(r.structured_content)  # SVG stays out of the model's view
    assert "Curved arrows: lp:1 → bond:1-2" in r.content[0].text


def test_fallback_redraws_identically_with_highlight():
    r = sn2_step(highlights=[StepHighlight(side="reactants", index=2, smarts="CBr", color="red")])
    d = r.structured_content
    again = render_reactant_scheme(
        [{"smiles": s["canonical_smiles"], "highlight": s.get("highlight")} for s in d["reactants"]],
        d["arrows"],
    )
    assert not again.is_error
    assert again.structured_content["svg"] == r.meta["reactant_scheme"]["svg"]


def test_bad_arrow_is_a_tool_error():
    r = show_reaction_step(
        reactants=[Species(smiles="[OH-:1]"), Species(smiles="CBr")],
        products=[Species(smiles="CO"), Species(smiles="[Br-]")],
        arrows=[Arrow(start="lp:1", end="atom:7")],
    )
    assert r.is_error and "arrow problem" in r.content[0].text and "not found" in r.content[0].text


def test_steps_without_arrows_are_unchanged():
    r = show_reaction_step(reactants=[Species(smiles="CC=O")], products=[Species(smiles="CC=O")])
    assert "arrows" not in r.structured_content and "reactant_scheme" not in r.meta
