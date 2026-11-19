"""Level 3: whole-mechanism check (lint + balance + continuity), and the tools using it."""

import pytest

from chemcard.chem.mechanism import Step, check_mechanism
from chemcard.server import Species, mcp, show_reaction_step
from chemcard.server import check_mechanism as check_tool
from chemcard.server import MechanismStep

# Robinson tropinone synthesis, first three steps
S1 = Step(["O=CCCC=O", "CN"], ["CN1C(O)CCC1O"], "H₂O")
S2 = Step(["CN1C(O)CCC1O", "[H+]"], ["C[N+]1=CCCC1O", "O"], "H⁺")
S3 = Step(["C[N+]1=CCCC1O", "OC(=O)CC(O)=CC(=O)O"], ["C[NH+]1C(CCC1O)C(C(=O)O)C(=O)CC(=O)O"], "pH 2.5–4.5")


@pytest.fixture
def anyio_backend():
    return "asyncio"


def test_correct_mechanism_passes():
    report = check_mechanism([S1, S2, S3])
    assert report.ok and report.issues == []


def test_regioisomer_breaks_the_chain():
    wrong = Step(["C[N+]1=CC(O)CC1", "OC(=O)CC(O)=CC(=O)O"], S3.products, S3.conditions)
    report = check_mechanism([S1, S2, wrong])
    assert not report.ok
    [err] = [i for i in report.issues if i.level == "error"]
    assert err.code == "broken_chain" and err.where == "step 2 → step 3"


def test_unrelated_step_warns():
    report = check_mechanism([S1, Step(["CCO"], ["C=C", "O"], "H⁺")])
    assert report.ok
    assert "no_continuity" in {i.code for i in report.issues}


def test_stereo_mismatch_warns():
    a = Step(["C[C@H](O)CC"], ["C[C@H]([OH2+])CC"], "H⁺")
    b = Step(["C[C@@H]([OH2+])CC"], ["C[CH+]CC", "O"], "H⁺")
    assert "stereo_mismatch" in {i.code for i in check_mechanism([a, b]).issues}


def test_invalid_smiles_stops_with_position():
    report = check_mechanism([S1, Step(["C1CC"], ["CC"])])
    assert not report.ok
    assert report.issues[-1].code == "invalid_smiles"
    assert "step 2" in report.issues[-1].where and "reactants[1]" in report.issues[-1].message


def test_errors_listed_first():
    bad = Step(["C[N+](C)C", "[OH-]"], ["CN(C)C", "O"], "H⁺")
    text = check_mechanism([bad]).text()
    assert text.splitlines()[1].startswith("[error]")


def test_reaction_step_refuses_radical_typo():
    r = show_reaction_step(reactants=[Species(smiles="C[N+](C)C")], products=[Species(smiles="CN(C)C"), Species(smiles="[H+]")])
    assert r.is_error and "C[NH+](C)C" in r.content[0].text


def test_reaction_step_draws_declared_radical_step():
    r = show_reaction_step(reactants=[Species(smiles="BrBr")], products=[Species(smiles="[Br]"), Species(smiles="[Br]")], radical=True)
    assert not r.is_error


def test_check_tool_wraps_report():
    r = check_tool([MechanismStep(**vars(s)) for s in (S1, S2, S3)])
    assert r.structured_content == {"ok": True, "issues": []}


@pytest.mark.anyio
async def test_check_tool_has_no_ui():
    tools = {t.name: t for t in await mcp.list_tools()}
    t = tools["check_mechanism"]
    assert not (t.meta or {}).get("ui")
    assert t.annotations.read_only_hint is True
