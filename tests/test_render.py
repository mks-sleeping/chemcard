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


# ---- scaffold templates (layout engines draw morphinans badly) ----------------

@pytest.mark.parametrize(
    "smiles, template",
    [
        ("COc1ccc2c(c1)[C@]13CCCC[C@@H]1[C@H](C2)N(C)CC3", "morphinan"),              # dextromethorphan
        ("CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5", "epoxymorphinan"),  # morphine
        ("O=C1CC[C@@]2(O)[C@H]3Cc4ccc(O)c5O[C@@H]1[C@]2(CCN3CC=C)c45", "epoxymorphinan"),  # naloxone
    ],
)
def test_morphinans_use_textbook_template(smiles, template):
    # CoordGen/RDKit put the N-methyl across the rings (2026-11-19 screenshot); PubChem-derived templates fix it
    from rdkit import Chem

    from chemcard.chem.render import _layout

    m = Chem.MolFromSmiles(smiles)
    assert _layout(m) == f"template:{template}"
    assert render_smiles(smiles).width > 0


BUPRENORPHINE = "CC(C)(C)[C@](C)([C@H]1C[C@@]23CC[C@]1([C@H]4[C@@]25CCN([C@@H]3CC6=C5C(=C(C=C6)O)O4)CC7CC7)OC)O"  # PubChem CID 644073


def test_template_that_adds_crossings_competes_instead_of_winning():
    # buprenorphine matches the epoxymorphinan scaffold, but its extra 6,14-bridge is not in the
    # template and got drawn across the rings (6 crossings vs CoordGen's 2): the template must not win
    from rdkit import Chem

    from chemcard.chem.render import _layout, _template_crossings, layout_quality

    m = Chem.MolFromSmiles(BUPRENORPHINE)
    assert not _layout(m).startswith("template:")
    assert layout_quality(m)["crossings"] < 6
    # the templates themselves have one crossing: the textbook bridge
    assert _template_crossings("epoxymorphinan") == 1 and _template_crossings("morphinan") == 1


@pytest.mark.parametrize("smiles", [
    "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",                 # morphine
    "O=C1CC[C@@]2(O)[C@H]3Cc4ccc(O)c5O[C@@H]1[C@]2(CCN3CC=C)c45",              # naloxone: crowded 14-OH label only
])
def test_template_still_wins_when_it_adds_no_crossing(smiles):
    from rdkit import Chem

    from chemcard.chem.render import _layout

    assert _layout(Chem.MolFromSmiles(smiles)) == "template:epoxymorphinan"


def test_ordinary_molecules_keep_coordgen():
    from rdkit import Chem

    from chemcard.chem.render import _layout

    for smi in (CAFFEINE, "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O", "C1COCCOCCOCCOCCOCCO1"):
        assert _layout(Chem.MolFromSmiles(smi)) == "coordgen"


@pytest.mark.parametrize("smiles", ["CN1C2CCC1CC(=O)C2", "COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C"])
def test_substituent_drawn_into_a_ring_loses(smiles):
    # CoordGen draws the N-methyl of tropinone/cocaine inside a ring; the score must notice and switch
    from rdkit import Chem
    from rdkit.Chem import rdCoordGen

    from chemcard.chem.render import _layout, layout_quality

    m = Chem.MolFromSmiles(smiles)
    rdCoordGen.AddCoords(m)
    assert layout_quality(m)["inside"] >= 1
    chosen = Chem.MolFromSmiles(smiles)
    assert _layout(chosen) != "coordgen"   # a template or projected perspective wins instead
    assert layout_quality(chosen)["inside"] == 0


def test_chosen_layout_is_never_clearly_worse_than_an_alternative():
    from rdkit import Chem

    from chemcard.chem.render import LAYOUT_ENGINES, LAYOUT_MIN_GAIN, _layout, layout_quality

    for smi in ("CC1(C)C2CCC1(C)C(=O)C2", "C1C2CC3CC1CC(C2)C3", "Clc1ccc(cn1)[C@@H]1C[C@H]2CC[C@@H]1N2", CAFFEINE):
        chosen = Chem.MolFromSmiles(smi)
        _layout(chosen)
        alternatives = []
        for _, engine in LAYOUT_ENGINES:
            m = Chem.MolFromSmiles(smi)
            engine(m)
            alternatives.append(layout_quality(m)["score"])
        assert layout_quality(chosen)["score"] <= min(alternatives) + LAYOUT_MIN_GAIN


def test_quality_counts_crossings_and_clashes():
    from rdkit import Chem
    from rdkit.Geometry import Point3D

    from chemcard.chem.render import layout_quality

    m = Chem.MolFromSmiles("CC.CC")  # two bonds drawn as an X
    conf = Chem.Conformer(4)
    for i, (x, y) in enumerate([(0, 0), (1, 1), (0, 1), (1, 0)]):
        conf.SetAtomPosition(i, Point3D(x, y, 0))
    m.AddConformer(conf)
    assert layout_quality(m)["crossings"] == 1


@pytest.mark.parametrize("smiles, same_side", [(r"C1/C=C\CCCCC1", True), ("C1/C=C/CCCCC1", False)])
def test_cyclooctene_geometry_on_the_final_card(smiles, same_side):
    from rdkit import Chem

    card = render_smiles(smiles)
    mol = Chem.MolFromSmiles(card.canonical_smiles)
    bond = next(b for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE)
    a, b = [card.atoms[i] for i in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())]
    signs = []
    for i in bond.GetStereoAtoms():
        p = card.atoms[i]
        signs.append((b["x"] - a["x"]) * (p["y"] - a["y"]) - (b["y"] - a["y"]) * (p["x"] - a["x"]))
    assert all(abs(v) > 1 for v in signs)
    assert (signs[0] * signs[1] > 0) == same_side


def test_rejects_invalid_geometry_even_in_the_last_resort(monkeypatch):
    from rdkit import Chem
    from rdkit.Chem import rdCoordGen
    from chemcard.chem import render

    smi = "C1/C=C/CCCCC1"

    def wrong_template(mol):
        trial = Chem.Mol(mol)
        rdCoordGen.AddCoords(trial)
        yield "template:test", trial

    monkeypatch.setattr(render, "layout_candidates", wrong_template)
    monkeypatch.setattr(render.rdDepictor, "Compute2DCoords", lambda m, **kw: rdCoordGen.AddCoords(m))
    with pytest.raises(render.UndrawableStructure, match="double-bond geometry"):
        render_smiles(smi)


def test_undrawable_is_not_reported_as_invalid_smiles(monkeypatch):
    # the SMILES is fine; the model must not be told to "fix" it
    from chemcard import server
    from chemcard.chem import render

    def refuse(mol):
        raise render.UndrawableStructure("no 2D layout reproduces the specified structure (test)")

    monkeypatch.setattr(render, "_layout", refuse)
    r = server.show_structure("C1/C=C/CCCCC1")
    assert r.is_error
    text = r.content[0].text
    assert "Invalid SMILES" not in text and "SMILES is valid" in text and "Do not change the structure" in text


# Pinned E/Z cases, checked independently: RDKit re-perceives the double-bond stereo
# from the final 2D coordinates and must get back the input. layout_problem relies on
# RDKit's stereo atoms for E/Z bonds; these catch a change in that across RDKit versions.
@pytest.mark.parametrize(
    "smiles",
    [
        "C/C=C/C", "C/C=C\\C", "OC(=O)/C=C/C(=O)O", "OC(=O)/C=C\\C(=O)O",
        "F/C(Cl)=C/F",                                   # priorities differ from the written neighbours
        "Cl/C(F)=C(\\Br)I",                             # tetrasubstituted, all halogens
        "CC/C(=C(\\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1",  # tamoxifen
        "CC1=C(/C=C/C(C)=C/C=C/C(C)=C/C=O)C(C)(C)CCC1",  # retinal: four E bonds
        "C1CCC/C=C/CC1", "C1CCC/C=C\\CC1", "C1CCCC/C=C/CCC1",  # E/Z in rings
    ],
)
def test_final_layout_reproduces_double_bond_geometry(smiles):
    from rdkit import Chem

    from chemcard.chem.render import _layout

    m = Chem.MolFromSmiles(smiles)
    _layout(m)
    back = Chem.MolFromMolBlock(Chem.MolToMolBlock(m))
    assert Chem.MolToSmiles(back) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_reference_molecule_keeps_its_substituent_coordinates():
    import math
    from rdkit import Chem
    from chemcard.chem.render import _layout, _scaffold_templates

    for name, ref, _ in _scaffold_templates():
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(ref))
        assert _layout(mol) == f"template:{name}"
        mapping = ref.GetSubstructMatch(mol, useChirality=True)
        source = ref.GetConformer().GetPositions()
        result = mol.GetConformer().GetPositions()
        assert len(mapping) == len(result)
        for i, j in enumerate(mapping):
            assert math.dist(result[i], source[j]) < 1e-5
