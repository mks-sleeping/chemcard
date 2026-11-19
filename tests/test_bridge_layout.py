"""Check final drawings, including the depth cue after reaction layout transforms."""

import math
import xml.etree.ElementTree as ET

import pytest
from rdkit import Chem

from chemcard.chem.bridged_layout import GAP_FRACTION, bridge_candidates
from chemcard.chem.render import _layout, layout_quality, render_smiles
from chemcard.chem.scheme import ArrowSpec, compose, resolve_arrows

CAMPHOR = "CC1(C)C2CCC1(C)C(=O)C2"
NORBORNANE = "C1CC2CCC1C2"
# PubChem CIDs 1201518 and 6973640; only the OH-bearing stereocenter differs.
BORNEOL = "C[C@]12CC[C@H](C1(C)C)C[C@H]2O"
ISOBORNEOL = "C[C@]12CC[C@H](C1(C)C)C[C@@H]2O"
SVG = "{http://www.w3.org/2000/svg}"


@pytest.mark.parametrize("smiles", [CAMPHOR, NORBORNANE, BORNEOL, ISOBORNEOL])
def test_clean_cage_view_with_one_depth_cue(smiles):
    mol = Chem.MolFromSmiles(smiles)
    assert _layout(mol).startswith("bridge-221-")
    q = layout_quality(mol)
    assert q["crossings"] == 1
    assert (q["clashes"], q["on_bond"], q["inside"], q["label_contacts"]) == (0, 0, 0, 0)
    assert sum(b.HasProp(GAP_FRACTION) for b in mol.GetBonds()) == 1
    card = render_smiles(smiles)
    assert len(ET.fromstring(card.svg).findall(f".//{SVG}mask")) == 1
    assert render_smiles(card.canonical_smiles).svg == card.svg


@pytest.mark.parametrize("smiles", ["C1CCCCCC1", "C1C2CC3CC1CC(C2)C3", "C1CC2CCC1N2",
                                    "C1CC2CCC1C2.C", "CCC1CC2CCC1C2", "C1=CC2CCC1C2"])
def test_other_ring_topologies_and_sidechains_use_existing_engines(smiles):
    assert list(bridge_candidates(Chem.MolFromSmiles(smiles))) == []
    assert "bridge-gap" not in render_smiles(smiles).svg


def test_epimers_keep_registration_and_show_opposite_oh_wedges():
    wedge_paths = []
    cards = [render_smiles(s) for s in (BORNEOL, ISOBORNEOL)]
    assert cards[0].atoms == cards[1].atoms
    for card in cards:
        mol = Chem.MolFromSmiles(card.canonical_smiles)
        oxygen = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 8)
        bond = oxygen.GetBonds()[0].GetIdx()
        paths = [p.attrib["d"] for p in ET.fromstring(card.svg).findall(f"{SVG}path")
                 if f"bond-{bond}" in p.get("class", "").split()]
        # RDKit emits polygons for a solid wedge and separate strokes for a hashed wedge.
        wedge_paths.append((any("Z" in d for d in paths), len(paths)))
    assert wedge_paths[0][0] != wedge_paths[1][0]
    assert max(count for _, count in wedge_paths) >= 3


def test_gaps_follow_each_species_into_rotated_highlighted_scheme():
    a = CAMPHOR.replace("=O", "=[O:1]")
    b = CAMPHOR.replace("=O", "=[O:2]")
    canonicals, arrows = resolve_arrows([a, b], [ArrowSpec("lp:1", "atom:2")])
    colors = []
    for smi in canonicals:
        mol = Chem.MolFromSmiles(smi)
        _layout(mol)
        rear = next(b for b in mol.GetBonds() if b.HasProp(GAP_FRACTION))
        colors.append(({}, {rear.GetIdx(): "green"}))
    sc = compose(canonicals, arrows, colors)
    root = ET.fromstring(sc.svg)
    masks = root.findall(f".//{SVG}mask")
    assert len(masks) == 2
    assert len({el.get("id") for el in masks}) == 2
    assert "mechanism-arrow" in sc.svg
    vx, vy, *_ = map(float, root.attrib["viewBox"].split())
    for species, mask in enumerate(masks):
        atoms = [a for a in sc.atoms if a["species"] == species]
        circle = mask.find(f"{SVG}circle")
        point = float(circle.get("cx")) - vx, float(circle.get("cy")) - vy
        # Independently find the two bonds whose interiors contain the cutout.
        mol = Chem.MolFromSmiles(canonicals[species])
        crossed = []
        for bond in mol.GetBonds():
            a, b = (atoms[i] for i in (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
            dx, dy = b["x"] - a["x"], b["y"] - a["y"]
            px, py = point[0] - a["x"], point[1] - a["y"]
            length = math.hypot(dx, dy)
            if abs(px * dy - py * dx) / length < 0.2 and 0 < px * dx + py * dy < length**2:
                crossed.append(bond.GetIdx())
        assert len(crossed) == 2
        paths = [p for p in root.findall(f"{SVG}path") if p.get("mask") == f"url(#{mask.get('id')})"]
        assert len(paths) >= 2  # both the rear bond's ink and its highlight are masked
