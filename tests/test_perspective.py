"""Perspective views of any bridged bicyclic core, projected from 3D."""

import math

import pytest
from rdkit import Chem

from chemcard.chem.bridged_layout import GAP_FRACTION
from chemcard.chem.perspective import STRAIGHT_DEG, _angle, bicyclic_cores, perspective_candidates
from chemcard.chem.render import _layout, layout_problem, layout_quality, render_smiles


@pytest.mark.parametrize(
    "smiles, bridges",
    [
        ("C1CC2CCC1C2", (2, 2, 1)),            # norbornane
        ("C1=CC2CCC1C2", (2, 2, 1)),           # norbornene: a double bond in a bridge
        ("C1CC2CCC1O2", (2, 2, 1)),            # 7-oxanorbornane: heteroatom bridge
        ("C1CN2CCN1CC2", (2, 2, 2)),           # DABCO: listed with three rings by RDKit
        ("CC1=CCC2CC1C2(C)C", (3, 1, 1)),      # alpha-pinene
        ("CN1C2CCC1CC(=O)C2", (3, 2, 1)),      # tropinone
        ("C1CC2CCCC(C1)C2", (3, 3, 1)),
    ],
)
def test_bicyclic_cores_found_with_bridge_lengths(smiles, bridges):
    [core] = bicyclic_cores(Chem.MolFromSmiles(smiles))
    assert tuple(len(b) for b in core["bridges"]) == bridges


@pytest.mark.parametrize("smiles", [
    "C1CCC2CCCCC2C1",                 # decalin: fused, bridgeheads bonded
    "C1C2CC3CC1CC(C2)C3",            # adamantane: tricyclic
    "C1CCC2(CC1)CCCC2",              # spiro
    "c1ccccc1",
])
def test_non_bicyclic_systems_are_left_alone(smiles):
    m = Chem.MolFromSmiles(smiles)
    assert bicyclic_cores(m) == [] and list(perspective_candidates(m)) == []


@pytest.mark.parametrize("smiles", ["C1=CC2CCC1C2", "C1CN2CCN1CC2", "CN1C2CCC1CC(=O)C2", "C1CC2CCCC(C1)C2"])
def test_views_have_one_marked_crossing_and_no_straight_angles(smiles):
    m = Chem.MolFromSmiles(smiles)
    cands = list(perspective_candidates(m))
    assert cands
    for name, trial in cands:
        assert name.startswith("bridge-")
        gaps = [b for b in trial.GetBonds() if b.HasProp(GAP_FRACTION)]
        assert len(gaps) == 1 and 0 < gaps[0].GetDoubleProp(GAP_FRACTION) < 1
        assert layout_quality(trial)["crossings"] == 1
        pos = trial.GetConformer().GetPositions()
        for atom in trial.GetAtoms():
            ns = [n.GetIdx() for n in atom.GetNeighbors()]
            for i in range(len(ns)):
                for j in range(i + 1, len(ns)):
                    assert _angle(pos[atom.GetIdx()], pos[ns[i]], pos[ns[j]]) <= STRAIGHT_DEG


@pytest.mark.parametrize(
    "smiles, expected",
    [
        ("C1CC2CCC1C2", "bridge-221-"),           # hand-tuned [2.2.1] keeps priority
        ("C1=CC2CCC1C2", "bridge-2.2.1-"),        # outside the hand-tuned family: projected
        ("C1CN2CCN1CC2", "bridge-2.2.2-"),
        ("CN1C2CCC1CC(=O)C2", "bridge-3.2.1-"),
    ],
)
def test_layout_choice(smiles, expected):
    assert _layout(Chem.MolFromSmiles(smiles)).startswith(expected)


@pytest.mark.parametrize("smiles", [
    "CC1(C)C2CCC(C)(C2)C1=O",                                     # fenchone: crowded views lose to CoordGen
    "COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C",     # cocaine
])
def test_crowded_projections_do_not_win(smiles):
    assert not _layout(Chem.MolFromSmiles(smiles)).startswith("bridge-")


@pytest.mark.parametrize("smiles", [
    "Clc1ccc(cn1)[C@@H]1C[C@H]2CC[C@@H]1N2",       # epibatidine
    "CC(=O)C1=CCC[C@@H]2CC[C@H]1N2",               # anatoxin-a
    "C[C@@]12CC[C@@H](C1)C(C)(C)C2=O",             # fenchone, one enantiomer
])
def test_every_candidate_keeps_the_specified_stereo(smiles):
    # the 3D embedding is of this stereoisomer, so the projection cannot contradict it
    m = Chem.MolFromSmiles(smiles)
    for _, trial in perspective_candidates(m):
        assert layout_problem(m, trial) is None
        back = Chem.MolFromMolBlock(Chem.MolToMolBlock(trial))
        assert Chem.MolToSmiles(back) == Chem.MolToSmiles(m)


def test_projection_is_deterministic():
    a = render_smiles("C1=CC2CCC1C2").svg
    assert a == render_smiles("C1=CC2CCC1C2").svg
    assert "bridge-gap" in a  # the rear bond is cut at the crossing


def test_projected_gap_survives_in_a_reaction_scheme():
    # reactants are rotated to face each other; the gap is stored as a fraction along the bond
    from chemcard.chem.scheme import compose

    sc = compose(["C1=CC2CCC1C2", "BrBr"], [{"from": ["bond", 1, 0, 1, 1], "to": ["atom", 0, 0], "style": "full"}])
    assert "bridge-gap" in sc.svg and math.isfinite(sc.width)


def test_big_molecules_skip_the_3d_embedding():
    from chemcard.chem.perspective import MAX_HEAVY_ATOMS

    big = Chem.MolFromSmiles("C1CC2CCC1C2" + "C" * (MAX_HEAVY_ATOMS - 6))
    assert bicyclic_cores(big) and list(perspective_candidates(big)) == []
