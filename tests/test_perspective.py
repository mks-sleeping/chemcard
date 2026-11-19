"""Drawings of bicyclo[a.b.c] systems: topology, perspective projections, constructed flat views."""

import math
import re

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from chemcard.chem import perspective
from chemcard.chem.bridged_layout import GAP_FRACTION
from chemcard.chem.perspective import (
    LINEAR_LABEL_DEG, MAX_CANDIDATES, MAX_PLACEMENTS, STRAIGHT_DEG, _angle, bicyclo_systems,
    fused_candidates, perspective_candidates,
)
from chemcard.chem.render import _layout, layout_problem, layout_quality, mol_from_smiles, render_smiles

DEFECTS = ("clashes", "on_bond", "inside", "label_contacts")


def card_mol(smiles):
    """The molecule in the card's own atom order, with the card's chosen layout."""
    card = render_smiles(smiles)
    m = mol_from_smiles(card.canonical_smiles)
    return card, m, _layout(m)


# ---- 1. topology, separate from the drawing ----------------------------------------------

@pytest.mark.parametrize(
    "smiles, abc",
    [
        ("C1CC2CCC1C2", (2, 2, 1)),            # norbornane
        ("C1=CC2CCC1C2", (2, 2, 1)),           # norbornene: a double bond in a bridge
        ("C1CC2CCC1O2", (2, 2, 1)),            # 7-oxanorbornane: heteroatom bridge
        ("C1CN2CCN1CC2", (2, 2, 2)),           # DABCO: listed with three rings by RDKit
        ("CC1=CCC2CC1C2(C)C", (3, 1, 1)),      # alpha-pinene
        ("CN1C2CCC1CC(=O)C2", (3, 2, 1)),      # tropinone
        ("C1CC2CCCC(C1)C2", (3, 3, 1)),
        ("C12CC(C1)C2", (1, 1, 1)),
        ("C12CCCC(CCC1)CCC2", (3, 3, 3)),
        ("C12CCCCC(CCCC1)CCCC2", (4, 4, 4)),
        ("C1CCC2CCCC(C1)CC2", (4, 3, 2)),
        ("C1CCC2CCCCC2C1", (4, 4, 0)),         # decalin: the bridgeheads share a bond
        ("c1ccc2ccccc2c1", (4, 4, 0)),         # naphthalene
        ("c1ccc2[nH]ccc2c1", (4, 3, 0)),       # indole
    ],
)
def test_bicyclo_systems_are_recognised_with_bridge_lengths(smiles, abc):
    [system] = bicyclo_systems(Chem.MolFromSmiles(smiles))
    assert system["abc"] == abc
    h1, h2 = system["heads"]
    for path in system["bridges"]:  # every path really joins the two heads
        chain = [h1] + path + [h2]
        m = Chem.MolFromSmiles(smiles)
        assert all(m.GetBondBetweenAtoms(chain[k], chain[k + 1]) for k in range(len(chain) - 1))


@pytest.mark.parametrize("smiles", [
    "C1C2CC3CC1CC(C2)C3",            # adamantane: tricyclic
    "C1CCC2(CC1)CCCC2",              # spiro
    "c1ccccc1",
    "C1CCC(CC1)C1CCCC1",             # two separate rings
])
def test_other_ring_systems_are_not_bicyclo(smiles):
    m = Chem.MolFromSmiles(smiles)
    assert bicyclo_systems(m) == []
    assert list(perspective_candidates(m)) == [] and list(fused_candidates(m)) == []


def test_decalin_is_drawn_flat_with_the_shared_bond():
    # c = 0 takes the fused route: no projection, no marked crossing, and the two rings
    # lie on opposite sides of the bridgehead bond, which is drawn once, in plain view.
    smiles = "C1CCC2CCCCC2C1"
    m = Chem.MolFromSmiles(smiles)
    assert list(perspective_candidates(m)) == []
    [(name, fused)] = list(fused_candidates(m))
    assert name == "bicyclo-4.4.0-fused1"

    card, cm, chosen = card_mol(smiles)
    for trial in (fused, cm):
        [system] = bicyclo_systems(trial)
        h1, h2 = system["heads"]
        assert trial.GetBondBetweenAtoms(h1, h2) is not None
        q = layout_quality(trial)
        assert q["crossings"] == 0 and not any(q[d] for d in DEFECTS)
        assert not any(b.HasProp(GAP_FRACTION) for b in trial.GetBonds())
        pos = trial.GetConformer().GetPositions()[:, :2]
        a, b = pos[h1], pos[h2]

        def side(i):
            return np.sign((b[0] - a[0]) * (pos[i][1] - a[1]) - (b[1] - a[1]) * (pos[i][0] - a[0]))

        left, right = ({side(i) for i in path} for path in system["bridges"][:2])
        assert len(left) == 1 and len(right) == 1 and left != right
    assert "bridge-gap" not in card.svg
    assert re.search(rf"class='bond-\d+ atom-{min(h1, h2)} atom-{max(h1, h2)}'", card.svg) or \
        re.search(rf"class='bond-\d+ atom-{max(h1, h2)} atom-{min(h1, h2)}'", card.svg)


# ---- 2. views are judged after the substituents are placed -------------------------------

def test_chloro_bicyclooctane_keeps_its_good_view():
    # The two best views of the bare [2.2.2] core put the Cl on a bond or into a ring once it
    # is placed. Judging only after placement finds the third view, which is clean.
    _, m, chosen = card_mol("ClC12CCC(CC1)CC2")
    assert chosen.startswith("bridge-2.2.2-view")
    assert int(chosen.rsplit("view", 1)[1]) > MAX_CANDIDATES  # beyond the old top-2 cut
    q = layout_quality(m)
    assert q["crossings"] == 1 and not any(q[d] for d in DEFECTS)
    assert sum(b.HasProp(GAP_FRACTION) for b in m.GetBonds()) == 1


def test_placement_budget(monkeypatch):
    calls = []
    real = perspective._place
    monkeypatch.setattr(perspective, "_place", lambda *a: calls.append(1) or real(*a))
    for smiles in ("ClC12CCC(CC1)CC2", "CC1(C)C2CCC(C)(C2)C1=O", "C12CC(C1)C2"):
        calls.clear()
        list(perspective_candidates(Chem.MolFromSmiles(smiles)))
        assert len(calls) <= 2 * MAX_PLACEMENTS  # projections, then at most as many flat views


@pytest.mark.parametrize("smiles", ["C1=CC2CCC1C2", "C1CN2CCN1CC2", "CN1C2CCC1CC(=O)C2", "C1CC2CCCC(C1)C2"])
def test_views_have_one_marked_crossing_and_no_straight_angles(smiles):
    cands = list(perspective_candidates(Chem.MolFromSmiles(smiles)))
    assert cands and all(name.startswith("bridge-") for name, _ in cands)
    for _, trial in cands:
        gaps = [b for b in trial.GetBonds() if b.HasProp(GAP_FRACTION)]
        assert len(gaps) == 1 and 0 < gaps[0].GetDoubleProp(GAP_FRACTION) < 1
        q = layout_quality(trial)
        assert q["crossings"] == 1 and not any(q[d] for d in DEFECTS)
        assert_no_straight_angles(trial)


def assert_no_straight_angles(trial):
    pos = trial.GetConformer().GetPositions()
    for atom in trial.GetAtoms():
        ns = [n.GetIdx() for n in atom.GetNeighbors()]
        limit = LINEAR_LABEL_DEG if len(ns) == 2 else STRAIGHT_DEG
        for i in range(len(ns)):
            for j in range(i + 1, len(ns)):
                assert _angle(pos[atom.GetIdx()], pos[ns[i]], pos[ns[j]]) <= limit


# ---- 3. constructed flat drawings when no projection works -------------------------------

@pytest.mark.parametrize("smiles, abc", [
    ("C12CC(C1)C2", "1.1.1"),
    ("C12CCCC(CCC1)CCC2", "3.3.3"),
    ("C12CCCCC(CCCC1)CCCC2", "4.4.4"),
])
def test_bicycles_without_a_usable_projection_get_a_constructed_view(smiles, abc):
    m = Chem.MolFromSmiles(smiles)
    cands = list(perspective_candidates(m))
    assert cands and all(name.startswith(f"bicyclo-{abc}-flat") for name, _ in cands)
    for _, trial in cands:
        q = layout_quality(trial)
        assert q["crossings"] == 0 and not any(q[d] for d in DEFECTS)
        assert_no_straight_angles(trial)
    card, cm, chosen = card_mol(smiles)
    assert chosen.startswith(f"bicyclo-{abc}-flat")
    # no carbon gets a stray "C" label: the card draws no atom labels at all
    assert not re.search(r"<path class='atom-\d+'", card.svg)


def test_existing_projection_successes_are_kept():
    for smiles, prefix in [("C1CCC2CCCC(C1)CC2", "bridge-4.3.2-view"), ("C1CN2CCN1CC2", "bridge-2.2.2-view"),
                           ("CN1C2CCC1CC(=O)C2", "bridge-3.2.1-view"), ("C1=CC2CCC1C2", "bridge-2.2.1-view"),
                           # quinine: one bridged system (quinuclidine) next to a fused one (quinoline)
                           ("COc1ccc2nccc(C(O)C3CC4CCN3CC4C=C)c2c1", "bridge-2.2.2-view")]:
        assert _layout(Chem.MolFromSmiles(smiles)).startswith(prefix)
    assert _layout(Chem.MolFromSmiles("C1CC2CCC1C2")).startswith("bridge-221-")  # hand-tuned [2.2.1] first


@pytest.mark.parametrize("smiles", [
    "CC1(C)C2CCC(C)(C2)C1=O",                                     # fenchone
    "COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C",     # cocaine
])
def test_crowded_projections_do_not_win(smiles):
    assert not _layout(Chem.MolFromSmiles(smiles)).startswith("bridge-")


# ---- 4. stereo: the final picture, checked against the known configuration ---------------

def rear_bond_from_3d(m, seeds=(7, 11)):
    """Which bond of the drawn crossing is farther from the viewer, from independent 3D conformers.

    Embed the same stereoisomer afresh (other seeds, ten conformers each) and fit the drawing's
    core coordinates as a linear projection of each. For every conformer the drawing is a
    view of, the projection's normal points at the viewer (a mirror-image drawing flips it)
    and depth along it orders the two crossing bonds. Returns the set of answers.
    """
    [system] = bicyclo_systems(m)
    core = sorted(system["atoms"])
    P = m.GetConformer().GetPositions()[core, :2]
    Pc = P - P.mean(axis=0)
    bond = np.median([np.linalg.norm(P[core.index(b.GetBeginAtomIdx())] - P[core.index(b.GetEndAtomIdx())])
                      for b in m.GetBonds() if b.GetBeginAtomIdx() in core and b.GetEndAtomIdx() in core])
    mh = Chem.AddHs(m)
    rears = set()
    for seed in seeds:
        rears |= _rears_for_seed(m, mh, core, P, Pc, bond, seed)
    return rears


def _rears_for_seed(m, mh, core, P, Pc, bond, seed):
    conf_ids = list(AllChem.EmbedMultipleConfs(mh, numConfs=10, randomSeed=seed))
    assert conf_ids
    rears = set()
    for cid in conf_ids:
        X = mh.GetConformer(cid).GetPositions()[core]
        Xc = X - X.mean(axis=0)
        A, *_ = np.linalg.lstsq(Xc, Pc, rcond=None)          # Pc ≈ Xc @ A
        if np.sqrt(((Xc @ A - Pc) ** 2).sum(axis=1).mean()) < 0.1 * bond:  # a view of this shape
            rears.add(_rear(m, core, P, Xc, A))
    return rears


def _rear(m, core, P, Xc, A):
    normal = np.cross(A[:, 0], A[:, 1])
    depth = dict(zip(core, Xc @ normal))
    pos = dict(zip(core, P))
    bonds = [b for b in m.GetBonds() if b.GetBeginAtomIdx() in core and b.GetEndAtomIdx() in core]
    crossings = []
    for i, b1 in enumerate(bonds):
        for b2 in bonds[i + 1:]:
            a, b = b1.GetBeginAtomIdx(), b1.GetEndAtomIdx()
            c, d = b2.GetBeginAtomIdx(), b2.GetEndAtomIdx()
            if len({a, b, c, d}) == 4 and perspective._segments_cross(pos[a], pos[b], pos[c], pos[d]):
                crossings.append((b1, b2))
    [(b1, b2)] = crossings

    def depth_at_crossing(bd, other):
        p, q = pos[bd.GetBeginAtomIdx()], pos[bd.GetEndAtomIdx()]
        t = perspective._intersection_param(p, q, pos[other.GetBeginAtomIdx()], pos[other.GetEndAtomIdx()])
        return depth[bd.GetBeginAtomIdx()] + t * (depth[bd.GetEndAtomIdx()] - depth[bd.GetBeginAtomIdx()])

    return b1.GetIdx() if depth_at_crossing(b1, b2) < depth_at_crossing(b2, b1) else b2.GetIdx()


# Stereocentres at the bridgeheads: the configuration fixes which bond is behind.
MUST_PROJECT = [
    ("O=C1N[C@H]2CC[C@@H]1CC2", "bridge-2.2.2-"),            # a chiral 2-azabicyclo[2.2.2]octan-3-one
    ("O=C1CC[C@@H]2CC[C@H]1C2", "bridge-3.2.1-"),            # a chiral bicyclo[3.2.1]octan-2-one
    ("OC(=O)[C@@H]1C[C@H]2C=C[C@@H]1C2", "bridge-2.2.1-"),   # endo-5-norbornene-2-carboxylic acid
]


@pytest.mark.parametrize("smiles, prefix", MUST_PROJECT)
def test_projected_view_puts_the_gap_on_the_bond_that_is_behind(smiles, prefix):
    card, m, chosen = card_mol(smiles)
    assert chosen.startswith(prefix)                    # a perspective was produced at all
    rears = rear_bond_from_3d(m)
    assert len(rears) == 1                              # every matching conformer agrees
    [rear] = rears
    gaps = [b.GetIdx() for b in m.GetBonds() if b.HasProp(GAP_FRACTION)]
    assert gaps == [rear]
    # and the final SVG cuts exactly that bond
    assert re.findall(r"<path mask='url\(#bridge-gap-(\d+)\)' class='bond-(\d+) ", card.svg) == [(str(rear),) * 2]
    assert len(re.findall(r"<mask id='bridge-gap-", card.svg)) == 1
    # supplement: the coordinates (with RDKit's wedges) still encode the input stereo
    assert_stereo_round_trips(m, smiles)


def assert_stereo_round_trips(m, smiles):
    assert layout_problem(m, m) is None
    back = Chem.MolFromMolBlock(Chem.MolToMolBlock(m))
    assert Chem.MolToSmiles(back) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


@pytest.mark.parametrize("smiles", ["O[C@@H]1CN2CCC1CC2", "O[C@H]1CN2CCC1CC2"])  # (R)/(S)-quinuclidin-3-ol
def test_stereo_off_the_bridgeheads_is_carried_by_the_wedge(smiles):
    # The only stereocentre is a bridge carbon; the bare core is achiral, and conformers of
    # this one enantiomer fit the drawing with either crossing bond behind. The depth cue
    # therefore says nothing about the configuration: the wedge at C3 does, and survives.
    card, m, chosen = card_mol(smiles)
    assert chosen.startswith("bridge-2.2.2-") and "bridge-gap" in card.svg
    assert len(rear_bond_from_3d(m)) == 2
    assert_stereo_round_trips(m, smiles)


def test_the_rear_bond_check_tells_enantiomers_apart():
    # Same drawing, mirror-image 3D shape: the other bond would have to be behind.
    _, m, _ = card_mol("O=C1N[C@H]2CC[C@@H]1CC2")
    [rear] = rear_bond_from_3d(m)
    mirror = Chem.Mol(m)
    for atom in mirror.GetAtoms():
        if atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
        elif atom.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
    mirror_rears = rear_bond_from_3d(mirror)
    assert mirror_rears and rear not in mirror_rears


@pytest.mark.parametrize("smiles", [
    "CC(=O)C1=CCC[C@@H]2CC[C@H]1N2",               # anatoxin-a: projections crowded, a flat view instead
    "C[C@@]12CC[C@@H](C1)C(C)(C)C2=O",             # fenchone, one enantiomer: crowded projections lose
    "Clc1ccc(cn1)[C@@H]1C[C@H]2CC[C@@H]1N2",       # epibatidine: its one projection has a straight angle
])
def test_expected_fallbacks_keep_the_stereo(smiles):
    card, m, chosen = card_mol(smiles)
    assert not chosen.startswith("bridge-") and "bridge-gap" not in card.svg
    for name, trial in perspective_candidates(Chem.Mol(m)):
        assert layout_problem(m, trial) is None
        if name.startswith("bridge-"):  # a projection offered only because none was clean
            q = layout_quality(trial)
            assert any(q[d] for d in DEFECTS) or q["crossings"] != 1
    assert_stereo_round_trips(m, smiles)


# ---- misc ----------------------------------------------------------------------------------

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
    assert bicyclo_systems(big) and list(perspective_candidates(big)) == []


@pytest.mark.parametrize("smiles", [
    "CC1(C)[C@@H]2CC[C@@]1(C)C(=O)C2",   # camphor: hand-tuned view wins, no 3D embedding needed
    "C1CCC2CCCCC2C1", "COc1ccc2cc(ccc2c1)[C@H](C)C(=O)O",  # fused: engines clean, construction skipped
    "ClC12CCC(CC1)CC2", "CN1C2CCC1CC(=O)C2", "C12CC(C1)C2",
])
def test_skipping_stages_never_changes_the_choice(smiles, monkeypatch):
    from chemcard.chem import render

    gated = Chem.MolFromSmiles(smiles)
    name = _layout(gated)
    full = render.layout_candidates
    monkeypatch.setattr(render, "layout_candidates", lambda mol, worth_trying=None: full(mol))
    ungated = Chem.MolFromSmiles(smiles)
    assert _layout(ungated) == name
    assert np.allclose(gated.GetConformer().GetPositions(), ungated.GetConformer().GetPositions())
