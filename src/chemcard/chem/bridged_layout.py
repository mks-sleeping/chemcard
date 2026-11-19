"""Oblique views of a saturated bicyclo[2.2.1]heptane ring system."""

from __future__ import annotations

import math
import re
from collections.abc import Iterator
from itertools import permutations

from rdkit import Chem

# Fraction measured along the rear bond, relative to its RDKit begin atom.
# Bond properties survive CombineMols, unlike molecule-local bond indices.
GAP_FRACTION = "_chemcardBridgeGap"


def bridge_candidates(mol: Chem.Mol) -> Iterator[tuple[str, Chem.Mol]]:
    """Offer four assignments of the same seven ring positions to the molecular graph.

    The first family supports an isolated carbon cage with terminal substituents.
    More complex side chains, fused rings and mixtures use the existing engines.
    """
    ring = {a.GetIdx() for a in mol.GetAtoms() if a.IsInRing()}
    if len(ring) != 7:
        return
    for atom in mol.GetAtoms():
        if atom.GetIdx() in ring:
            if atom.GetAtomicNum() != 6:
                return
        elif atom.GetDegree() != 1 or atom.GetNeighbors()[0].GetIdx() not in ring:
            return
    neighbors = {i: [] for i in ring}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if a in ring and b in ring:
            if bond.GetBondType() != Chem.BondType.SINGLE:
                return
            neighbors[a].append(b)
            neighbors[b].append(a)
    heads = sorted(i for i, adj in neighbors.items() if len(adj) == 3)
    if len(heads) != 2 or any(len(adj) != 2 for i, adj in neighbors.items() if i not in heads):
        return
    view = 0
    for front, back in permutations(heads):
        paths = []
        for first in neighbors[front]:
            prev, curr = front, first
            path = []
            while curr != back:
                if curr == front or curr in path:
                    return
                path.append(curr)
                prev, curr = curr, next(i for i in neighbors[curr] if i != prev)
            paths.append(path)
        paths.sort(key=len)
        if [len(path) for path in paths] != [1, 2, 2]:
            return
        apex = paths[0][0]
        for left, right in permutations(paths[1:]):
            view += 1
            trial = _view(mol, front, back, apex, left, right, neighbors)
            if trial is not None:
                yield f"bridge-221-{view}", trial


def _view(mol, front, back, apex, left, right, neighbors):
    # The two long paths outline the deck. The single-carbon bridge rises above it.
    # Depth is used only to check whether this view agrees with specified bridgeheads.
    points = {front: (-0.22, -0.62, 1), back: (0.44, 0.52, -1), apex: (-0.54, 1.68, 0),
              left[0]: (-1.62, -0.92, 1), left[1]: (-1.24, 0.24, -1),
              right[0]: (1.22, -0.96, 1), right[1]: (1.85, 0.13, -1)}
    for i, ring_neighbors in neighbors.items():
        x, y, z = points[i]
        branches = [a.GetIdx() for a in mol.GetAtomWithIdx(i).GetNeighbors() if a.GetIdx() not in neighbors]
        if len(branches) == 2 and i == apex:
            for j, (dx, dy) in zip(branches, [(-0.92, 0.43), (0.19, 1.04)]):
                points[j] = (x + dx, y + dy, z)
        elif len(branches) == 1:
            vx, vy = 0.0, 0.0
            for j in ring_neighbors:
                dx, dy = x - points[j][0], y - points[j][1]
                length = math.hypot(dx, dy)
                vx += dx / length
                vy += dy / length
            scale = 1.18 / math.hypot(vx, vy)
            points[branches[0]] = (x + vx * scale, y + vy * scale, 2 * z)
        elif branches:
            return None

    trial = Chem.Mol(mol)
    trial.RemoveAllConformers()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, xyz in points.items():
        conf.SetAtomPosition(i, xyz)
    conf.Set3D(True)
    trial.AddConformer(conf)
    oriented = Chem.Mol(trial)
    Chem.AssignAtomChiralTagsFromStructure(oriented, replaceExistingTags=True)
    for i in (front, back):
        tag = mol.GetAtomWithIdx(i).GetChiralTag()
        if tag != Chem.ChiralType.CHI_UNSPECIFIED and tag != oriented.GetAtomWithIdx(i).GetChiralTag():
            return None
    conf = trial.GetConformer()
    conf.Set3D(False)
    for i, (x, y, _) in points.items():
        conf.SetAtomPosition(i, (x, y, 0.0))

    # Intersect the rear edge and the front bridge in this projection. Store the
    # affine-invariant fraction so later rotations/translations require no remapping.
    rear = trial.GetBondBetweenAtoms(back, left[1])
    a, b = (points[j] for j in (rear.GetBeginAtomIdx(), rear.GetEndAtomIdx()))
    p, q = points[front], points[apex]
    dx, dy = q[0] - p[0], q[1] - p[1]
    fraction = ((p[0] - a[0]) * dy - (p[1] - a[1]) * dx) / ((b[0] - a[0]) * dy - (b[1] - a[1]) * dx)
    rear.SetDoubleProp(GAP_FRACTION, fraction)
    return trial


def bridge_gaps(svg: str, mol: Chem.Mol, coords: list[tuple[float, float]]) -> str:
    """Cut only ink belonging to a rear bond, including its highlight, on any background."""
    masks = []
    for bond in mol.GetBonds():
        if not bond.HasProp(GAP_FRACTION):
            continue
        f = bond.GetDoubleProp(GAP_FRACTION)
        a, b = coords[bond.GetBeginAtomIdx()], coords[bond.GetEndAtomIdx()]
        x, y = (a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1]))
        idx = bond.GetIdx()
        mask = f"bridge-gap-{idx}"
        # RDKit emits separate paths for hatches, colour splits and highlights.
        # The mask region starts at user-space (0, 0): fine while bond ink stays at
        # non-negative coordinates (RDKit's canvas; scheme._expand only grows the
        # viewBox outwards for arrows). Anything that moves bonds to negative
        # coordinates must widen this region, or the rear bond disappears entirely.
        svg = re.sub(rf"<path class='bond-{idx}(?=[ '])", f"<path mask='url(#{mask})' class='bond-{idx}", svg)
        masks.append(f"<mask id='{mask}' maskUnits='userSpaceOnUse' x='0' y='0' width='100%' height='100%'>"
                     f"<rect width='100%' height='100%' fill='white'/>"
                     f"<circle cx='{x:.2f}' cy='{y:.2f}' r='3.2' fill='black'/></mask>")
    return svg.replace("</svg>", "<defs>" + "".join(masks) + "</defs></svg>") if masks else svg
