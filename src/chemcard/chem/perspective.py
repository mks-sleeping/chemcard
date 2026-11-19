"""Perspective views of any bridged bicyclic core, bicyclo[a.b.c] with c ≥ 1, by projecting 3D.

A textbook perspective drawing of norbornane, a tropane or DABCO is what the 3D
molecule looks like from a suitable angle. So instead of hand-writing coordinates
for every bridge-length combination:

1. find a bridged bicyclic ring system (two non-adjacent bridgeheads joined by three
   bridges; heteroatoms, double bonds and substituents allowed);
2. embed the molecule in 3D with its specified stereochemistry (ETKDG, fixed seed);
3. look at the core from many directions and keep views whose projection has exactly
   one bond crossing (the depth cue), no near-coincident atoms, no atom on a bond and
   little foreshortening;
4. turn the view so the shortest bridge points up, fix those core coordinates, let
   RDKit lay out the rest of the molecule around them, and mark the rear bond at the
   crossing so it is drawn with a small gap (bridged_layout.bridge_gaps).

The bridgehead configuration is right by construction: the projection is of the
stereoisomer that was asked for. Ring systems with more than two rings (cages,
morphinans) and bicycles fused to further rings are left to other generators.

Pure chemistry: no `mcp`.
"""

from __future__ import annotations

import math
from collections.abc import Iterator

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Geometry import Point2D

from chemcard.chem.bridged_layout import GAP_FRACTION

N_VIEWS = 400            # directions sampled on the half sphere
MAX_CANDIDATES = 2       # best distinct views offered to the layout chooser
EMBED_SEED = 0xC0FFEE
# Above this the bicycle is a small part of a big molecule: a 3D embedding costs up to
# seconds (amphotericin B: 1.5 s) and the projected view practically never wins.
MAX_HEAVY_ATOMS = 40


def bicyclic_cores(mol: Chem.Mol) -> list[dict]:
    """Bridged bicyclic ring systems: {"atoms", "heads": (h1, h2), "bridges": [[...], [...], [...]]}.

    Bridges are listed longest first; each is the atom path strictly between the heads.
    """
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    systems: list[list[set[int]]] = []
    for ring in rings:  # rings sharing ≥ 2 atoms belong to one fused/bridged system
        joined = [s for s in systems if any(len(ring & r) >= 2 for r in s)]
        merged = [ring] + [r for s in joined for r in s]
        systems = [s for s in systems if s not in joined] + [merged]
    cores = []
    for system in systems:
        atoms = set().union(*system)
        ring_bonds = {b.GetIdx() for b in mol.GetBonds()
                      if b.IsInRing() and b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms}
        # bicyclic = cyclomatic number 2. (RDKit's ring list adds symmetric extra rings:
        # bicyclo[2.2.2]octane is listed with three, so counting rings would miss it.)
        if len(ring_bonds) - len(atoms) + 1 != 2:
            continue
        adj = {i: [n.GetIdx() for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetIdx() in atoms
                   and mol.GetBondBetweenAtoms(i, n.GetIdx()).IsInRing()] for i in atoms}
        heads = [i for i in atoms if len(adj[i]) == 3]
        if len(heads) != 2 or mol.GetBondBetweenAtoms(*heads) is not None:
            continue  # fused (decalin-like) or something else
        h1, h2 = heads
        bridges = []
        for first in adj[h1]:
            path, prev, curr = [], h1, first
            while curr != h2:
                path.append(curr)
                nxt = [n for n in adj[curr] if n != prev]
                if len(nxt) != 1:
                    break
                prev, curr = curr, nxt[0]
            else:
                bridges.append(path)
                continue
            break
        if len(bridges) != 3 or min(len(b) for b in bridges) < 1:
            continue
        bridges.sort(key=len, reverse=True)
        cores.append({"atoms": atoms, "heads": (h1, h2), "bridges": bridges})
    return cores


def perspective_candidates(mol: Chem.Mol) -> Iterator[tuple[str, Chem.Mol]]:
    if mol.GetNumHeavyAtoms() > MAX_HEAVY_ATOMS:
        return
    cores = bicyclic_cores(mol)
    if len(cores) != 1:  # several bicycles would need a joint arrangement; not attempted
        return
    core = cores[0]
    xyz = _embed(mol)
    if xyz is None:
        return
    a, b, c = (len(br) for br in core["bridges"])
    for k, (coords2d, gap) in enumerate(_best_views(mol, core, xyz), start=1):
        trial = Chem.Mol(mol)
        coord_map = {i: Point2D(float(x), float(y)) for i, (x, y) in coords2d.items()}
        try:
            rdDepictor.Compute2DCoords(trial, coordMap=coord_map, canonOrient=False)
        except (RuntimeError, ValueError):
            continue
        if _has_straight_angle(trial):
            continue  # RDKit then prints a "C" label on the carbon and the drawing reads wrong
        bond_idx, fraction = gap
        trial.GetBondWithIdx(bond_idx).SetDoubleProp(GAP_FRACTION, fraction)
        yield f"bridge-{a}.{b}.{c}-view{k}", trial


STRAIGHT_DEG = 165.0  # two bonds at one atom closer to a straight line than this look like one bond


def _has_straight_angle(m: Chem.Mol) -> bool:
    pos = m.GetConformer().GetPositions()
    for atom in m.GetAtoms():
        if atom.GetHybridization() == Chem.HybridizationType.SP:
            continue  # alkynes, nitriles, allenes are straight on purpose
        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(nbrs) < 2:
            continue
        i = atom.GetIdx()
        for k in range(len(nbrs)):
            for j in range(k + 1, len(nbrs)):
                if _angle(pos[i], pos[nbrs[k]], pos[nbrs[j]]) > STRAIGHT_DEG:
                    return True
    return False


def _angle(centre, a, b) -> float:
    va = np.asarray(a[:2]) - np.asarray(centre[:2])
    vb = np.asarray(b[:2]) - np.asarray(centre[:2])
    cos = float(va @ vb) / (float(np.linalg.norm(va)) * float(np.linalg.norm(vb)) or 1.0)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def _embed(mol: Chem.Mol) -> np.ndarray | None:
    mh = Chem.AddHs(mol, addCoords=False)
    params = AllChem.ETKDGv3()
    params.randomSeed = EMBED_SEED
    if AllChem.EmbedMolecule(mh, params) != 0:
        return None
    pos = mh.GetConformer().GetPositions()
    return pos[: mol.GetNumAtoms()]  # AddHs appends H atoms after the original ones


def _hemisphere(n: int) -> np.ndarray:
    """Roughly even directions on the z ≥ 0 half sphere (Fibonacci spiral)."""
    k = np.arange(n) + 0.5
    z = k / n
    r = np.sqrt(1 - z * z)
    phi = k * math.pi * (3 - math.sqrt(5))
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def _best_views(mol, core, xyz):
    atoms = sorted(core["atoms"])
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetIdx()) for b in mol.GetBonds()
             if b.GetBeginAtomIdx() in core["atoms"] and b.GetEndAtomIdx() in core["atoms"]]
    centre = xyz[atoms].mean(axis=0)
    scored = []
    for v in _hemisphere(N_VIEWS):
        u = np.cross(v, [0.0, 0.0, 1.0] if abs(v[2]) < 0.9 else [1.0, 0.0, 0.0])
        u /= np.linalg.norm(u)
        w = np.cross(v, u)
        rel = xyz - centre
        p2 = {i: (rel[i] @ u, rel[i] @ w) for i in atoms}
        depth = {i: rel[i] @ v for i in atoms}  # larger = nearer the viewer
        result = _score_view(p2, bonds, core)
        if result is None:
            continue
        score, crossing = result
        scored.append((score, p2, depth, crossing))
    scored.sort(key=lambda s: s[0])

    chosen = []
    for score, p2, depth, crossing in scored:
        coords = _orient(p2, core)
        if any(_same_drawing(coords, prev) for prev, _ in chosen):
            continue
        chosen.append((coords, _gap(coords, depth, crossing, bonds)))
        if len(chosen) == MAX_CANDIDATES:
            break
    return chosen


def _score_view(p2, bonds, core):
    """None if unusable, else (score, the crossing bond pair). Lower score is better."""
    lengths = [math.dist(p2[a], p2[b]) for a, b, _ in bonds]
    unit = float(np.median(lengths))
    if unit < 1e-6 or min(lengths) < 0.55 * unit:
        return None  # a bond seen nearly end-on
    nbrs = {i: [] for i in p2}
    for a, b, _ in bonds:
        nbrs[a].append(b)
        nbrs[b].append(a)
    for i, ns in nbrs.items():
        if any(_angle(p2[i], p2[ns[k]], p2[ns[j]]) > STRAIGHT_DEG
               for k in range(len(ns)) for j in range(k + 1, len(ns))):
            return None
    h1, h2 = core["heads"]
    head_sep = math.dist(p2[h1], p2[h2]) / unit
    if head_sep < 0.8:
        return None  # looking down the bridgehead axis: the bridges pile up in the middle
    crossings = []
    for i in range(len(bonds)):
        for j in range(i + 1, len(bonds)):
            a, b, _ = bonds[i]
            c, d, _ = bonds[j]
            if len({a, b, c, d}) == 4 and _segments_cross(p2[a], p2[b], p2[c], p2[d]):
                crossings.append((i, j))
    if len(crossings) != 1:
        return None  # 0 = the flat view the generic engines already draw; >1 = unreadable
    pts = list(p2)
    min_sep = min(math.dist(p2[i], p2[j]) for k, i in enumerate(pts) for j in pts[k + 1:]) / unit
    if min_sep < 0.5:
        return None
    min_on_bond = min(
        _dist_point_segment(p2[i], p2[a], p2[b]) for i in pts for a, b, _ in bonds if i not in (a, b)
    ) / unit
    if min_on_bond < 0.3:
        return None
    spread = float(np.std(lengths)) / unit
    # clearance matters most; then even bond lengths; textbook views keep the bridgeheads apart
    score = (-2.0 * min(min_sep, 1.0) - 2.0 * min(min_on_bond, 0.8) + 3.0 * spread
             - 1.0 * min(head_sep, 2.0))
    return score, crossings[0]


def _orient(p2, core) -> dict[int, tuple[float, float]]:
    """Rotate so the shortest bridge points up, centre on the bridgeheads, scale to RDKit's 1.5 bond."""
    h1, h2 = core["heads"]
    mid = np.array([(p2[h1][0] + p2[h2][0]) / 2, (p2[h1][1] + p2[h2][1]) / 2])
    top = np.mean([p2[i] for i in core["bridges"][-1]], axis=0) - mid
    ang = math.pi / 2 - math.atan2(top[1], top[0])
    c, s = math.cos(ang), math.sin(ang)
    out = {}
    for i, (x, y) in p2.items():
        x0, y0 = x - mid[0], y - mid[1]
        out[i] = (c * x0 - s * y0, s * x0 + c * y0)
    # scale so the median ring bond is RDKit's standard 1.5
    h1, h2 = core["heads"]
    paths = [[h1] + br + [h2] for br in core["bridges"]]
    ds = [math.dist(out[p[k]], out[p[k + 1]]) for p in paths for k in range(len(p) - 1)]
    k = 1.5 / float(np.median(ds))
    return {i: (x * k, y * k) for i, (x, y) in out.items()}


def _same_drawing(a, b) -> bool:
    """Two views count as the same if every atom moved less than a quarter bond."""
    return max(math.dist(a[i], b[i]) for i in a) < 0.25 * 1.5


def _gap(coords, depth, crossing, bonds):
    """The bond behind at the crossing and where along it (fraction from its begin atom)."""
    (a, b, i1), (c, d, i2) = bonds[crossing[0]], bonds[crossing[1]]
    t = _intersection_param(coords[a], coords[b], coords[c], coords[d])
    s = _intersection_param(coords[c], coords[d], coords[a], coords[b])
    za = depth[a] + t * (depth[b] - depth[a])
    zc = depth[c] + s * (depth[d] - depth[c])
    return (i1, t) if za < zc else (i2, s)


def _intersection_param(p, q, r, s) -> float:
    """t such that p + t (q - p) lies on line r–s."""
    px, py = p
    dx, dy = q[0] - px, q[1] - py
    ex, ey = s[0] - r[0], s[1] - r[1]
    den = dx * ey - dy * ex
    return ((r[0] - px) * ey - (r[1] - py) * ex) / den


def _segments_cross(p1, p2, p3, p4) -> bool:
    def side(a, b, c):
        return (c[0] - a[0]) * (b[1] - a[1]) - (b[0] - a[0]) * (c[1] - a[1])

    return side(p3, p4, p1) * side(p3, p4, p2) < 0 and side(p1, p2, p3) * side(p1, p2, p4) < 0


def _dist_point_segment(p, a, b) -> float:
    dx, dy = b[0] - a[0], b[1] - a[1]
    L2 = dx * dx + dy * dy
    if L2 == 0:
        return math.dist(p, a)
    t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2))
    return math.dist(p, (a[0] + t * dx, a[1] + t * dy))
