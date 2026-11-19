"""Drawings of bicyclo[a.b.c] ring systems, chosen from the topology.

Topology first (bicyclo_systems): a ring system with cyclomatic number 2 and exactly
two ring atoms of ring degree 3 is a bicyclo[a.b.c] system. The two bridgeheads are
joined by three paths of a ≥ b ≥ c atoms; heteroatoms, double bonds and substituents
are allowed. c = 0 means the bridgeheads share a bond (decalin is [4.4.0]).

Then the drawing strategy, from the bridge lengths:

- c = 0 (fused): two regular polygons sharing the bridgehead bond, flat, no crossing
  (fused_candidates). The general engines draw these the same way, so this candidate
  only wins when it is clearly better.
- c ≥ 1 (bridged): the textbook perspective, projected from 3D (perspective_candidates):

  1. embed the molecule in 3D with its specified stereochemistry (ETKDG, fixed seed);
  2. look at the core from many directions; keep views whose projection of the core has
     exactly one bond crossing (the depth cue), no near-coincident atoms, no atom on a
     bond, and rank them;
  3. in that order, turn each view so the shortest bridge points up, fix the core, let
     RDKit place the rest of the molecule, mark the rear bond at the crossing (drawn
     with a small gap, bridged_layout.bridge_gaps), and only then judge it: stereo,
     straight angles and every layout_quality count. The bare core cannot tell whether
     a substituent will land on a bond, so views are not discarded before this step.
     At most MAX_PLACEMENTS views are placed;
  4. if no projected view comes out clean, construct flat drawings directly from the
     three path lengths: two paths form the outline, the third crosses inside it as a
     bowed bridge. No crossing, no depth cue needed (bicyclo[1.1.1]pentane has no
     usable projection at all).

The bridgehead configuration of a projection is right by construction: it is a view
of the stereoisomer that was asked for, and the rear bond is the one farther from the
viewer. Ring systems with more than two rings (cages, morphinans) are left to other
generators.

Pure chemistry: no `mcp`.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Iterator

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Geometry import Point2D

from chemcard.chem.bridged_layout import GAP_FRACTION

N_VIEWS = 400            # directions sampled on the half sphere
MAX_CANDIDATES = 2       # clean drawings offered to the layout chooser
MAX_PLACEMENTS = 12      # views given a full substituent placement and check (~1 ms each)
EMBED_SEED = 0xC0FFEE
# Above this the bicycle is a small part of a big molecule: a 3D embedding costs up to
# seconds (amphotericin B: 1.5 s) and the projected view practically never wins.
MAX_HEAVY_ATOMS = 40
STRAIGHT_DEG = 165.0     # two bonds at one atom closer to a straight line than this look like one bond
# A carbon with just two bonds this close to straight gets a "C" label from RDKit
# (measured: from about 161°), so the chain reads as a stray atom.
LINEAR_LABEL_DEG = 160.0
DEFECTS = ("clashes", "on_bond", "inside", "label_contacts")


# ---- topology ------------------------------------------------------------------------

def bicyclo_systems(mol: Chem.Mol) -> list[dict]:
    """Every bicyclo[a.b.c] ring system: {"atoms", "heads", "bridges", "abc"}.

    "bridges" are the three bridgehead-to-bridgehead paths (atoms strictly between the
    heads, listed from heads[0]), longest first; the last one is empty when c = 0.
    """
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    systems: list[list[set[int]]] = []
    for ring in rings:  # rings sharing ≥ 2 atoms belong to one fused/bridged system
        joined = [s for s in systems if any(len(ring & r) >= 2 for r in s)]
        merged = [ring] + [r for s in joined for r in s]
        systems = [s for s in systems if s not in joined] + [merged]
    found = []
    for system in systems:
        atoms = set().union(*system)
        ring_bonds = [b for b in mol.GetBonds()
                      if b.IsInRing() and b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms]
        # bicyclic = cyclomatic number 2. (RDKit's ring list adds symmetric extra rings:
        # bicyclo[2.2.2]octane is listed with three, so counting rings would miss it.)
        if len(ring_bonds) - len(atoms) + 1 != 2:
            continue
        adj: dict[int, list[int]] = {i: [] for i in atoms}
        for b in ring_bonds:
            adj[b.GetBeginAtomIdx()].append(b.GetEndAtomIdx())
            adj[b.GetEndAtomIdx()].append(b.GetBeginAtomIdx())
        heads = sorted(i for i in atoms if len(adj[i]) == 3)
        if len(heads) != 2 or any(len(adj[i]) != 2 for i in atoms if i not in heads):
            continue  # e.g. a spiro centre (ring degree 4)
        h1, h2 = heads
        bridges = []
        for first in adj[h1]:
            path, prev, curr = [], h1, first
            while curr != h2 and len(path) < len(atoms):
                path.append(curr)
                prev, curr = curr, next(n for n in adj[curr] if n != prev)
            bridges.append(path)
        bridges.sort(key=len, reverse=True)
        found.append({"atoms": atoms, "heads": (h1, h2), "bridges": bridges,
                      "abc": tuple(len(b) for b in bridges)})
    return found


def _single_system(mol: Chem.Mol, fused: bool) -> dict | None:
    """The molecule's one bridged (c ≥ 1) or one fused (c = 0) system, if it has exactly one.

    Several of a kind would need a joint arrangement; not attempted. One of each is fine:
    quinine's quinuclidine is projected while RDKit places its quinoline around it.
    """
    if mol.GetNumHeavyAtoms() > MAX_HEAVY_ATOMS:
        return None
    systems = [s for s in bicyclo_systems(mol) if (s["abc"][2] == 0) == fused]
    return systems[0] if len(systems) == 1 else None


# ---- c ≥ 1: perspective, or a constructed flat drawing --------------------------------

def perspective_candidates(mol: Chem.Mol) -> Iterator[tuple[str, Chem.Mol]]:
    """Bridged systems (c ≥ 1): projected views, else constructed flat drawings."""
    core = _single_system(mol, fused=False)
    if core is None:
        return
    a, b, c = core["abc"]
    clean = 0
    xyz = _embed(mol)
    if xyz is not None:
        clean = yield from _judged(mol, _projected_views(mol, core, xyz), f"bridge-{a}.{b}.{c}-view", 1)
    if clean == 0:
        yield from _judged(mol, _flat_views(core), f"bicyclo-{a}.{b}.{c}-flat", 0)


def fused_candidates(mol: Chem.Mol) -> Iterator[tuple[str, Chem.Mol]]:
    """Fused systems (c = 0): two regular polygons sharing the bridgehead bond."""
    core = _single_system(mol, fused=True)
    if core is None:
        return
    a, b, _ = core["abc"]
    h1, h2 = core["heads"]
    coords = {h1: (0.0, -0.5), h2: (0.0, 0.5)}
    for side, path in ((-1, core["bridges"][0]), (1, core["bridges"][1])):
        n = len(path) + 2
        apothem, radius = 1 / (2 * math.tan(math.pi / n)), 1 / (2 * math.sin(math.pi / n))
        for k, atom in enumerate(path, start=1):  # the long way round from h1 to h2
            phi = -math.pi / n - 2 * math.pi * k / n
            coords[atom] = (side * (apothem - radius * math.cos(phi)), radius * math.sin(phi))
    trial = _place(mol, _scaled(coords), None)
    if trial is not None:
        yield f"bicyclo-{a}.{b}.0-fused1", trial  # one drawing: the layout chooser judges it


def _judged(mol: Chem.Mol, views: Iterable, prefix: str, crossings: int):
    """Place the rest of the molecule around each core view, in order, and judge the result.

    Clean drawings have only the expected `crossings` and no other layout defect. The
    search stops after MAX_CANDIDATES clean ones or MAX_PLACEMENTS views; the clean ones
    are offered best full score first. If none is clean, the best of the others still
    compete on score. The name carries the view's rank by core score. Returns the clean count.
    """
    from chemcard.chem.render import layout_problem, layout_quality  # render imports this module

    clean, others = [], []
    for k, (coords, gap) in enumerate(views, start=1):
        if k > MAX_PLACEMENTS:
            break
        trial = _place(mol, coords, gap)
        if trial is None or layout_problem(mol, trial):
            continue
        q = layout_quality(trial)
        if q["crossings"] == crossings and not any(q[d] for d in DEFECTS):
            clean.append((q["score"], k, trial))
            if len(clean) == MAX_CANDIDATES:
                break
        else:
            others.append((q["score"], k, trial))
    for _, k, trial in sorted(clean or others, key=lambda o: o[:2])[:MAX_CANDIDATES]:
        yield f"{prefix}{k}", trial
    return len(clean)


def _place(mol, coords, gap):
    trial = Chem.Mol(mol)
    coord_map = {i: Point2D(float(x), float(y)) for i, (x, y) in coords.items()}
    try:
        rdDepictor.Compute2DCoords(trial, coordMap=coord_map, canonOrient=False)
    except (RuntimeError, ValueError):
        return None
    if _has_straight_angle(trial):
        return None  # RDKit then prints a "C" label on the carbon and the drawing reads wrong
    if gap is not None:
        bond_idx, fraction = gap
        trial.GetBondWithIdx(bond_idx).SetDoubleProp(GAP_FRACTION, fraction)
    return trial


def _has_straight_angle(m: Chem.Mol) -> bool:
    pos = m.GetConformer().GetPositions()
    for atom in m.GetAtoms():
        if atom.GetHybridization() == Chem.HybridizationType.SP:
            continue  # alkynes, nitriles, allenes are straight on purpose
        nbrs = [n.GetIdx() for n in atom.GetNeighbors()]
        if len(nbrs) > 3:
            continue  # four bonds drawn as a cross (a sulfonyl S, a quaternary C) are straight in pairs
        i = atom.GetIdx()
        limit = LINEAR_LABEL_DEG if len(nbrs) == 2 else STRAIGHT_DEG
        for k in range(len(nbrs)):
            for j in range(k + 1, len(nbrs)):
                if _angle(pos[i], pos[nbrs[k]], pos[nbrs[j]]) > limit:
                    return True
    return False


def _angle(centre, a, b) -> float:
    va = np.asarray(a[:2]) - np.asarray(centre[:2])
    vb = np.asarray(b[:2]) - np.asarray(centre[:2])
    cos = float(va @ vb) / (float(np.linalg.norm(va)) * float(np.linalg.norm(vb)) or 1.0)
    return math.degrees(math.acos(max(-1.0, min(1.0, cos))))


def _core_bonds(mol, core):
    return [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetIdx()) for b in mol.GetBonds()
            if b.GetBeginAtomIdx() in core["atoms"] and b.GetEndAtomIdx() in core["atoms"]]


def _scaled(coords):
    """Scale so the median core bond is RDKit's standard 1.5."""
    return {i: (1.5 * x, 1.5 * y) for i, (x, y) in coords.items()}


# ---- projected views -----------------------------------------------------------------

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


def _projected_views(mol, core, xyz):
    """Distinct usable views of the core, best core score first: (coords, (rear bond, fraction))."""
    atoms = sorted(core["atoms"])
    bonds = _core_bonds(mol, core)
    centre = xyz[atoms].mean(axis=0)
    rel = xyz - centre
    scored = []
    for v in _hemisphere(N_VIEWS):
        u = np.cross(v, [0.0, 0.0, 1.0] if abs(v[2]) < 0.9 else [1.0, 0.0, 0.0])
        u /= np.linalg.norm(u)
        w = np.cross(v, u)
        p2 = {i: (rel[i] @ u, rel[i] @ w) for i in atoms}
        result = _score_view(p2, bonds, core, crossings=1)
        if result is not None:
            depth = {i: rel[i] @ v for i in atoms}  # larger = nearer the viewer
            scored.append((result[0], p2, depth, result[1]))
    scored.sort(key=lambda s: s[0])
    seen = []
    for _, p2, depth, crossing in scored:
        coords = _orient(p2, core)
        if any(_same_drawing(coords, prev) for prev in seen):
            continue
        seen.append(coords)
        yield coords, _gap(coords, depth, crossing, bonds)


def _score_view(p2, bonds, core, crossings):
    """None if unusable, else (score, the crossing bond pairs). Lower score is better."""
    lengths = [math.dist(p2[a], p2[b]) for a, b, _ in bonds]
    unit = float(np.median(lengths))
    if unit < 1e-6 or min(lengths) < 0.55 * unit:
        return None  # a bond seen nearly end-on
    nbrs = {i: [] for i in p2}
    for a, b, _ in bonds:
        nbrs[a].append(b)
        nbrs[b].append(a)
    angles = [_angle(p2[i], p2[ns[k]], p2[ns[j]])
              for i, ns in nbrs.items() for k in range(len(ns)) for j in range(k + 1, len(ns))]
    if max(angles) > STRAIGHT_DEG:
        return None
    h1, h2 = core["heads"]
    head_sep = math.dist(p2[h1], p2[h2]) / unit
    if head_sep < 0.8:
        return None  # looking down the bridgehead axis: the bridges pile up in the middle
    crossing = []
    for i in range(len(bonds)):
        for j in range(i + 1, len(bonds)):
            a, b, _ = bonds[i]
            c, d, _ = bonds[j]
            if len({a, b, c, d}) == 4 and _segments_cross(p2[a], p2[b], p2[c], p2[d]):
                crossing.append((i, j))
    if len(crossing) != crossings:
        return None  # projections: 0 = the flat view the generic engines draw; >1 = unreadable
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
    return score, crossing[0] if crossing else None, min(angles)


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
    paths = [[h1] + br + [h2] for br in core["bridges"]]
    ds = [math.dist(out[p[k]], out[p[k + 1]]) for p in paths for k in range(len(p) - 1)]
    k = 1.0 / float(np.median(ds))
    return _scaled({i: (x * k, y * k) for i, (x, y) in out.items()})


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


# ---- constructed flat drawings ---------------------------------------------------------

def _flat_views(core):
    """Flat drawings built from the three path lengths, best first: (coords, None).

    The bridgeheads sit on a horizontal axis `span` apart. Two paths are drawn as
    circular arcs with unit bonds above and below the axis (the outline); the third
    crosses inside as a zigzag chain, bowed toward one side, its bonds shorter or longer
    as needed. Every choice of inner path, head distance, bow and zigzag is scored like
    a projection, except that no crossing is allowed.

    Each of the three paths is tried as the inner one, even when two have the same
    length: a carbonyl or a methyl on one of them can need the outside. The core score
    cannot see substituents, so the views are offered round-robin over the inner-path
    choices (best core score first within each), and the placement step judges them
    on the whole molecule.
    """
    h1, h2 = core["heads"]
    bridges = core["bridges"]
    bonds = [(p[k], p[k + 1], None) for p in ([h1] + br + [h2] for br in bridges) for k in range(len(p) - 1)]
    groups = []
    for inner in range(3):
        scored = []
        up, down = sorted((i for i in range(3) if i != inner), key=lambda i: -len(bridges[i]))
        for f in (0.5, 0.65, 0.8):
            span = f * (len(bridges[down]) + 1)
            s_up, s_down = _fit_sagitta(len(bridges[up]), span), _fit_sagitta(len(bridges[down]), span)
            for bow in (0.0, 0.2, -0.2, 0.4, -0.4):
                for zig in (0.0, 0.2, 0.35):
                    coords = {h1: (-span / 2, 0.0), h2: (span / 2, 0.0)}
                    coords.update(zip(bridges[up], _arc(len(bridges[up]), span, s_up)))
                    coords.update(zip(bridges[down], _arc(len(bridges[down]), span, -s_down)))
                    coords.update(zip(bridges[inner], _zigzag(len(bridges[inner]), span,
                                                              bow * (s_up if bow > 0 else s_down), zig)))
                    result = _score_view(coords, bonds, core, crossings=0)
                    if result is None or result[2] < 30:
                        continue
                    score = result[0] + max(0.0, 70 - result[2]) / 20  # penalise needle-sharp angles
                    scored.append((score, coords))
        if scored:
            groups.append(sorted(scored, key=lambda s: s[0]))
    groups.sort(key=lambda g: g[0][0])
    ranked = [g[rank] for rank in range(max(map(len, groups), default=0)) for g in groups if rank < len(g)]
    seen = []
    for _, coords in ranked:
        unit = float(np.median([math.dist(coords[a], coords[b]) for a, b, _ in bonds]))
        coords = _scaled({i: (x / unit, y / unit) for i, (x, y) in coords.items()})
        if any(_same_drawing(coords, prev) for prev in seen):
            continue
        seen.append(coords)
        yield coords, None


def _zigzag(n: int, span: float, bow: float, zig: float) -> list[tuple[float, float]]:
    """n points evenly spaced in x between the heads: a parabolic bow plus alternating offsets."""
    pts = []
    for k in range(n):
        t = (k + 1) / (n + 1)
        pts.append((-span / 2 + span * t, 4 * bow * t * (1 - t) + zig * (1 if k % 2 == 0 else -1)))
    return pts


def _arc(n: int, span: float, sagitta: float) -> list[tuple[float, float]]:
    """n points evenly spaced on the circular arc from (-span/2, 0) to (span/2, 0).

    The arc bulges by |sagitta| toward +y (sagitta > 0) or -y; 0 is the straight chord.
    """
    t = [(k + 1) / (n + 1) for k in range(n)]
    if abs(sagitta) < 1e-9:
        return [(-span / 2 + span * tk, 0.0) for tk in t]
    s, h = math.copysign(1.0, sagitta), abs(sagitta)
    theta = 4 * math.atan(2 * h / span)          # central angle
    radius = span / (2 * math.sin(theta / 2))
    return [(radius * math.cos(math.pi / 2 + theta / 2 - theta * tk),
             s * (h - radius + radius * math.sin(math.pi / 2 + theta / 2 - theta * tk))) for tk in t]


def _fit_sagitta(n: int, span: float) -> float:
    """Sagitta of the arc on which n atoms between the heads are spaced one bond apart."""
    if span >= n + 1:
        return 0.0  # the path cannot reach: draw it straight (and stretched)
    lo, hi = 0.0, 4.0 * (n + 1)
    for _ in range(60):
        mid = (lo + hi) / 2
        pts = [(-span / 2, 0.0)] + _arc(n, span, mid)
        if math.dist(pts[0], pts[1]) < 1.0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ---- geometry helpers -----------------------------------------------------------------

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
