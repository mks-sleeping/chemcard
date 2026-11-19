"""Reactant side of a mechanism step drawn as ONE picture, with curved arrows.

Curved (electron-pushing) arrows often run from one molecule to another
(OH⁻ attacking CH₃Br), so all reactants of the step share one coordinate
system: each molecule is laid out on its own, rotated so its reacting atom
faces its partner, placed left to right, and drawn by RDKit as a single
molecule. Arrows are quadratic Béziers added to RDKit's SVG.

Arrow endpoints are written by the model with atom-map numbers:
  lp:N       lone pair on atom N (source only)
  atom:N     atom N
  bond:N-M   bond between N and M (as a target it may be a bond being formed)

Pure chemistry: no `mcp`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Geometry import Point2D

from chemcard.chem.highlight import canonical_maps
from chemcard.chem.render import (
    BOND_LENGTH_PX,
    FONT_SIZE_PX,
    _draw,
    _label_dummies,
    _layout,
    mol_from_smiles,
    normalize_layout,
)

MAX_ARROWS = 8
LONE_PAIR_ELEMENTS = {"N", "O", "S", "P", "F", "Cl", "Br", "I", "Se"}
ARROW_COLOR = "#222222"
GAP = 2.4  # space between molecules, in bond lengths (RDKit coordinates use 1.5 per bond)


class InvalidArrow(ValueError):
    """An arrow cannot be resolved. The message goes to the model verbatim."""


@dataclass(frozen=True)
class ArrowSpec:
    src: str            # "lp:1", "bond:1-2" or "atom:1"
    dst: str            # "atom:2" or "bond:1-2"
    style: str = "full"  # "full" (electron pair) or "half" (single electron, fishhook)


# A resolved endpoint, JSON-friendly: ["lp", s, a] / ["atom", s, a] / ["bond", s, a, t, b]
# where s/t are species indices and a/b canonical atom indices in that species.
Endpoint = list


@dataclass
class Scheme:
    svg: str
    width: int
    height: int
    atoms: list[dict] = field(default_factory=list)    # {"x","y","el","species","index"} in SVG px
    species: list[dict] = field(default_factory=list)  # {"cx","bottom"} per species, for labels


# ---- resolving the model's arrow specs ----------------------------------------

def resolve_arrows(smiles_list: list[str], specs: list[ArrowSpec]) -> tuple[list[str], list[dict]]:
    """Map arrow specs onto (species, canonical atom index). Returns canonical SMILES and arrows."""
    if len(specs) > MAX_ARROWS:
        raise InvalidArrow(f"{len(specs)} arrows; draw at most {MAX_ARROWS} per step")
    canonicals, where = [], {}
    for s, smi in enumerate(smiles_list):
        canonical, maps = canonical_maps(smi)
        canonicals.append(canonical)
        for num, idx in maps.items():
            if num in where:
                raise InvalidArrow(f"atom-map number {num} is used in more than one reactant; number atoms uniquely across the step")
            where[num] = (s, idx)
    mols = [mol_from_smiles(c) for c in canonicals]

    def atom(num: int, text: str) -> tuple[int, int]:
        if num not in where:
            known = ", ".join(str(n) for n in sorted(where)) or "none"
            raise InvalidArrow(f"{text}: atom-map number {num} not found in the reactants (numbers present: {known}); write it in a SMILES, e.g. [OH-:1]")
        return where[num]

    def parse(text: str, role: str) -> Endpoint:
        kind, _, rest = text.partition(":")
        try:
            if kind == "lp" and role == "src":
                s, a = atom(int(rest), text)
                at = mols[s].GetAtomWithIdx(a)
                if at.GetSymbol() not in LONE_PAIR_ELEMENTS and at.GetFormalCharge() >= 0:
                    raise InvalidArrow(f"{text}: {at.GetSymbol()} has no lone pair here; start from a bond (bond:N-M) instead")
                return ["lp", s, a]
            if kind == "atom":
                s, a = atom(int(rest), text)
                return ["atom", s, a]
            if kind == "bond":
                n, m = (int(x) for x in rest.split("-"))
                (s, a), (t, b) = atom(n, text), atom(m, text)
                if role == "src" and (s != t or mols[s].GetBondBetweenAtoms(a, b) is None):
                    raise InvalidArrow(f"{text}: an arrow can only start from an existing bond; atoms {n} and {m} are not bonded")
                return ["bond", s, a, t, b]
        except ValueError as e:
            if isinstance(e, InvalidArrow):
                raise
        allowed = "lp:N, bond:N-M or atom:N" if role == "src" else "atom:N or bond:N-M"
        raise InvalidArrow(f"cannot read arrow {'start' if role == 'src' else 'end'} {text!r}; use {allowed}")

    arrows = []
    for spec in specs:
        if spec.style not in ("full", "half"):
            raise InvalidArrow(f"arrow style {spec.style!r}; use 'full' (two electrons) or 'half' (one electron)")
        arrows.append({"from": parse(spec.src, "src"), "to": parse(spec.dst, "dst"), "style": spec.style})
    return canonicals, arrows


# ---- layout ------------------------------------------------------------------

def compose(
    canonicals: list[str],
    arrows: list[dict],
    colors: list[tuple[dict[int, str], dict[int, str]]] | None = None,
) -> Scheme:
    """Draw the reactants side by side with `+` between them and the arrows on top."""
    colors = colors or [({}, {}) for _ in canonicals]
    mols = []
    for c in canonicals:
        m = mol_from_smiles(c)
        _label_dummies(m)
        _layout(m)
        normalize_layout(m)
        mols.append(m)

    _face_partners(mols, arrows)
    plus_x = _place(mols)

    combined, offsets = mols[0], [0]
    for m in mols[1:]:
        offsets.append(combined.GetNumAtoms())
        combined = Chem.CombineMols(combined, m)
    atom_colors, bond_colors = {}, {}
    bond_offsets = [0]
    for m in mols[:-1]:
        bond_offsets.append(bond_offsets[-1] + m.GetNumBonds())
    for (ac, bc), ao, bo in zip(colors, offsets, bond_offsets):
        atom_colors.update({i + ao: c for i, c in ac.items()})
        bond_colors.update({i + bo: c for i, c in bc.items()})

    out: list = []
    svg, width, height, coords = _draw(combined, atom_colors, bond_colors, layout=False, drawer_out=out)
    drawer = out[0]

    def px(x: float, y: float) -> tuple[float, float]:
        p = drawer.GetDrawCoords(Point2D(x, y))
        return (p.x, p.y)

    species_of = [s for s, m in enumerate(mols) for _ in range(m.GetNumAtoms())]
    local_of = [i for m in mols for i in range(m.GetNumAtoms())]
    labelled = [_has_label(combined.GetAtomWithIdx(i)) for i in range(combined.GetNumAtoms())]
    pluses = [px(x, 0.0) for x in plus_x]
    geo = _Geometry(combined, coords, labelled, obstacles=pluses)

    shapes = [geo.arrow(_glob(a["from"], offsets), _glob(a["to"], offsets), a["style"]) for a in arrows]

    extra = "".join(s.svg for s in shapes)
    extra += "".join(
        f"<text x='{x:.1f}' y='{y:.1f}' font-family='system-ui,-apple-system,Segoe UI,sans-serif' font-size='18' "
        f"fill='#1a1a1a' text-anchor='middle' dominant-baseline='central'>+</text>"
        for x, y in pluses
    )
    pts = [p for s in shapes for p in s.points] + [(x, y - 10) for x, y in pluses] + [(x, y + 10) for x, y in pluses]
    svg, width, height, (dx, dy) = _expand(svg, width, height, pts, extra)

    atoms = [
        {"x": round(x + dx, 1), "y": round(y + dy, 1), "el": _symbol(combined.GetAtomWithIdx(i)),
         "species": species_of[i], "index": local_of[i]}
        for i, (x, y) in enumerate(coords)
    ]
    species = []
    for s, m in enumerate(mols):
        idx = [i for i in range(len(coords)) if species_of[i] == s]
        xs = [coords[i][0] for i in idx]
        ys = [coords[i][1] for i in idx]
        species.append({"cx": round((min(xs) + max(xs)) / 2 + dx, 1), "bottom": round(max(ys) + dy + 12, 1)})
    return Scheme(svg=svg, width=width, height=height, atoms=atoms, species=species)


def _glob(end: Endpoint, offsets: list[int]) -> tuple:
    if end[0] == "bond":
        return ("bond", offsets[end[1]] + end[2], offsets[end[3]] + end[4])
    return (end[0], offsets[end[1]] + end[2])


def _face_partners(mols: list[Chem.Mol], arrows: list[dict]) -> None:
    """Rotate each molecule so the atom an arrow links to another molecule faces that molecule."""
    def atoms_of(end):
        return [(end[1], end[2])] + ([(end[3], end[4])] if end[0] == "bond" else [])

    for s, m in enumerate(mols):
        if m.GetNumAtoms() < 2:
            continue
        target = None
        for a in arrows:
            ends = atoms_of(a["from"]) + atoms_of(a["to"])
            mine = [i for sp, i in ends if sp == s]
            others = [sp for sp, _ in ends if sp != s]
            if mine and others:
                target = (mine[0], 1.0 if others[0] > s else -1.0)
                break
        if target is None:
            continue
        idx, side = target
        conf = m.GetConformer()
        pos = conf.GetPositions()[:, :2]
        cx, cy = pos.mean(axis=0)
        vx, vy = pos[idx][0] - cx, pos[idx][1] - cy
        if math.hypot(vx, vy) < 1e-3:
            continue
        rot = math.atan2(0.0, side) - math.atan2(vy, vx)
        c, s_ = math.cos(rot), math.sin(rot)
        for i, (x, y) in enumerate(pos):
            x0, y0 = x - cx, y - cy
            conf.SetAtomPosition(i, (cx + c * x0 - s_ * y0, cy + s_ * x0 + c * y0, 0.0))


def _place(mols: list[Chem.Mol]) -> list[float]:
    """Left to right, vertically centred; returns x of the '+' between neighbours."""
    unit = 1.5  # RDKit bond length
    cursor, plus = 0.0, []
    for k, m in enumerate(mols):
        conf = m.GetConformer()
        pos = conf.GetPositions()[:, :2]
        pad = 0.45 * unit  # room for atom labels at the edges
        minx, maxx = pos[:, 0].min() - pad, pos[:, 0].max() + pad
        cy = (pos[:, 1].min() + pos[:, 1].max()) / 2
        shift = cursor - minx
        for i, (x, y) in enumerate(pos):
            conf.SetAtomPosition(i, (x + shift, y - cy, 0.0))
        cursor = maxx + shift
        if k < len(mols) - 1:
            plus.append(cursor + GAP * unit / 2)
            cursor += GAP * unit
    return plus


def _has_label(atom: Chem.Atom) -> bool:
    return atom.GetSymbol() != "C" or atom.GetFormalCharge() != 0 or atom.GetDegree() == 0 or atom.HasProp("atomLabel")


def _symbol(atom: Chem.Atom) -> str:
    return "R" if atom.GetAtomicNum() == 0 else atom.GetSymbol()


# ---- arrow geometry ----------------------------------------------------------

@dataclass
class _Shape:
    svg: str
    points: list[tuple[float, float]]


class _Geometry:
    def __init__(self, mol: Chem.Mol, coords: list[tuple[float, float]], labelled: list[bool],
                 obstacles: list[tuple[float, float]] = ()):
        self.mol, self.p, self.labelled = mol, coords, labelled
        self.obstacles = list(obstacles)  # e.g. the '+' signs between molecules
        self.B = BOND_LENGTH_PX
        self.segments = [(self.p[b.GetBeginAtomIdx()], self.p[b.GetEndAtomIdx()], {b.GetBeginAtomIdx(), b.GetEndAtomIdx()})
                         for b in mol.GetBonds()]


    def anchor(self, end: tuple, toward: tuple[float, float] | None) -> tuple[tuple[float, float], set[int]]:
        """Centre point of an endpoint and the atoms it belongs to."""
        if end[0] == "bond":
            a, b = end[1], end[2]
            if self.mol.GetBondBetweenAtoms(a, b) is None:  # bond being formed: aim at the far atom
                return self.p[b], {a, b}
            (x1, y1), (x2, y2) = self.p[a], self.p[b]
            return ((x1 + x2) / 2, (y1 + y2) / 2), {a, b}
        i = end[1]
        if end[0] == "lp":
            return self.lone_pair_centre(i, toward), {i}
        return self.p[i], {i}

    def lone_pair_centre(self, i: int, toward: tuple[float, float] | None) -> tuple[float, float]:
        x, y = self.p[i]
        nbrs = [n.GetIdx() for n in self.mol.GetAtomWithIdx(i).GetNeighbors()]
        if nbrs:
            vx = -sum(self.p[n][0] - x for n in nbrs)
            vy = -sum(self.p[n][1] - y for n in nbrs)
        else:
            vx, vy = (toward[0] - x, toward[1] - y) if toward else (0.0, -1.0)
        if math.hypot(vx, vy) < 1e-6:  # symmetric atom: go "up"
            vx, vy = 0.0, -1.0
        n = math.hypot(vx, vy)
        r = self.extent(i, vx, vy) + 4
        return (x + vx / n * r, y + vy / n * r)

    def arrow(self, src: tuple, dst: tuple, style: str) -> _Shape:
        dst_c, dst_atoms = self.anchor(dst, None)
        src_c, src_atoms = self.anchor(src, dst_c)
        best = None
        for side in (1.0, -1.0):
            for frac in (0.35, 0.6):
                for theta in (45, 70, 90, 115):
                    cand = self._curve(src, dst, src_c, dst_c, src_atoms | dst_atoms, side, frac, math.radians(theta))
                    if cand and (best is None or cand[0] < best[0]):
                        best = cand
        if best is None:
            return _Shape("", [])
        _, p0, q, p2 = best
        return self._render(src, src_c, p0, q, p2, style)

    def _curve(self, src, dst, src_c, dst_c, own: set[int], side: float, frac: float, theta: float):
        # 1. endpoints on the label edges, entering from the bulge side; 2. then the bulge.
        # (Bulging first and trimming afterwards overshoots on short arrows, e.g. C=O → O.)
        if math.hypot(dst_c[0] - src_c[0], dst_c[1] - src_c[1]) < 1e-3:
            return None
        p0 = self._endpoint(src, src_c, dst_c, side * theta)
        p2 = self._endpoint(dst, dst_c, src_c, -side * theta)
        (x0, y0), (x2, y2) = p0, p2
        L = math.hypot(x2 - x0, y2 - y0)
        if L < 2:
            return None
        nx, ny = -(y2 - y0) / L, (x2 - x0) / L
        h = min(max(frac * L, 0.45 * self.B), 1.6 * self.B)
        q = ((x0 + x2) / 2 + side * h * nx, (y0 + y2) / 2 + side * h * ny)
        # short arrows (bond → its own atom) read badly: prefer the variant with a longer chord
        score = 0.02 * h + 0.01 * math.degrees(theta) + max(0.0, 0.8 * self.B - L) * 0.5
        for k in range(1, 12):
            t = k / 12
            bx = (1 - t) ** 2 * x0 + 2 * (1 - t) * t * q[0] + t ** 2 * x2
            by = (1 - t) ** 2 * y0 + 2 * (1 - t) * t * q[1] + t ** 2 * y2
            for i, (ax, ay) in enumerate(self.p):
                if i in own:
                    continue
                d = math.hypot(bx - ax, by - ay)
                limit = self.extent(i, bx - ax, by - ay) + 3
                if d < limit:
                    score += (limit - d) * 3
            for a, b, ends in self.segments:
                if ends & own:
                    continue
                if _dist_to_segment((bx, by), a, b) < 3:
                    score += 6
            for ox, oy in self.obstacles:
                d = math.hypot(bx - ox, by - oy)
                if d < 12:
                    score += (12 - d) * 3
        return score, p0, q, p2

    def extent(self, i: int, dx: float, dy: float) -> float:
        """Distance from atom i's centre to the edge of its label in direction (dx, dy)."""
        if not self.labelled[i]:
            return 3.0
        F = FONT_SIZE_PX
        rx, ry = max(0.6 * F, 0.33 * F * len(self.label_text(i))), 0.62 * F
        n = math.hypot(dx, dy) or 1.0
        ux, uy = dx / n, dy / n
        return 1.0 / math.sqrt((ux / rx) ** 2 + (uy / ry) ** 2) + 2

    def label_text(self, i: int) -> str:
        a = self.mol.GetAtomWithIdx(i)
        h = a.GetTotalNumHs()
        q = a.GetFormalCharge()
        return a.GetSymbol() + ("H" + (str(h) if h > 1 else "") if h else "") + ("+" if q > 0 else "-" if q < 0 else "")

    def _endpoint(self, end: tuple, centre, other, angle: float):
        """Point where the arrow leaves/meets this endpoint, rotated `angle` off the straight line."""
        x, y = centre
        dx, dy = other[0] - x, other[1] - y
        n = math.hypot(dx, dy) or 1.0
        ux, uy = dx / n, dy / n
        c, s = math.cos(angle), math.sin(angle)
        rx, ry = ux * c - uy * s, ux * s + uy * c
        if end[0] == "atom" or (end[0] == "bond" and self.mol.GetBondBetweenAtoms(end[1], end[2]) is None):
            i = end[1] if end[0] == "atom" else end[2]
            r = self.extent(i, rx, ry)
        elif end[0] == "lp":
            r = 3.0
        else:  # existing bond: start just off the bond's midpoint, on the bulge side
            r = 2.0
        return (x + rx * r, y + ry * r)

    def _render(self, src, src_c, p0, q, p2, style) -> _Shape:
        parts = [f"<path class='arrow' d='M {p0[0]:.1f} {p0[1]:.1f} Q {q[0]:.1f} {q[1]:.1f} {p2[0]:.1f} {p2[1]:.1f}' "
                 f"fill='none' stroke='{ARROW_COLOR}' stroke-width='1.4' stroke-linecap='round'/>"]
        tx, ty = p2[0] - q[0], p2[1] - q[1]
        n = math.hypot(tx, ty) or 1.0
        ux, uy = tx / n, ty / n
        L, w = 7.5, 3.4
        bx, by = p2[0] - ux * L, p2[1] - uy * L
        left = (bx - uy * w, by + ux * w)
        right = (bx + uy * w, by - ux * w)
        if style == "half":  # fishhook: one barb
            parts.append(f"<path d='M {p2[0]:.1f} {p2[1]:.1f} L {left[0]:.1f} {left[1]:.1f}' stroke='{ARROW_COLOR}' "
                         f"stroke-width='1.4' stroke-linecap='round'/>")
        else:
            parts.append(f"<path d='M {p2[0]:.1f} {p2[1]:.1f} L {left[0]:.1f} {left[1]:.1f} L {right[0]:.1f} {right[1]:.1f} Z' "
                         f"fill='{ARROW_COLOR}'/>")
        pts = [p0, q, p2, left, right]
        if src[0] == "lp":  # two dots for the lone pair, across the arrow's start direction
            cx, cy = src_c
            dx, dy = q[0] - cx, q[1] - cy
            m = math.hypot(dx, dy) or 1.0
            ox, oy = -dy / m * 3.2, dx / m * 3.2
            for sgn in (1, -1):
                parts.append(f"<circle cx='{cx + sgn * ox:.1f}' cy='{cy + sgn * oy:.1f}' r='1.6' fill='{ARROW_COLOR}'/>")
                pts.append((cx + sgn * ox, cy + sgn * oy))
        return _Shape("<g class='mechanism-arrow'>" + "".join(parts) + "</g>", pts)


def _dist_to_segment(p, a, b) -> float:
    (px, py), (ax, ay), (bx, by) = p, a, b
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    if L2 == 0:
        return math.hypot(px - ax, py - ay)
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / L2))
    return math.hypot(px - (ax + t * dx), py - (ay + t * dy))


def _expand(svg: str, width: int, height: int, pts, extra: str):
    """Grow the canvas so arrows outside the molecules are not clipped; append `extra`."""
    margin = 4
    xs = [x for x, _ in pts] or [0]
    ys = [y for _, y in pts] or [0]
    left = max(0, math.ceil(margin - min(xs)))
    top = max(0, math.ceil(margin - min(ys)))
    right = max(0, math.ceil(max(xs) + margin - width))
    bottom = max(0, math.ceil(max(ys) + margin - height))
    W, H = width + left + right, height + top + bottom
    import re

    svg = re.sub(r"width='\d+(?:\.\d+)?px'", f"width='{W}px'", svg, count=1)
    svg = re.sub(r"height='\d+(?:\.\d+)?px'", f"height='{H}px'", svg, count=1)
    svg = re.sub(r"viewBox='[^']*'", f"viewBox='{-left} {-top} {W} {H}'", svg, count=1)
    svg = svg.replace("</svg>", extra + "</svg>")
    return svg, W, H, (left, top)
