"""SMILES → 2D structure SVG.

Pure chemistry core: must not import `mcp`. Everything here is testable
without starting a server.
"""

from __future__ import annotations

import contextlib
import io
import math
import re
from collections.abc import Iterator
from dataclasses import dataclass
from functools import partial

from rdkit import Chem, rdBase
from rdkit.Chem import rdCoordGen, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

from chemcard.chem.bridged_layout import GAP_FRACTION, bridge_candidates, bridge_gaps

MAX_SMILES_LEN = 500
MAX_HEAVY_ATOMS = 150
BOND_LENGTH_PX = 28  # fixed bond length: small molecules stay small, big ones are not squeezed
FONT_SIZE_PX = 15    # atom labels, close to body text
LINE_WIDTH_PX = 1.5
PADDING = 0.06       # fraction of the drawing

rdDepictor.SetPreferCoordGen(False)  # layout engine is chosen per molecule in _layout()


class InvalidSmiles(ValueError):
    """The SMILES cannot be parsed or exceeds limits.

    The message is returned to the model verbatim, so write it so the model
    can fix its input (say what is wrong and, if possible, where).
    """


# Soft tints: highlights read like a marker pen behind the drawing, strokes stay black.
HIGHLIGHT_RGB = {
    "red": (1.0, 0.74, 0.74),
    "blue": (0.70, 0.82, 1.0),
    "green": (0.70, 0.91, 0.70),
    "orange": (1.0, 0.84, 0.58),
    "purple": (0.85, 0.76, 1.0),
}


@dataclass(frozen=True)
class Structure:
    svg: str               # standalone <svg>; no user-supplied text inside it
    canonical_smiles: str  # RDKit canonical SMILES; atom i of the drawing is atom i of this SMILES
    formula: str           # e.g. "C8H10N4O2"; generic substituents appear as R
    width: int             # intrinsic size in px, from the fixed bond length
    height: int
    atoms: tuple[dict, ...] = ()  # per canonical atom: {"x", "y", "el"} in SVG pixels (for click targets)


def render_smiles(
    smiles: str,
    atom_colors: dict[int, str] | None = None,
    bond_colors: dict[int, str] | None = None,
) -> Structure:
    """Parse `smiles` and draw it as a compact 2D SVG.

    - Raise `InvalidSmiles` for empty / unparsable input, len > MAX_SMILES_LEN,
      or more than MAX_HEAVY_ATOMS heavy atoms.
    - Draw with a fixed bond length on a canvas sized to the molecule
      (`MolDraw2DSVG(-1, -1)` flexicanvas), transparent background.
    - Never put user text (labels, legends) into the SVG: RDKit declares
      iso-8859-1, so CJK text would be garbled. Captions belong in the widget HTML.
    - Dummy atoms (`*`) are drawn as "R" (R1, R2, … when there are several).
    - Deterministic: the same input gives byte-identical output.
    - The drawing is made from the canonical SMILES (atom maps dropped), so
      `atom_colors` / `bond_colors` and the returned `atoms` use canonical
      indices (see chem.highlight).
    """
    mol = parse_smiles(smiles)
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    canonical = Chem.MolToSmiles(mol)
    mol = mol_from_smiles(canonical)
    formula = _formula(mol)
    symbols = ["R" if a.GetAtomicNum() == 0 else a.GetSymbol() for a in mol.GetAtoms()]
    _label_dummies(mol)
    svg, width, height, coords = _draw(mol, atom_colors or {}, bond_colors or {})
    atoms = tuple({"x": x, "y": y, "el": el} for (x, y), el in zip(coords, symbols))
    return Structure(svg=svg, canonical_smiles=canonical, formula=formula, width=width, height=height, atoms=atoms)


def _smiles_params() -> Chem.SmilesParserParams:
    # Keep H atoms the SMILES writes explicitly ([H], [H:5]): mechanism steps point
    # arrows at them (proton transfers). Ordinary implicit Hs are unaffected.
    params = Chem.SmilesParserParams()
    params.removeHs = False
    return params


def mol_from_smiles(smiles: str) -> Chem.Mol | None:
    """`Chem.MolFromSmiles` without RDKit log noise (stdio servers log stderr to the host)."""
    with rdBase.BlockLogs():
        return Chem.MolFromSmiles(smiles, _smiles_params())


def parse_smiles(smiles: str) -> Chem.Mol:
    smiles = (smiles or "").strip()
    if not smiles:
        raise InvalidSmiles("empty SMILES")
    if len(smiles) > MAX_SMILES_LEN:
        raise InvalidSmiles(f"SMILES is {len(smiles)} characters; the limit is {MAX_SMILES_LEN}")

    # RDKit reports the reason on its log, not in the exception; capture it for the model.
    log = io.StringIO()
    rdBase.LogToPythonStderr()
    with contextlib.redirect_stderr(log):
        mol = Chem.MolFromSmiles(smiles, _smiles_params())
    if mol is None:
        reason = _rdkit_reason(log.getvalue())
        raise InvalidSmiles(f"RDKit could not parse {smiles!r}" + (f": {reason}" if reason else ""))

    heavy = mol.GetNumHeavyAtoms()
    if heavy > MAX_HEAVY_ATOMS:
        raise InvalidSmiles(f"{heavy} heavy atoms; the limit is {MAX_HEAVY_ATOMS}")
    return mol


def _rdkit_reason(log: str) -> str:
    lines = [re.sub(r"^\[[\d:]+\]\s*", "", ln).strip() for ln in log.splitlines()]
    lines = [ln for ln in lines if ln]
    return "; ".join(dict.fromkeys(lines))[:300]


def _formula(mol: Chem.Mol) -> str:
    n_dummy = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
    formula = rdMolDescriptors.CalcMolFormula(mol)
    if n_dummy:
        formula = re.sub(r"\*\d*", "", formula) + ("R" if n_dummy == 1 else f"R{n_dummy}")
    return formula


def _label_dummies(mol: Chem.Mol) -> None:
    dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    for i, atom in enumerate(dummies, start=1):
        atom.SetProp("atomLabel", "R" if len(dummies) == 1 else f"R<sub>{i}</sub>")


def _draw(
    mol: Chem.Mol,
    atom_colors: dict[int, str],
    bond_colors: dict[int, str],
    *,
    layout: bool = True,
    drawer_out: list | None = None,
) -> tuple[str, int, int, list[tuple[float, float]]]:
    # Flexicanvas ignores `fixedBondLength`: its pixel size follows the
    # coordinates. So lay out, draw once to measure, rescale the coordinates
    # to hit BOND_LENGTH_PX exactly, then draw for real.
    # layout=False keeps the caller's coordinates (reaction schemes place several
    # molecules themselves); `drawer_out` receives the final drawer for further
    # coordinate queries.
    if layout:
        _layout(mol)
        normalize_layout(mol)
    measured = _median_bond_px(_drawer_with(mol), mol)
    if measured:
        factor = BOND_LENGTH_PX / measured
        conf = mol.GetConformer()
        for i in range(mol.GetNumAtoms()):
            p = conf.GetAtomPosition(i)
            conf.SetAtomPosition(i, (p.x * factor, p.y * factor, 0.0))
    drawer = _drawer_with(mol, atom_colors, bond_colors)
    svg = drawer.GetDrawingText()
    # Drop the XML declaration (it claims iso-8859-1) so the SVG embeds cleanly.
    svg = re.sub(r"^<\?xml[^>]*\?>\s*", "", svg)
    coords = []
    for i in range(mol.GetNumAtoms()):
        p = drawer.GetDrawCoords(i)
        coords.append((round(p.x, 1), round(p.y, 1)))
    svg = bridge_gaps(svg, mol, coords)
    if drawer_out is not None:
        drawer_out.append(drawer)
    return svg, int(drawer.Width()), int(drawer.Height()), coords


def _scaffold_templates() -> list[tuple[str, Chem.Mol, Chem.Mol]]:
    """(reference molecule with good 2D coordinates, scaffold query) pairs, loaded once.

    For skeletons no layout engine draws well. Each .mol file in templates/ holds a
    molecule with textbook coordinates; its first line names it and its source.
    """
    global _TEMPLATES
    if _TEMPLATES is None:
        from importlib.resources import files

        _TEMPLATES = []
        for name, smarts in _TEMPLATE_SCAFFOLDS.items():
            ref = Chem.MolFromMolBlock((files("chemcard.chem") / "templates" / f"{name}.mol").read_text())
            _TEMPLATES.append((name, ref, Chem.MolFromSmarts(smarts)))
    return _TEMPLATES


_TEMPLATES: list | None = None  # [(name, reference, scaffold)]
_TEMPLATE_SCAFFOLDS = {  # most specific first
    # 4,5-epoxymorphinan: morphine, codeine, naloxone, oxycodone… (the O bridge needs its own template)
    "epoxymorphinan": "[#6]1~[#6]~[#6]~[#6]2~[#6]3~[#6]~1~[#8]~[#6]1~[#6]~[#6]~[#6]~[#6]4~[#6](~[#6]~2)~[#7]~[#6]~[#6]~[#6]~3~1~4",
    # morphinan: dextromethorphan, levorphanol… (any bonds, so unsaturation/aromaticity may differ)
    "morphinan": "[#6]1~[#6]~[#6]~[#6]2~[#6](~[#6]~1)~[#6]13~[#6]~[#6]~[#6]~[#6]~[#6]~1~[#6](~[#6]~2)~[#7]~[#6]~[#6]~3",
}


def layout_candidates(mol: Chem.Mol) -> Iterator[tuple[str, Chem.Mol]]:
    """Yield independent drawings in preference order. Atom indices stay fixed."""
    mol = Chem.Mol(mol)
    for bond in mol.GetBonds():
        bond.ClearProp(GAP_FRACTION)
    for name, ref, scaffold in _scaffold_templates():
        if not mol.HasSubstructMatch(scaffold):
            continue
        trial = Chem.Mol(mol)
        try:
            # On the reference molecule itself retain the substituent coordinates too.
            if Chem.MolToSmiles(ref) == Chem.MolToSmiles(mol):
                mapping = ref.GetSubstructMatch(mol, useChirality=True)
                rdDepictor.GenerateDepictionMatching2DStructure(trial, ref, list(zip(mapping, range(len(mapping)))))
            else:
                rdDepictor.GenerateDepictionMatching2DStructure(trial, ref, refPatt=scaffold, acceptFailure=False)
        except (RuntimeError, ValueError):
            continue
        yield f"template:{name}", trial

    yield from bridge_candidates(mol)

    for name, engine in LAYOUT_ENGINES:
        trial = Chem.Mol(mol)
        try:
            engine(trial)
        except (RuntimeError, ValueError):
            continue
        yield name, trial


def layout_problem(mol: Chem.Mol, trial: Chem.Mol) -> str | None:
    """Validate candidate coordinates against the stereo specified on the input."""
    if not trial.GetNumConformers():
        return "missing coordinates"
    pos = trial.GetConformer().GetPositions()
    if not all(math.isfinite(v) for point in pos for v in point):
        return "non-finite coordinates"
    stereo = Chem.BondStereo
    for bond in mol.GetBonds():
        tag = bond.GetStereo()
        if tag not in (stereo.STEREOE, stereo.STEREOZ, stereo.STEREOCIS, stereo.STEREOTRANS):
            continue
        left, right = bond.GetStereoAtoms()
        start, end = pos[bond.GetBeginAtomIdx()], pos[bond.GetEndAtomIdx()]
        dx, dy = (end - start)[:2]
        # Signed areas on each side of the double-bond axis; tolerate coordinate noise.
        def side(i):
            x, y = (pos[i] - start)[:2]
            return dx * y - dy * x
        a, b = side(left), side(right)
        eps = 1e-6 * (dx * dx + dy * dy)
        same_side = tag in (stereo.STEREOZ, stereo.STEREOCIS)
        if min(abs(a), abs(b)) <= eps or (a * b > 0) != same_side:
            return f"double-bond geometry at bond {bond.GetIdx()}"
    return None


def _layout(mol: Chem.Mol) -> str:
    """Prefer valid reference and clear cage views, then compare geometric quality."""
    best = None
    for name, trial in layout_candidates(mol):
        if layout_problem(mol, trial):
            continue
        if name.startswith("template:"):
            best = (0, 0, name, trial)
            break
        quality = layout_quality(trial)
        # A clear, depth-cued cage keeps the conventional view. If substituents
        # collide, it competes on geometry alongside the general engines.
        clear_bridge = (name.startswith("bridge-") and quality["crossings"] == 1
                        and not any(quality[k] for k in ("clashes", "on_bond", "inside", "label_contacts")))
        rank = 0 if clear_bridge else 1
        score = quality["score"]
        if best is None or rank < best[0] or (rank == best[0] and score < best[1] - LAYOUT_MIN_GAIN):
            best = (rank, score, name, trial)
    if best is None:
        trial = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(trial, forceRDKit=True)
        problem = layout_problem(mol, trial)
        if problem:
            raise InvalidSmiles(f"Could not draw the specified structure faithfully ({problem})")
        best = (1, 0, "rdkit", trial)
    _, _, name, trial = best
    for target, source in zip(mol.GetBonds(), trial.GetBonds()):
        target.ClearProp(GAP_FRACTION)
        if source.HasProp(GAP_FRACTION):
            target.SetDoubleProp(GAP_FRACTION, source.GetDoubleProp(GAP_FRACTION))
    mol.RemoveAllConformers()
    mol.AddConformer(Chem.Conformer(trial.GetConformer()), assignId=True)
    return name


def normalize_layout(mol: Chem.Mol) -> None:
    # Preserve the intentional upright view of a bridge while normalizing its scale.
    has_bridge = any(b.HasProp(GAP_FRACTION) for b in mol.GetBonds())
    rdDepictor.NormalizeDepiction(mol, canonicalize=0 if has_bridge else 1)


def _rdkit_templates(m: Chem.Mol) -> None:
    rdDepictor.Compute2DCoords(m, useRingTemplates=True)


# (RDKit's random flip sampling was tried too: on the benchmark it never beat the
# plain template layout, so it is not worth the time.)
LAYOUT_ENGINES = (  # order = preference on equal scores
    ("coordgen", lambda m: rdCoordGen.AddCoords(m)),
    ("rdkit-templates", _rdkit_templates),
)
LAYOUT_MIN_GAIN = 0.5  # score points; one bond-length spread of 10% costs 0.5


def layout_quality(mol: Chem.Mol) -> dict:
    """Score a 2D layout; lower is better. Scale-free (distances in median bond lengths).

    - crossings: pairs of bonds that cross (not sharing an atom)
    - clashes: atom pairs closer than 0.5 bond length
    - on_bond: atoms closer than 0.3 bond length to a bond they are not part of
    - inside: atoms that are in no ring but sit inside a ring's polygon (a substituent
      drawn into a ring, e.g. CoordGen's N-methyl of tropinone). Bridge atoms inside
      another ring are the normal way to draw bridged systems and are not counted.
    - label_contacts: estimated atom-label overlaps with a nonincident bond
    - spread: coefficient of variation of the bond lengths
    """
    import itertools

    pos = mol.GetConformer().GetPositions()[:, :2]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    if not bonds:
        return {"score": 0.0, "crossings": 0, "clashes": 0, "on_bond": 0, "inside": 0, "label_contacts": 0, "spread": 0.0}
    lengths = [math.dist(pos[a], pos[b]) for a, b in bonds]
    unit = sorted(lengths)[len(lengths) // 2] or 1.0

    def side(a, b, c):
        return (c[0] - a[0]) * (b[1] - a[1]) - (b[0] - a[0]) * (c[1] - a[1])

    crossings = 0
    for (a, b), (c, d) in itertools.combinations(bonds, 2):
        if len({a, b, c, d}) < 4:
            continue
        if side(pos[c], pos[d], pos[a]) * side(pos[c], pos[d], pos[b]) < 0 and \
           side(pos[a], pos[b], pos[c]) * side(pos[a], pos[b], pos[d]) < 0:
            crossings += 1
    n = len(pos)
    clashes = 0
    crowding = 0.0
    for i, j in itertools.combinations(range(n), 2):
        separation = math.dist(pos[i], pos[j]) / unit
        if separation < 0.5:
            clashes += 1
            crowding += 1 - 2 * separation
    # Approximate label clearance at the fixed font/bond scale used on the card.
    # This stays rotation invariant when a reactant turns to face its partner.
    label_radius = []
    for atom in mol.GetAtoms():
        visible = (atom.GetAtomicNum() != 6 or atom.GetFormalCharge() or atom.GetIsotope())
        chars = len(atom.GetSymbol()) + bool(atom.GetTotalNumHs()) + bool(atom.GetFormalCharge())
        label_radius.append((0.52 + 0.20 * (chars - 1)) * FONT_SIZE_PX / BOND_LENGTH_PX if visible else 0)
    label_contacts = 0
    on_bond = 0
    for i in range(n):
        for a, b in bonds:
            if i in (a, b):
                continue
            separation = _point_segment_distance(pos[i], pos[a], pos[b]) / unit
            if separation < label_radius[i]:
                label_contacts += 1
            if separation < 0.3:
                on_bond += 1
    ring_info = mol.GetRingInfo()
    in_ring = {i for ring in ring_info.AtomRings() for i in ring}
    inside = 0
    for i in range(n):
        if i in in_ring:
            continue
        if any(_inside_polygon(pos[i], [pos[j] for j in ring]) for ring in ring_info.AtomRings()):
            inside += 1
    mean = sum(lengths) / len(lengths)
    spread = (sum((x - mean) ** 2 for x in lengths) / len(lengths)) ** 0.5 / mean if mean else 0.0
    # The intended crossover has a visible gap in the final SVG. Other crossings
    # and near-contacts on either bond still count.
    marked = sum(b.HasProp(GAP_FRACTION) for b in mol.GetBonds())
    score = (10 * max(0, crossings - marked) + 8 * (clashes + crowding) + 8 * inside
             + 6 * on_bond + 4 * label_contacts + 5 * spread)
    return {"score": round(score, 4), "crossings": crossings, "clashes": clashes, "on_bond": on_bond,
            "inside": inside, "label_contacts": label_contacts, "spread": round(spread, 4)}


def _inside_polygon(p, poly) -> bool:
    x, y = p
    inside = False
    for k in range(len(poly)):
        (x1, y1), (x2, y2) = poly[k], poly[(k + 1) % len(poly)]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            inside = not inside
    return inside


def _point_segment_distance(p, a, b) -> float:
    dx, dy = b[0] - a[0], b[1] - a[1]
    L2 = dx * dx + dy * dy
    if L2 == 0:
        return math.dist(p, a)
    t = max(0.0, min(1.0, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / L2))
    return math.dist(p, (a[0] + t * dx, a[1] + t * dy))


def _drawer_with(
    mol: Chem.Mol,
    atom_colors: dict[int, str] | None = None,
    bond_colors: dict[int, str] | None = None,
) -> rdMolDraw2D.MolDraw2DSVG:
    drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)  # canvas sized to the molecule
    opts = drawer.drawOptions()
    opts.fixedFontSize = FONT_SIZE_PX
    opts.bondLineWidth = LINE_WIDTH_PX
    opts.scaleBondWidth = False
    opts.padding = PADDING
    opts.clearBackground = False
    # RDKit otherwise adds an H at a back bridgehead, directly over the cage.
    # Keep input-explicit Hs and let the existing bonds carry the stereo wedges.
    if any(b.HasProp(GAP_FRACTION) for b in mol.GetBonds()):
        mol = rdMolDraw2D.PrepareMolForDrawing(mol, addChiralHs=False)
        opts.prepareMolsBeforeDrawing = False
        draw = drawer.DrawMolecule
    else:
        draw = partial(rdMolDraw2D.PrepareAndDrawMolecule, drawer)
    if atom_colors or bond_colors:
        opts.highlightRadius = 0.38                   # in bond-length-ish units: a clear halo
        opts.highlightBondWidthMultiplier = 14        # wide marker stroke behind highlighted bonds
        opts.standardColoursForHighlightedAtoms = True  # keep O red, N blue inside a highlight
        opts.continuousHighlight = True
        draw(
            mol,
            highlightAtoms=list(atom_colors or {}),
            highlightBonds=list(bond_colors or {}),
            highlightAtomColors={i: HIGHLIGHT_RGB[c] for i, c in (atom_colors or {}).items()},
            highlightBondColors={i: HIGHLIGHT_RGB[c] for i, c in (bond_colors or {}).items()},
        )
    else:
        draw(mol)
    drawer.FinishDrawing()
    return drawer


def _median_bond_px(drawer: rdMolDraw2D.MolDraw2DSVG, mol: Chem.Mol) -> float | None:
    def xy(i: int) -> tuple[float, float]:
        p = drawer.GetDrawCoords(i)
        return (p.x, p.y)

    lengths = sorted(math.dist(xy(b.GetBeginAtomIdx()), xy(b.GetEndAtomIdx())) for b in mol.GetBonds())
    return lengths[len(lengths) // 2] if lengths else None
