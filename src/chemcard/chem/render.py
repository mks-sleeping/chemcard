"""SMILES → 2D structure SVG.

Pure chemistry core: must not import `mcp`. Everything here is testable
without starting a server.
"""

from __future__ import annotations

import contextlib
import io
import math
import re
from dataclasses import dataclass

from rdkit import Chem, rdBase
from rdkit.Chem import rdCoordGen, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

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


def mol_from_smiles(smiles: str) -> Chem.Mol | None:
    """`Chem.MolFromSmiles` without RDKit log noise (stdio servers log stderr to the host)."""
    with rdBase.BlockLogs():
        return Chem.MolFromSmiles(smiles)


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
        mol = Chem.MolFromSmiles(smiles)
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


def _draw(mol: Chem.Mol, atom_colors: dict[int, str], bond_colors: dict[int, str]) -> tuple[str, int, int, list[tuple[float, float]]]:
    # Flexicanvas ignores `fixedBondLength`: its pixel size follows the
    # coordinates. So lay out, draw once to measure, rescale the coordinates
    # to hit BOND_LENGTH_PX exactly, then draw for real.
    _layout(mol)
    rdDepictor.NormalizeDepiction(mol)  # canonical orientation, uniform bond length
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
    return svg, int(drawer.Width()), int(drawer.Height()), coords


def _layout(mol: Chem.Mol) -> None:
    """Compute 2D coordinates with the engine that draws this molecule best.

    - CoordGen: cleanest for ordinary rings and macrocycles (default).
    - RDKit's own depictor with ring-system templates: bridged bicyclics
      (tropanes, norbornanes, …), where CoordGen pushes substituents into
      the ring. Also the fallback when CoordGen fails (it returns NaN
      coordinates for a lone H+).
    """
    if rdMolDescriptors.CalcNumBridgeheadAtoms(mol) == 0:
        rdCoordGen.AddCoords(mol)
        if mol.GetNumConformers() and not _has_nan(mol):
            return
    rdDepictor.Compute2DCoords(mol, useRingTemplates=True)


def _has_nan(mol: Chem.Mol) -> bool:
    return any(math.isnan(c) for p in mol.GetConformer().GetPositions() for c in p)


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
    if atom_colors or bond_colors:
        opts.highlightRadius = 0.38                   # in bond-length-ish units: a clear halo
        opts.highlightBondWidthMultiplier = 14        # wide marker stroke behind highlighted bonds
        opts.standardColoursForHighlightedAtoms = True  # keep O red, N blue inside a highlight
        opts.continuousHighlight = True
        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            mol,
            highlightAtoms=list(atom_colors or {}),
            highlightBonds=list(bond_colors or {}),
            highlightAtomColors={i: HIGHLIGHT_RGB[c] for i, c in (atom_colors or {}).items()},
            highlightBondColors={i: HIGHLIGHT_RGB[c] for i, c in (bond_colors or {}).items()},
        )
    else:
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer


def _median_bond_px(drawer: rdMolDraw2D.MolDraw2DSVG, mol: Chem.Mol) -> float | None:
    def xy(i: int) -> tuple[float, float]:
        p = drawer.GetDrawCoords(i)
        return (p.x, p.y)

    lengths = sorted(math.dist(xy(b.GetBeginAtomIdx()), xy(b.GetEndAtomIdx())) for b in mol.GetBonds())
    return lengths[len(lengths) // 2] if lengths else None
