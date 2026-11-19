"""Compare 2D layout engines on molecules that are hard to draw.

    uv run python scripts/layout_benchmark.py [out_dir]

Prints, per molecule, each engine's quality score (crossings, clashes, atoms on
bonds, bond-length spread) and which layout chemcard picks. With out_dir, also
writes one SVG per molecule and engine for eyeballing.
"""

from __future__ import annotations

import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

from chemcard.chem.render import LAYOUT_ENGINES, _layout, layout_quality

MOLECULES = {
    "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "penicillin G": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
    "tropinone": "CN1C2CCC1CC(=O)C2",
    "cocaine": "COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C",
    "camphor": "CC1(C)C2CCC1(C)C(=O)C2",
    "adamantane": "C1C2CC3CC1CC(C2)C3",
    "cubane": "C12C3C4C1C5C2C3C45",
    "quinine": "COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
    "artemisinin": "C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]32OO4",
    "strychnine": "O=C1C[C@@H]2OCC=C3CN4CC[C@]56[C@H]4C[C@H]3[C@H]2[C@H]6N1c1ccccc15",
    "cholesterol": "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
    "dextromethorphan": "COc1ccc2c(c1)[C@]13CCCC[C@@H]1[C@H](C2)N(C)CC3",
    "epibatidine": "Clc1ccc(cn1)[C@@H]1C[C@H]2CC[C@@H]1N2",
    "18-crown-6": "C1COCCOCCOCCOCCOCCO1",
}


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if out:
        out.mkdir(parents=True, exist_ok=True)
    names = [n for n, _ in LAYOUT_ENGINES]
    print(f"{'molecule':18s}" + "".join(f"{n:>21s}" for n in names) + "   chosen")
    for label, smi in MOLECULES.items():
        cells = []
        for name, engine in LAYOUT_ENGINES:
            m = Chem.MolFromSmiles(smi)
            try:
                engine(m)
                q = layout_quality(m)
                cells.append(f"{q['score']:6.2f} ({q['crossings']}x/{q['clashes']}c/{q['on_bond']}b/{q['inside']}i)")
            except Exception:
                cells.append("failed")
                continue
            if out:
                rdDepictor.NormalizeDepiction(m)
                d = rdMolDraw2D.MolDraw2DSVG(320, 240)
                rdMolDraw2D.PrepareAndDrawMolecule(d, m)
                d.FinishDrawing()
                (out / f"{label}__{name}.svg").write_text(d.GetDrawingText())
        chosen = _layout(Chem.MolFromSmiles(smi))
        print(f"{label:18s}" + "".join(f"{c:>21s}" for c in cells) + f"   {chosen}")


if __name__ == "__main__":
    main()
