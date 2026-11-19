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

from chemcard.chem.render import _draw, _layout, layout_candidates, layout_problem, layout_quality, mol_from_smiles, render_smiles

MOLECULES = {
    "caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "penicillin G": "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
    "tropinone": "CN1C2CCC1CC(=O)C2",
    "cocaine": "COC(=O)[C@H]1[C@@H](OC(=O)c2ccccc2)C[C@@H]2CC[C@H]1N2C",
    "cis-cyclooctene": r"C1/C=C\CCCCC1",
    "trans-cyclooctene": "C1/C=C/CCCCC1",
    "camphor": "CC1(C)C2CCC1(C)C(=O)C2",
    "adamantane": "C1C2CC3CC1CC(C2)C3",
    "cubane": "C12C3C4C1C5C2C3C45",
    "quinine": "COc1ccc2nccc([C@@H](O)[C@@H]3C[C@@H]4CCN3C[C@@H]4C=C)c2c1",
    "artemisinin": "C[C@@H]1CC[C@H]2[C@@H](C)C(=O)O[C@@H]3O[C@@]4(C)CC[C@@H]1[C@]32OO4",
    "strychnine": "C1CN2CC3=CCO[C@H]4CC(=O)N5[C@H]6[C@H]4[C@H]3C[C@H]2[C@@]61C7=CC=CC=C75",
    "cholesterol": "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
    "dextromethorphan": "COc1ccc2c(c1)[C@]13CCCC[C@@H]1[C@H](C2)N(C)CC3",
    "epibatidine": "Clc1ccc(cn1)[C@@H]1C[C@H]2CC[C@@H]1N2",
    "18-crown-6": "C1COCCOCCOCCOCCOCCO1",
}


def main() -> None:
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else None
    if out:
        out.mkdir(parents=True, exist_ok=True)
    print("molecule           candidate                    score    crossings/clashes/on-bond/inside")
    for label, smi in MOLECULES.items():
        card = render_smiles(smi)
        mol = mol_from_smiles(card.canonical_smiles)
        chosen = _layout(Chem.Mol(mol))
        if out:
            (out / f"{label}__selected.svg").write_text(card.svg)
        for name, trial in layout_candidates(mol):
            problem = layout_problem(mol, trial)
            marker = "*" if name == chosen else " "
            if problem == "non-finite coordinates" or not trial.GetNumConformers():
                print(f"{label:18s} {name:28s}  {problem}")
                continue
            q = layout_quality(trial)
            counts = "/".join(str(q[key]) for key in ("crossings", "clashes", "on_bond", "inside"))
            print(f"{label:18s} {name:28s}{marker} {q['score']:6.2f}   {counts}   {problem or ''}")
            if out:
                rdDepictor.NormalizeDepiction(trial)
                svg = _draw(trial, {}, {}, layout=False)[0]
                (out / f"{label}__{name.replace(':', '-')}.svg").write_text(svg)


if __name__ == "__main__":
    main()
