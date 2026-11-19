"""Highlights (model → figure) and selections (user → model), in canonical atom indices.

Every figure is drawn from the molecule's canonical SMILES, so "atom i" means
the same atom in the drawing, in a highlight and in a user's selection.

- The model names what to highlight by SMARTS ("C=O") or by atom-map numbers
  written in its SMILES ("[NH2:1]", atoms=[1]).
- The card reports a user's selection as canonical indices; we turn it into
  a fragment SMILES and an atom-mapped SMILES the model can read and reuse.

Pure chemistry: no `mcp`.
"""

from __future__ import annotations

from dataclasses import dataclass

from rdkit import Chem

from chemcard.chem.render import InvalidSmiles, mol_from_smiles, parse_smiles

COLORS = ("red", "blue", "green", "orange", "purple")


class InvalidHighlight(ValueError):
    """A highlight cannot be resolved. The message goes to the model verbatim."""


@dataclass(frozen=True)
class HighlightSpec:
    color: str = "red"
    smarts: str | None = None
    atoms: tuple[int, ...] = ()  # atom-map numbers from the SMILES


@dataclass(frozen=True)
class Resolved:
    canonical_smiles: str
    atom_colors: dict[int, str]  # canonical atom index → colour name
    bond_colors: dict[int, str]  # canonical bond index → colour name


def resolve(smiles: str, specs: list[HighlightSpec]) -> Resolved:
    """Canonicalise `smiles` (dropping atom maps) and map highlights onto canonical indices."""
    mol = parse_smiles(smiles)
    by_map = {a.GetAtomMapNum(): a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum()}
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    canonical = Chem.MolToSmiles(mol)
    order = list(mol.GetPropsAsDict(True, True)["_smilesAtomOutputOrder"])
    to_canon = {orig: pos for pos, orig in enumerate(order)}
    cmol = mol_from_smiles(canonical)

    atom_colors: dict[int, str] = {}
    bond_colors: dict[int, str] = {}
    for spec in specs:
        if spec.color not in COLORS:
            raise InvalidHighlight(f"unknown colour {spec.color!r}; use one of {', '.join(COLORS)}")
        atoms, bonds = _match(mol, spec, by_map, smiles)
        for i in atoms:
            atom_colors[to_canon[i]] = spec.color
        for a, b in bonds:
            bond = cmol.GetBondBetweenAtoms(to_canon[a], to_canon[b])
            if bond is not None:
                bond_colors[bond.GetIdx()] = spec.color
    return Resolved(canonical, atom_colors, bond_colors)


def _match(mol: Chem.Mol, spec: HighlightSpec, by_map: dict[int, int], smiles: str) -> tuple[set[int], set[tuple[int, int]]]:
    if bool(spec.smarts) == bool(spec.atoms):
        raise InvalidHighlight("each highlight needs exactly one of `smarts` or `atoms`")
    if spec.smarts:
        patt = Chem.MolFromSmarts(spec.smarts)
        if patt is None:
            raise InvalidHighlight(f"invalid SMARTS {spec.smarts!r}")
        matches = mol.GetSubstructMatches(patt)
        if not matches:
            raise InvalidHighlight(f"SMARTS {spec.smarts!r} matches nothing in {Chem.MolToSmiles(mol)}")
        atoms = {i for m in matches for i in m}
        bonds = {(m[b.GetBeginAtomIdx()], m[b.GetEndAtomIdx()]) for m in matches for b in patt.GetBonds()}
        return atoms, bonds
    missing = [n for n in spec.atoms if n not in by_map]
    if missing:
        raise InvalidHighlight(
            f"atom-map number(s) {missing} not found in {smiles!r}; write them in the SMILES, e.g. [NH2:1]"
        )
    atoms = {by_map[n] for n in spec.atoms}
    bonds = {(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
             if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms}
    return atoms, bonds


def describe_selection(canonical_smiles: str, atoms: list[int]) -> dict:
    """Turn a user's selection (canonical atom indices) into text the model can use."""
    mol = parse_smiles(canonical_smiles)
    n = mol.GetNumAtoms()
    atoms = list(dict.fromkeys(atoms))  # keep click order, drop repeats
    bad = [i for i in atoms if not 0 <= i < n]
    if not atoms or bad:
        raise InvalidSmiles(f"selection {atoms} is not a set of atom indices 0..{n - 1}")
    fragment = Chem.MolFragmentToSmiles(mol, atomsToUse=atoms)
    mapped = Chem.Mol(mol)
    for num, i in enumerate(atoms, start=1):
        mapped.GetAtomWithIdx(i).SetAtomMapNum(num)
    return {
        "fragment_smiles": fragment,
        "mapped_smiles": Chem.MolToSmiles(mapped),
        "atoms": [
            {"index": i, "map": num, "element": _symbol(mol.GetAtomWithIdx(i))}
            for num, i in enumerate(atoms, start=1)
        ],
    }


def _symbol(atom: Chem.Atom) -> str:
    return "R" if atom.GetAtomicNum() == 0 else atom.GetSymbol()
