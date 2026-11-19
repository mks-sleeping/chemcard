"""Checks for one reaction / mechanism step (A + B → C). Pure chemistry: no `mcp`."""

from __future__ import annotations

from collections import Counter


from chemcard.chem.render import mol_from_smiles


def _composition(smiles: str) -> tuple[Counter[str], int]:
    mol = mol_from_smiles(smiles)
    atoms: Counter[str] = Counter()
    charge = 0
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 0:  # R groups are not counted
            continue
        atoms[a.GetSymbol()] += 1
        atoms["H"] += a.GetTotalNumHs()
        charge += a.GetFormalCharge()
    return atoms, charge


def balance_note(reactants: list[str], products: list[str]) -> str | None:
    """Describe any atom or charge imbalance between the two sides, or None if balanced.

    Mechanism steps often omit H2O / H+, so an imbalance is reported to the
    model as information, not rejected. Inputs must already be valid SMILES.
    """
    left, right = Counter(), Counter()
    q_left = q_right = 0
    for s in reactants:
        c, q = _composition(s)
        left.update(c)
        q_left += q
    for s in products:
        c, q = _composition(s)
        right.update(c)
        q_right += q

    extra_left = left - right
    extra_right = right - left
    parts = []
    if extra_left:
        parts.append("reactant side has extra " + _fmt(extra_left))
    if extra_right:
        parts.append("product side has extra " + _fmt(extra_right))
    if q_left != q_right:
        parts.append(f"total charge {q_left:+d} → {q_right:+d}")
    return "; ".join(parts) if parts else None


def _fmt(c: Counter[str]) -> str:
    order = sorted(c, key=lambda e: (e not in ("C", "H"), e != "C", e))  # Hill-like: C, H, then alphabetical
    return " ".join(f"{e}{c[e] if c[e] > 1 else ''}" for e in order)
