"""Chemistry lint: bookkeeping rules RDKit's sanitizer does not enforce.

RDKit accepts some inputs that are almost always typos in a mechanism, e.g.
`C[N+](C)C` (missing H) is silently read as a radical cation. These rules
catch such cases. Errors mean "almost certainly wrong"; warnings mean
"unusual, double-check". Pure chemistry: no `mcp`.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass
from typing import Literal

from rdkit import Chem

from chemcard.chem.render import mol_from_smiles

Level = Literal["error", "warning"]


@dataclass(frozen=True)
class Issue:
    level: Level
    code: str
    message: str
    where: str = ""  # e.g. "step 2 reactants[1]"

    def as_dict(self) -> dict:
        return asdict(self)

    def text(self) -> str:
        return f"[{self.level}] {self.where + ': ' if self.where else ''}{self.message}"


# ---- species rules ----------------------------------------------------------

_PATTERNS = {
    name: Chem.MolFromSmarts(smarts)
    for name, smarts in {
        # carbocations without obvious stabilisation (no aryl, no C=X, no heteroatom lone pair next to it)
        "methyl_cation": "[CH3+]",
        "primary_cation": "[CX3H2+;!$([CH2+]a);!$([CH2+]C=*);!$([CH2+]C#*);!$([CH2+][O,N,S])]",
        # strong bases (implausible as free species under acidic conditions)
        "hydroxide": "[OX1H1-]",
        "alkoxide": "[OX1-][CX4]",
        "enolate": "[OX1-][CX3]=[CX3]",
        "carbanion": "[C-;!$([C-]#*);!$([C-]=*)]",
        "amide_anion": "[NX2-,NX1-;!$([N-]=*);!$([N-]#*);!$([N-]S(=O)=O);!$([N-]C=O)]",
        # strong acids (implausible as free species under basic conditions)
        "oxonium": "[OX3H3+,OX3H2+]",
        "protonated_carbonyl": "[OX2H1+]=[#6]",
    }.items()
}
_STRONG_BASES = ("hydroxide", "alkoxide", "enolate", "carbanion", "amide_anion")
_STRONG_ACIDS = ("oxonium", "protonated_carbonyl")


def lint_species(smiles: str, *, where: str = "", radicals_expected: bool = False,
                 medium: str | None = None) -> list[Issue]:
    """Check one (already valid) species. `medium` is "acidic", "basic" or None."""
    mol = mol_from_smiles(smiles)
    issues: list[Issue] = []

    for atom in mol.GetAtoms():
        n_rad = atom.GetNumRadicalElectrons()
        if n_rad and not radicals_expected:
            fix = _with_added_h(mol, atom.GetIdx(), n_rad)
            hint = f" Did you mean {fix}?" if fix else ""
            issues.append(Issue(
                "error", "unexpected_radical",
                f"{atom.GetSymbol()} (atom {atom.GetIdx() + 1}) has {n_rad} unpaired electron(s) in {smiles}: "
                f"usually a bracket atom missing its H, e.g. write [NH+] not [N+].{hint} "
                "If this step really is radical chemistry, set radical=true.",
                where,
            ))

    if mol.HasSubstructMatch(_PATTERNS["methyl_cation"]):
        issues.append(Issue("warning", "methyl_cation", f"free methyl cation in {smiles}: essentially never formed in solution (think SN2).", where))
    if mol.HasSubstructMatch(_PATTERNS["primary_cation"]):
        issues.append(Issue("warning", "primary_cation", f"unstabilised primary carbocation in {smiles}: rarely a real intermediate; expect a concerted step or a rearrangement.", where))

    net = Chem.GetFormalCharge(mol)
    if abs(net) >= 2:
        issues.append(Issue("warning", "multiply_charged", f"{smiles} has net charge {net:+d}; doubly charged intermediates are unusual.", where))

    if medium == "acidic":
        for name in _STRONG_BASES:
            if mol.HasSubstructMatch(_PATTERNS[name]):
                issues.append(Issue("warning", "base_in_acid", f"{name.replace('_', ' ')} in {smiles} under acidic conditions: use the neutral/protonated form (e.g. enol instead of enolate).", where))
    elif medium == "basic":
        for name in _STRONG_ACIDS:
            if mol.HasSubstructMatch(_PATTERNS[name]):
                issues.append(Issue("warning", "acid_in_base", f"{name.replace('_', ' ')} in {smiles} under basic conditions: use the neutral/deprotonated form.", where))
    return issues


def _with_added_h(mol: Chem.Mol, idx: int, n: int) -> str | None:
    """The SMILES with `n` H added to atom `idx` instead of radical electrons, if valid."""
    rw = Chem.RWMol(mol)
    a = rw.GetAtomWithIdx(idx)
    a.SetNumRadicalElectrons(0)
    a.SetNumExplicitHs(a.GetNumExplicitHs() + n)
    a.SetNoImplicit(True)
    try:
        Chem.SanitizeMol(rw)
        return Chem.MolToSmiles(rw)
    except Exception:
        return None


# ---- reaction medium from the conditions text --------------------------------

_SUP = str.maketrans("⁺⁻₀₁₂₃₄₅₆₇₈₉⁰¹²³⁴⁵⁶⁷⁸⁹", "+-01234567890123456789")
_ACID = re.compile(r"H\+|H3O\+|\bHCl\b|\bHBr\b|\bHI\b|H2SO4|HNO3|H3PO4|TsOH|TfOH|\bTFA\b|\bacid\b|酸性|酸催化", re.I)
_BASE = re.compile(r"OH-|NaOH|KOH|LiOH|NaOEt|NaOMe|KOtBu|tBuOK|t-BuOK|\bLDA\b|\bNaH\b|NaNH2|BuLi|\bbase\b|碱性|碱催化", re.I)
_PH = re.compile(r"pH\s*[≈~=]?\s*(\d+(?:\.\d+)?)(?:\s*[-–~]\s*(\d+(?:\.\d+)?))?", re.I)


def medium_from_conditions(conditions: str | None) -> str | None:
    """'acidic', 'basic' or None (unknown / buffered / mixed)."""
    if not conditions:
        return None
    text = conditions.translate(_SUP)
    m = _PH.search(text)
    if m:
        lo = float(m.group(1))
        hi = float(m.group(2) or lo)
        if hi < 6:
            return "acidic"
        if lo > 8:
            return "basic"
        return None
    acid, base = bool(_ACID.search(text)), bool(_BASE.search(text))
    if acid and not base:
        return "acidic"
    if base and not acid:
        return "basic"
    return None
