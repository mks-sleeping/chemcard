"""Whole-mechanism check: per-species lint, per-step balance, step-to-step continuity.

Pure chemistry: no `mcp`. Used by the `check_mechanism` tool so the model can
validate all steps at once before drawing them one by one.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from chemcard.chem.lint import Issue, lint_species, medium_from_conditions
from chemcard.chem.reaction import balance_note
from chemcard.chem.render import InvalidSmiles, mol_from_smiles, parse_smiles


@dataclass
class Step:
    reactants: list[str]
    products: list[str]
    conditions: str | None = None
    radical: bool = False


@dataclass
class Report:
    issues: list[Issue] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not any(i.level == "error" for i in self.issues)

    def text(self) -> str:
        if not self.issues:
            return "Mechanism check passed: all species valid, steps balanced and connected."
        head = "Mechanism check: " + ("no errors" if self.ok else "ERRORS found — fix them before drawing")
        ordered = sorted(self.issues, key=lambda i: i.level != "error")  # errors first
        return head + "\n" + "\n".join(i.text() for i in ordered)

    def as_dict(self) -> dict:
        return {"ok": self.ok, "issues": [i.as_dict() for i in self.issues]}


def lint_step(step: Step, *, label: str = "") -> list[Issue]:
    """Lint every species of one step and report imbalance. Raises InvalidSmiles."""
    prefix = f"{label} " if label else ""
    medium = medium_from_conditions(step.conditions)
    issues: list[Issue] = []
    for side in ("reactants", "products"):
        for i, smi in enumerate(getattr(step, side), start=1):
            where = f"{prefix}{side}[{i}]"
            try:
                parse_smiles(smi)
            except InvalidSmiles as e:
                raise InvalidSmiles(f"{where}: {e}") from None
            issues += lint_species(smi, where=where, radicals_expected=step.radical, medium=medium)
    note = balance_note(step.reactants, step.products)
    if note:
        issues.append(Issue("warning", "unbalanced", f"{note} (fine if byproducts such as H2O / H+ are omitted on purpose)", label))
    return issues


def check_mechanism(steps: list[Step]) -> Report:
    """Check all steps; continuity means step n+1 consumes at least one product of step n."""
    report = Report()
    for n, step in enumerate(steps, start=1):
        try:
            report.issues += lint_step(step, label=f"step {n}")
        except InvalidSmiles as e:
            report.issues.append(Issue("error", "invalid_smiles", str(e), f"step {n}"))
            return report  # later checks need valid SMILES

    for n in range(1, len(steps)):
        report.issues += _continuity(steps[n - 1], steps[n], n)
    return report


def _continuity(prev: Step, cur: Step, n: int) -> list[Issue]:
    where = f"step {n} → step {n + 1}"
    made = {_key(s): s for s in prev.products}
    used = {_key(s): s for s in cur.reactants}
    if made.keys() & used.keys():
        return []

    made_flat = {_key(s, stereo=False): s for s in prev.products}
    for k, s in used.items():
        if _key(s, stereo=False) in made_flat:
            return [Issue("warning", "stereo_mismatch",
                          f"{s} in step {n + 1} matches step {n} product {made_flat[_key(s, stereo=False)]} "
                          "only when stereochemistry is ignored; check the stereo descriptors.", where)]

    by_formula = {_formula(s): s for s in prev.products}
    for s in cur.reactants:
        if _formula(s) in by_formula:
            return [Issue("error", "broken_chain",
                          f"step {n + 1} starts from {s}, but step {n} produced {by_formula[_formula(s)]} "
                          f"(same formula {_formula(s)}, different structure). The two steps do not connect; "
                          "draw the same intermediate in both.", where)]
    return [Issue("warning", "no_continuity",
                  f"no reactant of step {n + 1} ({', '.join(cur.reactants)}) is a product of step {n} "
                  f"({', '.join(prev.products)}). Each step should start from the previous step's intermediate.", where)]


def _key(smiles: str, *, stereo: bool = True) -> str:
    mol = mol_from_smiles(smiles)
    for a in mol.GetAtoms():  # atom maps are annotations, not structure
        a.SetAtomMapNum(0)
    return Chem.MolToSmiles(mol, isomericSmiles=stereo)


def _formula(smiles: str) -> str:  # atom maps do not affect the formula
    return rdMolDescriptors.CalcMolFormula(mol_from_smiles(smiles))
