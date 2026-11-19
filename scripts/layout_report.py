"""Measure how well chemcard lays out a broad set of molecules (benchmarks/molecules.csv).

    uv run python scripts/layout_report.py                 # summary in the terminal
    uv run python scripts/layout_report.py --md benchmarks/REPORT.md --sheet /tmp/sheet.html

For every molecule the full card pipeline runs (layout choice, drawing), then the
final layout is scored with render.layout_quality. A molecule is "clean" when it
has no unmarked bond crossing, no clashing atoms, no atom on a bond and no
substituent inside a ring. Label contacts (an atom label estimated to touch a
bond) are reported separately: they are frequent and usually cosmetic.

--sheet writes an HTML page with the final card of every flagged molecule, for
looking at them rather than trusting the numbers.
"""

from __future__ import annotations

import argparse
import csv
import html
import time
from collections import Counter, defaultdict
from pathlib import Path

from rdkit import Chem

from chemcard.chem.render import (
    _layout, layout_candidates, layout_problem, layout_quality, mol_from_smiles, render_smiles,
)

ROOT = Path(__file__).resolve().parent.parent
FLAGS = ("crossing", "clash", "on_bond", "inside")


def evaluate(row: dict) -> dict:
    out = {**row, "engine": "", "flags": [], "label_contacts": 0, "score": 0.0, "ms": 0.0, "error": ""}
    t0 = time.perf_counter()
    try:
        card = render_smiles(row["smiles"])
    except Exception as e:  # noqa: BLE001 - report every failure, whatever its type
        out["error"] = f"{type(e).__name__}: {e}"[:120]
        return out
    out["ms"] = (time.perf_counter() - t0) * 1000
    out["svg"] = card.svg

    m = mol_from_smiles(card.canonical_smiles)  # same atom order as the card
    out["engine"] = _layout(m)
    q = layout_quality(m)
    marked = sum(b.HasProp("_chemcardBridgeGap") for b in m.GetBonds())
    counts = {
        "crossing": max(0, q["crossings"] - marked),
        "clash": q["clashes"],
        "on_bond": q["on_bond"],
        "inside": q["inside"],
    }
    out["counts"] = counts
    out["flags"] = [f for f in FLAGS if counts[f]]
    if out["flags"]:
        out["best_candidate"] = best_candidate(m)
    out["label_contacts"] = q.get("label_contacts", 0)
    out["score"] = q["score"]
    return out


def problem_counts(m: Chem.Mol) -> dict:
    q = layout_quality(m)
    marked = sum(b.HasProp("_chemcardBridgeGap") for b in m.GetBonds())
    return {"crossing": max(0, q["crossings"] - marked), "clash": q["clashes"], "on_bond": q["on_bond"],
            "inside": q["inside"], "score": q["score"]}


def best_candidate(m: Chem.Mol) -> str:
    """Would any candidate the layout code already generates have been clean?

    Tells a selection mistake (a clean candidate existed) from a coverage gap
    (no candidate is clean: a new generator or template is needed).
    """
    best = None
    for name, trial in layout_candidates(Chem.Mol(m)):
        if layout_problem(m, trial):
            continue
        c = problem_counts(trial)
        problems = sum(c[f] for f in FLAGS)
        if best is None or (problems, c["score"]) < best[0]:
            best = ((problems, c["score"]), name, c)
    if best is None:
        return "no valid candidate"
    (problems, _), name, c = best
    if problems == 0:
        return f"CLEAN: {name}"
    return f"none clean (best {name}: " + ", ".join(f"{f}×{c[f]}" for f in FLAGS if c[f]) + ")"


def summarise(results: list[dict]) -> list[str]:
    lines = []
    by_cat = defaultdict(list)
    for r in results:
        by_cat[r["category"]].append(r)
    header = "| category | n | clean | crossing | clash | on bond | inside | label contacts | errors |"
    lines += [header, "|" + "---|" * 9]

    def row(name, rs):
        n = len(rs)
        clean = sum(1 for r in rs if not r["flags"] and not r["error"])
        c = Counter(f for r in rs for f in r["flags"])
        lab = sum(1 for r in rs if r["label_contacts"])
        err = sum(1 for r in rs if r["error"])
        return (f"| {name} | {n} | {clean} ({100 * clean / n:.0f}%) | {c['crossing']} | {c['clash']} | "
                f"{c['on_bond']} | {c['inside']} | {lab} | {err} |")

    for cat, rs in by_cat.items():
        lines.append(row(cat, rs))
    lines.append(row("**all**", results))

    engines = Counter(r["engine"].split(":")[0] if r["engine"].startswith("template") else r["engine"]
                      for r in results if not r["error"])
    lines += ["", "Layouts chosen: " + ", ".join(f"{k} {v}" for k, v in engines.most_common())]
    times = sorted(r["ms"] for r in results if not r["error"])
    if times:
        lines.append(f"Render time per molecule: median {times[len(times) // 2]:.0f} ms, "
                     f"slowest {times[-1]:.0f} ms")

    flagged = sorted((r for r in results if r["flags"] or r["error"]), key=lambda r: -r["score"])
    lines += ["", "### Flagged molecules (worst first)", "",
              "| molecule | category | heavy atoms | layout | problems | best existing candidate |",
              "|---|---|---|---|---|---|"]
    for r in flagged:
        probs = r["error"] or ", ".join(f"{f} ×{r['counts'][f]}" for f in r["flags"])
        lines.append(f"| {r['name']} (CID {r['cid']}) | {r['category']} | {r['heavy_atoms']} | {r['engine']} | "
                     f"{probs} | {r.get('best_candidate', '')} |")
    missed = [r for r in flagged if r.get("best_candidate", "").startswith("CLEAN")]
    lines += ["", f"Selection mistakes (a clean candidate existed but was not chosen): {len(missed)}"
              + (" — " + ", ".join(r["name"] for r in missed) if missed else "")]
    return lines


def write_sheet(results: list[dict], path: Path) -> None:
    cards = []
    for r in sorted((r for r in results if r["flags"] and r.get("svg")), key=lambda r: -r["score"]):
        probs = ", ".join(f"{f}×{r['counts'][f]}" for f in r["flags"])
        img = "data:image/svg+xml;charset=utf-8," + html.escape(_quote(r["svg"]))
        cards.append(f"<figure><img src=\"{img}\"><figcaption><b>{html.escape(r['name'])}</b> · "
                     f"{html.escape(r['engine'])}<br>{html.escape(probs)}</figcaption></figure>")
    path.write_text(
        "<!doctype html><meta charset='utf-8'><title>chemcard flagged layouts</title>"
        "<style>body{font:13px system-ui;margin:16px;background:#fafaf9}"
        "main{display:grid;grid-template-columns:repeat(auto-fill,minmax(260px,1fr));gap:12px}"
        "figure{margin:0;background:#fff;border-radius:8px;padding:8px;text-align:center}"
        "img{max-width:100%;max-height:220px}figcaption{color:#555;margin-top:4px}</style>"
        f"<main>{''.join(cards)}</main>",
        encoding="utf-8",
    )


def _quote(svg: str) -> str:
    from urllib.parse import quote

    return quote(svg)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=str(ROOT / "benchmarks" / "molecules.csv"))
    ap.add_argument("--md", help="also write the summary as Markdown to this file")
    ap.add_argument("--sheet", help="write an HTML page with the flagged molecules' cards")
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.csv, encoding="utf-8")))
    results = [evaluate(r) for r in rows]
    lines = summarise(results)
    print("\n".join(lines))
    if args.md:
        Path(args.md).write_text(
            "# Layout report\n\nGenerated by `scripts/layout_report.py` on `benchmarks/molecules.csv` "
            f"({len(rows)} molecules, structures from PubChem).\n\n" + "\n".join(lines) + "\n",
            encoding="utf-8",
        )
    if args.sheet:
        write_sheet(results, Path(args.sheet))


if __name__ == "__main__":
    main()
