"""Measure how well chemcard lays out a broad set of molecules (benchmarks/molecules.csv).

    uv run python scripts/layout_report.py                 # summary in the terminal
    uv run python scripts/layout_report.py --md benchmarks/REPORT.md --sheet /tmp/sheet.html

For every molecule the full card pipeline runs (layout choice, drawing). Two kinds
of problems are counted:

- layout problems, from the chosen coordinates (render.layout_quality): an unmarked
  bond crossing, clashing atoms, an atom on a bond, a substituent inside a ring;
- render overlaps, from the final SVG (render.svg_label_overlaps): labels on labels
  or a bond through a label. These include atoms RDKit adds while drawing (an
  explicit H on a stereocentre), which the coordinates never see.

A molecule is "clean" when it has neither. Label contacts (an atom label estimated
from the coordinates to touch a bond) are reported separately.

For each flagged molecule the report also lays out every other candidate the code
generates and shows the trade-off against the chosen one ("coordgen: crossing −4,
clash +1"). Its summary counts only what those numbers can support: candidates
with all four layout counts at zero that were not chosen, and choices for which
another candidate is no worse on any layout count and better on one. Neither says
the choice was right; a trade-off still needs a look (--sheet).

The CSV may list one compound under two names or categories (eucalyptol and
1,8-cineole are both CID 2758). Category rows count samples; the overall rate
counts each CID once.
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
    svg_label_overlaps,
)

ROOT = Path(__file__).resolve().parent.parent
FLAGS = ("crossing", "clash", "on_bond", "inside")  # from the coordinates
ALL_FLAGS = FLAGS + ("render_overlap",)            # + from the final SVG


def layout_counts(m: Chem.Mol) -> dict:
    q = layout_quality(m)
    marked = sum(b.HasProp("_chemcardBridgeGap") for b in m.GetBonds())
    return {"crossing": max(0, q["crossings"] - marked), "clash": q["clashes"], "on_bond": q["on_bond"],
            "inside": q["inside"], "score": q["score"], "label_contacts": q.get("label_contacts", 0)}


def evaluate(row: dict) -> dict:
    out = {**row, "engine": "", "flags": [], "label_contacts": 0, "score": 0.0, "ms": 0.0, "error": "",
           "counts": {}, "render_issues": [], "alternatives": []}
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
    counts = layout_counts(m)
    out["render_issues"] = svg_label_overlaps(card.svg)
    counts["render_overlap"] = len(out["render_issues"])
    out["counts"] = counts
    out["flags"] = [f for f in ALL_FLAGS if counts[f]]
    out["label_contacts"] = counts["label_contacts"]
    out["score"] = counts["score"]
    if out["flags"]:
        out["alternatives"] = alternatives(m, out["engine"], counts)
    return out


def alternatives(m: Chem.Mol, chosen: str, chosen_counts: dict) -> list[dict]:
    """Every other valid candidate with its layout counts and the difference to the chosen one."""
    alts = []
    for name, trial in layout_candidates(Chem.Mol(m)):
        if name == chosen or layout_problem(m, trial):
            continue
        c = layout_counts(trial)
        delta = {f: c[f] - chosen_counts[f] for f in FLAGS}
        alts.append({"name": name, "counts": {f: c[f] for f in FLAGS}, "delta": delta})
    alts.sort(key=lambda a: (sum(a["counts"].values()), a["name"]))
    return alts


def missed_clean(r: dict) -> bool:
    """A candidate with all four layout counts at zero existed but another was chosen."""
    return any(r["counts"].get(f) for f in FLAGS) and any(not any(a["counts"].values()) for a in r["alternatives"])


def pareto_better(r: dict) -> list[str]:
    """Candidates no worse on any layout count and better on at least one."""
    return [a["name"] for a in r["alternatives"]
            if all(v <= 0 for v in a["delta"].values()) and any(v < 0 for v in a["delta"].values())]


def trade_text(a: dict) -> str:
    parts = [f"{f} {'+' if v > 0 else '−'}{abs(v)}" for f, v in a["delta"].items() if v]
    return f"{a['name']}: " + (", ".join(parts) if parts else "same counts")


def problems_text(r: dict) -> str:
    if r["error"]:
        return r["error"]
    parts = [f"{f} ×{r['counts'][f]}" for f in FLAGS if r["counts"][f]]
    for issue in r["render_issues"]:
        if issue["kind"] == "label/label":
            parts.append(f"labels {issue['atoms'][0]}/{issue['atoms'][1]} overlap ({issue['distance']} px apart)")
        else:
            parts.append(f"bond {'-'.join(map(str, issue['bond_atoms']))} through label {issue['atoms'][0]}")
    return ", ".join(parts)


def summarise(results: list[dict]) -> list[str]:
    lines = []
    by_cat = defaultdict(list)
    for r in results:
        by_cat[r["category"]].append(r)
    header = "| category | samples | clean | crossing | clash | on bond | inside | render overlap | label contacts | errors |"
    lines += [header, "|" + "---|" * 10]

    def row(name, rs):
        n = len(rs)
        clean = sum(1 for r in rs if not r["flags"] and not r["error"])
        c = Counter(f for r in rs for f in r["flags"])
        lab = sum(1 for r in rs if r["label_contacts"])
        err = sum(1 for r in rs if r["error"])
        return (f"| {name} | {n} | {clean} ({100 * clean / n:.0f}%) | {c['crossing']} | {c['clash']} | "
                f"{c['on_bond']} | {c['inside']} | {c['render_overlap']} | {lab} | {err} |")

    for cat, rs in by_cat.items():
        lines.append(row(cat, rs))
    lines.append(row("all samples", results))
    unique = list({r["cid"]: r for r in results}.values())  # one compound, one vote
    lines.append(row("**all, unique compounds (by CID)**", unique))
    dups = [cid for cid, k in Counter(r["cid"] for r in results).items() if k > 1]
    if dups:
        names = {cid: " = ".join(r["name"] for r in results if r["cid"] == cid) for cid in dups}
        lines += ["", f"{len(results)} samples, {len(unique)} unique compounds. Listed more than once: "
                  + "; ".join(f"{names[c]} (CID {c})" for c in dups) + "."]

    engines = Counter(r["engine"].split(":")[0] if r["engine"].startswith("template") else r["engine"]
                      for r in results if not r["error"])
    lines += ["", "Layouts chosen: " + ", ".join(f"{k} {v}" for k, v in engines.most_common())]
    times = sorted(r["ms"] for r in results if not r["error"])
    if times:
        lines.append(f"Render time per molecule: median {times[len(times) // 2]:.0f} ms, "
                     f"slowest {times[-1]:.0f} ms")

    flagged = sorted((r for r in results if r["flags"] or r["error"]), key=lambda r: -r["score"])
    missed = [r for r in flagged if missed_clean(r)]
    dominated = [(r, pareto_better(r)) for r in flagged if pareto_better(r)]
    trades = [r for r in flagged if any(any(v < 0 for v in a["delta"].values()) for a in r["alternatives"])]
    lines += [
        "", "### What the candidate comparison can and cannot say", "",
        "Only the four layout counts are compared (render overlaps need each candidate drawn).",
        "",
        f"- Missed fully clean candidates (all four layout counts zero, not chosen): **{len(missed)}**"
        + (" — " + ", ".join(r["name"] for r in missed) if missed else ""),
        f"- Choices with a Pareto-better candidate (no count worse, one better): **{len(dominated)}**"
        + (" — " + ", ".join(f"{r['name']} ({', '.join(p)})" for r, p in dominated) if dominated else ""),
        f"- Choices where another candidate improves some count at the cost of another: **{len(trades)}**; "
        "whether the trade is worth it is a judgement the counts cannot make (see the table and --sheet).",
    ]
    lines += ["", "### Flagged molecules (worst first)", "",
              "| molecule | category | heavy atoms | layout | problems | other candidates (Δ vs chosen) |",
              "|---|---|---|---|---|---|"]
    for r in flagged:
        alts = "; ".join(trade_text(a) for a in r["alternatives"][:3]) or "—"
        lines.append(f"| {r['name']} (CID {r['cid']}) | {r['category']} | {r['heavy_atoms']} | {r['engine']} | "
                     f"{problems_text(r)} | {alts} |")
    return lines


def write_sheet(results: list[dict], path: Path) -> None:
    cards = []
    for r in sorted((r for r in results if r["flags"] and r.get("svg")), key=lambda r: -r["score"]):
        img = "data:image/svg+xml;charset=utf-8," + html.escape(_quote(r["svg"]))
        alts = "<br>".join(html.escape(trade_text(a)) for a in r["alternatives"][:2])
        cards.append(f"<figure><img src=\"{img}\"><figcaption><b>{html.escape(r['name'])}</b> · "
                     f"{html.escape(r['engine'])}<br>{html.escape(problems_text(r))}"
                     f"{'<br><i>' + alts + '</i>' if alts else ''}</figcaption></figure>")
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
        n_unique = len({r["cid"] for r in rows})
        Path(args.md).write_text(
            "# Layout report\n\nGenerated by `scripts/layout_report.py` on `benchmarks/molecules.csv` "
            f"({len(rows)} samples, {n_unique} unique compounds; structures from PubChem).\n\n"
            + "\n".join(lines) + "\n",
            encoding="utf-8",
        )
    if args.sheet:
        write_sheet(results, Path(args.sheet))


if __name__ == "__main__":
    main()
