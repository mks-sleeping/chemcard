"""scripts/layout_report.py: what the benchmark report can see and what it claims."""

import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def _load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


report = _load("layout_report")
ROWS = {r["name"]: r for r in csv.DictReader(open(ROOT / "benchmarks" / "molecules.csv", encoding="utf-8"))}


# ---- 1. overlaps that only exist in the final drawing ------------------------------

def test_render_overlap_in_artemisinin_is_flagged(tmp_path):
    # RDKit adds an explicit H on a stereocentre while drawing; it landed 3.7 px from the
    # peroxide O. The coordinates never contain that H, so layout_quality saw nothing.
    r = report.evaluate(ROWS["artemisinin"])
    assert all(r["counts"][f] == 0 for f in report.FLAGS)          # coordinates look clean...
    assert "render_overlap" in r["flags"]                          # ...the picture is not
    labels = [i for i in r["render_issues"] if i["kind"] == "label/label"]
    assert labels and min(i["distance"] for i in labels) < 5
    sheet = tmp_path / "sheet.html"
    report.write_sheet([r], sheet)
    assert "artemisinin" in sheet.read_text()


def test_svg_overlap_check_is_quiet_on_ordinary_cards():
    from chemcard.chem.render import render_smiles, svg_label_overlaps

    for name in ("caffeine", "aspirin", "ibuprofen", "glucose"):
        assert svg_label_overlaps(render_smiles(ROWS[name]["smiles"]).svg) == []


# ---- 2. the candidate comparison and its honest name --------------------------------

def test_pre_fix_buprenorphine_chose_the_template_over_fewer_crossings(monkeypatch):
    # Before the template rule, the epoxymorphinan template won outright for buprenorphine
    # (6 crossings) although CoordGen had 4 fewer. Only the crossing counts are asserted
    # here: they are the same on Ubuntu and macOS, while CoordGen's clash count is not.
    from chemcard.chem import render

    monkeypatch.setattr(render, "_template_crossings", lambda name: 99)  # the old "template always wins"
    r = report.evaluate(ROWS["buprenorphine"])
    assert r["engine"] == "template:epoxymorphinan" and r["counts"]["crossing"] == 6
    coordgen = next(a for a in r["alternatives"] if a["name"] == "coordgen")
    assert coordgen["delta"]["crossing"] == -4


def flagged(name, counts, alternatives):
    """A flagged result with fixed candidate data, as evaluate() would return it."""
    alts = []
    for alt, alt_counts in alternatives.items():
        alt_counts = {f: alt_counts.get(f, 0) for f in report.FLAGS}
        alts.append({"name": alt, "counts": alt_counts, "delta": {f: alt_counts[f] - counts.get(f, 0) for f in report.FLAGS}})
    counts = {f: counts.get(f, 0) for f in report.FLAGS} | {"render_overlap": 0}
    return {"name": name, "cid": name, "category": "x", "heavy_atoms": "30", "engine": "template:t", "error": "",
            "flags": [f for f in report.FLAGS if counts[f]], "label_contacts": 0, "ms": 1.0, "score": 10.0,
            "counts": counts, "render_issues": [], "alternatives": alts}


def test_trade_off_text_and_count_from_fixed_candidates():
    # the buprenorphine situation, with fixed numbers: CoordGen trades 4 crossings for 1 clash
    r = flagged("bup", {"crossing": 6}, {"coordgen": {"crossing": 2, "clash": 1}})
    assert not report.missed_clean(r) and report.pareto_better(r) == [] and report.trade_offs(r) == ["coordgen"]
    text = "\n".join(report.summarise([r]))
    assert "coordgen: crossing −4, clash +1" in text
    assert "Missed fully clean candidates" in text and "Selection mistakes" not in text
    assert "improves some count at the cost of another: **1**" in text


def test_trade_offs_are_judged_per_candidate():
    # one candidate better on one count and worse on none, another better on one and worse on
    # another: the molecule belongs to both lists
    both = flagged("both", {"crossing": 2, "clash": 1},
                   {"x": {"crossing": 1, "clash": 1}, "y": {"crossing": 0, "clash": 2}})
    assert report.pareto_better(both) == ["x"] and report.trade_offs(both) == ["y"]
    # only a Pareto-better candidate: no trade-off (the old rule counted any improvement)
    pareto_only = flagged("pareto", {"crossing": 2}, {"x": {"crossing": 1}})
    assert report.pareto_better(pareto_only) == ["x"] and report.trade_offs(pareto_only) == []
    text = "\n".join(report.summarise([both, pareto_only]))
    assert "Pareto-better candidate (no count worse, one better): **2**" in text
    assert "improves some count at the cost of another: **1**" in text


def test_pareto_better_is_reported():
    r = {"counts": {"crossing": 2, "clash": 1, "on_bond": 0, "inside": 0},
         "alternatives": [{"name": "x", "counts": {}, "delta": {"crossing": -1, "clash": 0, "on_bond": 0, "inside": 0}},
                          {"name": "y", "counts": {}, "delta": {"crossing": -2, "clash": 1, "on_bond": 0, "inside": 0}}]}
    assert report.pareto_better(r) == ["x"]


# ---- 4. one compound, one vote --------------------------------------------------------

def test_overall_rate_counts_each_cid_once():
    assert ROWS["eucalyptol"]["cid"] == ROWS["1,8-cineole"]["cid"]
    fake = [
        {"name": "a", "cid": "1", "category": "x", "heavy_atoms": "5", "flags": [], "error": "", "label_contacts": 0, "engine": "coordgen",
         "ms": 1.0, "score": 0.0, "counts": {}, "render_issues": [], "alternatives": []},
        {"name": "a again", "cid": "1", "category": "y", "heavy_atoms": "5", "flags": [], "error": "", "label_contacts": 0,
         "engine": "coordgen", "ms": 1.0, "score": 0.0, "counts": {}, "render_issues": [], "alternatives": []},
        {"name": "b", "cid": "2", "category": "x", "heavy_atoms": "5", "flags": ["clash"], "error": "", "label_contacts": 0,
         "engine": "coordgen", "ms": 1.0, "score": 8.0, "counts": {"crossing": 0, "clash": 1, "on_bond": 0, "inside": 0},
         "render_issues": [], "alternatives": []},
    ]
    text = "\n".join(report.summarise(fake))
    assert "| all samples | 3 | 2 (67%)" in text
    assert "| **all, unique compounds (by CID)** | 2 | 1 (50%)" in text
    assert "a = a again (CID 1)" in text


# ---- the per-[a.b.c] table applies the layout chooser's stereo check first ------------

def test_bicyclo_candidates_rejected_for_double_bond_geometry_are_reported_as_such():
    # trans-cyclooctene fused to cyclohexane: the flat fused construction draws the 8-ring
    # as a regular octagon, which makes the double bond cis. The chooser discards it before
    # scoring, so the report must say so rather than "CoordGen scored better".
    from chemcard.chem.render import layout_problem

    smiles = "C1/C=C/CCC2CCCCC2C1"
    r = report.evaluate({"name": "trans-bicyclo[6.4.0]dodec-2-ene", "category": "test", "cid": "0",
                         "heavy_atoms": "12", "smiles": smiles})
    [row] = r["bicyclo"]
    assert row["abc"] == "6.4.0" and not row["adopted"]
    assert row["rejected"] and all(p.startswith("double-bond geometry") for p in row["rejected"].values())
    assert row["reason"] == "rejected: double-bond geometry"
    text = "\n".join(report.summarise([r]))
    assert "rejected: double-bond geometry (1: trans-bicyclo[6.4.0]dodec-2-ene)" in text
    # the final drawing itself keeps the trans bond
    from chemcard.chem.render import _layout, mol_from_smiles, render_smiles
    m = mol_from_smiles(render_smiles(smiles).canonical_smiles)
    assert not _layout(m).endswith("-fused1") and layout_problem(m, m) is None
