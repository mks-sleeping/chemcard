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

def test_pre_fix_buprenorphine_shows_the_trade_off(monkeypatch):
    # Before the template rule, the epoxymorphinan template won outright for buprenorphine
    # (6 crossings) although CoordGen had 2. No candidate was fully clean, so "missed clean
    # candidates" stays 0; the report must still show what the alternative would trade.
    from chemcard.chem import render

    monkeypatch.setattr(render, "_template_crossings", lambda name: 99)  # the old "template always wins"
    r = report.evaluate(ROWS["buprenorphine"])
    assert r["engine"] == "template:epoxymorphinan" and r["counts"]["crossing"] == 6
    coordgen = next(a for a in r["alternatives"] if a["name"] == "coordgen")
    assert coordgen["delta"]["crossing"] == -4 and coordgen["delta"]["clash"] == 1
    assert not report.missed_clean(r)
    text = "\n".join(report.summarise([r]))
    assert "coordgen: crossing −4, clash +1" in text
    assert "Missed fully clean candidates" in text and "Selection mistakes" not in text
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
