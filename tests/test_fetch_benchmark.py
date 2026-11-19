"""scripts/fetch_benchmark_molecules.py must never replace the CSV with a partial one."""

import importlib.util
import io
import json
import urllib.error
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("fetch_bm", ROOT / "scripts" / "fetch_benchmark_molecules.py")
fetch_bm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fetch_bm)

NAMES = {"teaching": ["ethanol", "benzene"]}
PROPS = {"ethanol": {"CID": 702, "SMILES": "CCO", "HeavyAtomCount": 3},
         "benzene": {"CID": 241, "SMILES": "c1ccccc1", "HeavyAtomCount": 6}}
OLD = "name,category,cid,heavy_atoms,smiles\nold,teaching,1,1,C\n"


class FakePubChem:
    """urlopen stand-in: fails a name `fails[name]` times with `error`, then answers."""

    def __init__(self, fails=None, error=None):
        self.fails = dict(fails or {})
        self.error = error or TimeoutError("read timed out")
        self.calls = []

    def __call__(self, req, timeout=None):
        name = req.full_url.split("/name/")[1].split("/")[0]
        self.calls.append(name)
        if self.fails.get(name, 0) > 0:
            self.fails[name] -= 1
            raise self.error
        body = json.dumps({"PropertyTable": {"Properties": [PROPS[name]]}}).encode()
        return io.BytesIO(body)


def run(tmp_path, pubchem):
    out = tmp_path / "molecules.csv"
    out.write_text(OLD)

    def fetcher(name):
        return fetch_bm.fetch(name, urlopen=pubchem, sleep=lambda s: None)

    status = fetch_bm.main(out, NAMES, fetcher=fetcher, sleep=lambda s: None)
    return status, out.read_text()


def test_transient_errors_are_retried_then_the_file_is_replaced(tmp_path):
    pubchem = FakePubChem(fails={"ethanol": 2})  # two timeouts, then an answer
    status, text = run(tmp_path, pubchem)
    assert status == 0
    assert pubchem.calls.count("ethanol") == 3
    assert "ethanol,teaching,702,3,CCO" in text and "benzene" in text and "old" not in text
    assert not (tmp_path / "molecules.csv.tmp").exists()


def test_persistent_failure_keeps_the_old_file_and_exits_nonzero(tmp_path, capsys):
    pubchem = FakePubChem(fails={"benzene": 99})
    status, text = run(tmp_path, pubchem)
    assert status != 0
    assert text == OLD                                   # the snapshot is untouched
    assert pubchem.calls.count("benzene") == fetch_bm.ATTEMPTS
    err = capsys.readouterr().err
    assert "benzene" in err and "left unchanged" in err


def test_not_found_is_not_retried(tmp_path):
    error = urllib.error.HTTPError("u", 404, "PUGREST.NotFound", {}, None)
    pubchem = FakePubChem(fails={"ethanol": 99}, error=error)
    status, text = run(tmp_path, pubchem)
    assert status != 0 and text == OLD
    assert pubchem.calls.count("ethanol") == 1


def test_server_errors_are_transient(tmp_path):
    error = urllib.error.HTTPError("u", 503, "Service Unavailable", {}, None)
    pubchem = FakePubChem(fails={"benzene": 1}, error=error)
    status, _ = run(tmp_path, pubchem)
    assert status == 0 and pubchem.calls.count("benzene") == 2
