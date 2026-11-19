"""Build benchmarks/molecules.csv: names → exact structures from PubChem.

    uv run python scripts/fetch_benchmark_molecules.py

Run once (needs the internet); the CSV is committed so the layout report itself
works offline. Transient errors (timeouts, 429/5xx) are retried; if any name still
fails, the failures are listed, the existing CSV is left untouched and the script
exits non-zero. Only a complete result replaces the file, atomically. Only the public compound names below are sent to PubChem.
Structures come from PubChem (public domain) — isomeric SMILES and CID.
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

try:  # verify TLS against the OS trust store (as curl does), e.g. behind a corporate proxy
    import truststore

    truststore.inject_into_ssl()
except ImportError:
    pass

OUT = Path(__file__).resolve().parent.parent / "benchmarks" / "molecules.csv"

# category → names. Chosen to cover what a chemistry explanation actually draws:
# everyday drugs, natural products, teaching examples, bridged/cage skeletons, macrocycles.
NAMES = {
    "drug": [
        "acetaminophen", "ibuprofen", "naproxen", "aspirin", "caffeine", "metformin", "atorvastatin",
        "simvastatin", "rosuvastatin", "lisinopril", "amlodipine", "metoprolol", "losartan", "omeprazole",
        "sertraline", "fluoxetine", "citalopram", "bupropion", "gabapentin", "levothyroxine", "warfarin",
        "clopidogrel", "sildenafil", "tadalafil", "diazepam", "alprazolam", "lorazepam", "zolpidem",
        "fentanyl", "tramadol", "methadone", "ketamine", "lidocaine", "procaine", "propranolol",
        "salbutamol", "montelukast", "cetirizine", "loratadine", "diphenhydramine", "ciprofloxacin",
        "amoxicillin", "penicillin G", "cefalexin", "doxycycline", "tetracycline", "acyclovir",
        "oseltamivir", "imatinib", "tamoxifen", "methotrexate", "chlorpromazine", "haloperidol",
        "risperidone", "olanzapine", "quetiapine", "clozapine", "levodopa", "donepezil", "memantine",
        "varenicline", "sumatriptan", "ondansetron", "furosemide", "hydrochlorothiazide", "allopurinol",
        "sitagliptin", "apixaban", "rivaroxaban", "remdesivir",
    ],
    "opioid/morphinan": [
        "morphine", "codeine", "oxycodone", "hydromorphone", "naloxone", "naltrexone", "buprenorphine",
        "dextromethorphan", "levorphanol", "heroin",
    ],
    "steroid/terpenoid": [
        "cholesterol", "testosterone", "estradiol", "progesterone", "ethinylestradiol", "prednisone",
        "dexamethasone", "hydrocortisone", "cholecalciferol", "retinol", "beta-carotene", "squalene",
        "alpha-tocopherol", "menthol", "limonene", "alpha-pinene", "beta-pinene", "eucalyptol", "thujone",
        "carvone", "geraniol", "farnesol", "artemisinin", "forskolin", "santonin", "ginkgolide B",
        "picrotoxinin", "paclitaxel", "camphor", "borneol", "isoborneol", "fenchone",
    ],
    "alkaloid": [
        "quinine", "strychnine", "brucine", "cocaine", "atropine", "scopolamine", "tropinone", "ephedrine",
        "reserpine", "vinblastine", "camptothecin", "capsaicin", "piperine", "berberine", "papaverine",
        "colchicine", "galantamine", "huperzine A", "lysergic acid", "lysergide", "psilocybin", "mescaline",
        "muscimol", "nicotine", "cytisine", "sparteine", "coniine", "tubocurarine", "epibatidine",
        "physostigmine", "yohimbine", "ajmalicine", "lupinine", "anatoxin-a",
    ],
    "cage/bridged": [
        "adamantane", "cubane", "norbornane", "norbornene", "norbornadiene", "bicyclo[2.2.2]octane",
        "1,4-diazabicyclo[2.2.2]octane", "quinuclidine", "methenamine", "triptycene", "twistane",
        "barrelene", "1,8-cineole", "7-oxabicyclo[2.2.1]heptane", "tricyclo[3.3.1.1~3,7~]decan-2-one",
        "memantine", "amantadine", "rimantadine"
    ],
    "macrocycle/large": [
        "erythromycin", "azithromycin", "clarithromycin", "sirolimus", "cyclosporin A", "vancomycin",
        "18-crown-6", "alpha-cyclodextrin", "valinomycin", "amphotericin B", "porphine", "phthalocyanine",
        "tacrolimus", "rifampicin",
    ],
    "teaching": [
        "ethanol", "acetone", "benzene", "toluene", "phenol", "aniline", "pyridine", "pyrrole", "furan",
        "thiophene", "indole", "naphthalene", "anthracene", "phenanthrene", "pyrene", "coronene",
        "cyclohexane", "cyclohexene", "cyclooctatetraene", "azulene", "styrene", "(E)-stilbene",
        "(Z)-stilbene", "fumaric acid", "maleic acid", "oleic acid", "glucose", "fructose", "sucrose",
        "ribose", "2-deoxyribose", "adenine", "guanine", "cytosine", "thymine", "uracil",
        "adenosine triphosphate", "glycine", "alanine", "tryptophan", "phenylalanine", "glutathione",
        "acetyl coenzyme A", "nicotinamide adenine dinucleotide", "ethylenediaminetetraacetic acid",
        "sodium dodecyl sulfate", "triolein", "buckminsterfullerene", "trans-cyclooctene",
        "cis-cyclooctene", "1,3-butadiene", "cyclopentadiene", "dicyclopentadiene", "isoprene",
    ],
}


TRANSIENT_HTTP = {429, 500, 502, 503, 504}  # PubChem throttling / server hiccups
ATTEMPTS = 4


class FetchError(Exception):
    """A name that could not be resolved after retrying (or at all, e.g. 404)."""


def fetch(name: str, *, urlopen=urllib.request.urlopen, sleep=time.sleep, attempts: int = ATTEMPTS) -> dict:
    """PubChem properties for `name`. Retries timeouts, connection errors and 429/5xx
    with backoff; raises FetchError otherwise or once the attempts are used up."""
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(name, safe="") + "/property/SMILES,HeavyAtomCount/JSON")
    req = urllib.request.Request(url, headers={"User-Agent": "chemcard-benchmark"})
    last = ""
    for attempt in range(1, attempts + 1):
        try:
            with urlopen(req, timeout=20) as r:
                return json.load(r)["PropertyTable"]["Properties"][0]
        except urllib.error.HTTPError as e:
            last = f"HTTP {e.code}"
            if e.code not in TRANSIENT_HTTP:
                raise FetchError(f"{name}: {last}") from None  # e.g. 404: retrying will not help
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last = f"{type(e).__name__}: {getattr(e, 'reason', e)}"
        except (KeyError, IndexError, ValueError) as e:
            raise FetchError(f"{name}: unexpected response ({e})") from None
        if attempt < attempts:
            sleep(2 ** attempt)  # 2, 4, 8 s
    raise FetchError(f"{name}: {last} (after {attempts} attempts)")


def main(out: Path = OUT, names: dict[str, list[str]] = NAMES, *, fetcher=fetch, sleep=time.sleep) -> int:
    """Fetch everything; replace `out` only if every name resolved. Returns the exit status."""
    rows, failures, seen = [], [], set()
    for category, group in names.items():
        for name in group:
            if name in seen:
                continue
            seen.add(name)
            try:
                props = fetcher(name)
            except FetchError as e:
                failures.append(str(e))
                print(f"  ! {e}")
                continue
            finally:
                sleep(0.25)  # PubChem asks for at most 5 requests per second
            rows.append({"name": name, "category": category, "cid": props["CID"],
                         "heavy_atoms": props.get("HeavyAtomCount", ""), "smiles": props["SMILES"]})
            print(f"{category:18s} {name:40s} CID {props['CID']}")
    if failures:
        print(f"\n{len(failures)} name(s) failed; {out} was left unchanged:", file=sys.stderr)
        for f in failures:
            print(f"  - {f}", file=sys.stderr)
        return 1
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_name(out.name + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["name", "category", "cid", "heavy_atoms", "smiles"])
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, out)  # atomic: readers see the old file or the complete new one
    print(f"{len(rows)} molecules → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
