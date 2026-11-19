"""Build benchmarks/molecules.csv: names → exact structures from PubChem.

    uv run python scripts/fetch_benchmark_molecules.py

Run once (needs the internet); the CSV is committed so the layout report itself
works offline. Only the public compound names below are sent to PubChem.
Structures come from PubChem (public domain) — isomeric SMILES and CID.
"""

from __future__ import annotations

import csv
import json
import time
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


def fetch(name: str) -> dict | None:
    url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
           + urllib.parse.quote(name, safe="") + "/property/SMILES,HeavyAtomCount/JSON")
    req = urllib.request.Request(url, headers={"User-Agent": "chemcard-benchmark"})
    try:
        with urllib.request.urlopen(req, timeout=20) as r:
            props = json.load(r)["PropertyTable"]["Properties"][0]
    except Exception as e:  # noqa: BLE001
        print(f"  ! {name}: {e}")
        return None
    return props


def main() -> None:
    OUT.parent.mkdir(exist_ok=True)
    rows, seen = [], set()
    for category, names in NAMES.items():
        for name in names:
            if name in seen:
                continue
            seen.add(name)
            props = fetch(name)
            time.sleep(0.25)  # PubChem asks for at most 5 requests per second
            if props:
                rows.append({"name": name, "category": category, "cid": props["CID"],
                             "heavy_atoms": props.get("HeavyAtomCount", ""), "smiles": props["SMILES"]})
                print(f"{category:18s} {name:40s} CID {props['CID']}")
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "category", "cid", "heavy_atoms", "smiles"])
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} molecules → {OUT}")


if __name__ == "__main__":
    main()
