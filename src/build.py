"""Orchestrator: fetch both registries, group tissues, emit docs/data.json.

Run locally:   python src/build.py
In CI:         see .github/workflows/update.yml

Writes:
  docs/data.json         what the site reads
  data/seen.json         trial_id -> date first seen, so the site can say
                         "new this week" without storing weekly snapshots
  data/ctgov_drops.json  nct_id -> which filter excluded it, consumed by
                         recall.py (not published, regenerated each run)
"""
from __future__ import annotations

import json
import pathlib
import re
import sys
from datetime import date, timedelta

import yaml

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import articanz  # noqa: E402
import ctgov  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
CONFIG = ROOT / "config.yaml"
OUT = ROOT / "docs" / "data.json"
SEEN = ROOT / "data" / "seen.json"
DROPS = ROOT / "data" / "ctgov_drops.json"

# Free-text cancer type -> tissue group, for ClinicalTrials.gov conditions.
# Order matters: GIST must beat "gastro", neuroendocrine must beat organ
# terms, and the basket catch-all must come last or it swallows everything.
CTG_RULES = [
    ("Sarcoma", r"sarcoma|gist\b|gastrointestinal stromal|desmoid|fibromatosis|osteosarcom|chondrosarcom|rhabdomyosarcom|bone tumou?r|ewing"),
    ("Endocrine Cancer", r"neuroendocrine|carcinoid|thyroid|adrenal|ph(a)?eochromocytoma|paraganglioma|pituitary|islet cell"),
    ("Peripheral Nervous System", r"neuroblastoma|schwannoma|peripheral nerve"),
    ("Central Nervous System", r"glioma|glioblastoma|astrocytoma|oligodendroglioma|medulloblastoma|meningioma|craniopharyng|ependymoma|\bbrain\b|\bcns\b|central nervous|leptomening"),
    ("Breast Cancer", r"breast"),
    ("Pulmonary Cancer", r"\blung\b|nsclc|sclc|non.small.cell|small.cell lung|mesotheliom|pulmonary|bronch|pleural"),
    ("Female Reproductive Cancer", r"ovarian|\bovary\b|endometri|cervic|fallopian|peritoneal|uterine|uterus|vulva|vagina|gyn(a)?ecolog"),
    ("Male Reproductive Cancer", r"prostat|\bcrpc\b|\bmcrpc\b|testicular|germ cell|penile|seminoma"),
    ("Renal & Urinary Cancer", r"renal|kidney|\brcc\b|bladder|urothelial|ureter|urinary|wilms"),
    ("Lower Gastro-Intestinal Tract Cancer", r"colorect|\bcolon\b|rectal|rectum|\bcrc\b|\banal\b|appendice"),
    ("Upper Gastro-Intestinal Tract Cancer", r"gastric|stomach|o?esophag|gastro.?o?esophageal|\bgej\b|gastrointestinal"),
    ("Pancreatic Cancer", r"pancrea"),
    ("Hepatobiliary Cancer", r"hepatocellular|\bhcc\b|\bliver\b|hepatic|biliary|cholangio|gallbladder|ampullary"),
    ("Head and Neck Cancers", r"head and neck|\bhnscc\b|nasopharyn|oropharyn|hypopharyn|laryn|oral cavity|salivary|tongue|tonsil"),
    ("Melanoma and Skin Cancer", r"melanoma|\bskin\b|cutaneous|merkel|basal cell carcinoma"),
    ("Lymphoma", r"lymphoma|hodgkin"),
    ("Leukaemias & Myelomas", r"leuk(a)?emia|myeloma|myelodysplas|myeloprolifer"),
    ("Pan-cancer / basket", r"solid tumou?r|advanced cancer|advanced malignan|neoplasm|unknown primary|metastatic cancer|\bcancer\b|\bcarcinoma\b|\btumou?r\b"),
]
CTG_RULES = [(g, re.compile(p, re.I)) for g, p in CTG_RULES]

# A trial is "active" if it is open, about to open, or still running.
ACTIVE_STATUSES = {
    "RECRUITING", "NOT_YET_RECRUITING", "ENROLLING_BY_INVITATION",
    "ACTIVE_NOT_RECRUITING", "AVAILABLE",
}


def load_config() -> dict:
    return yaml.safe_load(CONFIG.read_text(encoding="utf-8"))


def build_lookup(cfg) -> dict[str, str]:
    lut = {}
    for group, types in cfg["tissue_groups"].items():
        for t in types:
            lut[t.casefold()] = group
    return lut


def group_tissues(names, lut) -> list[str]:
    """Exact match against the config mapping first, then free-text regex."""
    out = []
    for n in names:
        g = lut.get(n.casefold())
        if g is None:
            for grp, pat in CTG_RULES:
                if pat.search(n):
                    g = grp
                    break
        out.append(g or "Unmapped")
    return sorted(set(out))


def main() -> int:
    cfg = load_config()
    lut = build_lookup(cfg)
    records: list[dict] = []
    art: list[dict] = []
    drops: dict[str, str] = {}

    if cfg["sources"]["articanz"]["enabled"]:
        print("ARTiCANZ")
        a = cfg["sources"]["articanz"]
        art = articanz.build(a["datasets"], a["intent_labels"])
        records += art
        print(f"  {len(art)} trials")

    if cfg["sources"]["ctgov"]["enabled"]:
        print("ClinicalTrials.gov")
        # ARTiCANZ curators name the ADCs; inherit that list rather than
        # hand-maintaining company codes. Falls back to the built-in list
        # when ARTiCANZ is disabled.
        harvested = articanz.adc_drug_names(art) if art else []
        ctg, drops = ctgov.build(cfg["sources"]["ctgov"], extra_drugs=harvested)
        records += ctg

    for r in records:
        r["tissue_groups"] = group_tissues(r["cancer_types"], lut)
        r["active"] = r["status"].upper() in ACTIVE_STATUSES
        r.setdefault("title", "")

    # first-seen tracking -> "new this week" without storing weekly snapshots
    today = date.today().isoformat()
    try:
        seen = json.loads(SEEN.read_text()) if SEEN.exists() else {}
    except json.JSONDecodeError:
        seen = {}
    is_seed = not seen
    for r in records:
        key = f"{r['registry']}:{r['trial_id']}"
        seen.setdefault(key, today)
        r["first_seen"] = seen[key]
    # On the very first build everything would look "new", which is noise.
    if is_seed:
        seed_date = (date.today() - timedelta(days=365)).isoformat()
        seen = {k: seed_date for k in seen}
        for r in records:
            r["first_seen"] = seed_date
    cutoff = (date.today() - timedelta(days=8)).isoformat()
    n_new = sum(1 for r in records if r["first_seen"] > cutoff)

    SEEN.parent.mkdir(parents=True, exist_ok=True)
    SEEN.write_text(json.dumps(seen, indent=0, sort_keys=True))
    DROPS.write_text(json.dumps(drops, indent=0, sort_keys=True))

    payload = {
        "generated": today,
        "site": cfg["site"],
        "display": cfg["display"],
        "counts": {
            "total": len(records),
            "new_this_week": n_new,
            "by_registry": {k: sum(1 for r in records if r["registry"] == k)
                            for k in sorted({r["registry"] for r in records})},
        },
        "trials": records,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, separators=(",", ":")))
    print(f"\nwrote {OUT.relative_to(ROOT)}  "
          f"{len(records)} trials, {n_new} new in the last week, "
          f"{OUT.stat().st_size / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
