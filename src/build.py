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

SEED_KEY = "__seed__"      # records the date of the very first build
NEW_DAYS = 8               # "this week", with a day of slack
STALE_POST_DAYS = 120      # see is_genuinely_new()

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


def is_genuinely_new(rec: dict, cutoff: str, seed_date: str | None, today: str) -> bool:
    """Distinguish "a new trial opened" from "we only just started seeing it".

    Two things previously inflated this badge to meaninglessness:

    1. The seed build stamps every trial with the same date, so on the first
       run everything looks new. Trials carrying the seed date are excluded.

    2. Widening the classifier makes previously invisible trials appear. When
       the ARTiCANZ drug harvest went in, ~75 ClinicalTrials.gov studies were
       detected for the first time -- several registered years earlier. Those
       are new to *us*, not new to the world. If a trial reports a
       registration date well before we first saw it, it is a late detection
       rather than a new trial.
    """
    fs = rec.get("first_seen")
    if not fs or fs <= cutoff:
        return False
    if seed_date and fs == seed_date:
        return False
    post = rec.get("first_post") or (f"{rec['year']}-01-01" if rec.get("year") else None)
    if post:
        try:
            if date.fromisoformat(post[:10]) < date.fromisoformat(today) - timedelta(days=STALE_POST_DAYS):
                return False
        except ValueError:
            pass
    return True


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
        harvested = articanz.adc_drug_names(art) if art else []
        tmap = articanz.adc_target_map(art) if art else {}
        ctg, drops = ctgov.build(cfg["sources"]["ctgov"],
                                 extra_drugs=harvested, target_map=tmap)
        records += ctg

    for r in records:
        r["tissue_groups"] = group_tissues(r["cancer_types"], lut)
        r["active"] = r["status"].upper() in ACTIVE_STATUSES
        r.setdefault("title", "")

    today = date.today().isoformat()
    try:
        seen = json.loads(SEEN.read_text()) if SEEN.exists() else {}
    except json.JSONDecodeError:
        seen = {}
    seed_date = seen.pop(SEED_KEY, None)
    is_seed = not seen
    if is_seed:
        seed_date = today
    elif seed_date is None:
        # Legacy file written before the seed date was recorded: the earliest
        # timestamp in it IS the seed, because that build stamped everything
        # with one date.
        seed_date = min(seen.values())

    for r in records:
        key = f"{r['registry']}:{r['trial_id']}"
        seen.setdefault(key, today)
        r["first_seen"] = seen[key]

    cutoff = (date.today() - timedelta(days=NEW_DAYS)).isoformat()
    for r in records:
        r["is_new"] = is_genuinely_new(r, cutoff, seed_date, today)
    n_new = sum(1 for r in records if r["is_new"])
    n_first_seen = sum(1 for r in records
                       if r["first_seen"] > cutoff and r["first_seen"] != seed_date)

    SEEN.parent.mkdir(parents=True, exist_ok=True)
    SEEN.write_text(json.dumps({**seen, SEED_KEY: seed_date}, indent=0, sort_keys=True))
    DROPS.write_text(json.dumps(drops, indent=0, sort_keys=True))

    payload = {
        "generated": today,
        "seed_date": seed_date,
        "site": cfg["site"],
        "display": cfg["display"],
        "counts": {
            "total": len(records),
            "new_this_week": n_new,
            "first_seen_this_week": n_first_seen,
            "by_registry": {k: sum(1 for r in records if r["registry"] == k)
                            for k in sorted({r["registry"] for r in records})},
        },
        "trials": records,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, separators=(",", ":")))
    print(f"\nwrote {OUT.relative_to(ROOT)}  {len(records)} trials, "
          f"{n_new} genuinely new, {n_first_seen} first seen this week "
          f"(incl. late detections), {OUT.stat().st_size / 1024:.0f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
