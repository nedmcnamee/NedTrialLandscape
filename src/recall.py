"""Cross-registry recall audit for the ClinicalTrials.gov ADC classifier.

ARTiCANZ is specialist-curated and calls ADCs from a drug-class ontology.
About 96% of its ADC trials also carry an NCT id. So any trial ARTiCANZ calls
an ADC, whose NCT id is absent from our ClinicalTrials.gov ADC set, is a
miss -- and it comes with drug names attached, which is what you need to
widen the pattern.

Crucially the misses split into two very different populations, and reporting
them as one number is misleading:

  * POLICY exclusions -- the trial WAS retrieved and recognised, then dropped
    on purpose (not yet recruiting, terminated, haematological, out of phase
    range). Nothing is broken; these reflect deliberate scope choices.
  * CLASSIFIER misses -- the trial was never retrieved, or was retrieved and
    not recognised as an ADC. These are real recall failures.

Only the second kind should drive pattern changes. data/ctgov_drops.json,
written by build.py, records which filter removed each candidate.

Run:  python src/recall.py          (after build.py)
Writes docs/recall.json for display on the site.
"""
from __future__ import annotations

import collections
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "docs" / "data.json"
DROPS = ROOT / "data" / "ctgov_drops.json"
OUT = ROOT / "docs" / "recall.json"

ADC = "Antibody-drug conjugate"

# Reasons that represent a deliberate scope decision rather than a failure.
POLICY_REASONS = {"not_yet_recruiting", "haematological", "no_drug_intervention"}

# Combination partners and comparators. They ride along on ADC trials but are
# not conjugates, so they must never be suggested as pattern additions -- and
# adding one to the confirmation pattern would match every chemotherapy trial
# using it.
NOT_ADC = {
    "topotecan", "irinotecan", "liposomal irinotecan", "exatecan", "belotecan",
    "docetaxel", "paclitaxel", "nab-paclitaxel", "carboplatin", "cisplatin",
    "oxaliplatin", "gemcitabine", "capecitabine", "etoposide", "doxorubicin",
    "cyclophosphamide", "fluorouracil", "pemetrexed", "vinorelbine",
    "trifluridine/tipiracil", "eribulin",
    "pembrolizumab", "nivolumab", "atezolizumab", "durvalumab", "avelumab",
    "ipilimumab", "bevacizumab", "trastuzumab", "pertuzumab", "cetuximab",
    "rituximab", "ramucirumab", "zanidatamab", "vepugratinib", "rilvegostomig",
    "placebo", "best supportive care", "standard of care", "radiotherapy",
}


def main() -> int:
    if not DATA.exists():
        print("recall: docs/data.json missing -- run build.py first")
        return 0

    trials = json.loads(DATA.read_text())["trials"]
    drops = json.loads(DROPS.read_text()) if DROPS.exists() else {}

    art = [t for t in trials if t["registry"] == "ARTICANZ"
           and any(c["m"] == ADC for c in t["classes"])]
    ctg_ids = {t["trial_id"] for t in trials
               if t["registry"] == "ClinicalTrials.gov"}

    if not ctg_ids:
        print("recall: no ClinicalTrials.gov trials in the build -- skipping")
        return 0

    art_nct = [t for t in art if t["trial_id"].startswith("NCT")]
    missing = [t for t in art_nct if t["trial_id"] not in ctg_ids]

    for t in missing:
        t["_reason"] = drops.get(t["trial_id"], "never_retrieved")
    policy = [t for t in missing if t["_reason"] in POLICY_REASONS
              or t["_reason"].startswith("status:")]
    classifier = [t for t in missing if t not in policy]

    matched = len(art_nct) - len(missing)
    raw_recall = matched / len(art_nct) if art_nct else 1.0
    # in-scope recall excludes trials we chose not to collect
    in_scope = len(art_nct) - len(policy)
    adj_recall = matched / in_scope if in_scope else 1.0

    print("\nCross-registry ADC recall")
    print(f"  ARTiCANZ ADC trials:            {len(art)}")
    print(f"  ...with an NCT id:              {len(art_nct)}")
    print(f"  ...found in our CT.gov ADC set: {matched}")
    print(f"  raw recall:                     {raw_recall * 100:.1f}%")
    print(f"  recall excluding policy drops:  {adj_recall * 100:.1f}%"
          f"   <- the number that reflects classifier quality")

    if policy:
        print(f"\n  {len(policy)} excluded on purpose:")
        for reason, n in collections.Counter(t["_reason"] for t in policy).most_common():
            print(f"    {n:3d}  {reason}")

    drugs = collections.Counter()
    for t in classifier:
        for d in t.get("drugs", []):
            if d.strip().casefold() not in NOT_ADC:
                drugs[d.strip()] += 1

    if classifier:
        print(f"\n  {len(classifier)} genuine classifier misses:")
        for reason, n in collections.Counter(t["_reason"] for t in classifier).most_common():
            print(f"    {n:3d}  {reason}")
        print("\n  Candidate drug names (comparators removed):")
        for name, n in drugs.most_common(30):
            print(f"    {n:3d}  {name}")

    OUT.write_text(json.dumps({
        "articanz_adc": len(art),
        "articanz_adc_with_nct": len(art_nct),
        "matched": matched,
        "recall_pct": round(raw_recall * 100, 1),
        "recall_in_scope_pct": round(adj_recall * 100, 1),
        "policy_excluded": len(policy),
        "classifier_missed": len(classifier),
        "policy_reasons": dict(collections.Counter(t["_reason"] for t in policy)),
        "missing": [{"nct_id": t["trial_id"], "reason": t["_reason"],
                     "drugs": t.get("drugs", []),
                     "targets": sorted({c["t"] for c in t["classes"] if c["m"] == ADC}),
                     "cancer_types": t["cancer_types"]}
                    for t in classifier],
        "candidate_drugs": [{"name": k, "n": v} for k, v in drugs.most_common()],
    }, indent=1))
    print(f"\n  wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
