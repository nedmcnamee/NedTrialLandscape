"""Cross-registry recall audit for the ClinicalTrials.gov ADC classifier.

ARTiCANZ is specialist-curated and calls ADCs from a drug-class ontology.
About 96% of its ADC trials also carry an NCT id. So any trial ARTiCANZ calls
an ADC, whose NCT id is absent from our ClinicalTrials.gov ADC set, is a
false negative in the regex -- and it comes with the drug names attached,
which is exactly what you need to widen the pattern.

This is a stronger audit than v9's per-antigen queries, because the ground
truth is a human-curated call rather than another regex.

Run:  python src/recall.py          (after build.py)
Writes docs/recall.json for display on the site.
"""
from __future__ import annotations

import collections
import json
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "docs" / "data.json"
OUT = ROOT / "docs" / "recall.json"

ADC = "Antibody-drug conjugate"


def main() -> int:
    if not DATA.exists():
        print("recall: docs/data.json missing -- run build.py first")
        return 0

    trials = json.loads(DATA.read_text())["trials"]
    art = [t for t in trials if t["registry"] == "ARTICANZ"
           and any(c["m"] == ADC for c in t["classes"])]
    ctg_ids = {t["trial_id"] for t in trials
               if t["registry"] == "ClinicalTrials.gov"}

    if not ctg_ids:
        print("recall: no ClinicalTrials.gov trials in the build -- skipping")
        return 0

    art_nct = [t for t in art if t["trial_id"].startswith("NCT")]
    missing = [t for t in art_nct if t["trial_id"] not in ctg_ids]
    recall = 1 - len(missing) / len(art_nct) if art_nct else 1.0

    print(f"\nCross-registry ADC recall")
    print(f"  ARTiCANZ ADC trials:            {len(art)}")
    print(f"  ...with an NCT id:              {len(art_nct)}")
    print(f"  ...found in our CT.gov ADC set: {len(art_nct) - len(missing)}")
    print(f"  recall: {recall * 100:.1f}%")

    # Drug names on the missed trials are the regex candidates
    drugs = collections.Counter()
    for t in missing:
        for d in t.get("drugs", []):
            drugs[d] += 1

    if missing:
        print(f"\n  {len(missing)} ARTiCANZ ADC trials not matched by the "
              f"ClinicalTrials.gov classifier.")
        print("  Drug names on those trials (candidates for ADC_DRUGS / "
              "CONFIRM in src/ctgov.py):")
        for name, n in drugs.most_common(30):
            print(f"    {n:3d}  {name}")

    OUT.write_text(json.dumps({
        "articanz_adc": len(art),
        "articanz_adc_with_nct": len(art_nct),
        "matched": len(art_nct) - len(missing),
        "recall_pct": round(recall * 100, 1),
        "missing": [{"nct_id": t["trial_id"], "drugs": t.get("drugs", []),
                     "targets": sorted({c["t"] for c in t["classes"] if c["m"] == ADC}),
                     "cancer_types": t["cancer_types"]}
                    for t in missing],
        "candidate_drugs": [{"name": k, "n": v} for k, v in drugs.most_common()],
    }, indent=1))
    print(f"\n  wrote {OUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
