"""ARTiCANZ ingest.

Pulls the annotated TSV exports and flattens the "prefix:value; prefix:value"
tag strings into tidy records. ARTiCANZ is CC BY 4.0; attribution is rendered
in the page footer from config.yaml.
"""
from __future__ import annotations

import csv
import io
import pathlib
import re
import sys

import requests

TIMEOUT = 120
GENERIC_CATYPE = {"Cancer", "Solid tumour", "Carcinoma"}


def fetch_tsv(url: str) -> list[dict]:
    """Accepts an http(s) URL or a local path, so you can develop offline
    against a downloaded copy by editing the paths in config.yaml."""
    csv.field_size_limit(sys.maxsize)
    if url.startswith(("http://", "https://")):
        r = requests.get(url, timeout=TIMEOUT,
                         headers={"User-Agent": "NedTrialLandscape/1.0"})
        r.raise_for_status()
        text = r.text
    else:
        text = pathlib.Path(url).expanduser().read_text(encoding="utf-8")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def tag_vals(cell: str | None, prefix: str) -> list[str]:
    """Values may contain ':' (facility strings) and ',' (payload suffixes),
    so split on '; ' only and strip one leading 'prefix:'."""
    if not cell:
        return []
    out, pre = [], prefix + ":"
    for tok in cell.split("; "):
        if tok.startswith(pre):
            v = tok[len(pre):].strip()
            if v and v not in out:
                out.append(v)
    return out


# ---- drug class -> modality / target ---------------------------------------

def modality_of(cls: str) -> str:
    c = cls.lower()
    if "antibody-drug_conjugate" in c or "antibody_drug_conjugate" in c:
        return "Antibody-drug conjugate"
    if "t-cell_engager" in c:
        return "T-cell engager"
    if "radioligand" in c or "radioconjugate" in c:
        return "Radioligand / radioconjugate"
    if c.startswith("bispecific") and c.endswith("antibody"):
        return "Bispecific antibody"
    if "monoclonal_antibody" in c:
        return "Monoclonal antibody"
    if "vaccine" in c:
        return "Vaccine"
    if re.search(r"lymphocyte|car.?t|cell_therapy|nk_cell", c):
        return "Cell therapy"
    if c.startswith("radiotherapy") or "radiation" in c:
        return "Radiotherapy"
    if re.search(r"taxane|platinum|antimetabolite|alkylating_agent|anthracycline"
                 r"|fluoropyrimidine|topoisomerase_inhibitor$|cytotoxic", c):
        return "Cytotoxic chemotherapy"
    if re.search(r"_inhibitor|_degrader|protac|_agonist|_antagonist", c):
        return "Small molecule"
    if c == "placebo":
        return "Placebo"
    return "Other"


_T_RULES = [
    (r"^anti-(.+?)_antibody-drug_conjugate$", ""),
    (r"^bispecific[-_](.+?)_antibody-drug_conjugate$", "*"),
    (r"^anti-(.+?)_(?:monoclonal_antibody|antibody)$", ""),
    (r"^bispecific[-_](.+?)_antibody$", "*"),
]


def target_of(cls: str) -> str:
    base = cls.split(",")[0]
    rest = cls.split(",", 1)[1] if "," in cls else ""
    for pat, suffix in _T_RULES:
        m = re.match(pat, base)
        if m:
            t = m.group(1) + suffix
            break
    else:
        m = re.search(r"([A-Za-z0-9/.-]+)-targeting", rest) or \
            re.search(r"([A-Za-z0-9/.-]+)-targeting", base)
        if m:
            t = m.group(1)
        else:
            m = re.match(r"^(.+?)_(?:inhibitor|degrader|agonist|antagonist)$", base)
            t = m.group(1) if m else base.replace("_", " ")
    # same bispecific written two ways; unknown-target ADCs
    if t == "EGFR/cMET":
        t = "EGFR/cMET*"
    if re.match(r"^antibody.drug conjugate$", t):
        t = "unknown target"
    return t


def cancer_types(row: dict) -> list[str]:
    spec = [s for s in tag_vals(row.get("catype_expanded"), "catype_specific")
            if not s.startswith("NOT ")]
    if spec:
        return spec
    broad = [s for s in tag_vals(row.get("catype_expanded"), "catype")
             if not s.startswith("NOT ") and s not in GENERIC_CATYPE]
    return broad or ["Solid tumour"]


def build(datasets: dict[str, str], intent_labels: dict[str, str]) -> list[dict]:
    """Return one record per trial (targets/tissues as lists)."""
    seen: dict[str, dict] = {}
    for key, url in datasets.items():
        intent = intent_labels.get(key, key)
        for row in fetch_tsv(url):
            tid = (row.get("trial_id") or "").strip()
            if not tid:
                continue
            classes = tag_vals(row.get("drug_list_expanded"), "focused_therapy_class")
            # Keep modality and target PAIRED. Flattening them into two lists
            # leaks targets across modalities: a trial combining an ADC with
            # pembrolizumab would otherwise show PD-1 as an ADC target.
            pairs = sorted({(modality_of(c), target_of(c)) for c in classes})
            rec = {
                "registry": "ARTICANZ",
                "trial_id": tid,
                "intent": intent,
                "cancer_types": cancer_types(row),
                "classes": [{"m": m, "t": t} for m, t in pairs],
                "modalities": sorted({m for m, _ in pairs}),
                "drugs": tag_vals(row.get("drug_list_expanded"), "focused_drug"),
                "phase": " / ".join(sorted(tag_vals(row.get("phase_expanded"), "phase")))
                         or "Unspecified",
                "status": (tag_vals(row.get("recruitmentstatus_expanded"),
                                    "recruitmentstatus") or ["UNKNOWN"])[0],
                "states": sorted(tag_vals(row.get("states_expanded"), "recruitmentstate")),
                "n_sites": len(tag_vals(row.get("locations_expanded"), "facility")),
                "url": f"https://articanz.org/trial/{tid}",
            }
            if tid in seen:
                # a trial listed in both files spans both intents
                seen[tid]["intent"] = "Both"
                for f in ("cancer_types", "modalities", "drugs", "states"):
                    seen[tid][f] = sorted(set(seen[tid][f]) | set(rec[f]))
                merged = {(c["m"], c["t"]) for c in seen[tid]["classes"] + rec["classes"]}
                seen[tid]["classes"] = [{"m": m, "t": t} for m, t in sorted(merged)]
            else:
                seen[tid] = rec
    return list(seen.values())
