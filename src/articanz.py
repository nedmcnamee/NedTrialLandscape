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
ADC = "Antibody-drug conjugate"


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


def trial_url(tid: str) -> str:
    """ARTiCANZ aggregates both registries, so link each trial to its own.
    ANZCTR wants the number WITHOUT the ACTRN prefix in the query string."""
    if tid.startswith("ACTRN"):
        return ("https://anzctr.org.au/Trial/Registration/TrialReview.aspx"
                f"?ACTRN={tid[len('ACTRN'):]}")
    if tid.startswith("NCT"):
        return f"https://clinicaltrials.gov/study/{tid}"
    return f"https://articanz.org/trial/{tid}"


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
        return ADC
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
            # The "extra" exports carry no free-text title. Fall back to any
            # title-ish column if a future export adds one.
            title = next((row[k] for k in ("title", "brief_title", "public_title")
                          if row.get(k)), "")
            rec = {
                "registry": "ARTICANZ",
                "trial_id": tid,
                "intent": intent,
                "title": title,
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
                "url": trial_url(tid),
            }
            if tid in seen:
                seen[tid]["intent"] = "Both"
                for f in ("cancer_types", "modalities", "drugs", "states"):
                    seen[tid][f] = sorted(set(seen[tid][f]) | set(rec[f]))
                merged = {(c["m"], c["t"]) for c in seen[tid]["classes"] + rec["classes"]}
                seen[tid]["classes"] = [{"m": m, "t": t} for m, t in sorted(merged)]
            else:
                seen[tid] = rec
    return list(seen.values())


def adc_drug_names(records: list[dict]) -> list[str]:
    """Drug names that ARTiCANZ only ever associates with ADC trials.

    ARTiCANZ lists every focused_drug on a trial, so a combination study
    contributes both the conjugate and its partner. Taking the set difference
    -- names seen in ADC trials MINUS names seen in any non-ADC trial --
    removes comparators and backbones automatically: pembrolizumab and
    bevacizumab appear all over the non-ADC corpus and drop out, while
    DB-1311 or Trastuzumab Deruxtecan appear only alongside ADC classes and
    survive. Naked trastuzumab correctly drops out; the conjugate does not.

    This is what lets the ClinicalTrials.gov drug list maintain itself: a new
    agent becomes searchable the week ARTiCANZ's curators annotate it.
    """
    in_adc: set[str] = set()
    in_other: set[str] = set()
    for r in records:
        bucket = in_adc if any(c["m"] == ADC for c in r["classes"]) else in_other
        for d in r["drugs"]:
            bucket.add(d.strip())
    return sorted(n for n in (in_adc - in_other)
                  if len(n) > 2 and n.casefold() not in DRUG_STOPLIST)


# The set difference above removes any agent that also appears in a non-ADC
# trial, which catches most comparators. It cannot catch one that happens to
# appear ONLY alongside ADCs in this corpus -- topotecan does exactly that,
# as the control arm of a B7-H3 ADC study. Feeding it into the confirmation
# pattern would then match every topotecan chemotherapy trial on
# ClinicalTrials.gov. Free cytotoxics and naked antibodies are therefore
# excluded by name.
DRUG_STOPLIST = {
    # free camptothecins -- share the -tecan ending with conjugated payloads
    "topotecan", "irinotecan", "liposomal irinotecan", "exatecan", "belotecan",
    # common backbone chemotherapy
    "docetaxel", "paclitaxel", "nab-paclitaxel", "carboplatin", "cisplatin",
    "oxaliplatin", "gemcitabine", "capecitabine", "etoposide", "doxorubicin",
    "cyclophosphamide", "fluorouracil", "pemetrexed", "vinorelbine",
    "trifluridine/tipiracil", "eribulin",
    # checkpoint inhibitors and naked antibodies
    "pembrolizumab", "nivolumab", "atezolizumab", "durvalumab", "avelumab",
    "ipilimumab", "bevacizumab", "trastuzumab", "pertuzumab", "cetuximab",
    "rituximab", "ramucirumab", "zanidatamab",
    # not a drug
    "placebo", "best supportive care", "standard of care",
}

# ============================================================================
# ADD THIS to the END of src/articanz.py (after adc_drug_names and the
# DRUG_STOPLIST block). Nothing else in that file changes.
# ============================================================================


def adc_target_map(records: list[dict]) -> dict[str, str]:
    """Drug name -> ADC target, learned from ARTiCANZ's curated annotations.

    Companion to adc_drug_names(). That function tells ClinicalTrials.gov
    *which trials* are ADC trials; this one tells it *what they target*, so a
    study titled only "A Study of DM005 in Patients With Advanced Solid
    Tumors" resolves to EGFR/cMET instead of Other/unknown. No amount of
    regex on that title could have worked -- the antigen simply is not in the
    text. ARTiCANZ knows because a specialist annotated it.

    Two guards keep the mapping trustworthy:

      * Only trials carrying exactly ONE ADC target contribute. A study of
        two conjugates cannot tell us which drug hits which antigen.
      * A name claimed by more than one target is dropped rather than
        guessed at.

    Restricted to adc_drug_names() output, so comparators and backbone
    chemotherapy never enter the mapping.
    """
    allowed = set(adc_drug_names(records))
    claims: dict[str, set[str]] = {}
    for r in records:
        targets = {c["t"] for c in r["classes"] if c["m"] == ADC}
        if len(targets) != 1:
            continue
        target = next(iter(targets))
        for d in r["drugs"]:
            name = d.strip()
            if name in allowed:
                claims.setdefault(name, set()).add(target)
    return {n: next(iter(t)) for n, t in claims.items() if len(t) == 1}
