"""ClinicalTrials.gov API v2 ingest.

Port of the ADC registry analysis v9 R pipeline, extended so the global set
stays comparable with, and learns from, ARTiCANZ:

  1. DROP REASONS. Every candidate that does not survive is recorded with the
     filter that removed it, so recall.py can separate a pattern gap from a
     deliberate scope choice.

  2. HARVESTED DRUG NAMES. ARTiCANZ's curators name the ADCs; we inherit that
     list weekly instead of hand-maintaining ~90 company codes.

  3. HARVESTED TARGETS. ARTiCANZ also knows what each agent hits. A title
     like "A Study of DM005 in Patients With Advanced Solid Tumors" contains
     no antigen at all -- no regex can rescue it -- but ARTiCANZ says
     EGFR/cMET, so we use that.

  4. CANONICAL TARGET NAMES. The v9 rules said HER2, TROP-2, FRa, c-MET;
     ARTiCANZ says ERBB2, Trop2, FR-alpha, cMET. Left alone, the same antigen
     renders as two separate columns and the registries can never be compared
     on one axis. Everything is normalised to ARTiCANZ's vocabulary.

  5. GENERIC INN SUFFIX FAMILIES. ADC names are systematic (<stem>-<payload>),
     so matching the payload family catches agents nobody has listed yet.

Two deliberate deviations from v9, both bug fixes:
  * v9's HAEM_PATTERN is case-insensitive and contains \\bALL\\b, so it matches
    the ordinary English word "all". Acronyms are matched case-sensitively here.
  * v9's END_DATE was "2026-04-31", not a real date. Uses today instead.
"""
from __future__ import annotations

import re
import time
from datetime import date

import requests

BASE = "https://clinicaltrials.gov/api/v2/studies"
PAGE_SIZE = 1000
PAUSE = 0.4
TIMEOUT = 120
QUERY_CHUNK = 40

ADC_DRUGS = [
    "trastuzumab deruxtecan", "trastuzumab emtansine", "sacituzumab govitecan",
    "enfortumab vedotin", "mirvetuximab soravtansine", "disitamab vedotin",
    "datopotamab deruxtecan", "patritumab deruxtecan", "telisotuzumab vedotin",
    "ladiratuzumab vedotin", "tusamitamab ravtansine", "upifitamab rilsodotin",
    "ifinatamab deruxtecan", "raludotatug deruxtecan", "zilovertamab vedotin",
    "antibody drug conjugate", "antibody-drug conjugate",
    "HS-20089", "LNCB74", "BG-C9074", "puxitatug", "AZD8205",
    "sacituzumab tirumotecan", "BL-M07D1", "SHR-A1811", "ZW49",
    "sigvotatug vedotin", "sonesitatug vedotin", "telisotuzumab adizutecan",
    "HS-20117", "SYS6043", "MHB088C", "CMG901", "LM-302", "BGB-B455",
    "CX-2009", "IBI343", "anetumab", "becotatug vedotin", "MRG003",
]
CONCEPT_QUERY = '"antibody-drug conjugate" OR "ADC" OR conjugate'

I = re.IGNORECASE

FREE_CAMPTOTHECIN = (r"irinotecan|topotecan|exatecan|belotecan|rubitecan"
                     r"|gimatecan|silatecan|camptothecin")
GENERIC_INN = (
    rf"\b(?!(?:{FREE_CAMPTOTHECIN})\b)\w{{2,}}tecan\b"
    r"|\b\w{2,}dotin\b"
    r"|\b\w{2,}tansine\b"
)

CONFIRM_BASE = "|".join([
    r"antibody[- ]drug conjugate",
    r"deruxtecan|vedotin|emtansine|govitecan|soravtansine",
    r"mafodotin|tesirine|ravtansine|tirumotecan|rilsodotin",
    r"sunirine|tazevibulin|pasudotox|sarotalocan",
    r"maytansine|calicheamicin|auristatin", r"adizutecan",
    GENERIC_INN,
    r"(SN-38|exatecan|MMAE|MMAF|DM1|DM4).{0,30}conjugat",
    r"topoisomerase.{0,20}inhibitor.{0,20}conjugat",
    r"microtubule.{0,20}inhibitor.{0,20}conjugat",
    r"\bDXd\b",
    r"BL-B01D1|BL-M07D1|RC48|AZD0901|ARX788|MRG002|MRG003",
    r"HS-20089|LNCB74|BG-C9074|puxitatug|AZD8205",
    r"SHR-A1811|\bZW49\b|SYS6043|MHB088C|CMG901",
    r"HS-20117|LM-302|BGB-B455|CX-2009|IBI343",
])

PRIOR_ADC = re.compile("|".join([
    r"(prior|previous|progressed on|refractory to|intolerant to|receipt of)"
    r".{0,50}(ADC|antibody[- ]drug conjugate|immunoconjugate)",
    r"(prior|previous|progressed on|refractory to|intolerant to|receipt of)"
    r".{0,50}(deruxtecan|vedotin|emtansine|govitecan|soravtansine|MMAE|MMAF|DM1|DM4)",
]), I)

HAEM_WORDS = re.compile("|".join([
    "lymphoma", "hodgkin", "dlbcl", "follicular lymph", "mantle cell",
    "marginal zone", "burkitt", "t-cell lymph", "b-cell lymph",
    "anaplastic large", "peripheral t.cell", "leukaemia", "leukemia",
    "acute myeloid", "acute lympho", "chronic lympho", "chronic myeloid",
    "hairy cell", "myeloma", "plasma cell", "plasmacytoma", "myelodysplastic",
    "myeloproliferative", "myelofibrosis", "polycythaemia", "polycythemia",
    "essential thrombocyth", "amyloidosis", "amyloid light.chain",
    "waldenstrom", "macroglobulin", "aplastic anaemia", "aplastic anemia",
    "polatuzumab", "inotuzumab", "gemtuzumab", "loncastuximab",
    "camidanlumab", "coltuximab", "pinatuzumab",
]), I)
HAEM_ACRONYMS = re.compile(r"\bAML\b|\bALL\b|\bCLL\b|\bCML\b|\bMDS\b|\bMPN\b|\bAL amyloid")


def is_haem(txt: str) -> bool:
    return bool(HAEM_WORDS.search(txt) or HAEM_ACRONYMS.search(txt))


# ---- target assignment ------------------------------------------------------

# v9 named antigens its own way; ARTiCANZ names them differently. Since both
# feed one heatmap axis, normalise to the curated vocabulary. Without this,
# HER2 and ERBB2 sit in separate columns and no registry comparison is valid.
CANONICAL = {
    "HER2": "ERBB2", "HER3": "ERBB3", "TROP-2": "Trop2", "FRa": "FR-alpha",
    "Nectin-4": "NECTIN4", "c-MET": "cMET", "B7-H3": "B7H3", "B7-H4": "B7H4",
    "LIV-1": "LIV1", "TAG-72": "TAG72", "5T4": "5T4",
}

TARGET_RULES = [
    ("Trop2",    r"sacituzumab|datopotamab|TROP.?2|TACSTD2|\bSKB264\b|sac.tmt"),
    ("ERBB2",    r"trastuzumab|\bHER2\b|ERBB2|disitamab|\bRC48\b"
                 r"|\bBL-M07D1\b|\bSHR-A1811\b|\bZW49\b"),
    ("NECTIN4",  r"enfortumab|Nectin.?4|NECTIN4"),
    ("FR-alpha", r"mirvetuximab|folate receptor|FR.alpha|FOLR1"),
    ("ERBB3",    r"patritumab|\bHER3\b|ERBB3"),
    # antigens stated in the title of otherwise anonymous studies
    ("CD205",    r"\bCD205\b|\bLY75\b|\bOBT076\b"),
    ("TF",       r"tissue factor|\bADCE-T02\b|tisotumab"),
    ("KIT",      r"c.?Kit\b|\bNN3201\b"),
    ("FGFR2b",   r"\bFGFR2b\b|fibroblast growth factor receptor 2b|\bBG-C137\b"),
    ("DLK1",     r"\bDLK1\b|\bADCT-701\b"),
    ("ROR2",     r"\bROR2\b|\bBA3021\b"),
    ("EphA5",    r"\bEphA5\b|\bMBRC-101\b"),
    ("CDH6",     r"\bCDH6\b|cadherin.6|raludotatug|\bR-DXd\b"),
    ("cMET",     r"telisotuzumab|c-MET|cMET|\bMET\b|adizutecan|\bHS-20117\b|\bBYON3521\b"),
    ("EGFR",     r"becotatug|\bMRG003\b|EGFR.{0,20}(ADC|conjugate)|anti.EGFR.{0,20}(ADC|conjugate)"),
    ("ROR1",     r"zilovertamab|\bROR1\b"),
    ("B7H3",     r"ifinatamab|\bB7.H3\b|CD276|\bSYS6043\b|\bMHB088C\b"),
    ("B7H4",     r"\bB7.H4\b|VTCN1|GSK5733584|AZD8205|SGN.B7H4|HS-20089|LNCB74|BG-C9074|puxitatug"),
    ("CLDN18.2", r"\bCLDN18\.2\b|claudin.18|\bLM-302\b|\bCMG901\b"),
    ("5T4",      r"\b5T4\b|TBCA4|WAIF1|SYD1875"),
    ("DLL3",     r"\bDLL3\b|rovalpituzumab|SC16|BI 764532"),
    ("MSLN",     r"\bMSLN\b|mesothelin|anetumab"),
    ("CEACAM5",  r"\bCEACAM5\b|tusamitamab|\bIBI343\b"),
    ("PTK7",     r"\bPTK7\b|cofetuzumab"),
    ("LIV1",     r"LIV.1|SLC39A6|ladiratuzumab"),
    ("GPC3",     r"\bGPC3\b|glypican"),
    ("MUC16",    r"\bMUC16\b|CA125.ADC|sofituzumab"),
    ("MUC1",     r"\bMUC1\b|mucin.?1\b"),
    ("PSMA",     r"\bPSMA\b|prostate.specific membrane"),
    ("TAG72",    r"\bTAG.72\b|CC49"),
    ("EGFRvIII", r"\bEGFRvIII\b|depatux"),
    ("CLDN6",    r"\bCLDN6\b|claudin.6|\bBGB-B455\b"),
    ("CD166",    r"\bCD166\b|\bALCAM\b|\bCX-2009\b"),
]
TARGET_RULES = [(n, re.compile(p, I)) for n, p in TARGET_RULES]


def build_matchers(target_map: dict[str, str] | None):
    """Compile harvested drug names into word-bounded patterns.

    Longest first, so "AMT-1160" can never be shadowed by "AMT-116". The
    lookarounds stop a code matching inside a longer token.
    """
    if not target_map:
        return []
    items = sorted(target_map.items(), key=lambda kv: -len(kv[0]))
    return [(re.compile(r"(?<![A-Za-z0-9])" + re.escape(n) + r"(?![A-Za-z0-9])", I), t)
            for n, t in items]


def assign_target(text: str, matchers=()) -> str:
    # A curated ARTiCANZ annotation beats our own text matching.
    for pat, tgt in matchers:
        if pat.search(text):
            return CANONICAL.get(tgt, tgt)
    for name, pat in TARGET_RULES:
        if pat.search(text):
            return CANONICAL.get(name, name)
    return "Other/unknown"


def _advanced(cfg) -> str:
    phases = " OR ".join(cfg["phases"])
    end = date.today().isoformat()
    return (f"AREA[StudyType]INTERVENTIONAL AND AREA[Phase]({phases}) AND "
            f"AREA[StudyFirstPostDate]RANGE[{cfg['start_date']},{end}]")


def _pages(query: str, advanced: str, log=print) -> list[dict]:
    out, token, page = [], None, 1
    while True:
        params = {"query.term": query, "filter.advanced": advanced,
                  "countTotal": "true", "pageSize": PAGE_SIZE}
        if token:
            params["pageToken"] = token
        r = requests.get(BASE, params=params, timeout=TIMEOUT,
                         headers={"User-Agent": "NedTrialLandscape/1.0"})
        if r.status_code != 200:
            log(f"  CT.gov error page {page}: {r.status_code}")
            break
        body = r.json()
        studies = body.get("studies") or []
        if not studies:
            break
        out.extend(studies)
        token = body.get("nextPageToken")
        if not token:
            break
        page += 1
        time.sleep(PAUSE)
    return out


def _parse(s: dict) -> dict:
    ps = s.get("protocolSection", {})
    ids = ps.get("identificationModule", {})
    st = ps.get("statusModule", {})
    des = ps.get("designModule", {})
    arms = ps.get("armsInterventionsModule", {})
    cond = ps.get("conditionsModule", {})
    spon = ps.get("sponsorCollaboratorsModule", {})
    eli = ps.get("eligibilityModule", {})
    iv = arms.get("interventions") or []
    return {
        "nct_id": ids.get("nctId", ""),
        "brief_title": ids.get("briefTitle", ""),
        "phase": "; ".join(des.get("phases") or []) or "Unspecified",
        "status": st.get("overallStatus", "UNKNOWN"),
        "first_post_date": (st.get("studyFirstPostDateStruct") or {}).get("date", ""),
        "sponsor": (spon.get("leadSponsor") or {}).get("name", ""),
        "conditions": "; ".join(cond.get("conditions") or []),
        "interventions": "; ".join(i.get("name", "") for i in iv),
        "intervention_types": "; ".join(sorted({i.get("type", "") for i in iv})),
        "n_drug_interventions": sum(1 for i in iv
                                    if (i.get("type") or "").upper() in ("DRUG", "BIOLOGICAL")),
        "eligibility_text": eli.get("eligibilityCriteria") or "",
        "enrollment": (des.get("enrollmentInfo") or {}).get("count"),
    }


def build(cfg, extra_drugs=(), target_map=None, log=print):
    """Returns (records, drops). drops maps nct_id -> reason it was excluded."""
    adv = _advanced(cfg)

    extra = sorted({d for d in extra_drugs if d})
    confirm = re.compile(
        CONFIRM_BASE + ("|" + "|".join(re.escape(d) for d in extra) if extra else ""), I)
    matchers = build_matchers(target_map)
    log(f"  {len(extra)} drug names and {len(matchers)} target mappings "
        f"harvested from ARTiCANZ")

    raw = []
    terms = ADC_DRUGS + extra
    for i in range(0, len(terms), QUERY_CHUNK):
        chunk = terms[i:i + QUERY_CHUNK]
        log(f"  drug-name query {i // QUERY_CHUNK + 1} ({len(chunk)} terms)")
        raw += _pages(" OR ".join(chunk), adv, log)
    log("  concept query")
    raw += _pages(CONCEPT_QUERY, adv, log)

    by_id: dict[str, dict] = {}
    for s in raw:
        d = _parse(s)
        if d["nct_id"]:
            by_id.setdefault(d["nct_id"], d)
    rows = list(by_id.values())
    log(f"  {len(rows)} unique candidates")

    drops: dict[str, str] = {}
    keep = []
    for d in rows:
        nid = d["nct_id"]
        if d["status"].upper() in cfg["exclude_statuses"]:
            drops[nid] = "status:" + d["status"]
            continue
        if cfg["exclude_not_yet_recruiting"] and d["status"].upper() == "NOT_YET_RECRUITING":
            drops[nid] = "not_yet_recruiting"
            continue
        blob = f"{d['conditions']} {d['brief_title']} {d['interventions']}"
        if cfg["exclude_haematological"] and is_haem(blob):
            drops[nid] = "haematological"
            continue
        if not re.search(r"DRUG|BIOLOGICAL", d["intervention_types"], I):
            drops[nid] = "no_drug_intervention"
            continue
        arm = f"{d['interventions']} {d['brief_title']}"
        full = f"{arm} {d['eligibility_text']}"
        in_arm = bool(confirm.search(arm))
        if not confirm.search(full):
            drops[nid] = "no_adc_vocabulary"
            continue
        if not in_arm:
            drops[nid] = "adc_not_in_arm"
            continue
        if PRIOR_ADC.search(d["eligibility_text"]) and not in_arm:
            drops[nid] = "prior_adc_context_only"
            continue
        keep.append(d)

    log(f"  {len(keep)} confirmed ADC trials ({len(drops)} dropped)")

    out, unknown = [], 0
    for d in keep:
        tgt = assign_target(f"{d['interventions']} {d['brief_title']}", matchers)
        if tgt == "Other/unknown":
            unknown += 1
        out.append({
            "registry": "ClinicalTrials.gov",
            "trial_id": d["nct_id"],
            "intent": "Unspecified",
            "title": d["brief_title"],
            "cancer_types": [c.strip() for c in d["conditions"].split(";") if c.strip()],
            "classes": [{"m": "Antibody-drug conjugate", "t": tgt}],
            "modalities": ["Antibody-drug conjugate"],
            "drugs": [i.strip() for i in d["interventions"].split(";") if i.strip()],
            "phase": d["phase"].replace("PHASE", "Phase "),
            "status": d["status"],
            "states": [],
            "n_sites": 0,
            "sponsor": d["sponsor"],
            "year": int(d["first_post_date"][:4]) if d["first_post_date"][:4].isdigit() else None,
            "first_post": d["first_post_date"],
            "combination": d["n_drug_interventions"] > 1,
            "url": f"https://clinicaltrials.gov/study/{d['nct_id']}",
        })
    log(f"  {unknown} still Other/unknown "
        f"({100 * unknown / max(len(out), 1):.0f}%)")
    return out, drops
