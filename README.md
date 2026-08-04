# NedTrialLandscape

Cancer drug trial landscape for Australia / New Zealand and globally — drug
target × tissue, curative vs advanced intent. Rebuilds itself weekly from two
APIs and publishes a static dashboard to GitHub Pages. No server, no database,
nothing to keep running.

**Live site:** https://nedmcnamee.github.io/NedTrialLandscape/

## Data sources

| Source | What it gives | Licence |
|---|---|---|
| [ARTiCANZ](https://articanz.org) | Specialist-curated ANZ cancer drug trials with drug-class, target, cancer-type, eligibility and site annotations. Refreshed Fridays. | CC BY 4.0 |
| [ClinicalTrials.gov API v2](https://clinicaltrials.gov/data-api/api) | Global interventional Phase 1–3 trials, filtered to confirmed antibody-drug conjugates. | Public domain |

ARTiCANZ is the better data — it ships a curated drug-class ontology
(`focused_therapy_class`) that nothing else provides, which is what makes the
target × tissue view possible at all. ClinicalTrials.gov is included for
global context; its `conditions` field is free text and is resolved to tissue
groups by regex, which is necessarily less reliable.

## Setup

1. Push this repo to GitHub.
2. **Settings → Pages → Source: GitHub Actions.**
3. **Actions → Update trial landscape → Run workflow** to seed the first build.
4. Site appears at `https://<user>.github.io/NedTrialLandscape/`.

After that it runs itself: 21:00 UTC Friday (Saturday morning AEST, after
ARTiCANZ's Friday refresh).

## Local development

```bash
pip install -r requirements.txt
python src/build.py       # writes docs/data.json
python src/validate.py    # coverage + sanity checks
python -m http.server -d docs 8000
```

To work offline, point `sources.articanz.datasets` in `config.yaml` at
downloaded `.tsv` files — local paths work as well as URLs.

## Configuration

Everything tunable lives in `config.yaml`, and in practice it's the only file
you'll edit:

- `sources.*.enabled` — turn either registry off
- `sources.ctgov.*` — date window, phases, exclusion rules
- `display.*` — what the page shows on first load
- `tissue_groups` — the cancer type → tissue group mapping (164 types, 20 groups)

## Layout

```
config.yaml                 the only file you normally edit
src/articanz.py             ARTiCANZ TSV/API ingest + tag parsing
src/ctgov.py                ClinicalTrials.gov API v2, ADC classifier (v9 port)
src/build.py                orchestrates, groups tissues, writes docs/data.json
src/validate.py             fails CI on coverage regressions
data/seen.json              trial_id -> first-seen date (drives "new this week")
docs/index.html             the dashboard (static, no build step)
docs/data.json              generated
```

## Design notes

**Modality and target stay paired.** Each trial carries
`classes: [{m: modality, t: target}]` rather than two flat lists. Flattening
leaks targets across modalities — a trial combining an ADC with pembrolizumab
would otherwise show PD-1 as an ADC target.

**"New this week" without snapshots.** `data/seen.json` records the date each
trial was first observed. That's a few hundred KB that grows slowly, instead of
committing a 2 MB registry dump every week.

**Validation fails the build.** When upstream adds a cancer type that isn't in
`config.yaml`, it lands in `Unmapped`. `validate.py` fails the job once that
exceeds 10% of trials and prints the offending strings, so drift surfaces as a
red X rather than a quietly wrong figure.

## Relationship to the ADC registry analysis (v9)

`src/ctgov.py` is a port of the v9 R pipeline's fetch query, ADC confirmation
filter, and target assignment. It deliberately does **not** port the per-target
recall audit, LLM pattern drafting, or n=50 boundary validation sampling —
those are batch research steps and belong in v9, not in a weekly cron.

Two v9 bugs are fixed here and should be fixed there too:

1. `HAEM_PATTERN` is case-insensitive and contains `\bALL\b`, so it matches the
   ordinary English word "all" and silently drops trials from the denominator.
   Acronyms are matched case-sensitively in this port.
2. `END_DATE` was `"2026-04-31"`, which is not a real date (April has 30 days).

Biomarker classification (IHC / FISH / NGS / ctDNA) is not yet in this port.

## Licence

Code MIT. Data belongs to its sources — cite ARTiCANZ per CC BY 4.0 if you
reuse it.
