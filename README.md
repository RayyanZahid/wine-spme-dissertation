# Wine SPME Arrow GC-MS/MS — Method-Development Scaffold

**Title:** Method Development Scaffold and Pre-Registered Protocol for SPME Arrow Headspace GC-MS/MS of Aroma-Active Volatiles in Cabernet Sauvignon Wine

**Status:** Scaffold built 2026-05-12 in an AI-orchestrated literature-synthesis run. Adversarial peer review on 2026-05-13 (via [brutalist-mcp](https://github.com/ejmockler/brutalist-mcp)) surfaced four methodological weaknesses (catalogued in §9.2) and several factual errors (fixed in commits since). Live URL is passphrase-gated (Tier 1 share-gate).

**Honesty disclosure (load-bearing):**
- No instrument runs have been performed.
- Chapters 1–6 reflect the published primary literature.
- Chapter 7 is a validation protocol pending execution.
- Chapter 8 is a pre-registered application architecture pending data.
- The previous draft's §7 LOQ / R² / recovery / RSD tables and §8 OPLS-DA accuracies were illustrative templates, not measurements. They have been withdrawn from this version.
- The artifact is the **method scaffold**: literature integration, transition database, scheduled-MRM segment planner, DoE / RSM / Pawliszyn-isotherm / Setschenow toolkit, validation framework, application architecture, and acknowledged-weaknesses chapter.

**Final deliverable:** `index.html` — single self-contained HTML, live at https://wine-spme-dissertation.vercel.app (passphrase-gated).

## Build pipeline that produced this

```
┌────────────────────────────────────────────────────────────────┐
│  research/topics.yaml   wine-spme-gcms topic added             │
│                          (67 queries, 16 seed DOIs, 28 authors)│
└──────────────────┬─────────────────────────────────────────────┘
                   │ python main.py research seed wine-spme-gcms
                   ▼
┌────────────────────────────────────────────────────────────────┐
│  data/research/papers/wine-spme-gcms.jsonl                     │
│   250 ranked papers (from 2650 unique S2 + OpenAlex hits)      │
└──────────────────┬─────────────────────────────────────────────┘
                   │ auto-load + embed on pgvector farm
                   ▼
┌────────────────────────────────────────────────────────────────┐
│  pgvector cluster @ immersivecommons:5432/research             │
│   249 wine-spme-gcms papers embedded (nomic-embed-text 768d)   │
└──────────────────┬─────────────────────────────────────────────┘
                   │ 6 parallel Agent subagents (WebSearch fanout)
                   ▼
┌────────────────────────────────────────────────────────────────┐
│  _research/01_spme_arrow_tech.md            (48 KB)            │
│  _research/02_cabsauv_aroma_chemistry.md    (43 KB)            │
│  _research/03_hs_spme_doe.md                (30 KB)            │
│  _research/04_mrm_optimization.md           (42 KB)            │
│  _research/05_validation_framework.md       (45 KB)            │
│  _research/06_faults_and_application.md     (33 KB)            │
└──────────────────┬─────────────────────────────────────────────┘
                   │ research distill (24 LLM tldr cards)
                   │ direct extraction of MRM transitions
                   ▼
┌────────────────────────────────────────────────────────────────┐
│  MRM/transitions.yaml                                          │
│   30 analytes × 2 transitions = 72 unique transitions          │
│   With CE, RT, threshold, range, descriptor per analyte        │
└──────────────────┬─────────────────────────────────────────────┘
                   │ code/mrm_optimizer.py
                   ▼
┌────────────────────────────────────────────────────────────────┐
│  37 acquisition segments, 47.3 min runtime                     │
│   ≥12 PPP at FWHM 5s, min dwell 36 ms                          │
│   Scheduled MRM, ion-ratio acceptance SANTE/11312/2021 ±30%    │
└──────────────────┬─────────────────────────────────────────────┘
                   │ CC-PRO compose (Bierut × Rams, editorial-academic)
                   ▼
┌────────────────────────────────────────────────────────────────┐
│  drafts/dissertation.html                                      │
│   10 chapters, ~22 000 words, 7 inline tables, 1 inline SVG    │
│   chromatogram figure (discord), compound-spotlight tooltip    │
│   microfeature, mobile-first, 50 numbered references           │
└────────────────────────────────────────────────────────────────┘
```

## File inventory

```
projects/wine-spme-dissertation/
├── BRIEF.md                              project brief + research questions
├── OUTLINE.md                            chapter outline (TOC)
├── README.md                             this file
├── _research/                            specialist briefs (parallel team)
│   ├── 01_spme_arrow_tech.md
│   ├── 02_cabsauv_aroma_chemistry.md
│   ├── 03_hs_spme_doe.md
│   ├── 04_mrm_optimization.md
│   ├── 05_validation_framework.md
│   └── 06_faults_and_application.md
├── code/                                 dissertation-specific tooling
│   ├── mrm_optimizer.py                  scheduled-MRM planner + CE-ramp helper
│   ├── spme_method_dev.py                DoE planner + RSM fit + isotherm fit
│   └── rag_synthesis.py                  pgvector queries per section
├── MRM/
│   ├── transitions.yaml                  authoritative analyte panel database
│   └── plans/<timestamp>/                generated acquisition plans
│       ├── acquisition_plan.json
│       ├── segments.csv
│       ├── panel_summary.md
│       ├── design_plackett_burman.csv
│       └── design_bbd.csv
└── drafts/
    └── dissertation.html                 the deliverable
```

## How to regenerate

```bash
# 1. Re-seed (refreshes the paper corpus)
python main.py research seed wine-spme-gcms

# 2. Re-distill (refreshes LLM cards)
python main.py research distill wine-spme-gcms --top 40

# 3. Re-build the MRM acquisition plan
python projects/wine-spme-dissertation/code/mrm_optimizer.py \
    --peak-fwhm 5 --target-ppp 12

# 4. Re-generate a DoE design
python projects/wine-spme-dissertation/code/spme_method_dev.py \
    optimize --design bbd --factors 3

# 5. Re-run RAG synthesis
python projects/wine-spme-dissertation/code/rag_synthesis.py

# 6. Open the dissertation
start projects/wine-spme-dissertation/drafts/dissertation.html
```

## CC-PRO composition choices

- **Persona**: Bierut × Rams hybrid — adaptive reuse of the dissertation form + every element pays rent.
- **Register**: editorial-academic (essay-adjacent, not commercial; thesis voice not sales voice).
- **Aesthetic DNA**: cream paper + ink black + Bordeaux red + pyrazine green. Color is provenance.
- **Discord**: §8 full-bleed quantitative aroma fingerprint — the academic structure breaks once to show what the instrument actually measures.
- **Microfeature**: compound spotlight tooltip — every italicized compound name reveals formula + MW + threshold + MRM transition on hover/tap.
- **Forces**: structure 0.85, density 0.55, warmth 0.40, stillness 0.85, volume 0.30.
- **Typography**: Source Serif 4 body, Inter UI/tables, JetBrains Mono transitions.
- **Mobile**: single-column at <760 px; full-width tables in overflow-x; margin notes inline; tap-to-show compound tooltips.

## Verification of method math

The MRM optimizer was smoke-tested with the actual panel:
```
37 segments, 126 transition-instances, 72 unique transitions
total runtime 47.3 min
max N_concurrent = 8 (in 22-30 min RT region)
min dwell = 36 ms; min PPP = 12.4
```

The SPME DoE helper was smoke-tested with synthetic data:
```
Plackett-Burman 7 factors → 8 runs, screening matrix correct
Box-Behnken 3 factors → 15 runs, design matrix correct
Pawliszyn isotherm fit on synthetic τ=10 min, n_inf=1000 → recovered exactly (R²=1.0)
Setschenow fit on synthetic kS=0.3 → recovered exactly (R²=1.0)
Derringer-Suich composite of [0.85, 0.60, 0.92] = 0.7771 (geometric mean)
```
