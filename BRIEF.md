# Dissertation — SPME Arrow Headspace GC-MS/MS Analysis of Cabernet Sauvignon

**Status:** in_progress · scaffolded 2026-05-12
**Build mode:** agentic RAG over pgvector farm + WebSearch fanout + CC-PRO compose
**Owner:** Rayyan (Life automation)

## Title
Method Development, MRM Optimization, and Validation of SPME Arrow Headspace GC-MS/MS for the Quantitative Profiling of Aroma-Active Volatiles in Cabernet Sauvignon Wine

## Research questions

1. **SPME Arrow vs traditional fiber SPME.** What is the analytically-resolvable gain (LOD, dynamic range, reproducibility) when substituting an Arrow device (DVB/CWR/PDMS, ~120 µm × 20 mm) for a conventional 100 µm PDMS or 50/30 DVB/CAR/PDMS fiber across the Cabernet aroma-compound space?
2. **Extraction-isotherm map.** Which Cab Sauv markers (methoxypyrazines, C13 norisoprenoids, fermentation esters, varietal sulfur compounds, smoke-taint guaiacols, fault markers) have headspace partition behavior compatible with non-equilibrium short-time extraction vs which require equilibrium sampling?
3. **MRM transition portfolio.** What is the minimum number of MRM transitions (quan + qual, per analyte) needed to cover the panel with ≥10 cycles/s across the chromatographic window, accommodating IBMP / IPMP / β-damascenone / β-ionone / TDN / TCA / 4-ethylguaiacol / 4-ethylphenol / DMS / 3-isobutyl-2-methoxypyrazine / linalool / hexyl acetate / ethyl octanoate / 2-phenylethanol?
4. **Method validation.** Does the developed method meet ICH Q2(R2) / AOAC SMPR / EURACHEM linearity (R² ≥ 0.99 over 2 orders of magnitude), accuracy (recovery 80–120 %), precision (RSD < 15 % at LOQ), and matrix-effect (slope ratio) acceptance criteria across model-wine spike-recovery and authentic Cab Sauv matrix?
5. **Application.** Can the method differentiate vintage, appellation (Napa vs Bordeaux vs Coonawarra vs Maipo), and oak-program (American/French, new/used, 225 L/300 L) by OPLS-DA on the quantitative aroma profile?

## Deliverables

- **Single dissertation HTML** rendered via CC-PRO (mobile-first, Rams × editorial register) at `drafts/dissertation.html`
- **MRM optimizer toolkit** at `code/mrm_optimizer.py` — precursor/product picker, CE optimizer, scheduled-MRM segment planner with cycle-time math
- **SPME Arrow method-dev helper** at `code/spme_method_dev.py` — DoE planner (Box-Behnken / CCD), extraction isotherm fitter, RSM plotter
- **Validation pack** at `validation/` — linearity, LOD/LOQ calculator, recovery + matrix-effect calculator, intermediate precision tracker
- **Reference library** at `_research/papers.yaml` — extracted MRM tables, validated parameter ranges, citations
- **MRM database** at `MRM/transitions.yaml` — full panel with literature CE, dwell, retention-time windows

## Build pipeline

```
research seed wine-spme-gcms        →  auto-loads to pgvector
  ↓
research distill wine-spme-gcms     →  per-paper tldr cards
  ↓
TeamCreate "wine-spme-deepsearch"   →  6 parallel WebSearch specialists
  ↓
research search "..."               →  semantic queries over pgvector
  ↓
code/* tooling build                →  MRM optimizer + SPME DoE
  ↓
CC-PRO compose                      →  single-HTML dissertation
```

## Source priority

1. Local-first: pgvector farm on immersivecommons (8955 papers, will grow with seed)
2. Semantic Scholar abstracts + recommendations (auto-fetched by seeder)
3. OpenAlex parallel queries (auto-fetched by seeder)
4. WebSearch + WebFetch by parallel team agents for vendor app notes, theses, AWRI Tech Reviews
5. Direct PDF ingestion via `research oa-find` + `pdf_worker` (Phase 2)

## Reference instrument

PerkinElmer / Shimadzu / Agilent / Thermo TSQ-class triple quadrupole assumed. Where vendor-specific implementations differ (Bruker EVOQ, Sciex 7500), generalize to the IUPAC MRM definition. Column: 60 m × 0.25 mm × 1.4 µm polar wax (DB-WAX UI / Stabilwax / SUPELCOWAX-10) — chosen for the polar/non-polar aroma compound spread.

## Voice & register

Technical-academic, Anglo-American thesis register. Numbered chapter / section / subsection structure. Tables and equations rendered natively. Inline citations as `(Author Year)` mapped to a final References block. No marketing voice; no em-dashes per Ray's voice memory (use periods, commas, parens).

## Non-goals

- Not a comprehensive review of wine flavor chemistry — focused on the analytical method, with chemistry only where it justifies a transition choice or extraction parameter.
- Not a vendor comparison (CTC/PAL3 vs Agilent 7693 vs Shimadzu AOC-6000) — autosampler-agnostic.
- Not a sensory study — quantitative aroma profile only; sensory thresholds cited where they motivate LOD targets.
