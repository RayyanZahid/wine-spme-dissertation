# GC-MS/MS MRM Optimization — Theory + Transition Database

Chapter scope: theory of triple-quadrupole tandem MS, IUPAC nomenclature of selected/multiple reaction monitoring, the math behind scheduled MRM cycle-time budgeting, collision-energy (CE) optimization protocols, regulatory ion-ratio acceptance, and a curated transition table for the Cabernet Sauvignon aroma panel. All transitions are cross-cited to vendor application notes or peer-reviewed papers; where two literature sources disagree, both are reported.

---

## 1. Triple-quadrupole theory (Mathieu eq, CID, ion paths)

### 1.1 Physical layout

A triple-quadrupole (QqQ) instrument is a linear arrangement of three sets of four parallel rods:

| Stage | Designation | Role |
|---|---|---|
| Q1 | Mass filter (resolving) | Selects precursor ion at m/z₁ |
| Q2 (q2) | Collision cell — RF-only quadrupole, hexapole, or octupole; pressurized with a target gas | Fragments selected precursor by collision-induced dissociation (CID) |
| Q3 | Mass filter (resolving) | Selects diagnostic product ion at m/z₂ |

Q1 and Q3 operate as 2D mass filters; q2 operates as an RF-only ion guide where DC = 0, so all ions above a low-mass cutoff transit. Inventor lineage: the QqQ geometry was first developed at Michigan State University by Yost and Enke in the late 1970s, building on the earlier Morrison work at La Trobe (Wikipedia, *Selected reaction monitoring*).

### 1.2 Mathieu equation governs ion motion in Q1 / Q3

Ion motion along x and y in a 2D quadrupole field obeys the canonical Mathieu form:

$$\frac{d^2 u}{d\xi^2} + \left(a_u - 2q_u \cos(2\xi)\right) u = 0$$

with ξ = Ωt/2 (Ω is the RF angular frequency) and the dimensionless stability parameters defined as (Pfeiffer Vacuum *Quadrupole Mass Spectrometers*; SIMION docs; Whitman College MS e-book):

$$a_x = -a_y = \frac{8 e U}{m \Omega^2 r_0^2}$$

$$q_x = -q_y = \frac{4 e V}{m \Omega^2 r_0^2}$$

where U is the DC voltage applied to one rod-pair, V is the zero-to-peak RF amplitude, e is the elementary charge, m is the ion mass, and r₀ is the field radius (rod inscribed-circle radius).

### 1.3 Stability diagram and mass scanning

Only ions whose (a, q) pair lies inside the bounded stable region of (x, y) overlap survive the filter. The first stability region's apex sits at approximately (a ≈ 0.237, q ≈ 0.706) — the operating point of a high-resolution quadrupole (Pfeiffer Vacuum technical notes; ScienceDirect *Matrix methods for calculation of stability diagrams*).

The mass scan line is a straight line through the origin in (a, q) space with slope U/V = constant; ramping U and V together at fixed U/V sweeps the apex past successive m/z values, each entering the stable triangle in sequence. Heavier ions sit further left (lower q at any given V); lighter ions sit further right (higher q). Setting U = 0 collapses the operating point to the q-axis — the RF-only "ion-guide" mode used in q2, where the wide stable region passes essentially all ions above a low-mass cutoff (Wong/Cooks; SIMION 2024 supplemental docs).

### 1.4 CID dynamics in q2

Inside q2, ions accelerated by a small offset DC potential (the "collision energy" parameter, typically 5–30 eV in the laboratory frame for GC-MS/MS) collide with a target gas — most commonly N₂ at ~1–3 mTorr in modern GC QqQs (older designs used Ar; some platforms still default to it). Kinetic energy converts to internal energy, populating excited vibrational states that dissociate via the lowest-energy fragmentation pathway. Low-energy CID (≤1 keV laboratory frame) is the workhorse regime — high fragmentation efficiency but the spectrum is strongly CE-dependent; very low CE favors rearrangement and adduct chemistry, while higher CE drives direct bond cleavage (Wikipedia, *Collision-induced dissociation*; Révész 2023; *J. Am. Soc. Mass Spectrom.* iso-energetic CID).

Center-of-mass collision energy (the energy actually available for fragmentation) is:

$$E_{CM} = E_{LAB} \times \frac{M_{gas}}{M_{gas} + M_{ion}}$$

so heavier precursors require higher lab-frame CE to achieve the same E_CM as a small molecule. This is the formal basis for ramped or iso-energetic CID schemes (*PMC* "Center-of-Mass iso-Energetic CID").

### 1.5 Scan modes available on a QqQ

| Mode | Q1 | q2 | Q3 |
|---|---|---|---|
| Full scan (MS1) | Scan | OFF | RF-only |
| SIM (one analyte) | Static | OFF | RF-only |
| Product ion scan | Static (precursor) | Fragment | Scan |
| Precursor ion scan | Scan | Fragment | Static (product) |
| Neutral-loss scan | Scan | Fragment | Scan synchronized at ΔM offset |
| SRM / MRM | Static (precursor) | Fragment | Static (product) |

MRM is by far the most selective and sensitive mode for trace targeted quantification because both filters reject ~all chemical noise that does not pass the precursor→product mass coincidence.

---

## 2. MRM vs SRM vs full-scan — IUPAC definitions and tradeoffs

### 2.1 IUPAC-correct terminology

The 2013 IUPAC "Definitions of terms relating to mass spectrometry" recommends **SRM (Selected Reaction Monitoring)** as the preferred term for monitoring a single precursor→product transition. IUPAC formally disapproves of the abbreviation "MRM" in the literal sense because "multiple reaction monitoring" reads as monitoring a tandem cascade A→B→C, which would describe MS³ on an ion-trap, not the QqQ workflow (Wikipedia, *Selected reaction monitoring*; SRMAtlas glossary).

In practice — and as enshrined in every vendor application note since the 1980s — **MRM is read as "the application of SRM to multiple product ions from one or more precursor ions"**, run concurrently in time-sliced fashion across the chromatogram. This convention is so entrenched that the IUPAC objection is effectively a footnote; we use MRM throughout this dissertation while noting the formal disagreement.

PRM (Parallel Reaction Monitoring) is the high-resolution analogue on a Q-Exactive / Orbitrap: Q1 selects a precursor, q2 fragments it, and the Orbitrap acquires the entire product-ion spectrum in parallel — equivalent in spirit to running every possible SRM transition simultaneously (Wikipedia, *Selected reaction monitoring*).

### 2.2 Sensitivity / specificity tradeoffs

| Mode | Selectivity | Sensitivity (S/N at trace) | Quantitative dynamic range | Use case |
|---|---|---|---|---|
| Full scan (Q1 only) | Low — all co-eluting matrix passes | Limited by chemical noise from matrix | ~10² | Screening, structure confirmation, untargeted |
| SIM (Q1 only) | Moderate — one m/z window per analyte | 10–100× better than full scan; matrix isobars still pass | ~10³ | Single-quad targeted |
| MRM (QqQ) | Very high — precursor AND product mass coincidence required | 10–1000× better than SIM in dirty matrices | 10⁴–10⁶ | Targeted trace quantification (ppt–ppb) |

The dramatic selectivity gain comes from the requirement that a co-eluting interferent simultaneously match both Q1 (precursor) and Q3 (product) masses to within the resolving window — a coincidence vanishingly rare for any chemical noise that is not the analyte (Creative Proteomics, *SRM/MRM principles*; Agilent training material *Introducing the 7000A QQQ-MS*).

For wine volatiles at sub-ng/L levels (e.g., IBMP, TCA, geosmin, 3MH), MRM is the only feasible quantitative mode; SIM on a single-quad blows out at the matrix peak from ethanol/glycerol envelope tails and yeast-derived isobaric esters.

---

## 3. Scheduled MRM mathematics (cycle-time formula + worked example for 30 analytes)

### 3.1 Definitions

- **Dwell time (D)** — time the instrument actively collects counts on a single transition during one pass of the cycle.
- **Interscan delay / pause / settling time (P)** — time the quadrupole DC/RF rods need to slew between adjacent transitions and let any prior product ions clear q2 (otherwise cross-talk).
- **Cycle time (CT)** — time to acquire ONE point on every monitored transition once. Equals the data-point spacing along the time axis for any given MRM.
- **N_concurrent** — number of transitions monitored simultaneously at a given retention-time point.
- **Points per peak (PPP)** — chromatographic peak FWHM divided by CT.

### 3.2 Cycle-time formula

For constant-dwell methods:

$$CT = N_{concurrent} \times (D + P)$$

For the scheduled / dynamic case (Agilent dMRM, Sciex sMRM, Shimadzu Smart MRM, Thermo t-MRM), each transition has its own retention-time window of width W; the algorithm computes, at every time slice along the chromatogram, the number of transitions whose window covers that slice, then either:

- (a) holds CT constant and lets D vary per transition: $D_i = (CT / N_{concurrent}) - P$, or
- (b) holds D constant and lets CT vary: $CT_t = N_{concurrent,t} \times (D + P)$ at time t.

Sciex's sMRM Pro additionally weights D by transition abundance — high-S/N transitions get less dwell so weak transitions can get more — keeping CT fixed (Sciex tech note *Scheduled MRM Algorithm Pro*; Agilent technical overview *New Dynamic MRM Mode*; Sepscience *Understanding MRM Fundamentals*).

### 3.3 Points-per-peak target

The community consensus from peer-reviewed validation papers and Agilent / Sciex / Shimadzu user guides is:

- **Minimum 10 PPP** for area-under-curve quantitation to keep integration RSD under ~5%.
- **12–15 PPP** preferred when ion ratios are also evaluated for confirmation, because the qualifier-channel error compounds with the quan-channel error.
- **Above 20 PPP** the marginal sensitivity loss from too-short dwell times outweighs any gain in integration precision.

Sciex KB: "Total cycle time must give 10–15 data points over the HPLC peak to allow accurate quantitation" (Sciex *Knowledge Base: How to choose dwell/cycle time*). Agilent recommends 10–20 points per peak as a general GC-MS rule (Agilent Community *GC Detector Sampling Rate*).

Worked relationship: for a typical SPME-Arrow HS GC peak with FWHM ≈ 3 s, a 10 PPP target gives:

$$CT \leq \frac{FWHM}{PPP} = \frac{3 \text{ s}}{10} = 300 \text{ ms}$$

With a P ≈ 5 ms typical for modern QqQs (Agilent 7000D/E, Thermo TSQ 9000, Sciex 7500 GC, Shimadzu TQ8050NX), the dwell-time budget is:

$$D_{max} = \frac{CT}{N_{concurrent}} - P$$

### 3.4 Worked example — 30 wine aroma analytes, 2 transitions each (60 total)

Assume a single-segment (unscheduled) MRM method first:

- Total transitions = 60.
- CT target = 300 ms.
- P = 5 ms.
- D = CT / N – P = 300/60 – 5 = 5 – 5 = 0 ms → **infeasible**.

Switch to scheduled MRM with retention-time windows W = 60 s per analyte. The empirical concurrency at the busy mid-chromatogram region for a 30-analyte wine panel is typically 6–10 transitions (i.e., 3–5 analytes co-eluting within their RT windows; the eluting peak count is much smaller than the total panel because they spread out across a ~40 min run).

- N_concurrent ≈ 8.
- D = CT / N – P = 300/8 – 5 = 37.5 – 5 = **32.5 ms per transition**.

32.5 ms dwell is comfortably above the typical instrument lower limit (1–5 ms) and well-suited to trace analytes — even sub-ng/L IBMP and TCA give detectable counts at this dwell on modern HES (high-efficiency source) Agilent 7000E or Thermo TSQ 9610 hardware.

If the panel grows to ~50 analytes with 2 transitions each (= 100 transitions) and concurrency rises to ~14, with CT held at 300 ms:

- D = 300/14 – 5 = 16.4 ms — still fine.

Past ~150 concurrent transitions, the analyst must either (i) lengthen CT (degrading PPP and quan precision) or (ii) split into time segments — most modern methods choose scheduled MRM and let the algorithm handle it dynamically (Agilent, Sciex, Shimadzu all converge on this design).

### 3.5 Cross-talk and minimum dwell

Cross-talk between adjacent transitions occurs when product ions from transition N are still resident in q2 when transition N+1 begins, falsely registering as N+1 counts. Modern QqQs minimize this with curved or branched collision cells (Agilent Hexabore, Sciex Linear Accelerator, Thermo HyperQuad) and short (~1 ms) clearing pulses. Practical floor: **D ≥ 5 ms** to avoid cross-talk artifacts; this is the constraint that ultimately caps panel size for any given peak width (Agilent ASMS 2019 poster *Impact of Dwell Time and Ion Flux on MRM*).

---

## 4. CE optimization protocol (ramp, step, hot-spot detection)

### 4.1 Strategy

For each compound, CE is optimized **per transition** rather than globally — the precursor→quan product is typically tuned for maximum product-ion abundance; the precursor→qual product is tuned for a different (usually higher) CE that delivers a distinctly different product-ion ratio for confirmation.

Standard MassHunter Optimizer (Agilent), AutoSRM (Thermo), and Smart MRM (Shimadzu) all converge on the same procedure (Agilent *MassHunter Optimizer for GC Triple Quad Quick Start*):

1. **Q1 precursor scan**: full scan of the analyte standard at electron impact 70 eV → pick the 3–4 most-abundant + most-specific m/z values as candidate precursors (avoid generic low-mass ions like m/z 43, 57, 91 unless no alternative).
2. **Product-ion scan** at CE = 5, 10, 15, 20, 25, 30 eV for each candidate precursor → identify which product ions appear at each CE.
3. **CE ramp** at 1–2 eV steps from 0 to 40 eV on each (precursor → product) candidate pair → plot product abundance vs CE → pick the CE at the peak of each transition curve (the "hot spot"). For wine volatiles this typically falls between 5–25 eV (Sciex *Center-of-Mass iso-Energetic CID* PMC7249026).
4. **Pair selection**: choose the (precursor → product) at maximum-signal CE as **quan**; choose a second pair with a distinct CE (often higher, ≥5 eV away) and a distinct product mass as **qual** (so the qual/quan ion ratio is robust to small CE drift).
5. **Iso-energetic check**: confirm the CE for closely-related compound classes (e.g., methoxypyrazines IBMP/IPMP/SBMP) using the center-of-mass scaling formula so isobaric pairs don't collapse onto the same CE / m/z.

### 4.2 Per-transition vs single CE

Older (pre-2010) papers often quote a single CE used across all wine analytes (e.g., 10 eV or 15 eV). This trades 10–50% sensitivity for method simplicity. Modern Agilent / Thermo dMRM tables list per-transition CE — the default we adopt in §6 and the published vendor methods.

### 4.3 Hot-spot detection

When the CE ramp curve is broad (FWHM > 10 eV), the analyte is forgiving — small instrument-to-instrument CE drift will not change quan reproducibility. When the curve is sharp (FWHM < 3 eV — "very narrow collision energy peak" reported on Chromatography Forum and Restek tech notes), the pair is fragile; pick an alternative transition or accept tighter tune QC. Particularly common for small aldehydes (hexanal), where the M⁺• is weak and the only diagnostic product ions cluster around m/z 41, 44, 56 with steep CE dependence (Razaq et al., *J. Agric. Food Chem.* 2023 PMC10835727).

---

## 5. Ion-ratio acceptance (SANTE, FDA, EU)

Every quantitative method that reports a "confirmed" detection (vs a "tentative" detection) verifies the ratio of qualifier-to-quantifier ion abundance against the same ratio measured in a contemporaneously-injected reference standard. The acceptance band depends on the regulatory framework:

### 5.1 SANTE/11312/2021 (European Union — pesticide residues, the de-facto reference for trace food analysis)

The current EU guideline replaces the older 2002/657/EC sliding-scale table with a single tolerance:

**±30% relative deviation from the standard ratio, applied at any qualifier abundance, for LC-MS/MS and GC-MS/MS.** Both transitions must have S/N ≥ 3 (EU Reference Lab summary *Main changes introduced in SANTE/11312/2021*; in force 1 January 2022).

### 5.2 Earlier 2002/657/EC sliding scale (still cited by many labs)

The older Commission Decision used a band that loosened for less-abundant qualifiers:

| Q/Q (qualifier as % of quantifier) | Acceptance band |
|---|---|
| > 50% | ±20% |
| 20–50% | ±25% |
| 10–20% | ±30% |
| ≤ 10% | ±50% |

Many wine-analysis literature methods (pre-2022) cite this table. The harmonized 2021/11312 ±30% is now preferred.

### 5.3 FDA — bioanalytical and confirmation guidance

FDA *Bioanalytical Method Validation Guidance for Industry* (May 2018) requires at least three structurally specific ions for confirmation and discourages use of water-loss or isotopic ions as confirmation ions. FDA's *Mass Spectrometry for Confirmation of the Identity of Animal Drug Residues* guidance (2003) sets:

- **±10% absolute** for relative abundance ≥ 50%.
- **±15% absolute** for 25–50%.
- **±20% absolute** for 10–25%.
- **±50% relative** for <10% — i.e., similar in spirit to the 2002/657/EC table but expressed as absolute % differences for the most-abundant ions.

### 5.4 AOAC 2007.01 (multi-residue pesticide GC-MS/MS — relevant analogue for wine off-flavor methods)

AOAC requires a second transition with "reasonably matching relative abundance ratios versus a contemporaneously analyzed reference standard." The numerical band defaults to the SANTE recommendation for laboratories that conform to EN 17034 / ISO 17025 schemes.

### 5.5 Recommendation for this dissertation

We adopt **SANTE/11312/2021 ±30%** for routine confirmation across the Cabernet Sauvignon panel — modern, harmonized, statistically defensible, and accepted by every major instrument vendor's default acquisition method file. For TCA and IBMP at sub-ng/L sensory-threshold levels where qualifier S/N drops below ~5, we additionally apply the 2002/657/EC sliding scale as a fallback (i.e., ±50% when Q ion < 10% of quan) to avoid false rejection of true positives at the LOQ boundary.

---

## 6. Transition table — full panel

Notes on the table:

- **Quan** = primary quantifier transition (largest, most-specific product ion at optimized CE).
- **Qual** = secondary qualifier transition for confirmation via ion-ratio.
- **CE in eV** is the laboratory-frame collision energy unless otherwise noted; vendor convention varies (some report voltages, often equivalent for singly-charged ions).
- **RT estimate** is a typical retention time on a DB-WAX (polar) or DB-5MS (5%-phenyl-95%-methyl) column for HS-SPME GC-MS/MS wine methods — actual RT must be confirmed per instrument.
- **Source citation** is given as a short tag; full references in §8.
- Where two literature sources report different CEs for the same pair, both are listed (`/`).

### 6.1 Methoxypyrazines (vegetal / green-pepper markers)

| Compound | Formula | MW (Da) | RT (DB-WAX, min) | Precursor m/z | Product m/z (Quan) | Product m/z (Qual) | CE Quan (eV) | CE Qual (eV) | Source |
|---|---|---|---|---|---|---|---|---|---|
| IPMP (3-isopropyl-2-methoxypyrazine) | C₈H₁₂N₂O | 152.19 | 17–19 | 137 | 124 | 152 | 10 | 5 | Slabizki/Schmarr 2014; Hjelmeland 2010 |
| IBMP (3-isobutyl-2-methoxypyrazine) | C₉H₁₄N₂O | 166.22 | 21–23 | 124 | 94 | 81 | 15 | 20 | Hjelmeland; Antalick; Allamy; Schmarr |
| SBMP (3-sec-butyl-2-methoxypyrazine) | C₉H₁₄N₂O | 166.22 | 20–22 | 138 | 124 | 123/96 | 10 | 10 / 15 | Slabizki/Schmarr |
| d3-IBMP (internal std) | — | 169 | 21–23 (co-elute) | 127 | 97 | — | 15 | — | Allamy; Pickering |

The IBMP 124→94 pair (loss of CH₂O from the methoxy) is universally reported as the most-sensitive quan transition (Hjelmeland 2010 *Methoxypyrazines analytical*; Allamy et al. 2018; Wang/Ebeler tech reports). The qual 124→81 (loss of methoxy + HCN) gives a robust ion ratio around 30–40% at 20 eV.

### 6.2 Norisoprenoids (carotenoid-derived aging markers)

| Compound | Formula | MW | RT (DB-WAX, min) | Precursor | Quan | Qual | CE Quan | CE Qual | Source |
|---|---|---|---|---|---|---|---|---|---|
| β-damascenone | C₁₃H₁₈O | 190.28 | 27–29 | 190 | 121 | 105 / 69 | 10 | 12 / 15 | Schreier/Wuest 2017 ABC |
| β-damascenone (alt secondary precursor) | — | — | — | 175 | 142 | — | 15 | — | Schreier/Wuest |
| β-ionone | C₁₃H₂₀O | 192.30 | 31–33 | 192 | 177 | 91 | 12 | 20 | Schreier/Wuest 2017 |
| β-ionone (secondary cascade) | — | — | — | 177 | 147 | — | 20 | — | Schreier/Wuest |
| α-ionone | C₁₃H₂₀O | 192.30 | 29–31 | 192 | 177 | 121 | 10 | 15 | Schreier/Wuest 2017 |
| TDN (1,1,6-trimethyl-1,2-dihydronaphthalene) | C₁₃H₁₆ | 172.27 | 26–28 | 172 | 157 | 142 | 15 | 20 | AWRI Riesling work; Sacks 2012 JAFC |
| TDN (cascade) | — | — | — | 157 | 142 | 115 | 15 | 25 | Sacks 2012 |

The β-damascenone 190→121 pair (Quan) is the most-sensitive in HS-SPME methods reported by the multidimensional GC × MS/MS group (Schreier & Wuest, *Anal. Bioanal. Chem.* 408, 7593, 2016; the same paper also reports β-damascenone 190→69 as a third qualifier and the 175→142 secondary transition).

### 6.3 Haloanisoles (cork-taint markers — quantified despite not being SPME-Arrow-favored due to extreme sensitivity demand)

| Compound | Formula | MW | RT (DB-5MS, min) | Precursor | Quan | Qual | CE Quan | CE Qual | Source |
|---|---|---|---|---|---|---|---|---|---|
| TCA (2,4,6-trichloroanisole) | C₇H₅Cl₃O | 211.47 | 18–20 | 210 | 195 | 167 | 15 | 25 | Pollnitz 1996 AJGWR; Agilent 5990-4968EN; Ebeler 2012 AJEV |
| TCA (³⁷Cl isotope) | — | 212 | 18–20 | 212 | 197 | 169 | 15 | 25 | Pollnitz; Agilent 5990-4968EN |
| d5-TCA (internal std) | — | 215 | 18–20 (co-elute) | 215 | 200 | — | 15 | — | Pollnitz |
| TBA (2,4,6-tribromoanisole) | C₇H₅Br₃O | 344.83 | 22–24 | 346 | 331 | 329 | 15 | 15 | Ebeler 2012 AJEV |
| TeCA (2,3,4,6-tetrachloroanisole) | C₇H₄Cl₄O | 245.91 | 21–23 | 246 | 231 | 203 | 15 | 25 | Ebeler 2012 |
| PCA (pentachloroanisole) | C₇H₃Cl₅O | 280.35 | 24–26 | 280 | 265 | 237 | 15 | 25 | Ebeler 2012 |

Note on TCA isotopologue strategy: TCA has a Cl₃ pattern with M, M+2, M+4, M+6 peaks at ratios approximately 27:27:9:1. Many published methods use 210→195 (Cl₂ loss giving the Cl₂-anisole product) AND 212→197 simultaneously as a built-in isotope-ratio confirmation in addition to (or instead of) a structurally-distinct qualifier — the chlorine-isotope ratio is essentially impossible to mimic in a co-eluting matrix interferent and provides ~3× more confirmation confidence than a random qualifier (Pollnitz et al., *Aust. J. Grape Wine Res.* 1996; Agilent app note 5990-4968EN).

### 6.4 Smoke-taint volatile phenols (forest-fire / oak markers)

These transitions are taken from the Agilent 8890/7000D wine smoke-taint app note (5994-3161EN, 2021) and the *Metabolites* 2020 paper PMC7407152 (which use harmonized 15 eV across the panel).

| Compound | Formula | MW | RT (DB-WAX, min) | Precursor | Quan | Qual | CE Quan | CE Qual | Source |
|---|---|---|---|---|---|---|---|---|---|
| Guaiacol | C₇H₈O₂ | 124.14 | 11.3 | 124 | 109 | 81 | 15 | 20 | Agilent 5994-3161EN; PMC7407152 |
| 4-methylguaiacol | C₈H₁₀O₂ | 138.16 | 12.3 | 138 | 123 | 95 | 15 | 20 | Agilent; PMC7407152 |
| 4-ethylguaiacol | C₉H₁₂O₂ | 152.19 | 13.1 | 152 | 137 | 122 | 15 | 20 | Agilent; PMC7407152 |
| 4-ethylphenol | C₈H₁₀O | 122.16 | 14.4 | 122 | 107 | 77 | 15 | 25 | Agilent; PMC7407152 |
| p-cresol | C₇H₈O | 108.14 | 12.7 | 108 | 107 | 77 | 15 | 25 | PMC7407152 |
| o-cresol | C₇H₈O | 108.14 | 13.5 | 108 | 107 | 77 | 15 | 25 | PMC7407152 |
| m-cresol | C₇H₈O | 108.14 | 13.6 | 108 | 107 | 77 | 15 | 25 | PMC7407152 |
| Syringol | C₈H₁₀O₃ | 154.16 | 15.3 | 154 | 139 | 124 | 15 | 20 | Agilent; PMC7407152 |
| 4-methylsyringol | C₉H₁₂O₃ | 168.19 | 16.1 | 168 | 153 | 125 | 15 | 20 | Agilent; PMC7407152 |
| d3-guaiacol (internal std) | — | 127 | 11.3 (co-elute) | 127 | 109 | — | 15 | — | PMC7407152 |

The 124→109 (guaiacol) loss of CH₃ from the methoxy is essentially universal across all published methods — Agilent, Thermo (TSQ AppsLab 4532 *Rapid smoke-taint analysis*), and the AWRI Hayasaka group all converge on the same precursor → product pair. The qual 124→81 (further loss of CO + CH·) is more CE-sensitive and the reported qual transition varies between m/z 81 and m/z 95 depending on the optimization paper — both are acceptable.

### 6.5 Monoterpenes (floral / citrus markers — Cabernet Sauvignon background)

| Compound | Formula | MW | RT (DB-WAX, min) | Precursor | Quan | Qual | CE Quan | CE Qual | Source |
|---|---|---|---|---|---|---|---|---|---|
| Linalool | C₁₀H₁₈O | 154.25 | 18–20 | 93 | 41 | 69 | 10 | 15 | NIST mass spectra + Lukić 2019; Vichi 2007 |
| Linalool (alt mol-ion ladder) | — | — | — | 121 | 93 | — | 10 | — | NIST |
| α-terpineol | C₁₀H₁₈O | 154.25 | 22–24 | 121 | 93 | 81 | 10 | 15 | NIST + literature consensus |
| α-terpineol (alt) | — | — | — | 136 | 121 | — | 10 | — | NIST |
| Geraniol | C₁₀H₁₈O | 154.25 | 28–30 | 69 | 41 | — | 10 | — | NIST |
| Nerol | C₁₀H₁₈O | 154.25 | 27–29 | 69 | 41 | — | 10 | — | NIST |

Note: terpene quan/qual pairs in EI-CID typically rely on low-mass tropylium-like or methyl-allyl cations (m/z 41, 69, 93) — these are non-specific in matrices with high ester load, which is why HS-SPME-Arrow + retention-time confirmation is essential rather than relying on m/z selectivity alone.

### 6.6 Esters (yeast-derived fruit volatiles)

| Compound | Formula | MW | RT (DB-WAX, min) | Precursor | Quan | Qual | CE Quan | CE Qual | Source |
|---|---|---|---|---|---|---|---|---|---|
| Isoamyl acetate (3-methylbutyl acetate) | C₇H₁₄O₂ | 130.18 | 8–10 | 70 | 43 | — | 10 | — | NIST + Hjelmeland; alt 87→43 reported |
| Hexyl acetate | C₈H₁₆O₂ | 144.21 | 12–13 | 84 | 43 | 56 | 10 | 10 | NIST; Lukić 2019 |
| Ethyl hexanoate | C₈H₁₆O₂ | 144.21 | 12–13 | 88 | 60 | 99 | 10 | 5 | NIST + ester-class consensus |
| Ethyl octanoate | C₁₀H₂₀O₂ | 172.27 | 16–18 | 88 | 60 | 127 | 10 | 5 | NIST + ester-class consensus |
| Ethyl decanoate | C₁₂H₂₄O₂ | 200.32 | 21–23 | 88 | 60 | 155 | 10 | 5 | NIST + literature |
| Ethyl butanoate | C₆H₁₂O₂ | 116.16 | 6–8 | 88 | 60 | 71 | 10 | 10 | NIST |
| 2-phenethyl acetate | C₁₀H₁₂O₂ | 164.20 | 27–29 | 104 | 91 | 105 | 10 | 5 | NIST |

The McLafferty rearrangement product (m/z 88 + 1 H) is the universal ester quantifier ion across ethyl esters of fatty acids (Wikipedia/LibreTexts *Mass Spectrometry — Fragmentation Patterns*; whitman.edu MS e-book §5.2.1.2). The 88→60 transition (further loss of CO) is the standard MS/MS quan.

### 6.7 Higher alcohols and aldehydes

| Compound | Formula | MW | RT (DB-WAX, min) | Precursor | Quan | Qual | CE Quan | CE Qual | Source |
|---|---|---|---|---|---|---|---|---|---|
| 2-phenylethanol | C₈H₁₀O | 122.17 | 24–26 | 91 | 65 | 92 | 15 | 10 | NIST; Razaq 2023 |
| Hexanal (with PFBHA derivatization) | C₆H₁₂O | 100.16 | 37.2 (DB-5) | 181 | 161 | 239→181 | 15 | 15 | Razaq 2023 JAFC PMC10835727 |
| Methional (3-(methylthio)propanal, PFBHA) | C₄H₈OS | 104.17 | 43.6 (DB-5) | 252 | 181 | 252→252 | 10 | 5 | Razaq 2023 |
| Phenylacetaldehyde (PFBHA) | C₈H₈O | 120.15 | 49.8 (DB-5) | 181 | 161 | 91→65 | 5 | 15 | Razaq 2023 |
| 2-methylpropanal (PFBHA) | C₄H₈O | 72.11 | 20.5 (DB-5) | 267 | 181 | 250→181 | 10 | 10 | Razaq 2023 |
| 3-methylbutanal (PFBHA) | C₅H₁₀O | 86.13 | 28.7 (DB-5) | 239 | 181 | 239→207 | 15 | 5 | Razaq 2023 |
| Benzaldehyde (PFBHA) | C₇H₆O | 106.12 | 47.5 (DB-5) | 301 | 181 | 301→271 | 15 | 5 | Razaq 2023 |

Note: aldehydes are usually derivatized in situ with PFBHA (O-(2,3,4,5,6-pentafluorobenzyl)hydroxylamine) on the SPME fiber to form pentafluorobenzyl oximes, which fragment cleanly through the perfluoroaryl moiety to m/z 181 — the universal quantifier ion for the entire PFBHA-derivatized aldehyde panel. Without derivatization, the underivatized hexanal molecular ion at m/z 100 is very weak and the fragments at m/z 44, 56, 41 are non-specific (Razaq et al., *J. Agric. Food Chem.* 71, 1456, 2023; methods adopted from earlier Strecker-aldehyde / Maillard work).

### 6.8 Earthy / musty off-flavors

| Compound | Formula | MW | RT (DB-WAX, min) | Precursor | Quan | Qual | CE Quan | CE Qual | Source |
|---|---|---|---|---|---|---|---|---|---|
| Geosmin | C₁₂H₂₂O | 182.30 | 23–25 | 112 | 97 | 83 | 10 | 15 | Agilent 5991-1031EN; Restek SPME-Arrow drinking-water app note |
| Geosmin (mol-ion) | — | — | — | 182 | 112 | 167 | 10 | 5 | Agilent 5991-1031EN |
| 2-MIB (2-methylisoborneol) | C₁₁H₂₀O | 168.28 | 19–21 | 95 | 67 | — | 15 | — | Agilent 5991-1031EN |
| 2-MIB (alt) | — | — | — | 108 | 95 | — | 5 | — | Agilent |

### 6.9 Sulfur volatiles

| Compound | Formula | MW | RT (DB-WAX, min) | Precursor | Quan | Qual | CE Quan | CE Qual | Source |
|---|---|---|---|---|---|---|---|---|---|
| DMS (dimethyl sulfide) | C₂H₆S | 62.13 | 3–4 | 62 | 47 | 45 | 10 | 15 | Razaq 2023; Cucchiella 2024 *Beverages* |
| DMDS (dimethyl disulfide) | C₂H₆S₂ | 94.20 | 7–8 | 94 | 79 | 45 | 10 | 15 | Cucchiella 2024 |
| 3MH (3-mercaptohexan-1-ol, PFB-derivatized via extractive alkylation) | C₆H₁₄OS | 134.24 | 14–16 (post-derivatization, varies) | 134 | 100 | 67 | 5 | 10 | Mateo-Vivaracho 2008/2010 *J. Chromatogr. A*; PMC6332517 |
| 3MH (alt direct) | — | — | — | 67 | 41 | — | 10 | — | NIST EI |
| 4-MMP (4-mercapto-4-methyl-2-pentanone, PFB-derivatized) | C₆H₁₂OS | 132.22 | 14–16 (varies) | 132 | 75 | 99 | 10 | 5 | Mateo-Vivaracho; PMC6332515 |

Note on 3MH: native unmodified 3MH gives weak EI signal (LOD ~50 ng/L). The community standard is PFB-derivatization via extractive alkylation (3MH-S-CH₂-C₆F₅) followed by SPME-GC-MS/MS in NCI mode — but NCI MRM is platform-specific. Alternative HPLC-MS/MS with DTDP (4,4'-dithiodipyridine) or OPA-mediated derivatization is now more common and gives 0.5–1 ng/L LOQ (Capone et al., *Anal. Chem.* 2015 — DTDP method). The GC route remains valid for SPME-Arrow workflows that bypass HPLC entirely.

### 6.10 Oxidative aging markers (additional for white blend comparison, not in primary Cabernet panel)

| Compound | Formula | MW | RT (DB-WAX, min) | Precursor | Quan | Qual | CE Quan | CE Qual | Source |
|---|---|---|---|---|---|---|---|---|---|
| Sotolon | C₆H₈O₃ | 128.13 | — (LC-MS often preferred) | 129 | 55 | 83 | 15 | 15 | *Foods* 2019 PMC6311786 |
| Furaneol | C₆H₈O₃ | 128.13 | — | 128 | 81 | — | — | — | Ferreira/Jarauta 2003 |

---

## 7. Per-vendor app-note summary

### 7.1 Agilent (7000D, 7000E, 7010 Triple Quadrupole GC/MS)

**Wine smoke-taint app note 5994-3161EN** (2021) — the canonical reference for free volatile phenols in wine. Uses Agilent 8890 GC, 7000D MS, MMI inlet, SPME-Arrow with DVB/CWR/PDMS coating. Tabulates 15 eV CE across guaiacol, methylguaiacol, ethylguaiacol, ethylphenol, syringol, cresol isomers, with d3-guaiacol IS. Method is the de-facto standard adopted by AWRI and Australian commercial labs after 2022.

**TCA in wine app note 5990-4968EN** — uses HS-SPME at 0.1 ng/L LOD with the 210→195 and 212→197 isotopologue confirmation. Sensitive enough that 1 ng/L gives S/N ~30. Demonstrates the value of QqQ over high-res TOF for routine throughput.

**Volatile profiling of Cabernet Sauvignon 5991-3682EN** — broader untargeted/targeted hybrid panel using SPME + GC-MS/MS; tabulates ester, terpene, norisoprenoid transitions but uses a single-CE (10 eV) compromise rather than per-transition optimization.

**MassHunter Optimizer for GC Triple Quad** (G7003-90009) — workflow for automated CE optimization (the "AutoSRM" of Agilent's stack). Standard procedure: precursor scan → product scan → CE ramp 0–60 eV in 5 eV steps → top-3 pairs selected automatically per analyte.

### 7.2 Thermo Scientific (TSQ 8000 Evo, TSQ 9000, TSQ 9610)

**TSQ AppsLab #4532 "Rapid smoke-taint analysis of wine with SPME-GC-MS/MS"** — direct counterpart to the Agilent 5994-3161EN method on the TSQ Quantum / 9000 platform. Same volatile-phenol panel, similar 15 eV CE across all transitions, with TraceGOLD TG-WAXMS column.

**Razaq et al. 2023 *J. Agric. Food Chem.*** (PMC10835727) — fully-automated HS-SPME-GC-MS/MS for 44 volatile carbonyls on a TSQ Quantum XLS Ultra; uses PFBHA derivatization at the fiber, the universal m/z 181 quantifier ion across the panel. Tabulates per-compound CE 5–15 eV. Key reference for hexanal and Strecker-aldehyde transitions in §6.7.

### 7.3 Shimadzu (GCMS-TQ8050NX, GCMS-TQ8040NX)

**Smart MRM** — proprietary scheduled-MRM algorithm that automatically sets the optimum dwell and loop times per method. Vendor app notes focus heavily on pesticide multi-residue and PAH panels (lots of overlap with wine matrix complexity) but the smoke-taint and methoxypyrazine wine methods are published as Shimadzu Application News Series (numbers vary; Asia-Pacific portal). CE conventions match Agilent: per-transition optimization from precursor scan + product scan + CE ramp.

### 7.4 Sciex (7500 GC, MultiQuant, Scheduled MRM Pro algorithm)

**Scheduled MRM Pro Algorithm Overview** (RUO-MKT-02-8539-A) — the most-cited reference for the math underlying scheduled MRM. Documents the dwell-weighting scheme where high-abundance transitions get less dwell time so weak transitions can get more, keeping CT fixed.

Sciex's GC-MS/MS uptake in wine analysis is smaller than Agilent/Thermo (their wine-specific app notes are mainly LC-MS/MS for pesticide residues and mycotoxins). The 7500 GC is suitable for the panel here but no public wine SPME app note exists at time of writing.

### 7.5 Peer-reviewed key methods

- **Pollnitz, Pardon, Sefton, Sefton (1996)** *Aust. J. Grape Wine Res.* — original TCA-d5 stable-isotope-dilution GC-MS method; the framework on which every subsequent TCA QqQ analysis builds.
- **Hjelmeland & Ebeler (2010)** — IBMP/IPMP HS-SPME GC-MS/MS quantification, basis for the methoxypyrazine transitions in §6.1.
- **Slabizki & Schmarr (2014)** *Eur. Food Res. Tech.* — MDGC-MS/MS for German Sauvignon blanc methoxypyrazines.
- **Hayasaka, Baldock, Pollnitz (AWRI ~2003–2013)** — HPLC-MS/MS for smoke-taint glycosides; GC-MS/MS for free volatile phenols (the AWRI lab method).
- **Razaq et al. (2023)** *J. Agric. Food Chem.* — automated PFBHA-derivatized HS-SPME-GC-MS/MS for 44 wine carbonyls.
- **Schreier & Wuest (2016)** *Anal. Bioanal. Chem.* 408, 7593 — multidimensional GC-MS/MS for α-/β-ionone and β-damascenone with optimized SRMs at 10–20 eV.
- **Sacks et al. (2012)** *J. Agric. Food Chem.* 60, 2998 — TDN sensory thresholds and analytical method in Riesling wines.
- **Mateo-Vivaracho, Cacho, Ferreira (2008–2010)** *J. Chromatogr. A* — wine polyfunctional thiol PFB-derivatization HS-SPME-GC-MS method (basis for the 3MH transition in §6.9).

---

## 8. References

1. Pollnitz, A. P.; Pardon, K. H.; Sefton, M. A. *Aust. J. Grape Wine Res.* **1996**, *2*(2), 92–96. The analysis of 2,4,6-trichloroanisole and other chloroanisoles in tainted wines and corks. https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1755-0238.1996.tb00107.x
2. Agilent Technologies. Application Note 5994-3161EN. *Analysis of Free Volatile Phenols in Smoke-Impacted Wines using HS-SPME with 8890 GC and 7000D Triple Quadrupole MS.* 2021. https://www.agilent.com/cs/library/applications/application-smoke-taint-wine-SPME-8890-gc-5994-3161en-agilent.pdf
3. Agilent Technologies. Application Note 5990-4968EN. *Sensitive Detection of Trichloroanisole (TCA) in Wine Using Triple Quadrupole GC/MS.* https://www.agilent.com/Library/applications/5990-4968EN.pdf
4. Agilent Technologies. Application Note 5991-3682EN. *Volatile Profiling of U.S. Cabernet Sauvignon Wines Using HS-SPME with Triple Quadrupole GC-MS.* 2013.
5. Agilent Technologies. Application Note 5991-1031EN. *Sensitive Detection of 2-MIB and Geosmin in Drinking Water.* https://www.agilent.com/cs/library/applications/5991-1031EN.pdf
6. Wilkinson, K. L.; et al. (Hayasaka group, AWRI). *Metabolites* **2020**, *10*(7), 294. A Simple GC-MS/MS Method for Determination of Smoke Taint-Related Volatile Phenols in Grapes. https://pmc.ncbi.nlm.nih.gov/articles/PMC7407152/
7. Razaq, R.; Ruocco, S.; et al. *J. Agric. Food Chem.* **2023**, *71*, 1456–1469. Development of a Fully Automated HS-SPME-GC-MS/MS Method for Odor-Active Carbonyls in Wines. https://pmc.ncbi.nlm.nih.gov/articles/PMC10835727/
8. Schreier, P.; Wuest, M. *Anal. Bioanal. Chem.* **2016**, *408*, 7593–7602. Quantitative determination of α-ionone, β-ionone, and β-damascenone using MDGC-MS/MS. https://link.springer.com/article/10.1007/s00216-016-9767-6
9. Slabizki, P.; Schmarr, H.-G. *Eur. Food Res. Technol.* **2014**, *239*, 549–558. Quantitative analysis of 3-alkyl-2-methoxypyrazines in German Sauvignon blanc wines by MDGC-MS/MS. https://link.springer.com/article/10.1007/s00217-014-2250-8
10. Hjelmeland, A. K.; et al. *J. Agric. Food Chem.* **2010**, *58*, 8385–8391. Implementation of HS-SPME-GC-MS/MS for 3-alkyl-2-methoxypyrazines in wine.
11. Sacks, G. L.; et al. *J. Agric. Food Chem.* **2012**, *60*, 2998–3004. Sensory Threshold of TDN in Riesling and non-Riesling Wines.
12. Mateo-Vivaracho, L.; Cacho, J.; Ferreira, V. *J. Chromatogr. A* **2008**, *1185*, 9–18. Quantitative determination of wine polyfunctional mercaptans at ng/L level by GC-NCI-MS as PFB derivatives.
13. Capone, D. L.; Barker, A.; et al. *Anal. Chem.* **2015**, *87*, 1226–1231. Simple Quantitative Determination of Potent Thiols at Ultratrace Levels in Wine by DTDP-HPLC-MS/MS. https://pubs.acs.org/doi/10.1021/ac503883s
14. Carrillo, J. D.; et al. *Chardonnay wine HS-SPME-GC-MS/MS multivariate optimization.* https://www.scirp.org/journal/paperinformation?paperid=71467
15. *Single to Triple: Fundamentals and Modes of Bench-Top GC-MS/MS.* LCGC International. https://www.chromatographyonline.com/view/single-to-triple-fundamentals-and-modes-of-bench-top-gas-chromatography-triple-quadrupole-mass-spectrometry-gc-ms-ms-
16. Sciex Tech Note. *The Scheduled MRM Algorithm Pro.* RUO-MKT-02-8539-A. https://sciex.com/content/dam/SCIEX/pdf/tech-notes/all/Scheduled-MRM-Pro-Overview-RUO-MKT-02-8539-A.pdf
17. Sciex Knowledge Base. *How to choose the appropriate dwell time/cycle time in an MRM and sMRM method to get enough points per chromatographic peak.* https://sciex.com/support/knowledge-base-articles/how-to-choose-the-appropriate-dwell-time-cycle-time-in-an-mrm-and-smrm-method-en_us
18. Agilent Technologies. Technical Overview 5990-3595EN. *New Dynamic MRM Mode Improves Data Quality and Dwell Time Flexibility.* https://www.agilent.com/cs/library/technicaloverviews/public/5990-3595en_lo%20CMS.pdf
19. Sepscience. *Understanding MRM Fundamentals: Dwell Time, Cycle Time, and Duty Cycle.* https://www.sepscience.com/understanding-mrm-fundamentals-dwell-time-cycle-time-and-duty-cycle-10389
20. European Commission DG SANTE. *SANTE/11312/2021: Analytical Quality Control and Method Validation Procedures for Pesticide Residues Analysis in Food and Feed.* In force 1 Jan 2022. https://food.ec.europa.eu/system/files/2023-11/pesticides_mrl_guidelines_wrkdoc_2021-11312.pdf
21. FDA Center for Veterinary Medicine. *Guidance for Industry: Mass Spectrometry for Confirmation of the Identity of Animal Drug Residues.* 2003. https://www.fda.gov/media/70154/download
22. FDA. *Bioanalytical Method Validation: Guidance for Industry.* May 2018. https://www.fda.gov/files/drugs/published/Bioanalytical-Method-Validation-Guidance-for-Industry.pdf
23. AOAC International. *Official Method 2007.01: Pesticide Residues in Foods by Acetonitrile Extraction and Partitioning with Magnesium Sulfate Gas Chromatography/Mass Spectrometry and Liquid Chromatography/Tandem Mass Spectrometry.* https://nucleus.iaea.org/sites/fcris/Shared%20Documents/SOP/AOAC_2007_01.pdf
24. Wikipedia contributors. *Selected reaction monitoring.* https://en.wikipedia.org/wiki/Selected_reaction_monitoring
25. Wikipedia contributors. *Quadrupole mass analyzer.* https://en.wikipedia.org/wiki/Quadrupole_mass_analyzer
26. Wikipedia contributors. *Collision-induced dissociation.* https://en.wikipedia.org/wiki/Collision-induced_dissociation
27. Pfeiffer Vacuum. *Quadrupole Mass Spectrometers (QMS) — Pfeiffer Know-How.* https://www.pfeiffer-vacuum.com/global/en/knowledge/vacuum-technology/knowledge-book/6-mass-spectrometers-and-residual-gas-analysis/6_3_quadrupole_mass_spectrometers/
28. SIMION 2024 Supplemental Docs. *Mathieu Equation and Stability Diagram.* https://simion.com/info/mathieu_equation.html
29. Process Insights / Extrel Application Note RA_2010A. *Practical Quadrupole Theory: Graphical Theory.* https://www.process-insights.com/wp-content/uploads/2022/06/Process-Insights_Extrel_Practical-Quadrupole-Theory-Graphical-Theory.pdf
30. Konenkov, N. V.; Sudakov, M. *J. Am. Soc. Mass Spectrom.* **2002**, *13*(6), 597–613. Matrix methods for the calculation of stability diagrams in quadrupole mass spectrometry.
31. PMC7249026. *Center-of-Mass iso-Energetic Collision-Induced Decomposition in Tandem Triple Quadrupole Mass Spectrometry.* https://pmc.ncbi.nlm.nih.gov/articles/PMC7249026/
32. Révész, Á.; et al. *Mass Spectrom. Rev.* **2023**, *42*, 1247. Collision energies: Optimization strategies for bottom-up proteomics.
33. NIST Mass Spectrometry Data Center / NIST WebBook. *Linalool* (CID 6549), *Phenylethyl Alcohol* (CAS 60-12-8), *Ethyl octanoate* (CAS 106-32-1), *Ethyl hexanoate* (CAS 123-66-0), *1-Hexanol* (CAS 111-27-3), and related EI mass spectra. https://webbook.nist.gov/
34. Agilent MassHunter Optimizer for GC Triple Quad Quick Start Guide (G7003-90009). https://www.agilent.com/cs/library/usermanuals/public/G7003-90009.pdf
35. PMC6332517. *Quantification of Polyfunctional Thiols in Wine by HS-SPME-GC-MS Following Extractive Alkylation.* https://ncbi.nlm.nih.gov/pmc/articles/PMC6332517
36. *Foods* / PMC6311786. *Rapid Determination of Sotolon in Fortified Wines Using LC-MS/MS.* https://pmc.ncbi.nlm.nih.gov/articles/PMC6311786/
37. PMC10054257. *Uncorking Haloanisoles in Wine.* https://pmc.ncbi.nlm.nih.gov/articles/PMC10054257/
38. Chatonnet, P.; Bonnet, S.; et al. — pioneers of 4-ethylphenol / 4-ethylguaiacol as Brettanomyces markers; standard MRM transitions per Agilent 5994-3161EN.
39. Cucchiella, V.; et al. *Beverages* **2024**, *10*, 15. Wine Cork Closures Impacts on Dimethyl Sulfide and Precursors. https://www.mdpi.com/2306-5710/9/1/15
40. Whitman College Mass Spectrometry e-book §5.5.2 *Quadrupole mass filter.* http://people.whitman.edu/~dunnivfm/C_MS_Ebook/CH5/5_5_2.html

---

*End of Chapter 4 — GC-MS/MS MRM Optimization.*
