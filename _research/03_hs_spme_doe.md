# HS-SPME DoE and Optimization — Strategy Brief

**Specialist 3/6 — HS-SPME Extraction Optimization**
Dissertation: *Method Development, MRM Optimization, and Validation of SPME Arrow Headspace GC-MS/MS for the Quantitative Profiling of Aroma-Active Volatiles in Cabernet Sauvignon Wine.*
Compiled 2026-05-12 from primary literature (Saha 2018; Šikuten 2021; Welke 2024 OENO One; Rossi 2023; Burin 2013; Arcanjo 2015; Zhao 2019 Baijiu CCD; Siebert 2005 SIDA; Antalick 2014; Pizarro 2009; Bee-DiGregorio 2009; Câmara 2006; IOFI Working Group 2010; Stilo 2021 critical review).

---

## 1. Screening designs (Plackett-Burman, fractional factorial)

### 1.1 Why screen first

For an SPME Arrow method on Cabernet Sauvignon, the candidate factor set is large (≥7): extraction temperature, extraction time, NaCl concentration, sample volume / headspace ratio, ethanol dilution, pH, agitation speed, pre-incubation time. A full factorial at 2 levels = 2^7 = 128 runs; even k=5 → 32 runs. Plackett-Burman (PB) and 2-level fractional factorials (FrF) reduce this to 8-12 runs to identify significant main effects, then commit to a response-surface design (BBD / CCD) for the 3-4 surviving variables.

### 1.2 Plackett-Burman in wine SPME — literature exemplars

**Salinari et al. / Kallithraka group on Greek primary aromatics (Sereli et al. 2013).** A 2-level **Plackett-Burman** design screened the SPME method for free and bound primary aromatics. *Significant variables identified:* NaCl concentration, ethanol concentration in the sample, extraction time, pH, extraction temperature. Surviving variables were further optimised with a modified Simplex procedure. (Sereli, Kallithraka et al., *J. Chromatogr. A* 2013).

**Geometric properties.** PB designs are resolution-III, examine N-1 factors in N runs where N is a multiple of 4. The most commonly used in volatile analysis is **PB-12** (11 factors, 12 runs) or PB-8 (7 factors, 8 runs). PB cannot distinguish main effects from 2-factor interactions, so interactions are presumed negligible during screening (a defensible assumption for SPME because partition coefficients are largely additive).

### 1.3 Fractional factorial in wine SPME

**Two-level 2^4 full factorial** was used for oak-derived volatiles in wine to establish the influence of fiber type, temperature, pre-incubation time, and NaCl (cited in Stilo, Bicchi, Cordero 2021 review). For 4 factors this is 16 runs and resolves all 2-factor interactions.

**2-level screening for fiber × condition factors** was used by Rossi et al. 2023 for the Pecorino/Trebbiano study with CW/DVB and DVB/CAR/PDMS fibers, then BBD for PDMS.

### 1.4 Recommended screening factor list for a SPME Arrow / Cab Sauv project

| # | Factor | Low (-1) | High (+1) | Rationale |
|---|---|---|---|---|
| 1 | Extraction temperature (°C) | 30 | 60 | Equilibrium driver; pyrazine artifact risk >60 |
| 2 | Extraction time (min) | 15 | 45 | Pre-equilibrium vs equilibrium |
| 3 | NaCl (g per 10 mL) | 0 | 3.0 | Salting-out range from lit. |
| 4 | Ethanol of diluted sample (% v/v) | 5 | 13 | Matrix matching range |
| 5 | Pre-incubation time (min) | 5 | 20 | Headspace equilibration |
| 6 | Sample volume (mL in 20 mL vial) | 5 | 10 | Headspace/sample ratio |
| 7 | Agitation (rpm) | 250 | 500 | Mass transfer rate |

PB-8 (or equivalent FrF 2^(7-4)) suffices to identify the top 3-4 drivers before committing to RSM.

---

## 2. Optimization designs (BBD vs CCD)

### 2.1 Comparison at a glance

| Property | Box-Behnken (BBD) | Central Composite (CCD) |
|---|---|---|
| Levels per factor | 3 (-1, 0, +1) | 5 (-α, -1, 0, +1, +α) typical CCC; 3 if face-centered (α=1) |
| Total runs (k=3) | 12 + center reps (typ. 15) | 8 factorial + 6 axial + center reps (typ. 20) |
| Total runs (k=4) | 24 + 3-5 center = 27-29 | 16 factorial + 8 axial + 3-6 center = 27-30 |
| Rotatability | Nearly rotatable for k=3; not for k=4 | Rotatable when α = (n_factorial)^(1/4) |
| Operability | All points within ±1 hypercube — no extreme runs | Axial points exceed ±1 hypercube unless α=1 (face-centered) |
| Best use | When extreme runs (e.g., very high temperature → solvent loss) are infeasible | When operating window allows axial excursions; gives better edge prediction |

For SPME Arrow on wine, **BBD is the dominant choice in the literature** (Saha 2018; Rossi 2023; Šikuten 2021; Welke OENO One 2024) because:
1. The factor extremes are physically risky (60 °C × 60 min × 0% ethanol matrix would dry fiber / overshoot pyrazine artifact threshold).
2. BBD is statistically efficient — fewer runs than CCC for k=3 (15 vs 20).
3. Three levels per factor are sufficient to fit a second-order polynomial in a smooth response region.

**CCD/CCC is preferred when** the screening identifies linearity that needs axial verification, or when extrapolating to robustness ranges. Baijiu work by Zhao 2019 used CCD precisely because they wanted to verify behavior past the safe envelope.

### 2.2 Reference BBD designs in wine SPME literature

**Saha et al. 2018 — *Foods* 7(8):127. doi:10.3390/foods7080127.**
3-factor BBD, 15 runs.
- Extraction temperature: 30 / 50 / 70 °C
- Extraction time: 15 / 30 / 45 min
- Sample volume: 7 / 10 / 13 mL in 20 mL vial
- No NaCl in this design
- Response: peak area of 23 impact odorants across compound classes (esters, alcohols, acids, terpenes, sulfur, phenolic)
- *Optimum* (depends on ethanol level): 48 °C / 29 min / vol-adjusted (5% v/v ethanol); 46 °C / 43 min (8% v/v); 44 °C / 42 min (13% v/v).
- *Key finding:* extraction temperature was the dominant factor; the optimum shifts as ethanol rises (the depressive ethanol effect).

**Šikuten et al. 2021 — *Molecules* 26(23):7409. doi:10.3390/molecules26237409.**
SPME-Arrow + 4-factor BBD for free VOCs from Merlot grape skins (transferable principle).
- Temperature: 40 / 50 / 60 °C
- Incubation: 10 / 20 / 30 min
- Exposure: 30 / 45 / 60 min
- Desorption: 5 / 7.5 / 10 min
- *Optimum free VOCs:* 60 °C / 20 min incub / 49 min exp / 7 min desorb on DVB/CWR/PDMS Arrow.

**Welke et al. 2024 — OENO One v.58. doi:10.20870/oeno-one.2024.58.X (article view/7914).**
SPME Arrow + BBD for wine VOCs.
- Temperature: 40 / 50 / 60 °C
- Incubation: 10 / 20 / 30 min
- Exposure: 30 / 45 / 60 min
- *Optimum (white wine):* 50 °C / 10 min incub / 60 min exposure
- *Optimum (red wine):* 60 °C / 17 min incub / 53 min exposure
- *Direct relevance for Cabernet Sauvignon* — this is the closest published SPME Arrow / red wine optimization.

**Burin et al. 2013 — *Talanta* 117:87-93. doi:10.1016/j.talanta.2013.08.037.**
Heterocyclic compounds (furans, thiophenes, thiazoles, pyrazines) in wine with **CCD** + RSM (not BBD).
- Variables optimized: pH, NaCl, extraction time.
- Validated repeatability 2.7-12% RSD across 29 French wines.

**Rossi et al. 2023 — *Molecules* 28(4):1534. doi:10.3390/molecules28041534.**
PDMS fiber + BBD on white wine.
- Temperature: 30 / 40 / 50 °C
- NaCl: 10 / 20 / 30 % w/v
- Exposure: 10 / 20 / 30 min
- Two-level designs for CW/DVB and DVB/CAR/PDMS fibers (screening before commit).
- Internal standards: ethyl decanoate (apolar marker), 3-methyl-1-butanol (polar marker).
- *Two local maxima* on the response surface: (10% NaCl, 50 °C) and (30% NaCl, 30 °C) — illustrates the salt × temperature antagonism that classical RSM exposes.

**Arcanjo et al. 2015 — *Food Sci. Technol. Campinas* 35(4). doi:10.1590/1678-457X.6815.**
Rotatable CCD (RCCD), 17 runs (8 factorial + 6 axial + 3 center), 3 variables:
- Equilibrium time: 7 / 15 / 23 min (axial 5 / 25)
- Extraction time: 10 / 35 / 60 min
- Temperature: 13 / 30 / 47 °C
- *Optimum:* 15 min equilibrium / 35 min extraction / 30 °C on CAR/PDMS/DVB
- Model R² = 0.9785 for isoamyl acetate.

**Zhao et al. 2019 — *J. Chromatogr. A* 1610:460537 (PMID 31607446).**
SPME Arrow + CCD on Baijiu. 5 fiber chemistries screened, then CCD on the winner (DVB/CAR/PDMS).
- *Optimum:* 5 mL Baijiu diluted to 10% ethanol + 1.5 g NaCl + 45 °C + 45 min on DVB/CAR/PDMS.
- 82 aroma compounds identified (esters, alcohols, acids, aldehydes, furans, pyrazines, sulfur, phenols, terpenes, lactones).

### 2.3 BBD design table (k=3, 15 runs) — example for the dissertation

| Run | A (Temp °C) | B (Time min) | C (NaCl g/10 mL) |
|---|---|---|---|
| 1 | -1 (40) | -1 (30) | 0 (1.5) |
| 2 | +1 (60) | -1 (30) | 0 (1.5) |
| 3 | -1 (40) | +1 (60) | 0 (1.5) |
| 4 | +1 (60) | +1 (60) | 0 (1.5) |
| 5 | -1 (40) | 0 (45) | -1 (0) |
| 6 | +1 (60) | 0 (45) | -1 (0) |
| 7 | -1 (40) | 0 (45) | +1 (3) |
| 8 | +1 (60) | 0 (45) | +1 (3) |
| 9 | 0 (50) | -1 (30) | -1 (0) |
| 10 | 0 (50) | +1 (60) | -1 (0) |
| 11 | 0 (50) | -1 (30) | +1 (3) |
| 12 | 0 (50) | +1 (60) | +1 (3) |
| 13 | 0 (50) | 0 (45) | 0 (1.5) |
| 14 | 0 (50) | 0 (45) | 0 (1.5) |
| 15 | 0 (50) | 0 (45) | 0 (1.5) |

(Three center replicates give pure-error estimate for lack-of-fit testing.)

---

## 3. Response variables and desirability functions

### 3.1 Single-response is rare; multi-response is standard

A wine SPME method usually quantifies ≥20 compounds across ester, alcohol, acid, terpene, norisoprenoid, methoxypyrazine, and volatile-sulfur classes. Optima for each class diverge — Welke (OENO One 2024) shows the optimum extraction temperature is 40 °C for terpenes, 44 °C for alcohols, 54 °C for acids, and 60 °C for esters and norisoprenoids. A single "global" optimum requires reconciling these via a multi-response strategy.

### 3.2 Derringer-Suich desirability

**Definition (Derringer & Suich 1980, *J. Quality Tech.* 12:214-219).**
For each response Y_i(x), define a desirability d_i(Y_i) ∈ [0,1] where 0 = worst, 1 = best:
- Maximization: d_i = ((Y_i − L_i) / (T_i − L_i))^r_i for L ≤ Y ≤ T; 0 below L; 1 above T.
- Target: piecewise, peaks at the target value.
- Minimization: mirror image of the maximization form.

The overall desirability D is the geometric mean:
**D = (d_1 · d_2 · ... · d_n)^(1/n)**, optionally with weights w_i in the exponent.

D drops to 0 if *any* response is unacceptable — geometric mean is intentional, it punishes outlier-bad responses unlike arithmetic mean.

### 3.3 Application in wine SPME

Saha 2018 used a desirability function across 23 odorants, with the parameter r set to 1 (linear desirability) for most analytes and r > 1 for those with consistently low responses (down-weighted). The resulting "global optimum" found one set of T/t/V conditions that maximized the geometric mean of peak areas across all 23.

Rossi 2023 saw two desirability maxima (10% NaCl/50 °C vs 30% NaCl/30 °C), illustrating that the optimum can be multimodal — both deliver acceptable global D but emphasize different compound classes.

### 3.4 Weighting strategies

For a *quantitative* dissertation method, the desirability calc should weight:
- High weight: target/marker compounds with regulatory or sensory importance (e.g., IBMP for "green pepper" character in Cab Sauv; β-damascenone; rotundone if relevant).
- Medium weight: fermentation esters and aliphatic alcohols (background but quantitative).
- Low weight: high-abundance compounds that always saturate (e.g., ethyl octanoate at 50 °C × 45 min on DVB/CAR/PDMS).
- Skip: compounds with response factors so low that no condition will deliver adequate sensitivity (escalate to selective enrichment instead).

---

## 4. Reported factor ranges (table)

Compiled from the cited literature. "Typical optimum" is the cluster-center of the reported optima in red-wine / dark-fruit matrices, where applicable.

| Factor | Reported low | Reported high | Typical optimum (red wine) | Notes / source |
|---|---|---|---|---|
| **Extraction temperature (°C)** | 25 | 70 | 45-55 (compound-class dependent) | Saha 2018 (30-70); Welke 2024 (40-60); Šikuten 2021 (40-60); Burin 2013; >60 risks pyrazine artifact and may demolish heat-labile thiols |
| **Extraction time (min)** | 5 | 90 | 30-45 (pre-equilibrium common) | Saha 2018 (15-45); Welke 2024 (30-60); Arcanjo 2015 (10-60); Antalick 2014 = 30; Bee-DiGregorio = 40 |
| **Pre-incubation / equilibrium (min)** | 0 | 30 | 10-20 | Welke 2024 (10-30); Arcanjo 2015 (7-23, axial to 25); Bee-DiGregorio = 10 |
| **NaCl (g per 10 mL)** | 0 | 5 (saturation ≈ 3.6 at 25 °C) | 1.5-3 (15-30% w/v) | Rossi 2023 (10-30%); Zhao 2019 (1.5 g/5 mL = 30%); Antalick (saturated); Saha 2018 used no NaCl |
| **Sample volume (mL in 20 mL vial)** | 1 | 13 | 5-10 (1:1 to 1:3 sample:HS) | Saha 2018 (7-13); Welke 2024 (5); Zhao 2019 (5 mL diluted) |
| **Headspace/sample ratio** | 0.5 | 19 | 1-3 | Function of sample volume above |
| **Ethanol of analyzed sample (% v/v)** | 5 | 14 | 8-12 (matrix-matched or diluted) | Saha 2018 (5/8/13); Zhao 2019 diluted Baijiu from 50% to 10%; for Cab Sauv typically diluted from 13-14% to 10% for matrix consistency |
| **pH** | 2.5 | 4.0 | 3.2-3.5 (typical Cab Sauv) | Sereli 2013 PB found pH significant; usually left at native wine pH unless headspace acidity needed for specific analyte class |
| **Agitation (rpm)** | 0 | 1000 | 250-500 | Vortex preferred over magnetic stir (IOFI 2010); Antalick = 500; Arcanjo = "constant agitation" |
| **Desorption temperature (°C)** | 230 | 280 | 250-270 | Šikuten 2021; Zhao 2019 (250 °C, 1 min); manufacturer-fiber dependent |
| **Desorption time (min)** | 1 | 10 | 5-7 | Šikuten 2021 (5-10, optimum 7); Welke 2024 (5); typically tested for carryover on second blank desorption |

### Special case: pyrazine artifact

Methoxypyrazines (3-isobutyl-2-MP, 3-isopropyl-2-MP, 3-sec-butyl-2-MP) drive Cab Sauv's "green pepper" character at ng/L sensory threshold. The artifact concern at >60 °C arises because:
1. Maillard chemistry between amino acids and reducing sugars in wine matrix at extended heating can *generate* alkylpyrazines in-vial.
2. Heat-induced ester hydrolysis releases acids that depress fiber affinity for analytes.

Bindon and others (academia.edu paper on HS-SPME-GC-MS/MS implementation) use **40 °C / 40 min** with 10 min pre-incubation specifically to stay below this regime for IBMP/IPMP/SBMP quantification with d3-IBMP internal standard.

For a multi-analyte global method on Cab Sauv, **50-55 °C is the upper safe ceiling** if methoxypyrazines are quantified in the same run.

---

## 5. Extraction isotherms — kinetics per compound class

### 5.1 Compound-class extraction profiles (from Welke 2024 + Šikuten 2021)

| Class | Equilibration behavior | Optimum exposure |
|---|---|---|
| Ethyl esters (ethyl hexanoate, octanoate, decanoate) | Fast equilibrium, ≤20-30 min; long times *decrease* signal due to competitive displacement and fiber overloading | 25-30 min @ 50-60 °C |
| Acetate esters (isoamyl acetate, hexyl acetate) | Fast, ≤20 min | 20-30 min |
| Aliphatic alcohols (isobutanol, isoamyl alcohol, phenylethanol) | Slow, plateau ~50 min | 45-60 min @ 60 °C |
| Fatty acids (hexanoic, octanoic, decanoic) | Slow, plateau ~60 min | 60+ min @ 54 °C (Welke 2024 reports 54 °C optimum) |
| Terpenes (linalool, α-terpineol, geraniol) | Slow, longer fiber loading; benefit from thicker SPME-Arrow phase | 60 min @ 40-44 °C |
| Norisoprenoids (β-damascenone, β-ionone) | Slow, high optimum temp | 50-53 min @ 60 °C |
| Methoxypyrazines (IBMP, IPMP) | Fast at moderate T; concentrated at saturation NaCl + lower vol | 40 min @ 40 °C |
| Volatile sulfur (DMS, H2S, mercaptans) | Very volatile, fast; risk of breakdown >50 °C | 15-25 min @ 35-40 °C; lower T critical |

### 5.2 Pre-equilibrium sampling — tradeoffs

When operating *below* equilibrium time (e.g., 25 min when full equilibrium needs 60), the SPME signal is still proportional to bulk concentration *if* convection (agitation), temperature, and fiber-phase contact are held strictly constant. The pre-equilibrium SPME equation is:
**n = C₀ · A · D / δ · t** for short t,
where A = fiber surface area, D = analyte diffusivity, δ = diffusion-layer thickness (a function of stirring rate). Pre-equilibrium operation requires tight RSD on agitation and incubation timing (typical CV < 2% on the autosampler) — achievable with PAL-style auto SPME but a documented risk for manual SPME.

For the dissertation, **operating at or beyond equilibrium for the dominant analytes (esters, methoxypyrazines)** is the defensible choice; some slower compounds (long-chain acids) will not be fully equilibrated, but their CV remains acceptable when the autosampler timing is tight.

### 5.3 SPME Arrow specific kinetics

Arrow's thicker phase (15.3 µL vs 0.6 µL for classical 100 µm PDMS) means:
- Higher capacity → less risk of saturation at high analyte concentrations.
- Longer equilibration time for compounds whose partition into the polymer is mass-transfer-limited.
- A 60-min Arrow exposure is roughly equivalent in *recovery* to 30-min classical SPME for the same analyte class, with 6-20× higher signal-to-noise (Helin 2015; Restek/PAL bulletins).

---

## 6. Internal standards used in published wine SPME work

### 6.1 Surrogate internal standards (commonly available, single-analog)

| IS | Used for | Citation |
|---|---|---|
| **2-octanol** | General volatile profiling, often as single IS | Many; flagged for fenchone co-elution issue → 3-octanol substituted |
| **3-octanol** | General volatile profiling | Replacement for 2-octanol when terpenes co-elute |
| **4-methyl-2-pentanol** | Polar volatiles (alcohols, esters) | Rossi 2023; Welke 2024; typical concentration 2 mg/L |
| **Methyl isobutyl ketone (MIBK)** | Carbonyl tracker | IOFI 2010 |
| **n-Dodecane / n-Tridecane** | Apolar tracker (esters, hydrocarbons) | IOFI 2010 |
| **Ethyl nonanoate** | Mid-range ester surrogate | IOFI 2010; Antalick database 2014 |
| **Ethyl decanoate** | Apolar volatile marker | Rossi 2023 |
| **3-methyl-1-butanol** | Polar volatile marker | Rossi 2023 |

### 6.2 Deuterated (isotope-dilution) internal standards — SIDA

The gold standard for quantification is **stable isotope dilution analysis (SIDA)**: a deuterated analog of each target analyte is added at known concentration to the sample matrix. Because the deuterated analog has near-identical partition coefficient, ionization response, and chromatographic behavior (with sufficient m/z shift for MS resolution), the analyte-to-IS ratio cancels matrix effects, fiber-to-fiber variation, and instrument drift.

| Analyte class | Deuterated IS used | Reference |
|---|---|---|
| 3-isobutyl-2-methoxypyrazine (IBMP) | d₃-IBMP | Bindon et al. wine MS/MS protocol; Bee-DiGregorio 2009 (whole-berry HS-SPME → GCxGC-TOFMS predicting wine levels) |
| Linalool | d₃-linalool or d₂-linalool | Bee-DiGregorio 2009 |
| β-damascenone | d₄-damascenone | Siebert 2005 ABC 381:937-947 — synthesized for SIDA panel |
| Volatile sulfur compounds (3-MH, 4-MMP, 3-MHA) | d₂- and d₃-deuterated analogs | Mateo-Vivaracho 2006 (referenced in PubMed 12744649); Siebert 2005 |
| Esters (ethyl hexanoate, etc.) | d₅-ethyl esters | Siebert 2005 — 31 compounds; 9 commercial, 22 synthesized |
| Acetate esters | d₃-isoamyl acetate, etc. | Siebert 2005 |
| C13-norisoprenoids (β-ionone) | d₃-β-ionone | Slaghenaufi 2012 PubMed 22980832 — Piedmont wines |

**Siebert et al. 2005 — *Anal. Bioanal. Chem.* 381(4):937-947. doi:10.1007/s00216-004-2992-4.**
The most comprehensive published SIDA protocol for wine SPME, quantifying 31 fermentation products simultaneously. Calibration linearity R² = 0.995-1.000 across full matrix range (ethanol, ionic strength, pH). This is the *canonical* citation for SIDA in wine HS-SPME-GC-MS.

### 6.3 Addition timing

- **For SIDA**: add deuterated IS to the wine sample *before* any sample handling (before NaCl addition, before dilution). The IS goes through every step the analyte goes through.
- **For surrogate IS**: add at the *same step* the analyte is being measured against. If quantitation is by external calibration with IS correction, the IS is added at the sample loading step.
- **Concentration**: typical IS concentration is at mid-range of the analyte calibration curve. Saha 2018, Rossi 2023, Welke 2024 use IS at 0.5-2 mg/L final in sample.
- **Never add the IS to the headspace** — must be in the liquid phase to equilibrate properly.

### 6.4 Recommendation for Cab Sauv SPME-Arrow dissertation

A **tiered IS strategy** is defensible and aligns with current SOTA:

| Priority analyte | IS approach | Justification |
|---|---|---|
| IBMP, IPMP, SBMP (methoxypyrazines) | d₃-IBMP (SIDA) | Sub-ng/L sensory thresholds; matrix-sensitive |
| β-damascenone, β-ionone | d₄-damascenone (SIDA) | Sensory-active norisoprenoids; ethanol-sensitive partition |
| Polyfunctional thiols (3-MH, 4-MMP, 3-MHA) | d₂-3-MH (SIDA) | Trace level + reactivity |
| Major fermentation esters (8-12 compounds) | d₅-ethyl hexanoate as group IS, OR ethyl nonanoate as surrogate | Cost-vs-rigor tradeoff |
| Higher alcohols, acids | 4-methyl-2-pentanol or 2-octanol (surrogate) | High abundance; matrix effects manageable |
| Terpenes (linalool, α-terpineol) | d₃-linalool (SIDA preferred) or 3-octanol (surrogate fallback) | Class-specific |

For 4-5 SIDA standards + 2-3 surrogates, the budget impact is meaningful but defensible for a quantitative dissertation. The full 31-IS panel from Siebert 2005 is gold-standard but probably overkill.

---

## 7. Recommended DoE for the dissertation experimental design

### 7.1 Two-phase strategy

**Phase A — Plackett-Burman screening (12 runs).** Identify which 3-4 of the 7 candidate factors actually drive response. Factor list per §1.4. Response = sum of normalized peak areas of ~20 representative analytes (1-2 per class) on synthetic Cab Sauv matrix (model wine: 12% ethanol, 4 g/L tartaric acid, pH 3.4) spiked with known concentrations of target analytes.

**Phase B — Box-Behnken (15 runs) on surviving 3 factors.** Build a second-order response surface model:
**Y = β₀ + Σβ_i·X_i + Σβ_ii·X_i² + Σβ_ij·X_i·X_j + ε**
Fit and check lack-of-fit using the 3 center replicates. Use Derringer-Suich desirability across analyte classes to identify global optimum.

If the optimum lands at a BBD vertex (edge of design space), escalate to **face-centered CCD** to explore axial extensions in that direction.

### 7.2 Why not a one-shot CCD?

- 17-30 runs is high cost when most factors will end up insignificant.
- The PB screen costs 12 runs but saves the BBD/CCD from carrying 4 useless factors that bloat the design 4-8×.
- Industry consensus (Stilo 2021; Lord & Pawliszyn 1997 SPME guidance) is screen first, optimize second.

### 7.3 Specific Phase B design table

Assuming Phase A surfaces (extraction temperature, extraction time, NaCl) as the three drivers (most likely outcome based on literature consensus):

| Factor | Low (-1) | Center (0) | High (+1) |
|---|---|---|---|
| Temp (°C) | 40 | 50 | 60 |
| Time (min) | 30 | 45 | 60 |
| NaCl (g per 10 mL) | 0 | 1.5 | 3.0 (saturating-ish) |

Run 15-experiment BBD with 3 center replicates. Analyze each on the SPME Arrow DVB/CAR/PDMS (1.10 mm × 120 µm × 20 mm) fiber. Hold pre-incubation at 10 min, sample volume at 5 mL in 20 mL vial (giving 3:1 HS:sample), ethanol level at 10% (diluted from 13% Cab Sauv).

### 7.4 Confirmation runs

Run the model-predicted optimum in triplicate. Acceptable if predicted vs. observed peak area is within ±10% for ≥80% of analytes and within ±20% for all analytes. Then perform validation (linearity, LOD, LOQ, accuracy at 3 levels, intra/inter-day precision) under those locked conditions.

### 7.5 Robustness check (Phase C)

After the optimum is locked, run a **2^(k-p) fractional factorial** on tight ±5% perturbations of the optimum factors plus 2 nuisance factors (vial-to-vial variability, fiber age in injection count). This is the ICH-style robustness assessment — typically 8-16 runs.

---

## 8. References (DOIs + URLs)

Retrieved 2026-05-12.

1. **Saha, B., Longo, R., Torley, P., Saliba, A., & Schmidtke, L. (2018).** SPME Method Optimized by Box-Behnken Design for Impact Odorants in Reduced Alcohol Wines. *Foods* 7(8):127. doi:10.3390/foods7080127. https://pmc.ncbi.nlm.nih.gov/articles/PMC6112000/

2. **Šikuten, I., Štambuk, P., Karoglan Kontić, J., Maletić, E., Tomaz, I., & Preiner, D. (2021).** Optimization of SPME-Arrow-GC/MS method for determination of free and bound volatile organic compounds from grape skins. *Molecules* 26(23):7409. doi:10.3390/molecules26237409. https://pmc.ncbi.nlm.nih.gov/articles/PMC8659239/

3. **Welke, J. E. et al. (2024).** Optimisation of SPME Arrow GC/MS method for determination of wine volatile organic compounds. *OENO One* (article view 7914). https://oeno-one.eu/article/view/7914

4. **Rossi, L., Foschi, M., Biancolillo, A., et al. (2023).** Optimization of HS-SPME-GC/MS Analysis of Wine Volatiles Supported by Chemometrics for the Aroma Profiling of Trebbiano d'Abruzzo and Pecorino White Wines Produced in Abruzzo (Italy). *Molecules* 28(4):1534. doi:10.3390/molecules28041534. https://pmc.ncbi.nlm.nih.gov/articles/PMC9962864/

5. **Burin, V. M., Marchand, S., de Revel, G., & Bordignon-Luiz, M. T. (2013).** Development and validation of method for heterocyclic compounds in wine: optimization of HS-SPME conditions applying a response surface methodology. *Talanta* 117:87-93. doi:10.1016/j.talanta.2013.08.037. https://pubmed.ncbi.nlm.nih.gov/24209315/

6. **Arcanjo, N. M. O. et al. (2015).** Optimization of the HS-SPME-GC/MS technique for determining volatile compounds in red wines made from Isabel grapes (*Vitis labrusca*). *Food Sci. Technol. Campinas* 35(4). doi:10.1590/1678-457X.6815. https://www.scielo.br/j/cta/a/QrYYPhg6JhfTPF5KtCxLjxP/?lang=en

7. **Zhao, Y. et al. (2019).** Optimization and validation of a head space solid-phase microextraction-arrow gas chromatography-mass spectrometry method using central composite design for determination of aroma compounds in Chinese liquor (Baijiu). *J. Chromatogr. A* 1610:460537. PMID:31607446. https://pubmed.ncbi.nlm.nih.gov/31607446/

8. **Siebert, T. E., Smyth, H. E., Capone, D. L., Neuwöhner, C., Pardon, K. H., Skouroumounis, G. K., Herderich, M. J., Sefton, M. A., & Pollnitz, A. P. (2005).** Stable isotope dilution analysis of wine fermentation products by HS-SPME-GC-MS. *Anal. Bioanal. Chem.* 381:937-947. doi:10.1007/s00216-004-2992-4. https://link.springer.com/article/10.1007/s00216-004-2992-4

9. **Antalick, G., Perello, M.-C., & de Revel, G. (2014).** Esters in Wines: New Insight through the Establishment of a Database of French Wines. *Am. J. Enol. Vitic.* 65(3) early 2014. https://www.ajevonline.org/content/early/2014/05/13/ajev.2014.13133

10. **Pizarro, C., Pérez-del-Notario, N., & González-Sáiz, J. M. (2009).** Headspace solid-phase microextraction for direct determination of volatile phenols in cider [Doehlert experimental design example]. *J. Sep. Sci.* 32(18):3271-3279. doi:10.1002/jssc.200900347. https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/full/10.1002/jssc.200900347

11. **Bee-DiGregorio, M. Y., Feng, H., Pan, B. S., & Dokoozlian, N. K. (2009 / later).** Rapid measurement of 3-alkyl-2-methoxypyrazine content of winegrapes to predict levels in resultant wines. *J. Agric. Food Chem.* 57(20):9367-9374. doi:10.1021/jf9019695. https://pubs.acs.org/doi/10.1021/jf9019695

12. **Sereli, A., Kallithraka, S. et al. (2013).** Chemometrical development and comprehensive validation of a solid phase microextraction/gas chromatography–mass spectrometry methodology for the determination of important free and bound primary aromatics in Greek wines [Plackett-Burman + Simplex]. *J. Chromatogr. A* 1311:11-23. https://www.sciencedirect.com/science/article/abs/pii/S002196731301025X

13. **Câmara, J. S., Alves, M. A., & Marques, J. C. (2006).** Development of headspace solid-phase microextraction-gas chromatography-mass spectrometry methodology for analysis of terpenoids in Madeira wines. *J. Chromatogr. A* (canonical wine SPME reference). (frequently cited foundation work).

14. **Derringer, G., & Suich, R. (1980).** Simultaneous Optimization of Several Response Variables. *J. Quality Technology* 12(4):214-219. doi:10.1080/00224065.1980.11980968. https://www.tandfonline.com/doi/abs/10.1080/00224065.1980.11980968

15. **Stilo, F., Bicchi, C., Cordero, C. (2021).** Quantification of Volatile Compounds in Wines by HS-SPME-GC/MS: Critical Issues and Use of Multivariate Statistics in Method Optimization. *Processes* 9(4):662. doi:10.3390/pr9040662. https://www.mdpi.com/2227-9717/9/4/662 [authoritative review of DoE in wine SPME]

16. **Tufariello, M., Anaya, J. A., et al. (2018).** Development and optimization of a HS-SPME-GC-MS methodology to quantify volatile carbonyl compounds in Port wines [4-factor BBD example for Port]. *Food Chemistry*. doi: see ScienceDirect. https://www.sciencedirect.com/science/article/abs/pii/S030881461831238X

17. **Mateo-Vivaracho, L., Cacho, J., & Ferreira, V. (2006-2007).** Quantitative determination of wine highly volatile sulfur compounds by HS-SPME and GC-PFPD. *J. Chromatogr. A*. https://pubmed.ncbi.nlm.nih.gov/17207804/

18. **Slaghenaufi, D., Indorato, C., Troilo, M., et al. (2012).** Quantification by solid phase micro extraction and stable isotope dilution assay of norisoprenoid compounds in red wines obtained from Piedmont rare varieties. https://pubmed.ncbi.nlm.nih.gov/22980832/

19. **IOFI Working Group on Methods of Analysis. (2010).** Guidelines for solid-phase micro-extraction (SPME) of volatile flavour compounds for gas-chromatographic analysis. *Flavour Fragr. J.* 25:404-406. doi:10.1002/ffj.1991. https://onlinelibrary.wiley.com/doi/full/10.1002/ffj.1991

20. **NIST/SEMATECH e-Handbook of Statistical Methods**, §5.3.3 — Response Surface Designs (CCD, BBD comparisons, alpha selection, rotatability). https://www.itl.nist.gov/div898/handbook/pri/section3/pri3361.htm

21. **Helin, A. et al. (2015).** Comparison of SPME using classical fibers versus mini-Arrows applying multiple headspace extraction and various agitation techniques. *Chromatographia*. doi:10.1007/s10337-018-3659-1. https://link.springer.com/article/10.1007/s10337-018-3659-1

22. **Agilent Application 5994-3159EN.** Use of Salt to Increase Analyte Concentration in SPME — Smoke-Taint Wine Application Note. (PDF; access intermittent — referenced for salt-out protocol context.)

23. **PAL System Smart SPME Arrow product documentation.** DVB/Carbon-WR/PDMS 120 µm × 20 mm Arrow specification. https://www.palsystem.com/fileadmin/user_upload/content_hub/Files/Brochures/PAL_Smart_SPME_Arrow_Brochure_screen.pdf

---

*End of brief — Specialist 3/6.*
