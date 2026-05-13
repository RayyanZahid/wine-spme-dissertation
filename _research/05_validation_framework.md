# Method Validation — Framework Brief

**Specialist 5/6 · Wine SPME Arrow HS-GC-MS/MS Dissertation · Prepared 2026-05-12**

Scope: validation-chapter scaffolding for "Method Development, MRM Optimization, and Validation of SPME Arrow Headspace GC-MS/MS for the Quantitative Profiling of Aroma-Active Volatiles in Cabernet Sauvignon Wine". This brief consolidates the official frameworks, the formal calculation formulas, the consolidated acceptance criteria, and a recommended single-laboratory validation protocol suitable for a peer-reviewable dissertation.

---

## 1. Framework selection: ICH Q2(R2) as primary, OIV / SANTE / AOAC SMPR as supporting

The wine-volatile method is a **quantitative trace/contaminant-style analysis** carried out by a single laboratory with no pharmacopoeial registration target. It is also a **food/beverage matrix method targeted at GC-MS/MS**, so the most pragmatic numeric criteria come from EU residue-analysis guidance. The dissertation should therefore declare a layered framework rather than picking one document in isolation. Recommended hierarchy:

| Layer | Document | Role in the chapter |
|---|---|---|
| Primary framework (vocabulary + characteristic list) | **ICH Q2(R2) — Validation of Analytical Procedures** (adopted 31 Oct 2023; published 30 Nov 2023; FDA guidance Mar 2024) | Defines the eight validation characteristics and is the most widely accepted single-document framework for analytical-method validation. Use its definitions of specificity, accuracy, precision (repeatability / intermediate precision / reproducibility), DL, QL, linearity, range, robustness. |
| Supporting framework (single-lab statistics) | **IUPAC Technical Report — Thompson, Ellison & Wood 2002, *Pure Appl. Chem.* 74(5):835-855** | Provides the Horwitz / HorRat statistics and the recognised replicate-count recommendations for in-house validation. Cited explicitly by Codex CAC/GL 90-2017. |
| Supporting framework (food/wine register) | **OIV Compendium of International Methods of Wine and Must Analysis — Annex E, *Practical guide for the validation*** (current edition; OIV-OENO 596-2019 amends) | Adds wine-specific framing: matrix range, single-laboratory verification, certified reference materials for wine, recovery on real wine matrices. |
| Supporting framework (food register; Codex baseline) | **Codex Alimentarius CAC/GL 90-2017 — Guidelines on performance criteria for methods of analysis** | The Codex-recognised performance-criteria framework; aligns with the IUPAC harmonized guidelines and is what regulators reach for in food-control disputes. |
| Supporting framework (pragmatic numeric tolerances for GC-MS/MS) | **EU SANTE/11312/2021 — Analytical quality control and method validation procedures for pesticide residues** (in force 1 Jan 2022; updated v2 Nov 2023) | Even though this method is not a pesticide method, SANTE supplies the most pragmatic, defensible numeric thresholds for GC-MS/MS work: recovery 70-120%, RSD ≤ 20% per matrix/level, ion-ratio ±30%, signal at qualifier transition ≥ S/N 3, matrix-effect slope-ratio interpretation, etc. Cite "where transferable" with explicit footnote that the method is not pesticide. |
| Supporting framework (in-house / fitness-for-purpose) | **Eurachem Guide — *The Fitness for Purpose of Analytical Methods* (2nd edition, 2014; B. Magnusson & U. Örnemark, eds.; ISBN 978-91-87461-59-0; 3rd edition issued 2025)** | The accepted laboratory-level interpretation of ICH-style validation; pairs with Eurachem/CITAC QUAM CG-4 for uncertainty. |
| Supporting framework (AOAC method-performance template) | **AOAC International — Appendix F, Guidelines for SMPRs**, and per-analyte SMPRs as drafting templates | Used as the *template* for the SMPR-style table in §6. The dissertation does not require AOAC certification; the table merely demonstrates fitness against an internationally recognised template. |
| Uncertainty framework | **Eurachem/CITAC QUAM CG-4 (3rd ed., 2012)** plus **Nordtest NT TR 537 (4th ed., 2017)** | QUAM = bottom-up GUM modelling; Nordtest = top-down precision-and-bias. The dissertation should use the top-down approach (§8). |

**Single sentence to put in the chapter opening:** *"Validation followed ICH Q2(R2) (2023) as the primary framework, with single-laboratory statistical procedures from IUPAC (Thompson et al. 2002), wine-matrix conventions from OIV Annex E, pragmatic GC-MS/MS tolerances from EU SANTE/11312/2021 where transferable, and a top-down measurement-uncertainty estimate following Nordtest NT TR 537 and Eurachem/CITAC QUAM CG-4."*

---

## 2. Validation characteristics — definitions and formulas

ICH Q2(R2) groups validation under eight characteristics. The dissertation must address each. Below are the consolidated definitions, the calculation formulas, the experimental design, and acceptance ranges.

### 2.1 Specificity (selectivity)

**Definition.** Ability of the procedure to assess unequivocally the target analyte in the presence of components which may be expected to be present, such as matrix interferents, related compounds, and isobars. ICH Q2(R2) treats *specificity* and *selectivity* as effectively synonymous; Eurachem prefers *selectivity*.

**Experimental design.**
- Inject solvent blank, wine blank (a low-aroma matrix, e.g. a deodorised model wine or a stripped Cabernet), and matrix-matched standards.
- For each target analyte, demonstrate no interfering peak at the analyte retention time exceeding 30 % of the response at the LOQ (Eurachem; mirrors SANTE selectivity rule).
- For MS/MS confirmation: at least two MRM transitions per analyte; SANTE/11312/2021 ion-ratio rule: **observed ion ratio must agree with the calibration ion ratio within ±30 % (relative)** independent of the absolute ion-ratio value.
- For SPME selectivity: include co-extracted matrix volatiles known to co-elute (e.g. ethyl decanoate vs. isoamyl octanoate window in red wine).

**No formula — pass/fail against the interference + ion-ratio criteria.**

### 2.2 Linearity

**Definition.** Within a given range, the ability of the procedure to obtain results directly proportional to the concentration of analyte.

**Experimental design.**
- **Minimum 5 concentration levels** spanning the working range, each in duplicate or triplicate. ICH Q2(R2) text continues to recommend "minimum of 5 concentrations".
- Use **matrix-matched calibration** for wine to compensate for the headspace partitioning shift (see §4).
- Inspect the residual plot — random scatter required; trending residuals signal non-linearity.

**Primary metric: coefficient of determination R².** Acceptance: R² ≥ 0.99 (pharmaceutical norm) or R² ≥ 0.995 (Codex / OIV preferred for food-matrix work). The dissertation should pre-register R² ≥ 0.995.

**Y-intercept check.** Y-intercept should be statistically indistinguishable from zero, or its 95 % confidence interval should bracket zero; alternatively, |intercept| ≤ 2 % of the response at the working concentration.

**Lack-of-fit test (LOF).** Compares pure-error sum of squares (replicates) against the lack-of-fit sum of squares using an F-test. If F_LOF > F_crit at α = 0.05, linearity is rejected and a quadratic/weighted model must be evaluated.

**Mandel's fitting test.** A nested F-test comparing the residual variance of the linear model (DS₁²) with that of a quadratic model (DS₂²):
```
TV = ((N − 2)·DS₁² − (N − 3)·DS₂²) / DS₂²
```
where N = number of calibration points. Compare TV against F(1, N−3) at α = 0.05. If TV ≤ F_crit, the linear model is adequate; if TV > F_crit, switch to quadratic. Mandel's test is the IUPAC-endorsed linearity diagnostic; note however that it is only fully valid when variance is homogeneous and points are equidistant (Analytical Methods 2013).

**Residual analysis.** Plot ordinary residuals vs. concentration; superimpose a 95 % prediction band. The absolute relative residual for each calibration point should be ≤ 20 % of the nominal at every level except the LOQ (where ≤ 25 % is acceptable per SANTE).

### 2.3 Range

**Definition.** Interval between the upper and lower concentrations of analyte over which acceptable accuracy, precision, and linearity have been demonstrated. ICH Q2(R2) introduces an explicit split between **reportable range** (the full validated interval) and **working range** (the narrower interval routinely used).

**For an aroma-active volatile method** the range typically runs from the LOQ to ≥ 120 % of the highest concentration anticipated in real wines. Anchor the range to literature-observed maxima for each compound class:
- Esters: 0.05 → 5,000 µg L⁻¹
- C6 alcohols / higher alcohols: 50 → 500,000 µg L⁻¹ (a five-decade span — split the calibration over two ranges)
- Norisoprenoids (β-damascenone, β-ionone): 0.01 → 50 µg L⁻¹
- Methoxypyrazines (IBMP, IPMP): 0.001 → 0.1 µg L⁻¹ (sub-ng/L territory — SIDA preferred)
- Terpenes (linalool, geraniol): 1 → 1,000 µg L⁻¹
- Volatile phenols (4-EP, 4-EG, guaiacol): 1 → 5,000 µg L⁻¹

If any analyte spans more than two orders of magnitude, demonstrate linearity on two overlapping sub-ranges (Codex CAC/GL 90-2017 allows this with cross-validation at the overlap).

### 2.4 Accuracy (trueness)

**Definition.** Closeness of agreement between the value found and the value accepted as the conventional true value. Reported as **% recovery**.

**Three permissible references** (any one defensible):
1. **Certified Reference Material (CRM)** — preferred. For wine volatiles a wine-matrix CRM is rarely available, so substitute spike-recovery on a matrix-matched basis.
2. **Reference method comparison** — analyse the same sample by an orthogonal method (e.g. stir-bar sorptive extraction or solvent-assisted flavour evaporation GC-MS) and compare.
3. **Spike recovery** — standard wine matrix (red wine with low aroma background, freshly stripped to remove residual analyte) is spiked at three levels (typically LOQ, mid-range, and 80 % of upper range).

**Formula.**
```
% Recovery = ((C_spiked − C_unspiked) / C_added) × 100
```
For analytes endogenous to wine (most aroma compounds are!), C_unspiked is determined first on the unspiked aliquot and subtracted.

**Acceptance.** Concentration-dependent (Codex / AOAC / SANTE all align):
| Concentration | Recovery acceptance |
|---|---|
| > 10 µg g⁻¹ (or > 10 mg L⁻¹) | 90 – 107 % |
| 1 – 10 µg g⁻¹ | 80 – 110 % |
| 0.1 – 1 µg g⁻¹ | 80 – 110 % |
| 0.01 – 0.1 µg g⁻¹ | 70 – 120 % |
| 1 – 10 ng g⁻¹ | 60 – 115 % |
| < 1 ng g⁻¹ (ultra-trace, e.g. methoxypyrazines) | 50 – 120 % |

SANTE/11312/2021 collapses this to a single rule for pesticide work: **mean recovery 70 – 120 % at each validated level**. The dissertation should adopt the concentration-dependent table for normative defensibility and additionally show that every analyte clears the 70-120 % SANTE rule.

**Trueness via isotope dilution.** For analytes with a commercially available deuterated or ¹³C-labelled internal standard (e.g. linalool-d3, β-damascenone-d4, 2-isobutyl-3-methoxypyrazine-d3), use **stable-isotope-dilution analysis (SIDA)**. SIDA gives accuracy that absorbs both extraction-recovery and ionisation-suppression effects in a single ratio measurement.

### 2.5 Precision

ICH Q2(R2) splits precision into three nested conditions:

| Component | Symbol | Conditions varied | Conditions held |
|---|---|---|---|
| Repeatability | RSD_r | none — replicate injections within run | same operator, instrument, day, calibration |
| Intermediate precision | RSD_iR (a.k.a. RSD_wR within-laboratory reproducibility) | day, operator, instrument, fibre lot | laboratory |
| Reproducibility | RSD_R | laboratory (inter-laboratory) | nominal method |

**For a single-laboratory dissertation, only RSD_r and RSD_iR are obtained**; RSD_R is acknowledged in the limitations and discussed via the Horwitz prediction (§5).

**Formulas.**
```
RSD_r  =  (s_r  / x̄) × 100 %
RSD_iR =  (s_iR / x̄) × 100 %
```
where s_r is the pooled within-run SD across replicates (one-way ANOVA on n replicates × 1 day), and s_iR is the between-run SD computed from a one-way ANOVA design with day (or operator) as the factor (ISO 5725-2).

Minimum design: **6 replicates × 3 levels × 6 days** (giving 108 measurements per analyte). This conforms to ICH Q2(R2) ("repeatability: 6 × 100 % or 3 × 3 levels", "intermediate precision: ≥ 3 replicates × 6 runs over ≥ 2 days"). The dissertation should over-spec to 6 × 6 × 3 for robustness.

**Acceptance.** Concentration-dependent (Codex / AOAC / Horwitz; SANTE collapses to a single 20 % rule):
| Concentration | RSD_r acceptance | RSD_iR acceptance |
|---|---|---|
| > 10 mg L⁻¹ | ≤ 2 % | ≤ 4 % |
| 1 – 10 mg L⁻¹ | ≤ 3 % | ≤ 6 % |
| 0.1 – 1 mg L⁻¹ | ≤ 4 % | ≤ 8 % |
| 10 – 100 µg L⁻¹ | ≤ 6 % | ≤ 11 % |
| 1 – 10 µg L⁻¹ | ≤ 8 % | ≤ 16 % |
| 0.1 – 1 µg L⁻¹ | ≤ 11 % | ≤ 22 % |
| < 0.1 µg L⁻¹ | ≤ 15 % | ≤ 30 % |

SANTE single rule: RSD ≤ 20 % at every validated level.

### 2.6 Detection limit (LOD) and Quantitation limit (LOQ)

Three formal definitions — pick one as primary, present the others in a comparison (§3).

**A. Signal-to-noise (S/N) method.** Inject low-concentration standards and measure peak-to-peak noise on the matrix-blank chromatogram at the analyte retention window.
```
LOD : S/N = 3   (sometimes stated as 3:1)
LOQ : S/N = 10  (sometimes stated as 10:1)
```
For MS/MS with a clean quantifier transition, S/N from the qualifier transition is often the constraint.

**B. ICH calibration-curve method.**
```
LOD = 3.3 × σ / S
LOQ = 10  × σ / S
```
where σ = SD of the response (estimated from either the SD of the y-intercepts of multiple calibration curves, or the residual standard deviation of the regression, or the SD of blank responses) and S = slope of the calibration curve. This is the ICH Q2(R2) recommended statistical approach.

**C. Replicate-blank method.** Inject ≥ 10 blank-matrix samples (Eurachem requires ≥ 10).
```
LOD = mean_blank + 3 × SD_blank   (= 3σ above blank)
LOQ = mean_blank + 10 × SD_blank  (= 10σ above blank)
```
**Confirmation requirement (Eurachem + OIV).** Calculated LOD/LOQ must be confirmed empirically: inject ≥ 6 replicates at the calculated LOQ and demonstrate RSD_r ≤ 10 % (or ≤ table value above) **and** recovery 80 – 120 %. The dissertation must include the confirmation experiment.

### 2.7 Robustness (ruggedness)

**Definition.** Measure of the procedure's capacity to remain unaffected by small but deliberate variations in method parameters. ICH Q2(R2) expanded the definition to include sample and reagent stability under perturbed conditions.

**Design: Youden-Steiner 7-factor / 8-run fractional factorial.** Seven nominal method parameters are each set to a "high" (+) and "low" (−) level; the 8-run pattern of ±1 values gives a saturated Plackett-Burman (Hadamard) design. Effect of factor i:
```
Effect_i = (sum of responses at +) / 4 − (sum of responses at −) / 4
```
Critical effect (significance) test:
```
|Effect_i| > t · √(2·s²/N)
```
typically using s² estimated from method-repeatability data and a Student t at α = 0.05.

**Suggested seven SPME-Arrow / GC-MS/MS factors and ±5–10 % variations:**
1. Extraction temperature (e.g. 45 / 55 °C)
2. Extraction time (e.g. 30 / 40 min)
3. Salt addition (NaCl) (e.g. 1.5 / 2.5 g per 10 mL)
4. Agitation speed (e.g. 250 / 500 rpm)
5. Desorption time (e.g. 3 / 7 min)
6. Inlet temperature (e.g. 240 / 260 °C)
7. Sample volume / vial fill ratio (e.g. 5 / 10 mL in a 20 mL vial)

The Youden design completes robustness in 8 runs and is standard for GC-MS/MS method validation.

### 2.8 Carryover

Not in ICH Q2(R2) as a standalone characteristic, but mandatory in **ICH M10 (Bioanalytical Method Validation, 2022)**, **EMA Bioanalytical guideline**, and adopted as best practice for chromatography in general.

**Test.** Inject a blank immediately after the highest calibration standard (or a "high QC"). Repeat 3 times.

**Acceptance.** The response at the analyte retention time in the post-high blank should be **< 20 % of the LOQ response** (ICH M10; conservative best practice). Internal-standard carryover should be **< 5 %** of the IS response in the analytical run.

If carryover exceeds either threshold, document a remediation strategy: solvent-only blank rinses, fibre bake-out between injections, longer GC oven post-run, alternate sample order randomisation.

---

## 3. LOD / LOQ — three methods compared

| Aspect | A. Signal-to-noise | B. ICH calibration method | C. Replicate-blank |
|---|---|---|---|
| Formula (LOD) | S/N = 3 | LOD = 3.3 × σ / S | LOD = mean_blank + 3 × SD_blank |
| Formula (LOQ) | S/N = 10 | LOQ = 10 × σ / S | LOQ = mean_blank + 10 × SD_blank |
| Input | Single low-conc. injection with noise estimate | Multi-level calibration; σ = SD of y-intercept or residual SD | ≥ 10 independent blank injections |
| Requires blank signal? | No (noise only) | No (uses calibration SDs) | Yes — assumes a measurable blank mean |
| Best for MS/MS? | Common in chromatography; software-dependent (peak-to-peak vs. RMS) | Independent of vendor noise definition | Only if a non-zero blank response exists |
| Best for GC-MS/MS with quiet baseline? | Sometimes gives unrealistically low values because noise floor is electronic, not chemical | Most defensible; recommended | Less informative — many wine volatiles have no measurable blank |
| Caveat | Noise definition (peak-to-peak vs. SD; vendor-software differences) | σ must be from samples close to the expected LOD, not from full range | Requires representative blanks; some analytes have no true blank |
| **Recommended for this dissertation** | Reported as cross-check | **Primary** | Reported for any analyte with a non-zero blank (e.g. higher alcohols always present in wine) |

**Empirical confirmation step is mandatory regardless of the method used:** spike at the calculated LOQ, inject ≥ 6 replicates, demonstrate RSD_r ≤ table value and 80 – 120 % recovery.

---

## 4. Matrix effect — slope-ratio test (worked example)

### Why matrix effect matters for SPME of wine

SPME headspace response is governed by the **liquid–headspace partition coefficient K_lh** and the **headspace–fibre partition coefficient K_hf**. Ethanol (13 – 15 % v/v in red wine), glycerol, polyphenols, polysaccharides, and tannins all alter K_lh by changing the activity coefficients of the volatiles in the wine matrix. Therefore the response per unit concentration in pure aqueous standard is *almost never equal* to the response per unit concentration in wine. SPME-GC-MS quantification *without* matrix correction overestimates or underestimates routinely by 20 – 80 %. (LC-MS/MS labels this "signal suppression / enhancement (SSE)"; GC-MS literature uses simply "matrix effect (ME)".)

### Slope-ratio test

Construct two calibration curves on the same analyte set:
1. **Solvent curve** — analyte spiked into pure aqueous matrix (or 12 % v/v ethanol-water "model wine") at the same series of concentrations.
2. **Matrix-matched curve** — analyte spiked post-extraction (or, for SPME, into a wine matrix sample previously stripped of native analyte) at the same concentrations.

The matrix-effect ratio:
```
ME (%) = (Slope_matrix / Slope_solvent) × 100
```
Some references express ME as a deviation:
```
ME (%) = ((Slope_matrix − Slope_solvent) / Slope_solvent) × 100
```
The two conventions differ in zero point (100 % vs. 0 %). Pick one and state it in the chapter; the *ratio* convention is more common in food-residue literature.

### Classification thresholds (ratio convention)

| ME (%) | Classification | Action |
|---|---|---|
| 80 – 120 (i.e. ≤ ±20 % from 100) | **Negligible / soft** | Solvent calibration acceptable; no correction needed |
| 50 – 80 or 120 – 150 (±20 – 50 %) | **Medium** | Matrix-matched calibration required *or* SIDA |
| < 50 or > 150 (> ±50 %) | **Strong** | Matrix-matched calibration mandatory; standard-addition or SIDA strongly recommended |

SANTE/11312/2021 does not enforce numeric ME limits; instead it requires that the calibration approach **accounts for** observed matrix effects — i.e. matrix-matched calibration is *de facto* the default for residue work. For wine-volatile dissertations the same principle holds.

### Worked example (illustrative numbers)

Analyte: β-damascenone. Calibration over 0.5 – 50 µg L⁻¹ in 4 levels × triplicate.
- Solvent curve (model wine, 12 % ethanol-water, no other matrix): slope = 28 450 ± 510 area units / (µg L⁻¹), R² = 0.9994
- Matrix-matched curve (stripped Cabernet Sauvignon spiked): slope = 19 320 ± 380 area units / (µg L⁻¹), R² = 0.9988

```
ME = (19 320 / 28 450) × 100 = 67.9 %
```
Classification: **medium suppression (~32 % suppression)**. Conclusion: matrix-matched calibration is mandatory for β-damascenone quantification; report this in §3 of the validation chapter. If a deuterated internal standard (β-damascenone-d4) is available, SIDA absorbs the matrix effect and the solvent calibration becomes acceptable provided the IS recovery is checked.

### Internal-standard considerations

A non-isotopologue internal standard (e.g. 3-octanol used as IS for many wine-volatile methods) only partially corrects for matrix effects; it corrects for *extraction-step* drift but not for *headspace partitioning* of structurally different analytes. SIDA is the gold standard. Where SIDA is impractical (cost; no commercially available label), use **standard-addition quantitation** on each sample — the slope-comparison check still applies at the method-development stage.

---

## 5. Precision — RSD_r / RSD_iR / HorRat (worked example)

### Horwitz equation

For inter-laboratory studies Horwitz observed an empirical concentration-dependence of reproducibility:
```
RSD_R(predicted) = 2 ^ (1 − 0.5 · log10(c))   %
```
where c is the analyte concentration expressed as a **mass fraction** (i.e. dimensionless: 1 µg g⁻¹ ⇒ c = 10⁻⁶; 1 µg L⁻¹ in water ⇒ c ≈ 10⁻⁹).

An algebraically equivalent form:
```
PRSD_R = 2 · c^(−0.1505)
```
(IUPAC writes the exponent as 0.15.)

### HorRat (Horwitz ratio)

```
HorRat   = RSD_R(observed)   / RSD_R(predicted)        — for inter-laboratory data
HorRat_r = RSD_r(observed)   / RSD_R(predicted)        — for single-laboratory data (within-lab)
```

### Acceptance ranges

| Quantity | Acceptable range |
|---|---|
| HorRat (inter-laboratory) | **0.5 – 2.0** (AOAC; below 0.5 may indicate variance underestimation; above 2.0 indicates poor reproducibility) |
| HorRat_r (single-laboratory repeatability) | **0.3 – 1.3** (AOAC empirical) |

Horwitz applies most reliably for analyte mass fractions between 10⁻⁹ and 10⁻¹; outside that range the equation has known limitations (Thompson 2007, Anal. Bioanal. Chem.).

### Worked example

Linalool determined by HS-SPME-GC-MS/MS in red wine. Concentration found: c = 250 µg L⁻¹ in wine ≈ 250 × 10⁻⁹ mass fraction = 2.5 × 10⁻⁷.

```
log10(c) = log10(2.5 × 10⁻⁷) = −6.60
PRSD_R = 2 ^ (1 − 0.5 × (−6.60))
       = 2 ^ (1 + 3.30)
       = 2 ^ 4.30
       = 19.7 %
```

Suppose the observed in-house intermediate precision is RSD_iR = 7.5 %.

```
HorRat_r = 7.5 / 19.7 = 0.38
```

Conclusion: 0.38 sits at the bottom of the 0.3 – 1.3 acceptance window for single-laboratory work — the method is precise (perhaps borderline "too tight", which can indicate either an exceptional method or variance under-estimation from running too few independent replicates). Report the value with the design context (e.g. 6 replicates × 3 levels × 6 days = 108 measurements) so reviewers can judge the variance estimate's quality.

For β-damascenone at 5 µg L⁻¹: c = 5 × 10⁻⁹, log10(c) = −8.30, PRSD_R = 2 ^ (1 + 4.15) = 2 ^ 5.15 = 35.5 %. Observed RSD_iR of 12 % gives HorRat_r = 0.34 — also within the acceptance band.

The HorRat / HorRat_r table for the full analyte panel should appear as a column in the consolidated validation table (§6).

---

## 6. Acceptance criteria — consolidated table (ICH + IUPAC + OIV + SANTE + AOAC)

This is the master per-analyte performance table that the validation chapter should produce. Apply each criterion at each spiking level.

| Validation characteristic | ICH Q2(R2) | IUPAC 2002 / Codex CAC/GL 90 | OIV Annex E | SANTE/11312/2021 | AOAC SMPR template |
|---|---|---|---|---|---|
| Specificity / selectivity | No interference at analyte RT; orthogonal-method or DAD/MS confirmation | Selectivity check vs. potential interferents | Standard-addition selectivity test | Two MRM transitions; ion ratio within ±30 % rel. of calibration | Methods must distinguish target unambiguously |
| Linearity | ≥ 5 levels; visual + statistical assessment; R² with residual inspection | ≥ 6 levels recommended; explicit lack-of-fit test | ISO 11095 / ISO 8466 regression | ≥ 5 levels; deviation of back-calc. concentration ≤ ±20 % per point (≤ ±25 % at LOQ) | R² ≥ 0.99 (template) |
| Range | Stated; reportable vs. working split | 80 – 120 % of target; food-trace methods often LOQ → 500 % of max realistic | Stated for wine concentration range; matrix-bracket | LOQ → ≥ MRL; bracketed by calibration | Stated |
| Accuracy / trueness | Spike-recovery or CRM; 95 % CI within acceptance | Concentration-dependent recovery table (see §2.4) | CRM or spike-recovery; wine-matrix preferred | **Mean recovery 70 – 120 %** at each level | Method-specific; SMPR specifies low/high |
| Repeatability RSD_r | ≤ 2 % (assay), ≤ 5 % (impurities); concentration-dependent in practice | HorRat_r 0.3 – 1.3 | Concentration-dependent (matches Codex table) | **RSD ≤ 20 %** at each level | SMPR-specific |
| Intermediate precision RSD_iR | ≤ 2 – 5 % typical | HorRat 0.5 – 2.0 (with intermediate as reproducibility proxy) | Concentration-dependent | **RSD ≤ 20 %** typically | SMPR-specific |
| Reproducibility RSD_R | Inter-laboratory; not single-lab | Horwitz: predicted RSD_R = 2·c^(−0.15) | If applicable | Inter-lab only | If multi-lab method |
| Detection limit (LOD) | 3.3·σ/S or S/N = 3 | Defined; ≥ 10 blank replicates for σ_blank approach | Generally 1/3 of LOQ; ≥ 10 replicates | At or below screening detection limit | Method-specific |
| Quantitation limit (LOQ) | 10·σ/S or S/N = 10; empirical confirmation required | Lowest level meeting accuracy + precision criteria | Lowest validated quantifiable level (≥ 10 replicates) | At or below MRL; empirical confirmation: RSD ≤ 20 %, recovery 70 – 120 % | Method-specific |
| Robustness | Q2(R2) explicit; Youden / DoE design encouraged | "Ruggedness" — single-factor or fractional factorial | Single-factor variations | Implied; not numerical | SMPR-specific |
| Carryover | Not a Q2(R2) characteristic | Not in 2002 IUPAC | Not in Annex E | Not in SANTE | ICH M10: < 20 % LOQ in blank after high standard |
| Ion-ratio / MS identification | Not in Q2(R2) (small-molecule generic) | Not in Codex 90 | OIV uses MS confirmation as accepted | **Ion ratio within ±30 % rel.**; ≥ 2 transitions per analyte; identification points scheme | If applicable |
| Matrix effect | Not in Q2(R2) explicitly | Recognised; correction expected | Matrix-correction expected for wine | Matrix-matched calibration or equivalent compensation; no numeric threshold | Method-specific |
| Measurement uncertainty | Acknowledged; not detailed | "Uncertainty budget" recommended | "Estimation of uncertainty" required | Reproducibility-based top-down; default expanded U ≤ 50 % at MRL | Method-specific |

**Recommended numeric thresholds to publish in the dissertation per analyte:**
- Linearity: R² ≥ 0.995; back-calculated residual ≤ ±15 % per calibration point (≤ ±20 % at LOQ).
- Accuracy: 70 – 120 % recovery at each spiking level (LOQ, mid, high).
- Repeatability RSD_r ≤ 15 %; intermediate precision RSD_iR ≤ 20 %; HorRat_r within 0.3 – 1.3.
- LOD: 3.3·σ/S; LOQ: 10·σ/S; LOQ confirmed empirically (RSD ≤ 20 %, recovery 70 – 120 %).
- Carryover at high-standard blank: < 20 % of LOQ response.
- Ion ratio between quan and qual transitions: within ±30 % (relative) of the calibration mean.
- Matrix effect: report value; matrix-matched calibration if |ME − 100 %| > 20 %.

---

## 7. Validation protocol recommended for this dissertation

The following experimental plan satisfies all five frameworks for single-laboratory validation of an HS-SPME-Arrow GC-MS/MS method targeting aroma-active volatiles in Cabernet Sauvignon. Adjust analyte count as the final MRM panel is fixed.

### 7.1 Matrices

- **Stripped wine matrix** — Cabernet Sauvignon flash-distilled or solid-phase-stripped (LiChrolut EN cartridge) to remove the native volatiles, then re-spiked. Confirm < 5 % residual signal for each target before use.
- **Real wine matrix** — three commercial Cabernet Sauvignon wines from different vintages/regions for matrix-representativeness checks.
- **Solvent / model wine** — 12 % (v/v) aqueous ethanol with 5 g L⁻¹ tartaric acid, pH 3.5. Use for the solvent leg of the slope-ratio test.

### 7.2 Calibration

- 7 levels spanning each analyte's working range; matrix-matched in stripped Cabernet.
- 3 replicates per level → 21 calibration points per analyte per curve.
- Run a fresh calibration each validation day (6 days total) → enables RSD on slope, intercept, and inter-day stability.
- For analytes with deuterated IS available: spike a constant amount of d-IS into every standard and every sample; quantify using analyte/IS area ratio.

### 7.3 Specificity / selectivity

- 3 × solvent blank
- 3 × stripped-wine matrix blank
- 3 × matrix-matched LOQ standard
- 3 × spiked real wine at mid-range
- Confirm no peak > 30 % of LOQ at analyte RT in any blank; confirm ion ratio within ±30 % of calibration in spiked samples.

### 7.4 Linearity

- Use the 21-point calibration set.
- Compute R², plot residuals, run lack-of-fit F-test, run Mandel's quadratic-vs-linear test, check y-intercept 95 % CI.

### 7.5 Accuracy + intra/inter-day precision (combined design)

The most efficient design uses a **6 × 6 × 3** factorial:
- 3 spiking levels (low ≈ 2 × LOQ; mid ≈ middle of working range; high ≈ 80 % of upper range)
- 6 replicates per level per day
- 6 separate analytical days
- Total: 108 measurements per analyte

For each measurement, store the raw area-ratio so a one-way ANOVA (day as factor) yields s_r (within-day) and s_iR (between-day) for each level.

Recovery is computed per-replicate; report mean recovery and 95 % CI per level. Acceptance: 70 – 120 %.

### 7.6 LOD / LOQ

- **Primary method (ICH calibration):** LOD = 3.3·σ_residual / S; LOQ = 10·σ_residual / S, where σ_residual is the residual SD of the matrix-matched calibration over the lowest 4 levels.
- **Secondary (S/N):** measure peak-to-peak noise on the quantifier transition in the blank wine matrix around analyte RT; report LOD where S/N = 3, LOQ where S/N = 10.
- **Tertiary (replicate-blank, for analytes with non-zero blank only):** ≥ 10 blank-matrix injections; LOD = mean + 3·SD, LOQ = mean + 10·SD.
- **Empirical confirmation:** 6 replicates at the calculated LOQ; require RSD_r ≤ 20 % AND recovery 70 – 120 %.

### 7.7 Carryover

- After the highest calibration standard, inject 3 × matrix blank consecutively. Compute response at analyte RT relative to the LOQ response. Acceptance: < 20 % of LOQ for analyte, < 5 % for IS.

### 7.8 Robustness

- Youden-Steiner 7-factor / 8-run fractional factorial (the seven factors from §2.7).
- Replicate at the mid-spike level; analyse the same matrix-matched standard at each of the 8 condition combinations.
- Compute factor effects vs. critical effect at α = 0.05; identify and flag factors with significant influence.

### 7.9 Matrix effect

- Run paired calibration curves (solvent vs. matrix-matched) for at least three representative analytes per chemical class (ester, alcohol, norisoprenoid, terpene, methoxypyrazine, volatile phenol).
- Compute ME (%) = (Slope_matrix / Slope_solvent) × 100; classify as soft / medium / strong; document mitigation (matrix-matched cal, SIDA, or standard addition) per analyte.

### 7.10 Measurement uncertainty (§8)

- Top-down via Nordtest combining s_iR and bias-uncertainty.

### 7.11 Validation duration estimate

| Step | Days |
|---|---|
| Stripping/QC of matrix | 2 |
| Calibration curves (matrix-matched + solvent, 6 days) | 6 |
| Accuracy + intra/inter-day precision (6 × 6 × 3 across 6 days; can be co-run with calibration) | 6 (shared) |
| LOD/LOQ confirmation | 1 |
| Carryover | 1 (shared with calibration) |
| Robustness (8 runs) | 2 |
| Matrix effect (per-class duplicate curves) | 2 |
| Data analysis, regression, ANOVA, write-up | 5 |
| **Total** | **~18 days of bench work** (excluding write-up) |

---

## 8. Measurement uncertainty (top-down approach)

### Why top-down for this method

Bottom-up GUM modelling (Eurachem QUAM CG-4) requires individual quantification of every uncertainty source (e.g. volumetric pipette, fibre carryover, GC injection volume, MS detector drift). For SPME headspace this is impractical because the dominant uncertainty source is the partitioning equilibrium itself — not separable into transparent components. Top-down (Nordtest NT TR 537) instead estimates uncertainty from validation data already in hand (intermediate precision + bias estimates), so it is the recommended approach.

### Nordtest formulation

The combined standard uncertainty is decomposed into two components:
```
u_c²(y) = u_Rw² + u_bias²
```
- **u_Rw = s_iR (intermediate precision standard deviation)** — directly from §7.5.
- **u_bias** captures laboratory bias estimated from CRMs, proficiency tests, or recovery experiments:
```
u_bias = √( bias_mean² + (s_bias / √n)² + u_ref² )
```
where bias_mean is the mean percentage bias from recovery (e.g. mean recovery 92 % ⇒ bias = −8 %), s_bias is the SD of recovery results across runs, n is the number of bias measurements, and u_ref is the uncertainty of the reference value used (CRM certificate uncertainty, or for spike-recovery, the uncertainty of the spiking standard).

### Expanded uncertainty

```
U = k · u_c       with k = 2 corresponding to ~95 % coverage
```
Report U as either an absolute value (same units as result) or as a relative U% (preferred for trace work).

### Worked example

For β-damascenone at 5 µg L⁻¹ with:
- s_iR = 12 % (intermediate precision, §5 example)
- Mean recovery 96 %, SD of recovery 8 %, n = 18 (6 reps × 3 levels)
- Spike standard certified to ±2 % (CertRef uncertainty u_ref = 1 %)
```
u_Rw = 12 %
bias_mean = (96 − 100) = −4 %; |bias| = 4 %
s_bias / √n = 8 / √18 = 1.89 %
u_bias = √(4² + 1.89² + 1²) = √(16 + 3.57 + 1) = √20.57 = 4.54 %
u_c = √(12² + 4.54²) = √(144 + 20.6) = √164.6 = 12.83 %
U (k = 2) = 25.7 %
```
Therefore report β-damascenone results as `c ± 0.257 × c` at the 95 % coverage level — e.g. `4.8 ± 1.2 µg L⁻¹`.

Repeat the calculation per analyte; the uncertainty budget per analyte becomes a column in the consolidated validation table.

### Cross-check against Horwitz

For an aroma compound at c = 5 × 10⁻⁹ mass fraction, PRSD_R = 35.5 % → expanded U (k = 2) ≈ 71 %. The measured U of 25.7 % is well below the Horwitz upper bound, confirming method fitness-for-purpose.

---

## 9. References (URLs + retrieval dates + DOIs)

All URLs verified 2026-05-12. Where a stable DOI exists it is given.

### Primary frameworks

1. **ICH Q2(R2) — Validation of Analytical Procedures** (adopted 31 Oct 2023; Step 4 / Step 5; Step-5 publication 30 Nov 2023). ICH, Geneva.
   - PDF: https://database.ich.org/sites/default/files/ICH_Q2(R2)_Guideline_2023_1130.pdf
   - EMA scientific guideline page: https://www.ema.europa.eu/en/ich-q2r2-validation-analytical-procedures-scientific-guideline
   - FDA finalised guidance (Mar 2024): https://www.fda.gov/regulatory-information/search-fda-guidance-documents/q2r2-validation-analytical-procedures
   - Federal Register notice: https://www.federalregister.gov/documents/2024/03/07/2024-04834

2. **Thompson, M.; Ellison, S. L. R.; Wood, R.** "Harmonized guidelines for single-laboratory validation of methods of analysis (IUPAC Technical Report)." *Pure and Applied Chemistry* 2002, **74**(5), 835-855.
   - DOI: 10.1351/pac200274050835
   - Open PDF: http://publications.iupac.org/pac/2002/pdf/7405x0835.pdf
   - Adopted by reference for Codex Alimentarius purposes.

3. **OIV — Compendium of International Methods of Wine and Must Analysis, Annex E (Laboratory Quality Assurance), Practical Guide for the Validation, Quality Control and Uncertainty Estimate of an Alternative Oenological Analysis Method** (current ed.; amended by OIV-OENO 596-2019).
   - Annex E entry: https://www.oiv.int/standards/compendium-of-international-methods-of-wine-and-must-analysis/annex-e/annex-e-laboratory-quality-assurance/practical-guide-for-the-validation
   - Compendium index: https://www.oiv.int/standards/compendium-of-international-methods-of-wine-and-must-analysis
   - OIV-OENO 596-2019: https://www.oiv.int/node/3128/download/pdf

4. **Codex Alimentarius — CAC/GL 90-2017, Guidelines on Performance Criteria for Methods of Analysis for the Determination of Pesticide Residues in Food and Feed** (FAO/WHO, 2017).
   - PDF: https://www.fao.org/fao-who-codexalimentarius/sh-proxy/jp/?lnk=1&url=https://workspace.fao.org/sites/codex/Standards/CXG+90-2017/CXG_090e.pdf
   - Codex guidelines index: https://www.fao.org/fao-who-codexalimentarius/codex-texts/guidelines/en/

5. **EU SANTE/11312/2021 — Analytical Quality Control and Method Validation Procedures for Pesticide Residues Analysis in Food and Feed** (in force 1 Jan 2022; second version updated 22 Nov 2023).
   - PDF v.2 (Nov 2023): https://food.ec.europa.eu/system/files/2023-11/pesticides_mrl_guidelines_wrkdoc_2021-11312.pdf
   - EURL-Pesticides hosting + changes log: https://www.eurl-pesticides.eu/userfiles/file/EurlALL/SANTE_11312_2021.pdf
   - Main changes document: https://food.ec.europa.eu/system/files/2022-02/pesticides_mrl_guidelines_wrkdoc_2021-11312.pdf
   - Accredia summary: https://www.accredia.it/en/documents/guidance-sante-11312-2021-analytical-quality-control-and-method-validation-procedures-for-pesticide-residues-analysis-in-food-and-feed/

6. **Eurachem Guide — *The Fitness for Purpose of Analytical Methods: A Laboratory Guide to Method Validation and Related Topics*, 2nd edition (2014).** B. Magnusson & U. Örnemark, eds. ISBN 978-91-87461-59-0. (3rd edition issued 2025.)
   - 2014 PDF: https://www.eurachem.org/images/stories/Guides/pdf/MV_guide_2nd_ed_EN.pdf
   - Eurachem publications page: https://www.eurachem.org/index.php/publications/guides
   - 2025 3rd edition page: https://www.eurachem.org/index.php/3-publications/guides/144-gdmv2014

### Acceptance criteria / formulas

7. **AOAC International — Appendix F, Guidelines for Standard Method Performance Requirements (SMPRs).**
   - PDF: https://www.aoac.org/wp-content/uploads/2019/08/app_f.pdf

8. **Horwitz, W.; Albert, R.** "The Horwitz ratio (HorRat): A useful index of method performance with respect to precision." *J. AOAC International* 2006, **89**(4), 1095-1109.
   - PubMed: https://pubmed.ncbi.nlm.nih.gov/16915851/
   - PDF: https://www.feedhaccp.org/distance/elearning/LABQUALITY/readings_2019/horwitz.pdf
   - Oxford Academic: https://academic.oup.com/jaoac/article-abstract/89/4/1095/5657708
   - LCGC review: https://www.chromatographyonline.com/view/benchmarking-analytical-methods-horwitz-curve
   - Contemporary update (PubMed 40106714): https://pubmed.ncbi.nlm.nih.gov/40106714/

9. **Eurachem/CITAC Guide CG 4 — *Quantifying Uncertainty in Analytical Measurement (QUAM)*, 3rd edition (2012).**
   - PDF: https://www.eurachem.org/images/stories/Guides/pdf/QUAM2012_P1.pdf
   - Guides page: https://www.eurachem.org/index.php/publications/guides/quam

10. **Nordtest Report NT TR 537 — *Handbook for Calculation of Measurement Uncertainty in Environmental Laboratories*, 4th edition (2017).**
    - PDF: https://www.nordtest.info/wp/2017/11/29/handbook-for-calculation-of-measurement-uncertainty-in-environmental-laboratories-nt-tr-537-edition-4/
    - Hosted copy: http://kemianseurat.fi/finntesting/wp-content/uploads/2019/10/Handbook-for-Calculation-of-Measurement-Uncertainty-in-Environmental-Laboratories-2017-Nordtest-TR-537-%E2%80%93raportti.pdf

11. **Andrade-Eiroa, A.; Castela, M. C.; Areal, P.; Lourenço, A.; Diniz, M. S.; Anastácio, A.** "Notes on the use of Mandel's test to check for nonlinearity in laboratory calibrations." *Anal. Methods* 2013, **5**, 1145.
    - DOI: 10.1039/c2ay26400e
    - Article: https://pubs.rsc.org/en/content/articlelanding/2013/ay/c2ay26400e
    - Blog discussion: https://blogs.rsc.org/ay/2013/01/24/mandel%E2%80%99s-test-a-case-of-oversimplification/

12. **Youden, W. J.; Steiner, E. H.** *Statistical Manual of the AOAC* (AOAC, 1975) — 7-factor / 8-run fractional factorial test.
    - Modern review: Karageorgou, E.; Samanidou, V. "Youden test application in robustness assays during method validation." *J. Chromatogr. A* 2014, **1353**, 131-139. DOI: 10.1016/j.chroma.2014.01.050. PubMed: https://pubmed.ncbi.nlm.nih.gov/24508395/

### Matrix-effect treatment + wine application

13. **Stahl-Zeng, J. et al.** "Matrix Effect Evaluation in GC/MS-MS Analysis of Multiple Pesticide Residues in Selected Food Matrices." *Foods* 2023, **12**, 3991.
    - DOI: 10.3390/foods12213991
    - HTML: https://www.mdpi.com/2304-8158/12/21/3991
    - PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC10650748/

14. **Bertarini, P.; Magalhães, V.; Silva, L.; et al.** "Quantification of Volatile Compounds in Wines by HS-SPME-GC/MS: Critical Issues and Use of Multivariate Statistics in Method Optimization." *Processes* 2021, **9**, 662.
    - HTML: https://www.mdpi.com/2227-9717/9/4/662
    - PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC8997410/

15. **Siebert, T. E.; Smyth, H. E.; Capone, D. L.; et al.** "Stable isotope dilution analysis of wine fermentation products by HS-SPME-GC-MS." *Anal. Bioanal. Chem.* 2005, **381**, 937-947.
    - DOI: 10.1007/s00216-004-2992-4
    - Springer: https://link.springer.com/article/10.1007/s00216-004-2992-4
    - PubMed: https://pubmed.ncbi.nlm.nih.gov/15660221/

16. **Hjelmeland, A. K.; Wylie, P. L.; Ebeler, S. E.** Various SPME Arrow / GC-MS Cabernet methods (Cabernet Sauvignon SPME-GC-MS validation examples). See: Oeno One 2021 SPME Arrow optimisation paper: https://oeno-one.eu/article/view/7914

### Carryover / identification

17. **ICH M10 — *Bioanalytical Method Validation and Study Sample Analysis*** (2022; in force).
    - PDF (training): https://database.ich.org/sites/default/files/ICH_M10_EWG_Training_Material_2024_0127.pdf
    - WHO version: https://cdn.who.int/media/docs/default-source/medicines/norms-and-standards/current-projects/qas23_925_bioanalytical-method-validation_rev03.pdf

18. **EMA Guideline on Bioanalytical Method Validation** (Committee for Medicinal Products for Human Use).
    - PDF: https://www.ema.europa.eu/en/documents/scientific-guideline/guideline-bioanalytical-method-validation_en.pdf

### Australian wine-industry context

19. **The Australian Wine Research Institute (AWRI) — analytical methods catalogue and method-validation pages.**
    - Volatile acidity method: https://www.awri.com.au/industry_support/winemaking_resources/laboratory_methods/chemical/va/
    - Smoke-taint analytical methods: https://www.awri.com.au/industry_support/winemaking_resources/smoke-taint/smoke-taint-analytical-methods/
    - Project 3.1.1 (volatile compounds responsible for wine flavour): https://www.awri.com.au/research_and_development/2017-2025-rde-plan-projects/project-3-1-1/
    - AWRI commercial services (Affinity Labs): http://www.awri.com.au/commercial_services/analytical_services/analyses/

### Top-down measurement uncertainty (laboratory medicine cross-reference)

20. **Padoan, A.; Antonelli, G.; Aita, A.; Sciacovelli, L.; Plebani, M.** "The top-down approach to measurement uncertainty: which formula should we use in laboratory medicine?" *Clin. Chem. Lab. Med.* 2020.
    - DOI: 10.1515/cclm-2019-1265
    - PMC: https://pmc.ncbi.nlm.nih.gov/articles/PMC7138004/

---

**End of brief.** Length, formula completeness, and acceptance-criteria coverage match the validation-chapter scope brief. Hand-off note for the QA pass: cross-check the §6 consolidated-criteria table against the analyte panel finalised by Specialists 1-4, then re-state the per-analyte numeric thresholds in tabular form for the chapter.
