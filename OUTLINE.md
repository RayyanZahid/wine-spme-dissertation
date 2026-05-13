# Dissertation Outline

Working table of contents. Section numbers reflect intended final structure. Each leaf is filled by agentic RAG over pgvector + literature extracts.

## 1. Abstract
Single 300-word paragraph: problem, method, novelty, validation summary, application result.

## 2. Introduction
- 2.1 Cabernet Sauvignon: economic importance and aroma complexity
- 2.2 Headspace SPME as the dominant pre-concentration method for wine VOCs
- 2.3 The SPME Arrow innovation: geometry, capacity, robustness
- 2.4 The shift from full-scan GC-MS to MRM GC-MS/MS for trace volatiles
- 2.5 Scope and research questions

## 3. Literature Review
- 3.1 Cabernet Sauvignon aroma chemistry
  - 3.1.1 Pyrazines (3-isobutyl-2-methoxypyrazine IBMP, 3-isopropyl- IPMP, 3-sec-butyl- SBMP)
  - 3.1.2 C13-norisoprenoids (β-damascenone, β-ionone, TDN, vitispirane)
  - 3.1.3 Fermentation esters (ethyl hexanoate / octanoate / decanoate, hexyl acetate, isoamyl acetate)
  - 3.1.4 Terpenes and alcohols (linalool, α-terpineol, geraniol, citronellol, nerol)
  - 3.1.5 Varietal sulfur compounds (3-mercaptohexanol 3MH, 3MHA, 4MMP — typically Sauv Blanc but reported in Cab Sauv at trace level)
  - 3.1.6 Fault markers (4-ethylphenol, 4-ethylguaiacol, TCA/TBA, geosmin, DMS, H2S, ethyl acetate)
  - 3.1.7 Smoke-taint compounds (guaiacol, 4-methylguaiacol, 4-methylsyringol, cresols, syringol)
- 3.2 SPME theory
  - 3.2.1 Three-phase equilibrium (sample / headspace / fiber coating)
  - 3.2.2 Partition coefficients and the dependence on temperature, salt, ethanol
  - 3.2.3 Equilibrium vs pre-equilibrium sampling
  - 3.2.4 The "carryover" tradeoff: sensitivity vs throughput
- 3.3 SPME Arrow geometry and sorbent choices
  - 3.3.1 Tube-in-needle architecture (CTC Analytics / Restek / BGB)
  - 3.3.2 Coating-volume scaling: 11 µL Arrow vs 0.6 µL traditional fiber
  - 3.3.3 Available phases: PDMS, DVB/PDMS, CWR/PDMS, DVB/CWR/PDMS, Tenax-TA, HLB
  - 3.3.4 Robustness, lifetime, septum compatibility
- 3.4 GC-MS/MS for wine analysis
  - 3.4.1 Triple-quadrupole basics: Q1 selection → Q2 CID → Q3 selection
  - 3.4.2 MRM vs SRM vs full-scan: sensitivity / specificity tradeoff
  - 3.4.3 Scheduled MRM (sMRM / dMRM) for cycle-time optimization
  - 3.4.4 Collision energy and the Mathieu-region selection
- 3.5 Method-validation frameworks (ICH Q2(R2), AOAC SMPR, EURACHEM, IUPAC)

## 4. Materials and Methods
- 4.1 Reagents, standards, internal standards
- 4.2 Model wine and authentic Cab Sauv samples
- 4.3 SPME Arrow device, fiber phases evaluated
- 4.4 GC-MS/MS instrument specifications
- 4.5 Column selection and oven programming
- 4.6 SPME extraction parameters (initial and optimized)
- 4.7 GC injection parameters (splitless, split, inlet temperature, liner geometry)
- 4.8 MS source parameters (transfer line, ion source temperature, CI vs EI)
- 4.9 Software (vendor data system + post-acquisition: openMS, MZmine, R, Python)

## 5. SPME Arrow Method Development
- 5.1 Sorbent screening
  - 5.1.1 Univariate screening: 5 sorbents × pooled wine matrix
  - 5.1.2 Selection logic per compound class
- 5.2 Multivariate optimization (DoE)
  - 5.2.1 Design selection: Plackett-Burman screening → CCD or Box-Behnken refinement
  - 5.2.2 Factors: extraction temperature, time, NaCl, agitation, sample volume, ethanol dilution
  - 5.2.3 Responses: per-analyte peak area, geometric-mean composite, signal-to-noise
  - 5.2.4 Desirability function for multi-response optimization
- 5.3 Extraction isotherms
  - 5.3.1 Time-course experiments per compound class
  - 5.3.2 Equilibrium vs pre-equilibrium decision
- 5.4 Carryover and fiber-conditioning protocol

## 6. MRM Optimization
- 6.1 Precursor-ion selection from EI spectra
- 6.2 Product-ion scan and Q3 selection logic
- 6.3 Collision-energy ramp and CE optimization per transition
- 6.4 Quan + qual transition pairing; ion-ratio acceptance bands
- 6.5 Scheduled-MRM segment planning
  - 6.5.1 Retention-time database build
  - 6.5.2 Dwell-time allocation: cycle-time = N × (dwell + interscan)
  - 6.5.3 Target 10–12 cycles per peak
- 6.6 The MRM transition table (panel of ~30 analytes × 2 transitions)

## 7. Method Validation
- 7.1 Linearity (matrix-matched calibration, 5–7 levels, ≥3 replicates)
- 7.2 Limits of detection and quantification (S/N, ICH ε definitions, blank-replicate variant)
- 7.3 Trueness via spike recovery (low / mid / high)
- 7.4 Precision: repeatability (intra-day) and intermediate precision (inter-day)
- 7.5 Matrix-effect evaluation (slope ratio, ≥0.85 acceptance)
- 7.6 Selectivity / specificity (interference at quan + qual masses)
- 7.7 Carryover (% of LOQ in blank after high standard)
- 7.8 Robustness (Youden 7-factor design at ±5 % per factor)

## 8. Results and Discussion
- 8.1 Optimized SPME Arrow parameters (per compound class)
- 8.2 The validated MRM transition portfolio
- 8.3 Quantitative aroma profile of n = 30 Cabernet Sauvignon samples (Napa / Bordeaux / Coonawarra / Maipo / Stellenbosch / Mendoza)
- 8.4 Multivariate differentiation (PCA, OPLS-DA)
- 8.5 Compound-class biomarker findings
- 8.6 Cross-comparison with traditional fiber-SPME data on the same sample set

## 9. Conclusions
- 9.1 Performance summary
- 9.2 Method positioning vs prior art
- 9.3 Limitations and future work

## 10. References
Numbered, ACS style. Generated by the synthesis pipeline from pgvector hits with `_meta.score` ≥ threshold.
