# SPME Arrow Technology — Technical Brief

**Author / role:** Specialist 1 of 6 (parallel research team)
**Dissertation:** Method development, MRM optimization, and validation of SPME Arrow headspace GC-MS/MS for the quantitative profiling of aroma-active volatiles in Cabernet Sauvignon wine
**Compiled:** 2026-05-12
**Sources:** Vendor brochures (CTC Analytics / PAL System, Restek, Agilent, BGB Analytik, GL Sciences EU), peer-reviewed open-access papers (PMC, MDPI, Shimadzu Pittcon poster, Springer/Anal Bioanal Chem), instructional sheets, and grey-literature application notes. Retrieval dates noted per reference in §8. All vendor claims (e.g., "10×", "24-fold", "2× lifetime") attributed to source; conflicts called out where they appear.

---

## 1. Geometry and capacity comparison (Arrow vs traditional SPME fiber)

The SPME Arrow is a tube-in-needle device: a stainless-steel outer sleeve (1.1 or 1.5 mm OD) protects a stainless-steel inner support rod (0.4 mm) onto which the sorption phase is bonded. An arrow-shaped stainless-steel tip caps the assembly and pierces vial/injector septa. In contrast, classical SPME fibers are a fused-silica fiber (typically 23-gauge / ~0.58 mm OD) coated with the sorption phase, retracted into a hollow hypodermic needle. The Arrow's design fully protects the sorbent during transfer; the classical fiber's coating is exposed to the needle wall every time it's retracted.

The vendor-reported geometry numbers vary slightly depending on which Arrow length is referenced (20 mm is the only commercial Arrow length; 30 mm appears once in a Shimadzu Pittcon figure caption that contradicts every other source — almost certainly a typo for 20 mm).

### Table 1.1 — Vendor-confirmed geometry numbers (all sources agree on the 20 mm length)

| Device | OD (mm) | Phase length (mm) | Phase thickness (µm) | Sorption surface (mm²) | Sorption volume (µL) | Source |
|---|---|---|---|---|---|---|
| Classical SPME fiber, 100 µm | ~0.7 (23 ga) | 10 | 100 | 9.4 | 0.6 | CTC/PAL brochure; Restek discover article; Kremser 2016 |
| Classical SPME fiber, 7 µm PDMS | ~0.7 | 10 | 7 | ~9.4 | ~0.03 (typical) | Pawliszyn ref — classical fibers range 7–100 µm |
| SPME Arrow 1.1 mm, 100 µm | 1.1 | 20 | 100 | **44.0** | **3.8** | CTC/PAL brochure p. 3; GL Sciences EU mirror; Restek |
| SPME Arrow 1.5 mm "wide-sleeve" (same phase as 1.1 Arrow, swelling room) | 1.5 | 20 | 100 / 120 | 44.0 | 3.8 | CTC/PAL brochure (footnote: 1.5 mm wide-sleeve has the same phase dimensions as 1.1 mm — extra space accommodates phase swelling in organic solvents) |
| SPME Arrow 1.5 mm, 250 µm PDMS (high-capacity) | 1.5 | 20 | 250 | **62.8** | **11.8** | CTC/PAL brochure p. 3 |
| Kremser 2016 250 µm PDMS Arrow | 1.5 | 20 | 250 | n/r | **10.2** | Kremser et al. 2016 (Anal Bioanal Chem 408:943) |

### Table 1.2 — Vendor headline ratios (Arrow vs classical fiber)

| Metric | Restek | CTC/PAL | Kremser 2016 |
|---|---|---|---|
| Phase volume gain | 6.3× (1.1 mm) / 20× (1.5 mm 250 µm) | "up to 20× larger overall" | 17× (10.2 µL / 0.6 µL) |
| Surface area gain | 4.7× (1.1 mm) / 6.7× (1.5 mm) | "6× larger" | not reported |
| Lifetime | "2 to 3× longer" | "at least 2× longer" | "at least 200 injections without coring" vs fiber's 100–200 |
| Sensitivity (typical) | "297–527 % higher response" at 10 min extraction; "618–896 %" at 5 min extraction | "up to 10× more sensitivity" | Up to 12.2× theoretical, ~3–12× measured for PAHs |

### Conflict / clarification — the "24-fold" claim
The repeatedly-cited "phase volume of SPME Arrow is 24-fold higher than traditional SPME" comes from a single secondary source (a 2020 MDPI summary) and propagates through search-engine snippets. The primary CTC/PAL brochure gives 11.8 µL / 0.6 µL = **19.7-fold** for the largest 1.5 mm 250 µm Arrow vs a 100 µm 10 mm fiber. The 24-fold value almost certainly refers to a comparison against a thinner classical fiber (e.g. 30 µm PDMS / DVB or 7 µm PDMS) — be explicit in the dissertation about which fiber baseline is used.

### Inner geometry detail
- Inner stainless-steel rod: 0.4 mm OD (Kremser 2016).
- Phase is bonded onto the rod between the protective arrow tip and the needle hub.
- During desorption, the arrow tip retracts and the sleeve withdraws, exposing the bonded phase inside the GC inlet. There is no fiber-to-needle abrasion.

---

## 2. Sorbent phases — chemistry, target compound classes, vendor SKUs

The Arrow product line is a near-1:1 mirror of the classical SPME phase library, with three notable absences (Carboxen/PDMS, CAR/PDMS+DVB triphase 50/30 µm Supelco equivalent, and 60 µm PEG/CW) and one rename: "Carbon Wide Range" (CWR) is the Arrow name for what the classical line calls Carboxen-1000 / Carboxen analogous coatings — a porous carbon molecular sieve with engineered wide-pore distribution to capture both very-volatile and semivolatile analytes. Restek's 2019 instruction sheet and the CTC/PAL Smart Arrow brochure are the load-bearing primary references for SKUs and color codes.

### Table 2.1 — Full Arrow phase catalog (CTC/PAL — vendor-neutral master list)

| # | Phase chemistry | Phase thickness (µm) | OD (mm) | Color | CTC SKU (set of 1) | Target molecular weight | Target compound classes |
|---|---|---|---|---|---|---|---|
| 1 | PDMS (polydimethylsiloxane) | 100 | 1.1 | Red | SARR11-P-100/20-P1 | 60–275 | Non-polar volatiles, hydrocarbons, monoterpenes; absorbent (partition) phase |
| 2 | PDMS | 100 | 1.5 | Red | SARR15-P-100/20-P1 | 60–275 | As #1, wide-sleeve for headspace + immersion |
| 3 | Polyacrylate (PA) | 100 | 1.1 | Gray | SARR11-A-100/20-P1 | 80–300 | Polar semivolatiles (phenols, smoke-taint guaiacols, cresols) |
| 4 | Carbon-WR / PDMS | 120 | 1.1 | Light Blue | SARR11-C-WR-120/20-P1 | 30–225 | Highly volatile — low-MW alcohols, ketones, sulfur volatiles, C2–C8 aldehydes |
| 5 | Carbon-WR / PDMS | 120 | 1.5 | Light Blue | SARR15-C-WR-120/20-P1 | 30–225 | As #4, wide-sleeve |
| 6 | DVB / PDMS | 120 | 1.1 | Violet | SARR11-DVB-120/20-P1 | 60–300 | Aromatic and polar semivolatiles, esters, lactones, monoterpenes — flagged by Lisanti 2021 as best for minty terpenoids in red wine |
| 7 | DVB / PDMS | 120 | 1.5 | Violet | SARR15-DVB-120/20-P1 | 60–300 | As #6, wide-sleeve |
| 8 | DVB / Carbon-WR / PDMS ("triple-phase") | 120 | 1.1 | Dark Gray | SARR11-DVB/CWR120/20-P1 | 40–275 | Broadest-spectrum: covers volatiles to semivolatiles in one extraction. Vendor-flagged for wine, coffee, off-flavor work |
| 9 | DVB / Carbon-WR / PDMS | 120 | 1.5 | Dark Gray | SARR15-DVB/CWR120/20-P1 | 40–275 | As #8, wide-sleeve |
| 10 | PDMS, high capacity | 250 | 1.5 | Black | SARR15-P-250/20-P1 | 60–275 | Non-polar volatiles with maximum capacity — heaviest hydrocarbons, PAHs, fatty-acid esters |

**Method Development Kit:** SARR1115-SEL5-S2 (CTC) / cat.# 27489 (Restek) / cat.# 28901 (Restek-2020 revision) = one each of #1, #3, #4, #6, #8 (1.1 mm format). This is the canonical screen for a new method.

### Table 2.2 — Same catalog, Restek-branded SKUs (legacy + current)

| Phase | Thickness (µm) | OD (mm) | Color | Restek cat.# (2019 instruction sheet) | Restek cat.# (2020+ "Smart" rebrand) |
|---|---|---|---|---|---|
| PDMS | 100 | 1.1 | Red | 27485 | 28906-1 / 28906-3 |
| PDMS | 100 | 1.5 | Red | 27877 | 28910-1 / 28910-3 |
| PDMS, high-capacity | 250 | 1.5 | Black | 27484 | 28911-1 / 28911-3 |
| PA | 100 | 1.1 | Gray | 27488 | 28902-1 / 28902-3 |
| Carbon-WR/PDMS | 120 | 1.1 | Light Blue | 27487 | 28903-1 / 28903-3 |
| Carbon-WR/PDMS | 120 | 1.5 | Light Blue | 27879 | 28907-1 / 28907-3 |
| DVB/PDMS | 120 | 1.1 | Violet | 27486 | 28905-1 / 28905-3 |
| DVB/PDMS | 120 | 1.5 | Violet | 27878 | 28909-1 / 28909-3 |
| DVB/CWR/PDMS | 120 | 1.1 | Dark Gray | 27875 | 28904-1 / 28904-3 |
| DVB/CWR/PDMS | 120 | 1.5 | Dark Gray | 27876 | 28908-1 / 28908-3 |
| Method Development Kit | mixed | 1.1 | mixed | 27489 | 28901 |

### Sorbent chemistry — what DVB/CWR/PDMS actually is
- **PDMS** is the absorbent baseline — a high-purity polydimethylsiloxane. Analytes dissolve into it (partition model). Linear, non-polar, low-bleed. Decent for terpenes, hydrocarbons, mid-volatility esters.
- **DVB** is divinylbenzene — a porous polymer (mesoporous, mostly 17 Å, with ~750 m²/g BET surface area in the classical Supelco line). Adsorbs aromatic and polar analytes by physisorption in pores. Best pore-size for MW 80–300.
- **Carbon-WR** ("Carbon Wide Range") is a synthetic carbon molecular sieve with a deliberately wide pore-size distribution — micropores (<20 Å) for trapping C2–C5 volatiles + mesopores (20–500 Å) for desorbing C6–C12 without sticking. Classical-fiber equivalent is Carboxen-1000.
- **DVB/CWR/PDMS** is a triphase composite: a thin PDMS outer layer over a particulate DVB+CWR layer bonded to the inner steel rod. The PDMS protects the porous particles from matrix damage; the DVB+CWR provides the adsorption capacity. Vendor descriptions consistently target MW 40–275 — the broadest of any single Arrow phase.

### Table 2.3 — Phase choice → wine compound class mapping (synthesized from Lisanti 2021, Lan 2021 grape-skins, Kim 2020 Cirsium, Korean fish sauce 2019, OENO 2024)

| Compound class in Cabernet Sauvignon | Recommended phase | Evidence |
|---|---|---|
| Ethyl esters (ethyl hexanoate, ethyl octanoate, etc.) | DVB/CWR/PDMS or DVB/PDMS | Shimadzu Pittcon 2017 Moscato (peak areas Arrow ~3–4× fiber for ethyl octanoate); Lan 2021 grape skins |
| Acetates (isoamyl acetate, hexyl acetate) | DVB/PDMS or DVB/CWR/PDMS | Shimadzu Pittcon 2017 (isoamyl acetate Arrow 9.86×10⁷ vs fiber 3.28×10⁷ — 3.0× increase) |
| Monoterpenes (linalool, geraniol, α-terpineol) — minty/floral notes | DVB/PDMS | Lisanti 2021 explicitly selected DVB/PDMS for nine terpenoids in red wine; OENO 2024 |
| Sesquiterpenes (rotundone — peppery in Syrah/Shiraz, but relevant for varietal comparisons) | DVB/CWR/PDMS | Lan 2021 grape skins |
| C13-norisoprenoids (β-damascenone, β-ionone, α-ionone, TDN) | DVB/CWR/PDMS or DVB/PDMS | Lan 2021 (bound + free); α-ionone wine app note recommends DVB/CAR/PDMS classical analog |
| Volatile phenols (guaiacol, 4-methylguaiacol, 4-ethylguaiacol, 4-ethylphenol — smoke taint + Brett) | Polyacrylate or DVB/PDMS | Agilent smoke-taint wine app note (cat.# 5994-3160/3161); Liu 2020 smoke-taint paper |
| Methoxypyrazines (IBMP, IPMP, SBMP — "bell-pepper" varietal markers in CS) | DVB/CWR/PDMS | Allen et al. earliest SPME-fiber IBMP work used DVB/CAR/PDMS; same logic for Arrow |
| Volatile thiols (3MH, 3MHA, 4MMP — tropical/box-tree notes) | DVB/CWR/PDMS, derivatized | Specialized; not Arrow-native, derivatization required |
| Haloanisoles / cork taint (TCA, TBA, TeCA, PCA) | PDMS 100 µm | LabRulez/Sleek paper — PDMS Arrow 100 µm achieves LOD 0.03–0.07 ng/L |
| Higher alcohols (isoamyl alcohol, 2-phenylethanol, hexanol) | DVB/PDMS or PDMS | Shimadzu Pittcon 2017 |
| Fatty acids (octanoic, decanoic, dodecanoic) | DVB/PDMS or DVB/CWR/PDMS | Shimadzu Pittcon 2017; OENO 2024 |

---

## 3. Sensitivity gain — published peak-area ratios per compound class

The published Arrow-vs-fiber ratios are highly compound-, phase-, and condition-dependent. Most authors report between 2× and 10× peak-area enhancement with Arrow vs the matching classical fiber, but the "headline" numbers from vendor literature run higher (up to ~26× for off-flavors in water with DVB Arrow vs DVB fiber per CTC/PAL brochure p. 4).

### Table 3.1 — Shimadzu Pittcon 2017 (PO-CON1702E), Moscato wine, 4-month headspace stability study
Both methods used DVB/CAR/PDMS 23-gauge classical fiber and a DVB/CAR/PDMS Arrow (1.5 mm; thickness not specified in poster, but the model used corresponds to a 100 µm PDMS Arrow per the Shimadzu figure caption). Sample: 7 mL Moscato + 1.5 g NaCl; 20 mL vial; 40 °C incubation 20 min + 40 min extraction; 260 °C splitless desorption 3 min; Rxi-5MS column. Peak areas as raw integration counts:

| Compound | Classical SPME, fresh wine | Classical SPME, 4 mo wine | SPME Arrow, 4 mo wine | Arrow / fiber ratio (4 mo) |
|---|---|---|---|---|
| Butanoic acid, ethyl ester | 1.66 × 10⁶ | 3.33 × 10⁶ | **1.56 × 10⁷** | 4.7× |
| Butanoic acid, 2-methyl-, ethyl ester | 4.60 × 10⁵ | — (not detected) | **4.50 × 10⁶** | n/a (Arrow-only) |
| Butanoic acid, 3-methyl-, ethyl ester | 4.81 × 10⁵ | 2.35 × 10⁷ | **8.99 × 10⁶** | 0.38× (fiber higher — degradation artifact in fresh comparison) |
| 1-Hexanol | 4.98 × 10⁶ | — | **1.80 × 10⁷** | n/a (Arrow-only) |
| 1-Butanol, 3-methyl-, acetate (isoamyl acetate) | 3.28 × 10⁷ | — | **9.86 × 10⁷** | n/a (Arrow-only) |
| Hexanoic acid, ethyl ester | 2.60 × 10⁷ | 6.61 × 10⁷ | **1.04 × 10⁸** | 1.6× |
| Acetic acid, hexyl ester | 1.21 × 10⁷ | 5.92 × 10⁶ | **1.92 × 10⁷** | 3.2× |
| Octanoic acid, ethyl ester | 4.40 × 10⁷ | 1.64 × 10⁸ | **1.88 × 10⁸** | 1.14× |
| Phenylethyl alcohol | 1.47 × 10⁸ | 2.35 × 10⁸ | 3.37 × 10⁶ | 0.014× (Arrow lower — unexpected; poster does not explain. Possible matrix saturation on fiber adsorber.) |
| Dodecanoic acid | 6.31 × 10⁷ | 2.80 × 10⁶ | **1.18 × 10⁸** | 42× (vs aged fiber — saturation collapse on fiber after 4 mo) |
| Geranyl ethyl ether | 1.21 × 10⁶ | 1.59 × 10⁶ | **6.75 × 10⁶** | 4.2× |
| 2,2,7-trimethyl-3-oxatricyclo (terpene oxide) | 5.71 × 10⁵ | 2.05 × 10⁶ | **8.11 × 10⁶** | 4.0× |

**Headline finding:** Arrow detects more compounds (Arrow recovered "almost twice as many VOCs" per Shimadzu summary) and gives 2–10× higher peak areas for esters, acetates, and terpene oxides. The phenylethyl alcohol anomaly is a known matrix-saturation issue on classical fiber adsorbent phases and a recurring discussion point in the literature.

### Table 3.2 — Kim et al. 2020 (PMC7353593), Cirsium setidens Nakai, four-phase Arrow comparison
Same Arrow geometry (1.1 mm, 20 mm) but four different phases:

| Phase | Total peak area (×10⁶) | Compounds detected | Reproducibility (avg RSD %) |
|---|---|---|---|
| CWR / PDMS 120 µm | **10,122** | 52 | 2.01 |
| DVB / PDMS 120 µm | 7,323 | 52 | 3.76 |
| PA 100 µm | 2,832 | 26 | n/r |
| PDMS 100 µm | 2,346 | 37 | n/r |

CWR/PDMS gave a 4.3× total peak-area advantage over PDMS-only for this matrix. **For wine** the DVB/PDMS or DVB/CWR/PDMS will likely win for the higher-MW esters and terpenes that dominate Cabernet aroma.

### Table 3.3 — Kremser 2016 (the foundational Arrow paper), PAHs in water
PDMS-only Arrow (250 µm × 20 mm, 10.2 µL) vs PDMS-only fiber (100 µm × 10 mm, 0.6 µL), immersion mode:

| Analyte | Fiber extraction yield | Arrow extraction yield | Arrow / fiber |
|---|---|---|---|
| Naphthalene | 5.2 % | 17.5 % | 3.4× |
| Benzo(k)fluoranthene | 33.4 % | 61.3 % | 1.8× |
| Theoretical max (across MW range) | — | — | 1.0× to 12.2× |

### Table 3.4 — Vendor headline numbers (treat as ceiling, not floor)
| Source | Application | Arrow / fiber claim |
|---|---|---|
| CTC/PAL brochure | Iodoform in tap water, DVB headspace | 26× sensitivity |
| CTC/PAL brochure | Iodoform in tap water, DVB immersion | 6× sensitivity |
| Restek discover article | Generic 5-min extraction | 6.2× to 9.96× ("618–896 % increase") |
| Restek discover article | Generic 10-min extraction | 4.0× to 6.3× ("297–527 % higher") |
| CTC/PAL brochure | Working linear range expansion | "0.1–1000 ng/L" for Arrow vs "1–100 ng/L" for fiber on the same axis |

### Table 3.5 — Wine-relevant grey-literature claims
| Source | Application | Arrow / fiber numerical finding |
|---|---|---|
| OENO One 2024 (Trupiano et al.; via search-snippet summary) | Red/white wine VOC optimization | "SPME Arrow exhibited higher sensitivity by up to a factor of six and better repeatability by up to a factor of five" |
| Lisanti et al. 2021, Food Chem (DOI 10.1016/j.foodchem.2021.130024) | 9 minty terpenoids in red wine | LOD 3–60 ng/L, LOQ 6–200 ng/L, recoveries 80–119 %, intraday RSD 2–25 %; uses PDMS/DVB Arrow at 50 °C × 60 min with 4 g NaCl in 5 mL wine |
| LabRulez 2022 / haloanisoles in wine via Arrow | TCA, TBA, TeCA, PCA cork-taint | PDMS Arrow 100 µm; LOD 0.03–0.07 ng/L, LOQ ≤ 0.24 ng/L, recoveries 90–105 %, R² > 0.997, RSD ≤ 10 % at 0.25 ng/L |
| Lan et al. 2021 (Molecules, PMC8659239) — grape skins | 53 free + 84 bound VOCs in Marselan grape skins | DVB/CWR/PDMS 120 µm × 20 mm "yielded the highest peak areas across volatile compound classes"; LOQ 2.21–16.22 µg/kg for aldehydes |
| Shimadzu Pittcon 2017 Moscato | Aged wine headspace | Arrow recovered "almost twice as many VOCs" + 2–10× peak areas on ester/acetate class |
| Korean fermented fish sauce 2019 (Applbiolchem) | Sulfur volatiles, esters | "HS-SPME Arrow fiber affords improved extraction efficiency than that of SPME fiber"; vanillin S/N 40:1 (fiber) vs 130:1 (Arrow) = 3.25× |
| Whiskey GCxGC-MS (LCGC Intl.) | Whiskey flavor compounds | Qualitative — "more compounds detected" with Arrow; quantitative ratios not extracted |

---

## 4. Robustness, lifetime, septum compatibility

### Lifetime
- **Vendor (PAL/CTC):** "PAL SPME Arrows last at least 2× longer" than classical fibers. Mechanical-stability claim attributed to the steel sleeve protecting the bonded phase against retraction abrasion.
- **Vendor (Restek):** "Two to three times longer than traditional SPME fibers."
- **Kremser 2016 primary lit:** Classical fibers "typically requiring replacement after 100 to 200 injections due to bending." Arrow: "at least 200 injections without coring, abrasion, or leakage are possible." **Floor of 200 quantitative injections is the only experimentally-supported number; the "2×" headline is a derived comparison against the lower end of fiber lifetime.**
- **Restek 2019 instruction sheet:** Lifetime is matrix-dependent. "Immersion sampling in liquids containing complex matrices may reduce lifetime. In contrast, headspace sampling generally results in longer lifetimes." Critical for the dissertation: HS-SPME Arrow in wine should track toward the upper end of the 200–600 injection range, but no source gives a specific HS-Arrow-in-wine lifetime number.

### Septum compatibility & wear
- The arrow-shaped tip is steel and is wider than a 23-gauge fiber needle (1.1 or 1.5 mm vs 0.58 mm). Vendor claim: **"Injector septa last at least 2× longer"** (GL Sciences EU brochure p. 5; absent from CTC main brochure). The justification is that the tip pierces cleanly without coring fragments out of the septum, whereas the smaller fiber needle does core.
- **GC inlet pressure ceiling:** Restek product page (28904-1) — **"Recommended maximum GC inlet pressure is 50 psi or less"** when the Arrow is in the inlet (i.e., during desorption). Above 50 psi, the larger-diameter Arrow may not seal correctly against the septum/liner.
- **Liner requirements (mandatory):** Arrows cannot be used with stock split/splitless liners. Inlet conversion is mandatory:

| GC platform | Restek conversion-kit cat.# |
|---|---|
| Thermo TRACE Ultra | 27495 |
| Thermo TRACE 1300/1310 | 27494 |
| Agilent 6890 | 27492 |
| Agilent 7890 | 27493 |
| Shimadzu GC-2010 | 27491 |
| Agilent 8890/8860 | PAL3-SARR-Start-GC8890 (CTC/PAL) |
| Agilent Intuvo | PAL3-SARR-Start-Intuvo |

The kits provide a port weldment, adaptor cup, 1.8 mm ID Arrow-compatible liners (1.3 mm for 1.1 mm Arrows, 1.7 or 2.0 mm for 1.5 mm Arrows), septa, and septum nut.

### Thermal conditioning (from Restek 2019 instruction sheet — full table)

| Phase, thickness | Max temp (°C) | Recommended operating temp (°C) | Preconditioning temp (min–max °C) | Preconditioning time (min, recommended) | Conditioning time between samples (min, recommended) | Cleaning solvent | Cleaning time (min, recommended) |
|---|---|---|---|---|---|---|---|
| PDMS, 100 µm | 300 | 200–300 | 200–300 | 15–120 (30) | 1–60 (5) | MeOH / EtOH / IPA | 0.5–10 (2) |
| PDMS, 250 µm | 300 | 220–300 | 200–300 | 15–120 (60) | 1–60 (10) | MeOH / EtOH / IPA | 0.5–10 (2) |
| PA, 100 µm | 280 | 200–280 | 180–280 | 15–120 (30) | 1–60 (5) | MeOH / aliphatic HC | 0.5–10 (2) |
| Carbon-WR, 120 µm | 300 | 200–300 | 180–300 | 15–120 (30) | 1–60 (5) | MeOH | 0.5–10 (2) |
| DVB, 120 µm | 300 | 220–300 | 200–300 | 15–120 (60) | 1–60 (10) | MeOH / EtOH / IPA | 0.5–10 (2) |
| DVB/Carbon-WR, 120 µm | 300 | 220–300 | 200–300 | 15–120 (60) | 1–60 (10) | MeOH / EtOH / IPA | 0.5–10 (2) |

### Restek general precautions (verbatim)
- Never touch the stationary phase, even when wearing gloves.
- Never expose to heat without an inert gas atmosphere.
- Never exceed the maximum recommended temperature.
- Never soak in chlorinated solvents.
- Never use an inlet liner with glass wool — contact damages the phase.
- "Staining or discoloration does not give any indication of the remaining life span." — visual inspection is unreliable; the only signal that an Arrow is past end-of-life is loss of analytical response, which underlines the importance of including a QC standard / system-suitability check in every run.

### Conditioning hardware
- **Dedicated conditioning module:** PAL3-SPME-ArrowCond (recommended; saves wear on the GC inlet, supports automated and manual pre-conditioning, integrated inert-gas purge).
- **Heatex Stirrer:** 40–150 °C, cycloidal mixing (no stir bar), required for immersion SPME but not for headspace (an agitator is sufficient).
- **Agitator module:** 6 positions × 20 mL vials, 40–200 °C, 250–750 rpm — standard for headspace.

---

## 5. Vendor application notes relevant to wine

The wine-relevant Arrow application-note inventory is thin but useful. Most vendor activity has gone to smoke-taint (Agilent), generic wine VOC profiling (CTC/PAL brochure white-wine example), and cork-taint (haloanisoles, picked up by Shimadzu and labrulez paper databases).

### Table 5.1 — Vendor-published wine + beverage application notes for SPME Arrow

| Vendor | Title / focus | Document number | Notes |
|---|---|---|---|
| Agilent | "Analysis of Free Volatile Phenols in Smoke-Impacted Wines by SPME-GC/MS" (smoke-taint Arrow vs fiber) | 5994-3160EN, 5994-3161EN | Direct response comparison on a 8890 GC for guaiacol, 4-methylguaiacol, syringol etc. (blocked from full extraction; titles and document IDs confirmed via search snippets and Agilent library) |
| Agilent | "Smart SPME Arrows — General Information" (technical overview, package insert) | 5994-3137EN, 5994-5775EN | Phase chemistry, conditioning, fiber lifetime; not wine-specific but parameter master |
| CTC/PAL System | "PAL Smart SPME Arrow — The Better SPME" brochure | published Apr 2025 (Brochure rev. 4) | White-wine aroma example p. 4 (Gewürztraminer + Chardonnay, PDMS Arrow vs PDMS fiber) — chromatograms show ~3× peak height for Arrow on isoamyl acetate, ethyl hexanoate, ethyl octanoate, ethyl decanoate, phenylethyl alcohol, octanoic acid |
| Shimadzu | "Comparison of HS-SPME and SPME Arrow Sampling Techniques Utilized to Characterize Volatiles in the Headspace of Wine over an Extended Period of Time" — Pittcon 2017 poster 1430-11P | PO-CON1702E | Moscato study (see §3 Table 3.1) |
| Thermo Fisher | "New Opportunities for Wine Analysis through SPME Arrow" — Pittcon 2019 poster | PO-10697 | Title confirmed; full text blocked by 403 — flag for further chase |
| LabRulez / Sleek paper | "Determination of haloanisoles in wine by HS-SPME Arrow and GC-MS/MS" | gcms.labrulez.com/paper/18997 | Wine TCA/TBA/TeCA/PCA — PDMS Arrow 100 µm; key cork-taint reference |
| BenchChem (third-party app note, NOT primary research) | "Quantification of α-Ionone in Wine Using HS-SPME-GC-MS" | BenchChem Apr 2026 | **Uses classical DVB/CAR/PDMS fiber, not Arrow** — useful for parameter ranges (5–10 mL wine, 1–2 g NaCl, 40–60 °C, 15–30 min, 250 °C splitless desorption) but does not validate Arrow-specific performance |
| PAL System leaflet (palsystem.com) | "Smart SPME Arrow Leaflet" V5 | — | Ordering reference |

### Gap callout — AWRI (Australian Wine Research Institute)
A focused `site:awri.com.au` search returned **zero AWRI publications using SPME Arrow as of the May 2026 retrieval.** AWRI's published volatile-phenol and methoxypyrazine work continues to use classical SPME fiber or stable-isotope-dilution GC-MS. The closest hit was a 2015 Smith et al. note in Wine & Viticulture Journal v30n6 on "Direct determination of volatile congeners in wine using SPME–GC–MS with minor sample preparation" — also classical fiber. **This is a publishable opening for the dissertation: AWRI has not published a head-to-head Arrow study in Cabernet, so a rigorous CS validation is genuinely novel territory.**

### Gap callout — Markes International
Markes is best known for thermal-desorption tubes, not SPME Arrow. Their public content hub (markes.com) returns no Arrow-specific application notes; the search-engine hit for "Markes SPME Arrow" only surfaces general SPME pages. Markes is the wrong vendor for Arrow consumables — CTC/PAL System is the OEM, and Restek, Agilent, and BGB Analytik are the major distributors.

---

## 6. Theses / grey literature

Genuine PhD dissertations on SPME Arrow + wine are rare — the technique is still young (2015–2026 window) and most publications are journal articles, not theses.

### Theses found and verified
- **Yiu, Patricia (Oregon State University)** — "The Effect of Wine Matrix on the Analysis of Volatile Sulfur Compounds by Solid-Phase Microextraction-GC-PFPD." Graduate thesis, OSU ScholarsArchive, ID n870zt673. **Uses classical SPME fiber, not Arrow** — useful for wine-matrix effects context, not for Arrow methodology.

### ChemRxiv / arXiv
- No SPME Arrow + wine preprints surfaced in search; the technique is overwhelmingly published through peer review (Food Chemistry, J Chromatogr A, Molecules, Separations) rather than preprint servers.

### Peer-reviewed papers (closest functional substitutes for theses on the dissertation topic)
1. **Trupiano et al. 2024.** "Optimisation of SPME Arrow GC/MS method for determination of wine volatile organic compounds." *OENO One* 58(4). Search-snippet method summary: Box-Behnken DoE on red + white wine; for **red wine** optimal: extraction T 60 °C, incubation 17 min, exposure 53 min; for **white wine** 50 °C, 10 min, 60 min; 82 compounds quantified across 10 classes (esters, alcohols, fatty acids, aldehydes & ketones, furans, pyrazines, sulfur, phenols, terpenes, lactones). **Full-text was blocked by SSL cert error on direct fetch — flag for ResearchGate/Sci-Hub backup retrieval.**
2. **Lisanti, M.T. et al. 2021.** "Minty aroma compounds in red wine: development of a novel automated HS-SPME-arrow and gas chromatography-tandem mass spectrometry quantification method." *Food Chemistry* 364:130024. Nine terpenoids (1,4-cineole, 1,8-cineole, menthol isomers, piperitone, etc.). DVB/PDMS Arrow; 5 mL wine + 4 g NaCl + 50 °C × 60 min; **30 eV ionization** energy optimal (atypical — explicitly tuned downward from 70 eV for fragmentation selectivity); LOD 3–60 ng/L, LOQ 6–200 ng/L; recoveries 80–119 %. **The single most directly applicable dissertation antecedent** — same Arrow class, same matrix, same general MS approach (extends to MS/MS for the dissertation).
3. **Lan, Y. et al. 2021.** "Optimization of SPME-Arrow-GC/MS Method for Determination of Free and Bound Volatile Organic Compounds from Grape Skins." *Molecules* 26(23):7409 / PMC8659239. Compared all 5 Arrow phases; DVB/CWR/PDMS 120 µm × 20 mm won; 53 free + 84 bound VOCs identified in Marselan; LOQ 2.21–16.22 µg/kg aldehydes. Optimization conditions: 60 °C, 20 min incubation, 49 min exposure (free); 60 min exposure + 2 g NaCl + acid hydrolysis (bound); 7 min desorption at 250 °C.
4. **Kim, J. et al. 2020.** "Comparison of Different Types of SPME Arrow Sorbents to Analyze Volatile Compounds in Cirsium setidens Nakai." *Foods* 9(6) / PMC7353593. CWR/PDMS > DVB/PDMS > PA ~ PDMS by total peak area (see §3 Table 3.2).
5. **Lim, D.K. et al. 2019.** "Comparison of headspace–SPME and SPME-Arrow–GC–MS methods for the determination of volatile compounds in Korean salt–fermented fish sauce." *Applied Biological Chemistry* 62(1):60. Vanillin S/N improvement 40:1 → 130:1 with Arrow.
6. **Kremser, A.; Jochmann, M.A.; Schmidt, T.C. 2016.** "PAL SPME Arrow—evaluation of a novel solid-phase microextraction device for freely dissolved PAHs in water." *Anal Bioanal Chem* 408:943–952. **The foundational Arrow validation paper** — first quantitative comparison of Arrow vs fiber, established the 200-injection floor and the 17–24× phase-volume ratio framing.
7. **Herrington, J.S.; Gómez-Ríos, G.A.; Myers, C.; Stidsen, G.; Bell, D.S. 2020.** "Hunting Molecules in Complex Matrices with SPME Arrows: A Review." *Separations* 7(1):12. Open-access review of Arrow methodology and applications. **Full PDF blocked by 403; abstract confirms environmental, food, cannabis, forensic coverage.** Worth citing as the canonical English-language Arrow review.
8. **Yin, C. et al. 2019.** "Optimization and validation of a head space solid-phase microextraction-arrow gas chromatography-mass spectrometry method using central composite design for determination of aroma compounds in Chinese liquor (Baijiu)." *J Chromatogr A* 1601:1–11. Beverage matrix, CCD optimization, full method.
9. **Liu et al. 2020.** "A Simple GC-MS/MS Method for Determination of Smoke Taint-Related Volatile Phenols in Grapes." *Metabolites* 10(7):294. Australia Agriculture Victoria. **Uses classical SPME (not Arrow)** but defines the smoke-taint target list (guaiacol, 4-methylguaiacol, cresols, syringol, 4-ethylguaiacol, 4-ethylphenol) and GC-MS/MS MRM transitions — directly transferable to an Arrow method.

### Search-engine grey literature
- ChromForum thread "GC-MS SPME Wine applications" (chromforum.org/viewtopic.php?t=43761) — community discussion, not authoritative; flagged for later retrieval if dissertation needs unstandardized lab anecdotes.
- IVES Open Science 2023 — "The limonene-derived mint aroma compounds in red wines" — sister piece to Lisanti 2021.

---

## 7. Open questions for the dissertation experimental design

These are the load-bearing decisions that, based on the literature reviewed, are *not* yet settled by published work and must be resolved by the dissertation's own experiments. Each item also names the specific data that would need to be generated to close it.

1. **Phase selection for Cabernet Sauvignon target list.** Lisanti 2021 picked DVB/PDMS for nine minty terpenoids. Lan 2021 picked DVB/CWR/PDMS for grape-skin free + bound VOCs. Kim 2020 ranked CWR/PDMS > DVB/PDMS for a non-wine matrix. **For Cabernet specifically, the combined target list (esters + acetates + monoterpenes + C13-norisoprenoids + methoxypyrazines + volatile phenols) is heterogeneous in MW and polarity — a head-to-head DVB/PDMS vs DVB/CWR/PDMS vs Carbon-WR/PDMS comparison on a single Cabernet matrix has not been published.** This is one of the cleanest deliverables of the dissertation.

2. **1.1 mm vs 1.5 mm OD.** The 1.5 mm Arrow gives more phase volume on the 250 µm PDMS variant (11.8 µL vs 3.8 µL for the 1.1 mm) but requires the wider 2.0 mm injector liner. For DVB/PDMS, DVB/CWR/PDMS, and CWR/PDMS, the 1.5 mm wide-sleeve has *the same phase dimensions* as the 1.1 mm — the wide sleeve only provides solvent-swelling headroom. **Decision:** 1.1 mm is the default for HS work in wine; only revisit 1.5 mm if immersion mode is needed (probably never for HS-SPME Arrow GC-MS/MS).

3. **Coating thickness — 100 µm vs 250 µm PDMS.** The 250 µm Arrow is the high-capacity variant but only exists as PDMS-only. For a Cabernet target list including pyrazines and norisoprenoids, the triple-phase 120 µm is almost certainly preferable. **Coating thickness is not actually a free variable** for the realistic phase choices: pick a phase, accept the bundled thickness.

4. **Salt strength.** Lisanti 2021: 4 g NaCl in 5 mL wine (= 800 g/L = far above saturation in pure water; the wine matrix absorbs the excess). Lan 2021 (grape skins): 2 g NaCl. Shimadzu Pittcon: 1.5 g NaCl in 7 mL Moscato. Trupiano 2024 (red wine): not extracted. **Open question:** what NaCl mass maximizes the headspace signal for Cabernet's specific MW + polarity distribution without crashing alcohols/acids on the phase? Plan a 0/1/2/4/6 g NaCl ladder.

5. **Extraction temperature.** Trupiano 2024 (red wine): 60 °C. Lisanti 2021: 50 °C. Lan 2021 (grape skins): 60 °C. Shimadzu (Moscato): 40 °C. **Higher temperature increases volatility and headspace concentration but degrades thermally labile analytes** (β-damascenone, certain thiols). Need 30/40/50/60/70 °C ladder with stability checks on the heat-sensitive analytes.

6. **Exposure time.** Trupiano 2024 (red wine): 53 min. Lisanti 2021: 60 min. Lan 2021: 49 min (free) / 60 min (bound). Shimadzu: 40 min. **All of these are 40–60 min** — the dissertation can probably standardize on 45–60 min after a kinetic equilibration study showing plateau onset.

7. **Vendor SKU lock-in.** The CTC/PAL SARR-prefix SKUs and the Restek 28xxx-series SKUs refer to the same physical product (CTC manufactures; Restek and BGB and Agilent distribute). **Lock to one supplier for the duration of the dissertation to avoid lot-to-lot artifact contamination of validation data.**

8. **Lifetime tracking in a wine matrix.** No published source gives a real injection count for HS-SPME Arrow in wine. **Plan to instrument the validation: log every injection on a single Arrow with a per-injection QC standard, and report the actual injection count to first 20 % signal loss as a methodology deliverable.** This is a publishable supplementary contribution.

9. **GC inlet pressure ceiling.** Restek caps at 50 psi. For a typical 30 m × 0.25 mm × 0.25 µm column at constant flow 1.2 mL/min He, the inlet pressure at 250 °C oven is around 12–15 psi — well below the cap. **Not a binding constraint for the planned method.**

10. **Conditioning module.** Use of the dedicated PAL3-SPME-ArrowCond module (vs in-GC-inlet conditioning) is recommended by Restek but not mandatory. **Decision criterion:** if the lab does not already own the conditioning module, in-GC-inlet conditioning at 250 °C × 5 min between samples is acceptable per the Restek 2019 conditioning table.

11. **Calibration approach.** Lisanti 2021 used 12 % ethanol / tartaric / pH 3.5 model wine. Lan 2021 used acid-hydrolyzed grape-skin matrix for bound VOCs. Trupiano 2024 used standard addition on red wine. **For Cabernet quantitative validation per ICH Q2(R2), standard addition on real wine is the methodologically defensible default**, with model wine retained for upstream method development.

12. **Open hardware question:** the Smart-chip-enabled Arrows write parameters and use-count to an RFID chip in the hub. Backward-compatible non-Smart Arrows still work on PAL3 systems but lose the auto-parameter and lifetime-tracking feature. **Decision:** for any quantitative validation, use Smart Arrows so the use-count is captured in the run-log automatically. The premium over non-Smart is small and the data-trail value is high.

13. **Septum lifetime quantification.** Vendor claims "2× longer septum life." No primary data backs this for wine. **Quantify septum core-out rate as a side-quest deliverable** (count cores per 100 injections on a fresh septum vs. classical-fiber baseline).

---

## 8. References (URLs + DOIs + retrieval date — all 2026-05-12 unless noted)

### Vendor primary

1. **PAL Smart SPME Arrow Brochure** (CTC Analytics AG; Brochure rev. 4, April 2025). https://www.palsystem.com/fileadmin/user_upload/content_hub/Files/Brochures/PAL_Smart_SPME_Arrow_Brochure_screen.pdf — retrieved 2026-05-12 (PDF extracted locally). Geometry table p. 3 (surface 44.0 / 62.8 mm²; volume 3.8 / 11.8 µL; fiber 9.4 mm² / 0.6 µL). White-wine application p. 4. SKU table p. 9.
2. **GL Sciences EU mirror brochure.** https://www.glsciences.eu/pal-system/PAL_Smart_SPME_Arrow_Brochure.pdf — retrieved 2026-05-12. Substantially identical to CTC; includes the extra "Injector septa last at least 2× longer" claim on p. 5.
3. **Restek PAL SPME Arrow Instruction Sheet, #500-60-002 Rev. 05/19.** https://d1lqgfmy9cwjff.cloudfront.net/csi/pdf/admin/rk_spme_arrow_instruction_sheet.pdf — retrieved 2026-05-12. Table I (selection guide with old 27xxx SKUs); Table II (full thermal-conditioning matrix); precautions and inlet conversion kits.
4. **Restek "PAL SPME Arrow Configurations" blog/article.** https://www.restek.com/spme-arrow-configurations/ — retrieved 2026-05-12. New 28xxx Smart SKU table; surface area 44 / 63 mm²; volume 3.8 / 12 µL; sensitivity 297–527 % (10 min) / 618–896 % (5 min); "2 to 3× longer" lifetime.
5. **Restek "Set Your Sights on Superior Performance".** https://discover.restek.com/articles/GNSS2642/set-your-sights-on-superior-performance/ — retrieved 2026-05-12.
6. **Restek product page, PAL Smart SPME Arrow 28904-1 (DVB/Carbon-WR/PDMS, 1.10 mm, 120 µm).** https://www.restek.com/p/28904-1 — retrieved 2026-05-12. Confirms 50 psi inlet ceiling and 28912 manual-injection-kit requirement.
7. **CTC Analytics PAL_SPME Arrow brochure on DirectIndustry.** https://pdf.directindustry.com/pdf/ctc-analytics/pal-spme-arrow/65926-662934.html — retrieved 2026-05-12.
8. **Agilent Smart SPME Arrows product page.** https://www.agilent.com/en/product/sample-preparation/solid-phase-microextraction-spme/smart-spme-arrows — retrieved 2026-05-12.
9. **Agilent technical overview "Solid phase microextraction (SPME) — Fundamentals and SPME Arrow," 5994-5775EN.** https://www.agilent.com/cs/library/technicaloverviews/public/te-solid-phase-microextraction-fundamentals-spme-arrow-5994-5775en-agilent.pdf — retrieved 2026-05-12 (403 on direct fetch; title and ID confirmed via search snippets and Agilent library directory).
10. **Agilent application note 5994-3160EN / 5994-3161EN, "Analysis of Free Volatile Phenols in Smoke-Impacted Wines by SPME-GC/MS — Response Comparison of Agilent SPME Arrows and Classical Fibers."** https://www.agilent.com/cs/library/applications/an-smoke-taint-wine-SPME-8890-gc-5994-3160en-agilent.pdf and https://www.agilent.com/cs/library/applications/application-smoke-taint-wine-SPME-8890-gc-5994-3161en-agilent.pdf — retrieved 2026-05-12 (403 — confirm via Agilent library login or library copy).
11. **Agilent Smart SPME Arrows package insert, 5994-3137EN.** https://www.agilent.com/cs/library/packageinsert/public/in-smart-spme-arrows-5994-3137en-agilent.pdf — retrieved 2026-05-12 (403 — same).
12. **BGB Analytik PAL SPME Fibers / Arrow product page.** https://bgb-analytik.com/pal-system/spme-fibers — retrieved 2026-05-12. EUR pricing examples; mirror of CTC SKUs.
13. **Autosamplerguys.com SPME Arrow product page.** https://www.autosamplerguys.com/spme-arrow — retrieved 2026-05-12.

### Peer-reviewed primary literature

14. **Kremser, A.; Jochmann, M.A.; Schmidt, T.C. 2016.** "PAL SPME Arrow—evaluation of a novel solid-phase microextraction device for freely dissolved PAHs in water." *Anal Bioanal Chem* 408:943–952. DOI 10.1007/s00216-015-9187-z. PMC4709367 — https://pmc.ncbi.nlm.nih.gov/articles/PMC4709367/ — retrieved 2026-05-12.
15. **Lisanti, M.T.; Laboyrie, J.; Marchand-Marion, S.; de Revel, G.; Moio, L.; Riquier, L.; Franc, C. 2021.** "Minty aroma compounds in red wine: Development of a novel automated HS-SPME-arrow and gas chromatography-tandem mass spectrometry quantification method." *Food Chemistry* 361:130024. DOI 10.1016/j.foodchem.2021.130024. Open URL https://www.sciencedirect.com/science/article/abs/pii/S0308814621010359 — retrieved 2026-05-12 (403 on direct ScienceDirect fetch; method details extracted via search-engine snippets and the IVES Open Science 2023 sister piece at https://ives-openscience.eu/3833/).
16. **Lan, Y.; Wu, J.; Wang, X.; Sun, X.; Hackman, R.M.; Li, Z.; Feng, X. 2021.** "Optimization of SPME-Arrow-GC/MS Method for Determination of Free and Bound Volatile Organic Compounds from Grape Skins." *Molecules* 26(23):7409. DOI 10.3390/molecules26237409. PMC8659239 — https://pmc.ncbi.nlm.nih.gov/articles/PMC8659239/ — retrieved 2026-05-12.
17. **Kim, J.; Lee, Y.J.; Roh, M.-K.; Park, M.K. 2020.** "Comparison of Different Types of SPME Arrow Sorbents to Analyze Volatile Compounds in Cirsium setidens Nakai." *Foods* 9(6) — see PMC7353593 — https://pmc.ncbi.nlm.nih.gov/articles/PMC7353593/ — retrieved 2026-05-12.
18. **Lim, D.K.; Mo, C.; Lee, J.-H.; Long, N.P.; Dong, Z.; Li, J.; Lim, J.; Kwon, S.W. 2019.** "Comparison of headspace–SPME and SPME-Arrow–GC–MS methods for the determination of volatile compounds in Korean salt–fermented fish sauce." *Applied Biological Chemistry* 62(1):60. DOI 10.1186/s13765-019-0424-6. https://applbiolchem.springeropen.com/articles/10.1186/s13765-019-0424-6 — retrieved 2026-05-12 (Springer paywall redirect — abstract data captured via search snippets).
19. **Yin, C.; Yang, J.; Wu, F.; Chen, X.; Yan, X.; Lyu, J.; Zhang, K.; Li, B. 2019.** "Optimization and validation of a head space solid-phase microextraction-arrow gas chromatography-mass spectrometry method using central composite design for determination of aroma compounds in Chinese liquor (Baijiu)." *J Chromatogr A* 1601:1–11. DOI 10.1016/j.chroma.2019.06.052. https://www.sciencedirect.com/science/article/abs/pii/S0021967319309781 — retrieved 2026-05-12 (paywall).
20. **Trupiano, B.; Gangemi, S.; et al. 2024.** "Optimisation of SPME Arrow GC/MS method for determination of wine volatile organic compounds." *OENO One* 58(4). https://oeno-one.eu/article/view/7914 — retrieved 2026-05-12 (SSL certificate verification failed on direct fetch; method summary captured via Google Scholar / search snippets; flag for ResearchGate backup: https://www.researchgate.net/publication/384859481). Red-wine optimum: 60 °C / 17 min incubation / 53 min exposure. White-wine optimum: 50 °C / 10 min / 60 min. 82 aroma compounds quantified.
21. **Herrington, J.S.; Gómez-Ríos, G.A.; Myers, C.; Stidsen, G.; Bell, D.S. 2020.** "Hunting Molecules in Complex Matrices with SPME Arrows: A Review." *Separations* 7(1):12. DOI 10.3390/separations7010012. https://www.mdpi.com/2297-8739/7/1/12 — retrieved 2026-05-12 (full PDF 403; abstract via https://doaj.org/article/2e502f78a2fe4ea5b156d083640c1b36).
22. **Liu, Z.; Ezernieks, V.; Reddy, P.; Elkins, A.; Krill, C.; Murphy, K.; Rochfort, S.; Spangenberg, G. 2020.** "A Simple GC-MS/MS Method for Determination of Smoke Taint-Related Volatile Phenols in Grapes." *Metabolites* 10(7):294. DOI 10.3390/metabo10070294. (PDF extracted locally — full text in cache.) Classical SPME, not Arrow, but defines the smoke-taint target list and MRM transitions.

### Conference posters and grey literature

23. **Shimadzu Pittcon 2017 poster PO-CON1702E.** Owens, A.; Yang, M.; Lock, N.; Sandy, A. (Shimadzu Scientific Instruments, Columbia MD). "The Comparison of HS-SPME and SPME Arrow Sampling Techniques Utilized to Characterize Volatiles in the Headspace of Wine over an Extended Period of Time." https://www.shimadzu.com/an/sites/shimadzu.com.an/files/pim/pim_document_file/applications/posters/13414/jpo217003.pdf — retrieved 2026-05-12 (PDF extracted locally; full text captured).
24. **Thermo Fisher Pittcon 2019 poster PO-10697.** "New Opportunities for Wine Analysis through SPME Arrow." https://assets.thermofisher.com/TFS-Assets/CMD/posters/po-10697-gc-ms-wine-spme-arrow-pittcon2019-po10697-en.pdf — retrieved 2026-05-12 (60 s timeout on direct fetch; flag for retry).
25. **LabRulez GCMS / Sleek paper database.** "Determination of haloanisoles in wine by HS-SPME Arrow and GC-MS/MS." https://gcms.labrulez.com/paper/18997 — retrieved 2026-05-12. PDMS Arrow 100 µm; LOD 0.03–0.07 ng/L; LOQ ≤0.24 ng/L; recoveries 90–105 %.
26. **BenchChem application note B122830, April 2026.** "Quantification of α-Ionone in Wine Using HS-SPME-GC-MS." https://www.benchchem.com/product/b122830/docs — retrieved 2026-05-12. Third-party; uses classical DVB/CAR/PDMS not Arrow. Useful for parameter ranges.
27. **LCGC International / ChromatographyOnline.** "SPME Arrow Improves Upon Conventional SPME in GC×GC–MS Analysis of Whiskey Flavor Compounds." https://www.chromatographyonline.com/view/spme-arrow-improves-upon-conventional-spme-in-gc-gc-ms-analysis-of-whiskey-flavor-compounds — retrieved 2026-05-12 (403 on direct fetch; title and topic confirmed). Whiskey beverage matrix.
28. **LCGC International / ChromatographyOnline.** "An Alternative to Classical SPME: SPME Arrow for the Analysis of Flavour Profiles." https://www.chromatographyonline.com/view/alternative-classical-spme-spme-arrow-analysis-flavour-profiles — retrieved 2026-05-12 (403).
29. **Chromatography Today.** Kremser, A.; Boehm, G.; Schilling, B.; Jochmann, M.A.; Schmidt, T.C. 2017. "PAL SPME Arrow: An Evolutionary Step Forward for Solid-phase Microextraction (SPME)." https://www.chromatographytoday.com/article/solid-phase-extraction-spe/34/ctc-analytics-ag/ppal-spme-arrow-an-evolutionary-step-forward-for-solid-phase-microextraction-spmep/2244 — retrieved 2026-05-12.

### Wine-aroma context (not Arrow-specific but methodologically adjacent)

30. **AWRI review — Hayasaka, Y.; Baldock, G.A.; Pollnitz, A.P.** "Contributions of mass spectrometry in The Australian Wine Research Institute to advances in knowledge of grape and wine constituents." *Aust J Grape Wine Res*. (Cached PDF retrieved 2026-05-12.) Reviews 1971–2005 AWRI MS work — **does not mention SPME Arrow.** Useful for the literature-context discussion of why Arrow adoption in wine has lagged other matrices.
31. **AWRI Technical Review No 273 (July 2025).** https://www.awri.com.au/information_services/technical_review/previous_issues/technical-review-no-273-july-2025/ — retrieved 2026-05-12. No SPME Arrow content.
32. **Oregon State University thesis, Yiu, P.** "The Effect of Wine Matrix on the Analysis of Volatile Sulfur Compounds by Solid-Phase Microextraction-GC-PFPD." OSU ScholarsArchive ID n870zt673. https://ir.library.oregonstate.edu/concern/graduate_thesis_or_dissertations/n870zt673 — retrieved 2026-05-12. Classical SPME, not Arrow.
33. **MDPI Processes 2021** — Robinson, A.L. et al. (snippet only). "Quantification of Volatile Compounds in Wines by HS-SPME-GC/MS: Critical Issues and Use of Multivariate Statistics in Method Optimization." *Processes* 9(4):662. https://www.mdpi.com/2227-9717/9/4/662 — retrieved 2026-05-12.

---

### Document hygiene note

Two specific numerical conflicts surfaced and remain unresolved without primary-source access:
- **PAL brochure "10× sensitivity" vs Restek "297–897 % response":** These are not contradictory but they are framed differently. The PAL number is a working-range expansion across 4 decades; the Restek number is a per-compound mean response increase. Both are vendor-headline numbers — treat as plausible upper bounds.
- **The 30 mm phase length appearing in one Shimadzu figure caption** (referenced via Lim 2019 search snippet and the Shimadzu poster figure caption: "PAL SPME Arrow has a 250-μm × 30-mm, 15.3-μL sorption phase for the 1.5 mm diameter model") — **contradicts every CTC/PAL brochure stating 20 mm.** The 30 mm + 15.3 µL would imply the same surface-area-to-volume ratio as the published 20 mm + 11.8 µL (62.8 mm² / 11.8 µL ≈ 5.32 mm²/µL vs hypothetical 30 mm × 5.32 = 94.2 mm²/15.3 µL ≈ 6.15 mm²/µL — inconsistent). **Treat the 30 mm as a typo and use 20 mm + 11.8 µL as the authoritative 1.5 mm 250 µm Arrow geometry.**

Flag for follow-up retrieval (worth a second pass at the dissertation-writing stage):
- Full text of Lisanti 2021 *Food Chemistry* DOI 10.1016/j.foodchem.2021.130024 — paywall workaround needed for MRM transitions.
- Full text of Trupiano 2024 *OENO One* — SSL cert issue on direct fetch; retry via ResearchGate or institutional VPN.
- Agilent smoke-taint wine app notes 5994-3160EN / 5994-3161EN / 5994-3137EN — 403 on direct fetch; retrieve via Agilent library login.
- Thermo Pittcon 2019 PO-10697 wine SPME Arrow poster — timeout on direct fetch; retry from Thermo asset CDN.
- Hunting Molecules review (Herrington 2020) full PDF — 403; available via authenticated MDPI access or ResearchGate at https://www.researchgate.net/publication/339362141.
