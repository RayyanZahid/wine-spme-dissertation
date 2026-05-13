# Wine Faults, Taint, and Application — Brief

**Specialist 6/6** — application chapter source material for *Method Development, MRM Optimization, and Validation of SPME Arrow Headspace GC-MS/MS for the Quantitative Profiling of Aroma-Active Volatiles in Cabernet Sauvignon Wine*.

This chapter answers: once the validated SPME Arrow HS-GC-MS/MS method exists, what does it diagnose? Five fault classes plus the differentiation endpoint (appellation / vintage / oak program). All concentration ranges and thresholds below come from peer-reviewed literature, AWRI fact sheets, and AJEV studies surveyed 2026-05-12.

---

## 1. Smoke taint (volatile phenols and their glycosides)

Wildfire smoke — and to a lesser extent prescribed burns and stubble burns — pyrolyses lignin and releases free volatile phenols (guaiacol, the methyl- and ethyl- guaiacols, syringol, 4-methylsyringol, phenol, the three cresols) that are absorbed through grape skins, then glycosylated by the berry into largely odourless rutinoside / gentiobioside / pentose-glucoside conjugates. After fermentation and over months-to-years of bottle age the acidic wine matrix slowly hydrolyses these glycosides back to free phenols, giving smoke-affected wines an "ashtray / cold smoke / band-aid" finish that intensifies with age. The application chapter must quantify **both** the free phenols and a basket of glycoside markers to predict eventual sensory load. (Pollnitz/Sefton/Hayasaka AWRI lineage 2004-2013 established the canonical analyte set.)

### Free volatile phenols — descriptor / sensory threshold / range in clean vs smoke-affected Cab Sauv

| Compound | Descriptor | Threshold (red wine) | Clean Cab Sauv (typical) | Smoke-affected Cab Sauv (typical) |
|---|---|---|---|---|
| Guaiacol | Smoky, medicinal | 23 µg/L (red), 20 µg/L (white) | <5 µg/L (99th-pct background); 2-24 µg/L median in *oaked* commercial wine, max ~47 µg/L | 6 µg/L is the smoke-flavour off-flavour detection level for most tasters; affected wines commonly 10-50 µg/L, can exceed 100 µg/L in heavily exposed fruit |
| 4-Methylguaiacol | Smoky, spicy | ~65 µg/L | ≤5 µg/L typical; oaked wine max ~35 µg/L | up to 25 µg/L in smoke-affected, can climb higher with severe exposure |
| Syringol | Smoky, woody | ~570 µg/L (weak direct sensory, but reliable exposure marker) | Oaked Shiraz: median 3-47 µg/L, max 187 µg/L; Cab Sauv lower | Median 6-12 µg/L, max ~65 µg/L in affected wines |
| 4-Methylsyringol | Smoky | high (~570 µg/L equivalent) | Oaked Shiraz median up to 47 µg/L, max 96 µg/L; lower in Cab Sauv | Median <1-4 µg/L, max ~25 µg/L |
| o-Cresol | Tar, medicinal | ~62 µg/L (red wine, mixed reports) | <5 µg/L | >3 µg/L is suggestive of smoke exposure; affected wines 5-25 µg/L |
| m-Cresol | Tar, leather | ~20 µg/L | <5 µg/L | 3-15 µg/L typical in affected wines |
| p-Cresol | Tar, barnyard | ~64 µg/L | <5 µg/L | 3-15 µg/L; can be elevated by *Brett* independently of smoke, important confounder |
| Phenol | Smoky, tarry | ~ high | trace | a few µg/L in affected wines |

Aroma detection thresholds compiled in Parker et al. and reproduced widely; ResearchGate compilation of guaiacol/4-MG/4-EG/4-EP thresholds is the most-cited single table.

**Style-dependent thresholds.** Guaiacol off-flavour is flagged in delicate sparkling whites at 6-10 µg/L, medium-bodied reds at 15-25 µg/L, and full-bodied Shiraz at 30-40 µg/L. Cabernet Sauvignon sits at the upper end of this red-wine range because tannin and oak match-mask some guaiacol character. Guaiacol's *glucoside* has a threshold around 69 µg/L in water (≈10× free guaiacol).

### Glycoside-bound markers — the "total smoke" basket

AWRI / Hayasaka HPLC-MS/MS measures the bound forms directly; an alternative is direct acid hydrolysis followed by GC-MS quantitation of liberated free phenols (Singh et al. SPE-based adaptation). The two approaches do **not** produce numerically comparable values — the AWRI background database is glycoside-specific and HPLC-MS/MS-anchored.

| Glycoside marker | Clean Cab Sauv (oaked baseline) | Smoke-affected Cab Sauv |
|---|---|---|
| Syringol gentiobioside (SyGG) | Oaked Shiraz median 8 µg/L, max 18 µg/L; Cab Sauv similar low end | Commonly 13-123 µg/L, max 690 µg/L; the single most diagnostic marker |
| 4-Methylsyringol gentiobioside | Mean ~20 µg/L in just-affected wine sets the bottom of the affected range | 20-200 µg/L in fully-exposed Cab Sauv |
| Guaiacol rutinoside (GuRG) | Oaked Shiraz median 7 µg/L, max 20 µg/L | Range 11-85 µg/L in affected wines |
| 4-Methylguaiacol rutinoside | <20 µg/L oaked baseline (all cultivars) | tens of µg/L when smoke-exposed |
| Cresol rutinoside | <20 µg/L baseline | tens of µg/L when exposed |
| Phenol rutinoside | <20 µg/L baseline | tens of µg/L when exposed |

**AWRI risk thresholds (grape µg/kg, not wine — but cited because the dissertation method has to read across).** Wines made from grapes with guaiacol >7-10 µg/kg or SyGG >100-160 µg/kg fall into the AWRI "high risk" basket for perceivable smoke flavour. Extraction efficiency from grape to wine for the bound forms is ~67% for Cab Sauv (Caffrey et al. 2019), so a Cab Sauv smoke-risk model should multiply grape glycoside by ~0.67 to forecast wine load.

### Canonical AWRI / UC Davis literature

- Pollnitz, Pardon, Sykes & Sefton (AWRI, 2004) — stable-isotope-dilution GC-MS of guaiacol and 4-methylguaiacol in wine and oak; the canonical analytical paper for free phenols, foundation for all later SPME-GC-MS/MS work.
- Hayasaka et al. (AWRI, 2010 / 2013) — discovery of glycosylation of smoke-derived volatile phenols in grapes (J. Agric. Food Chem. 2010 and follow-ons); HPLC-MS/MS measurement of glycosides.
- Krstic, Johnson & Herderich (AJGWR review, 2015) — comprehensive review of smoke-derived volatile phenols and their glycosidic metabolites as smoke-exposure biomarkers.
- Wilkinson (Adelaide group) — passive sampling, machine-learning prediction, ozonation amelioration, stubble-burn studies (multiple papers 2013-2024).
- Caffrey, Lerno, Rumbaugh, Girardello, Zweigenbaum, Oberholster, Ebeler (UC Davis, AJEV 2019) — *Changes in Smoke-Taint Volatile-Phenol Glycosides in Wildfire Smoke-Exposed Cabernet Sauvignon Grapes throughout Winemaking* — landmark Napa 2017-vintage Cab Sauv study.
- Noestheden, Thiessen, Dennis, Tiet, Zandberg (2024 *AJEV*) — *Prevalence of Wildfire Smoke Exposure Markers in Oaked Commercial Wine* — establishes the oak vs smoke baseline that the Cab Sauv differentiation chapter must respect.
- Fudge / Boss / Wilkinson (Adelaide) — predictive ML on volatile phenol profiles for smoke classification.

---

## 2. Brettanomyces / "Brett character"

*Brettanomyces bruxellensis / Dekkera bruxellensis* decarboxylates hydroxycinnamic acids (p-coumaric → 4-vinylphenol → 4-ethylphenol; ferulic → 4-vinylguaiacol → 4-ethylguaiacol; caffeic → 4-vinylcatechol → 4-ethylcatechol). The ratio 4-EP:4-EG ≈ 8:1 to 10:1 is the *Brett* signature in red wine and the easiest discriminator from oak-derived volatile phenols (which carry a lower 4-EP-to-guaiacol ratio).

| Compound | Descriptor | Threshold in red wine | Typical clean red Cab Sauv | Typical Brett-affected red Cab Sauv |
|---|---|---|---|---|
| 4-Ethylphenol (4-EP) | Horse, sweaty leather, band-aid, medicinal | 605 µg/L alone; 369 µg/L when paired with 37 µg/L 4-EG (Chatonnet 1992) | <20 µg/L | 400-3000+ µg/L; some heavily affected reds 5000-6000 µg/L |
| 4-Ethylguaiacol (4-EG) | Smoky, clove, spicy, bacon | 140 µg/L alone; 110 µg/L in mixture | <5 µg/L | 50-500 µg/L typically |
| 4-Ethylcatechol (4-EC) | Horsey, animal | 60-400 µg/L (low study), 774 µg/L (Aust. neutral red), Merlot DT 823 µg/L / consumer rejection 1323 µg/L | <100 µg/L | 100-1000 µg/L in heavily affected wines; similar magnitude to 4-EG in European Cab Sauv / Pinot survey |
| 4-Vinylphenol / 4-Vinylguaiacol | Pharmaceutical (precursors, low in reds) | ~770 / ~440 µg/L | <40 µg/L | rarely elevated in reds (decarboxylated to ethyl forms) |

**Sum-of-evidence rule (Chatonnet).** When 4-EP and 4-EG are simultaneously present in a 10:1 ratio at >400 µg/L 4-EP + ~40 µg/L 4-EG, "Brett character" is overwhelmingly likely; the joint perception threshold collapses below the individual threshold because of additivity. Some Cab Sauv styles (e.g., older Australian) historically *welcome* this character at low levels — the application chapter should treat 4-EP and 4-EG as *style markers*, not as binary faults, and report both raw µg/L and the 4-EP:4-EG ratio.

Statistical re-evaluation by Pinto et al. of the 4-EP/4-EG threshold proposed lowering the working Brett-detection threshold below 605 µg/L in modern wine styles.

**Canonical literature.**
- Chatonnet, Dubourdieu, Boidron & Pons (J. Sci. Food Agric., 1992) — *The origin of ethylphenols in wines* — the founding paper, perception thresholds for 4-EP, 4-EG, and their 10:1 mixture.
- Curtin / Coulter / Cowey AWRI Brettanomyces FAQ — operational thresholds and intervention.
- Romano et al. — Brazilian survey correlating Brett with 4-EP.
- Pollnitz, Pardon & Sefton — analytical SIDA GC-MS.

---

## 3. Cork taint and "musty / earthy" taints (haloanisoles, geosmin, MIB)

| Compound | Descriptor | Threshold | Range in wine |
|---|---|---|---|
| 2,4,6-Trichloroanisole (TCA) | Mouldy, wet cardboard, damp basement | 1.4-4 ng/L (water/wine, trained); olfactometry 4-10 ng/L white wine, 50 ng/L red wine; consumer rejection 3.1 ng/L, detection 2.1 ng/L | Tainted bottles: 4-50 ng/L; SPME-GC-MS LOD ~0.3 ng/L; HS-SPME GC-MS/MS routine sub-ng/L; SPME Arrow + MS/MS competitive |
| 2,4,6-Tribromoanisole (TBA) | Musty, mushroom, mouldy | 0.5-3.4 ng/L (white wine, 4 ng/L red); even lower than TCA in some reports | Tainted bottles from cellar/cellar-paint contamination; ubiquitous fire-retardant phenol precursor (>>cork-only TCA in modern problem incidents) |
| 2,3,4,6-Tetrachloroanisole (TeCA) | Musty | 4-15 ng/L | Lower prevalence than TCA in cork; bleach + chlorophenol building source |
| Pentachloroanisole (PCA) | Musty | 4000 ng/L (high — rarely sensory limiting) | Higher levels possible; sub-sensory in most tainted wines |
| Geosmin | Earthy, beetroot, damp soil | Wine 50-90 ng/L (water 1-10 ng/L individually variable) | Affected wines tens to hundreds of ng/L; sometimes Botrytis-co-occurring |
| 2-Methylisoborneol (MIB) | Musty, earthy, camphor | Wine 30-55 ng/L (water 15 ng/L) | Tens of ng/L in affected wines |

**Method note.** HS-SPME (and SPME Arrow) GC-MS/MS in MRM mode achieves sub-ng/L LOQ for haloanisoles in wine — the Arrow's larger sorbent volume vs classical SPME fibres particularly suits this low-µg analyte class. Published methods reach 0.05 ng/L LOQ for TCA. The dissertation method should include TCA, TBA, geosmin, MIB in the same MRM acquisition window if the matrix permits — cork taint is a one-bottle decision and the same Arrow extract supports it.

**Canonical literature.**
- Pollnitz, Pardon, Liacopoulos, Skouroumounis & Sefton (AWRI, 1996) — *The analysis of 2,4,6-trichloroanisole and other chloroanisoles in tainted wines and corks* — Aust. J. Grape Wine Res. — founding paper.
- Chatonnet et al. (2004) — identification of TBA as a major contributor to musty taint independent of cork.
- Cravero (2015) *J. Inst. Brew.* — sensory evaluation of TCA in wines.
- Hjelmeland & Ebeler (UC Davis) — geosmin in wines.

---

## 4. Reductive faults (volatile sulfur compounds)

H2S and the low-MW mercaptans push wines through a four-stage progression (egg → onion → rubber → garlic) as the SH-bearing precursor pool oxidises to disulfides during ageing under low O2 / under screw cap. Cu addition removes mercaptans (and H2S) but *creates* a dissolved copper-sulfide reservoir that can re-release reductive compounds under further reductive conditions.

| Compound | Descriptor | Threshold (µg/L, wine) | Clean Cab Sauv | Reductive-fault Cab Sauv |
|---|---|---|---|---|
| Hydrogen sulfide (H2S) | Rotten egg, sewer | 1.1-1.6 (AWRI 1-2) | <1 | 5-50, rotten-egg perception above ~2 |
| Methanethiol (MeSH) | Cooked cabbage, putrefaction, rubber | 1.8-3.1 | <2 | 5-30 in cabbage-character wines |
| Ethanethiol (EtSH) | Onion, rubber, faecal, earthy | 1.1 | <0.5 | 2-15 |
| Dimethyl sulfide (DMS) | Asparagus, corn, blackcurrant (low), olive/cabbage (high) | 25 (model wine) — actual red-wine perception 25-60; Cab Sauv panels report 27-60 | 1-30 in young Cab Sauv | 42-910 in aged Australian Cab Sauv (Anabella Geneva data); high end is faulty |
| Diethyl sulfide | Cooked vegetable, garlic, onion | 0.92 | <1 | 1-10 |
| Dimethyl disulfide (DMDS) | Cooked cabbage, intense onion | 29 | <5 | 10-50 in faulty wines |
| Diethyl disulfide (DEDS) | Garlic, burnt rubber | 4.3 | <2 | 5-30 |
| Methionol (3-methylthio-1-propanol) | Cauliflower, cabbage, potato | 500 (in wine) | 1000-3000 typical in red wine; *suppresses fruitiness* below threshold | >3000 in some reductive wines |
| Methional | Boiled potato, cooked meat | 0.5-5 (very low; ~2 µg/L common) | <1 | Both a *reductive precursor* and a Strecker-aldehyde *oxidation* product; appears in both fault tables |

**Dual-role compounds.** DMS at 1-30 µg/L is an aroma *enhancer* in aged red Cab Sauv (lifts blackcurrant, blackberry); at 50-100 µg/L it tips into olive/truffle/cabbage, and >100 µg/L is generally faulty. The dissertation should treat DMS as a *style* descriptor across most of the analytical range and an *over-evolution / closure transfer* marker at the high end.

**Canonical literature.**
- Siebert, Solomon, Pollnitz & Jeffery (AWRI) — selective GC-SCD (sulfur chemiluminescence) determination of VSCs in wine; AWRI standard.
- Goode & Jamieson reductive-fault primer; Kreitman / Waterhouse copper-sulfide chemistry.
- Ugliano et al. — Cu, post-bottling reduction, screw cap.
- Franco-Luesma et al. — DMS in Cab Sauv from grape S-methylmethionine precursor.

---

## 5. Oxidative faults

| Compound | Descriptor | Threshold in wine | Clean red Cab Sauv | Oxidised Cab Sauv |
|---|---|---|---|---|
| Acetaldehyde | Bruised apple, sherry, nutty, "stuck ferment" | 100-125 mg/L sensory in red; 0.5 mg/L is a "minor descriptor" lower bound, not perception of fault | 30-75 mg/L typical young red | >125 mg/L in oxidised wines; sherry ≈300 mg/L (intentional) |
| Sotolon | Curry, fenugreek, walnut, *rancio* | 7-23 µg/L in fortified wines; ~112 µg/L in Madeira with endogenous 6 µg/L (matrix-driven); ~5-10 µg/L typical perception | <2 µg/L in young Cab Sauv | 10-100 µg/L in oxidised/premox dry wines; 2-516 µg/L across Madeira blends |
| Methional | Boiled potato, cooked-meat, masking fruit | ~0.5-5 µg/L | <1 µg/L young | Accumulates in oxidised reds well before acetaldehyde; useful early marker |
| Phenylacetaldehyde | Honey, rose, soapy | ~1-5 µg/L | <2 µg/L | 10-30+ µg/L in oxidised wines |
| 3-Methylbutanal | Malty, cocoa, fermented | ~5 µg/L | <5 µg/L | Tens of µg/L in oxidised wines |
| 2-Methylpropanal | Banana, malty | similar low single-digit µg/L | <5 µg/L | Tens of µg/L in oxidised wines |
| Furfural | Almond, caramel | 4-15 µg/L LOQ; far below sensory in red wine but caramelisation indicator | 1-15 µg/L from oak | Elevated with thermal stress |
| 5-Methylfurfural | Caramel, butter | 5-25 µg/L | <5 µg/L | Elevated with thermal/oxidative stress |
| 5-Hydroxymethylfurfural | Burnt sugar | ~10 µg/L LOQ | Trace | Increases under thermal abuse |

**Strecker aldehydes accumulate before acetaldehyde** — measuring methional and phenylacetaldehyde provides an *early* oxidative-stress flag in dry red Cab Sauv, well before the wine smells of bruised apple. The "aged red ARP" index (Ferreira, Bueno et al.) ties acetaldehyde accumulation potential to polyphenol reactivity.

**Canonical literature.**
- Bueno, Carrascón & Ferreira (J. Agric. Food Chem. 2018) — *Formation and Accumulation of Acetaldehyde and Strecker Aldehydes during Red Wine Oxidation*.
- Bueno, Marrufo-Curtido, Carrascón, Fernández-Zurbano, Escudero, Ferreira (Front. Chem. / J. Agric. Food Chem. 2015) — *Quantitative Analysis by GC-MS/MS of 18 Aroma Compounds Related to Oxidative Off-Flavor in Wines*. This is the canonical method paper for the application chapter — 18 oxidation-related compounds in one GC-MS/MS MRM run.
- Pereira et al. — Madeira sotolon evolution.
- Zea et al. (2015) — *Acetaldehyde as Key Compound for the Authenticity of Sherry Wines: A Study Covering 5 Decades*.

---

## 6. Application endpoint — appellation, vintage, oak

Once the SPME Arrow HS-GC-MS/MS method is validated for ~50-80 aroma-active analytes (esters, higher alcohols, monoterpenes, C13-norisoprenoids, volatile fatty acids, volatile phenols, methoxypyrazines, lactones, vanillin/eugenol, plus fault compounds 1-5 above), the natural application chapter is **multivariate discrimination of Cabernet Sauvignon by appellation × vintage × oak program**.

### Published differentiation studies on Cab Sauv

| Study | Regions / vintages | Method | Discriminant compounds | Performance |
|---|---|---|---|---|
| Niu et al. (PMC9801081) | Ningxia (China), California, Bordeaux | HS-SPME/GC-MS + PCA, OPLS-DA | 12 compounds with significant inter-region differences; citronellol + acetoin unique to Ningxia; 12 characteristic phenolics separate Ningxia from Bordeaux/CA | Region > variety as a driver of diversity; clear cluster separation |
| Pan et al. (Foods 2022, PMC8750599) | Helan Mountains 2013-2018 | HS-SPME-GC-MS + PCA, HCA | 17 VIP>1 markers in OPLS-DA across regions; vintage effects subordinate to region | OPLS-DA + HCA combined separates 6 sub-regions |
| Hu et al. (PMC8398669) | Cab Sauv smoke-tainted vs control | E-nose + ANN regression vs GC-MS volatiles | 8 volatile aromatic compounds; e-nose ANN R=0.97 (sensory) and 0.99 (chemical) | Demonstrates ML/sensory cross-validation for smoke-affected Cab Sauv |
| Souza Gonzaga et al. (Food Res. Int. 2020 / IVES) | Coonawarra, Margaret River, Yarra Valley, Bordeaux — 52 Cab Sauv 2015 commercial wines | Quantitative chemical panel + trained DA + multivariate stats | Coonawarra: mint, mallee leaf, dark fruit; Margaret River: floral, green pepper, violets, red fruit; Yarra: cool-climate green; Bordeaux: oak/spice axis | Regional chemical fingerprint validated against sensory typicity |
| Casassa et al. (2024 / 2026) | Mendoza Cab Sauv — geographical origin × vintage × plant material across 2018/2019/2022 | Phenolic + sensory + chemometric | Phenolic chemistry differentiates GIs | GI separation across 3 vintages |
| Castro et al. (Brazilian Cab Sauv) | Brazilian regions × altitude × temperature | GC-MS quantitative; IBMP, C13-norisoprenoids | IBMP scales inversely with bunch sun exposure / temperature; β-damascenone, β-ionone scale with vintage warmth | Vintage × altitude effect visible in norisoprenoid axis |

The chapter design should follow Pan/Souza Gonzaga — quantify, build PCA → OPLS-DA → random forest, report VIP scores, lock in the **chemical typicity fingerprint** for the source appellations of the sample set.

### Multivariate / chemometric methods to deploy

1. **Unsupervised: PCA + HCA.** First-pass cluster discovery; biplot of loading vectors against region/vintage/oak.
2. **Supervised: PLS-DA / OPLS-DA.** VIP score >1 retains the discriminant compounds; OPLS-DA used by Pan et al. and many recent Cab Sauv papers.
3. **Ensemble: Random Forest.** Bueno/Marrufo-Curtido used RF for oxidation markers; MDA score selects features robustly when collinearity is high.
4. **Validation: nested cross-validation** with stratification by vintage/region; permutation test for OPLS-DA Q²; ROC-AUC for binary tainted/clean classifiers; bootstrap CIs for VIP scores.

---

## 7. Oak compounds for oak-program differentiation

Oak introduces ~20-30 volatile compounds; these double as discriminators of *cooperage species, toast level, and barrel age*. Note guaiacol and 4-methylguaiacol are *legitimately* oak-derived in unmoke-affected wine — the same analytes that flag smoke at >10 µg/L are normal oak markers at 1-5 µg/L. The application chapter must separate **oak-baseline guaiacol** from **smoke-driven guaiacol** using the SyGG / GuRG glycoside ratio (only smoke loads the glycoside pool — see §1).

| Compound | Descriptor | Threshold (red wine) | Clean unoaked Cab Sauv | Oak-aged Cab Sauv |
|---|---|---|---|---|
| cis-Whisky lactone (cis-3-methyl-4-octanolide) | Coconut, sweet, oak | 74 µg/L red wine | <5 µg/L | 30-300 µg/L; markedly higher in *Q. alba* (American) than *Q. petraea/robur* (French) |
| trans-Whisky lactone | Grassy, celery, oak | 320 µg/L red wine | <5 µg/L | 20-150 µg/L; ratio cis/trans is much higher for American oak |
| Vanillin | Vanilla, pastry | ~60-200 µg/L red | <5 µg/L | 100-1500 µg/L; higher in *French* oak finished wines despite lower initial wood content, because of fermentation reduction balance |
| Syringaldehyde | Vanilla, smoky | ~50 mg/L (very high; rarely sensory) | <500 µg/L | >1 mg/L oaked |
| Acetovanillone | Vanilla | high | <100 µg/L | 100-1000 µg/L |
| Eugenol | Clove, spice | 6 µg/L (low!) — generally at or just above sensory in oaked wine | <2 µg/L | 5-50 µg/L; spicy intensity |
| Isoeugenol | Clove, spice | ~6-50 µg/L | <2 µg/L | 5-25 µg/L |
| 4-Allyl-2,6-dimethoxyphenol | Spice | ~ | trace | indicator of high toast |
| Furfural | Almond, caramel | 4-15 mg/L (>>sensory in red) | 1-15 µg/L (background) | 100-1000 µg/L oaked |
| 5-Methylfurfural | Caramel, butter | 16-45 mg/L | <5 µg/L | 50-200 µg/L oaked |
| 5-Hydroxymethylfurfural | Burnt sugar | high | trace | 100-500 µg/L oaked |
| Guaiacol (oak-derived) | Smoky, sweet | 23 µg/L | <2 µg/L | 1-10 µg/L from oak alone; higher = suspect smoke |
| 4-Methylguaiacol (oak-derived) | Smoky | ~65 µg/L | <2 µg/L | 1-5 µg/L from oak alone |
| Syringol (oak-derived) | Smoky, woody | high | <5 µg/L | 5-50 µg/L oaked Shiraz/Cab Sauv |
| Maltol | Caramel, sweet | high | trace | tens of µg/L oaked |
| 2-Furanmethanol (furfuryl alcohol) | Burnt sugar | high | <5 µg/L | tens of µg/L oaked |

**Cooperage × toast effects.**
- American oak (*Q. alba*): higher total whisky-lactone, higher cis/trans ratio (coconut), more vanillin extracted up-front but less retained after fermentation; expressive vanilla-coconut profile.
- French oak (*Q. petraea / Q. robur*): lower lactone load, higher fine-grain extraction, vanillin retained better post-fermentation, more eugenol/clove; finer "tobacco / cedar" character.
- Low toast: lactone-vanillin dominant.
- Medium toast: lactone + eugenol + vanillin + furfural balance.
- Heavy toast: guaiacol + 4-methylguaiacol + syringol + furanic compounds dominant; spicy/smoky/caramel.
- New vs second-use vs old: total volatile load roughly halves per use; volatile-phenol pool decays faster than vanillin/lactone.

The application chapter's oak-program discrimination axes are: (a) **cis-WL / trans-WL ratio** (American vs French), (b) **eugenol / vanillin ratio** (toast level), (c) **furanic + guaiacol load** (toast intensity), (d) **lactone + vanillin absolute** (new vs neutral barrel).

---

## 8. Multivariate analyses published on Cabernet Sauvignon aroma data

Beyond §6, three published categories of multivariate analysis are directly relevant to the dissertation's application chapter:

1. **Regional discrimination of Cab Sauv** — see Pan et al. (Foods 2022) on Ningxia 2013-2018, Souza Gonzaga et al. on Coonawarra/Margaret River/Yarra Valley/Bordeaux 2015 vintage (n=52), Niu et al. on Ningxia/California/Bordeaux. All use PCA + OPLS-DA with VIP>1 filtering. Citronellol, acetoin, IBMP, β-damascenone, β-ionone, hexanol, esters (ethyl hexanoate, ethyl octanoate, ethyl 3-methylbutanoate, isoamyl acetate) are the recurring discriminant compounds.

2. **Vintage discrimination within a single appellation** — Bueno, Carrascón & Ferreira (2018) on red wine oxidation; ARP (acetaldehyde reactive potential) as a vintage age-index. Pan et al.'s Ningxia panel showed vintage effects but subordinate to region. β-damascenone, β-ionone, ethyl esters (which hydrolyse with age), and the TDN class (1,1,6-trimethyl-1,2-dihydronaphthalene — sensory threshold 2 µg/L in Riesling, but found in Cab Sauv close to its threshold at 6-7 µg/L) are vintage markers.

3. **Smoke / fault classification** — Hu, Fuentes, Gonzalez Viejo et al. used E-nose + ANN regression to predict sensory smoke intensity in Cab Sauv; Noestheden et al. (AJEV 2024) classified oaked vs smoke-affected commercial wines using free volatile phenols + glycosides. The Wilkinson group's classifier validates that volatile phenol + glycoside features alone reach >90% smoke/no-smoke separation on Cab Sauv.

**Recommended workflow for the application chapter.**
1. Collect Cab Sauv panel (n ≥ 30 per appellation; ≥3 vintages; ≥3 oak programs per producer where possible) — at minimum two Napa/Sonoma + one Bordeaux + one Coonawarra/Margaret River subset.
2. Quantify ~60-80 analytes using validated SPME Arrow HS-GC-MS/MS MRM method (Chapters 3-5).
3. PCA → identify clusters and outliers.
4. OPLS-DA on appellation labels → VIPs → report top-20 discriminant compounds.
5. Random forest with permutation importance → cross-validate ranking.
6. Replicate analysis on **vintage** within best-represented appellation.
7. Build **oak program classifier** using the four ratios in §7.
8. Build **smoke classifier** using §1 free + bound markers; evaluate AUROC vs binary expert-panel ratings.
9. Validate all classifiers with stratified nested CV; report VIP, MDA, balanced accuracy, ROC-AUC, and (for OPLS-DA) permutation Q² p-value.

---

## 9. References

(References below are the canonical / accessible sources surveyed for this brief. Full bibliographic detail to be expanded into the dissertation reference manager.)

### Smoke taint
- Pollnitz, Pardon, Sykes, Sefton (2004) — *J. Agric. Food Chem.* — Stable-isotope-dilution GC-MS of guaiacol/4-MG; Anal. Bioanal. Chem. variant on sample prep and injection effects.
- Hayasaka, Baldock, Parker, Pardon, Black, Herderich, Jeffery (2010, 2013) — *J. Agric. Food Chem.* — Glycosylation of smoke-derived volatile phenols; HPLC-MS/MS quantitation.
- Krstic, Johnson, Herderich (2015) — *Aust. J. Grape Wine Res.* — Review.
- AWRI Smoke taint resource page; AWRI Sensory impact of smoke exposure fact sheet (Nov 2025 update); AWRI Smoke taint analysis FAQ; AWRI Smoke taint analytical methods page. — https://www.awri.com.au/industry_support/winemaking_resources/smoke-taint/
- Caffrey, Lerno, Rumbaugh, Girardello, Zweigenbaum, Oberholster, Ebeler (2019) — *AJEV* 70:4, 373 — *Changes in Smoke-Taint Volatile-Phenol Glycosides in Wildfire Smoke-Exposed Cabernet Sauvignon Grapes throughout Winemaking*.
- Noestheden, Thiessen, Dennis, Tiet, Zandberg (2024) — *AJEV* 75:1, 0750017 — *Prevalence of Wildfire Smoke Exposure Markers in Oaked Commercial Wine*.
- Fudge, Boss, Wilkinson — predictive ML / passive sampler studies.
- Yang, Xie, Hu, Fuentes, Gonzalez Viejo et al. (2021) — *Beverages / Foods* — E-nose ANN smoke prediction for Cab Sauv; PMC8398669.
- Caffrey et al. — Uptake and Glycosylation of Smoke-Derived Volatile Phenols by Cabernet Sauvignon Grapes (PMC7464031).
- Williams et al. (1995) — glycosyl-glucose (GG) assay foundation for acid hydrolysis approach.

### Brettanomyces
- Chatonnet, Dubourdieu, Boidron, Pons (1992) — *J. Sci. Food Agric.* — Origin of ethylphenols in wines; 4-EP, 4-EG perception thresholds.
- Chatonnet sequel work and AWRI Brettanomyces FAQ.
- Pinto et al. — Statistical Evaluation of 4-EP and 4-EG Concentrations to Support Sensory Evaluation of "Brett Character" — proposed revised threshold.
- Romano et al. (2009) — Brazilian survey of *Brettanomyces* and ethylphenols.
- Schumaker / Curtin / Coulter — wine style alters Brett impact (multiple).
- Lima et al. — 4-EC determination by HPLC-coulometric.

### Cork taint / haloanisoles / earthy taints
- Buser, Zanier, Tanner (1982) — original identification of TCA as the principal cork-taint compound.
- Pollnitz, Pardon, Liacopoulos, Skouroumounis, Sefton (1996) — *Aust. J. Grape Wine Res.* — TCA and other chloroanisoles in tainted wines and corks.
- Chatonnet et al. (2004) — TBA identification and responsibility.
- Cravero (2015) *J. Inst. Brew.* — Sensory evaluation of TCA in wines.
- Hjelmeland & Ebeler — geosmin in wine (J. Agric. Food Chem. 2000 baseline).
- Liger-Belair et al. — *Uncorking Haloanisoles in Wine* review (PMC10054257).
- HS-SPME-GC-QQQ-MS sub-ng/L methods (AJEV 2012 Coelho et al.).
- Thermal-desorption Vocus method (PMC8008374).

### Reductive faults
- Siebert, Solomon, Pollnitz, Jeffery (AWRI) — Selective GC-SCD determination of volatile sulfur compounds in wine.
- Ugliano, Kreitman, Waterhouse — copper / sulfide chemistry of post-bottling reduction.
- Franco-Luesma et al. — *J. Agric. Food Chem.* — DMS in Cab Sauv and grape S-methylmethionine.
- Goldberg & co — *Wine Faults: State of Knowledge in Reductive Aromas* (Molecules 2022, PMC9182507).
- Rauhut — wine sulfur chemistry chapter (in *Wine Microbiology* / *Wine Chemistry* texts).
- Ferreira et al. — Synthesis and SIDA SPME-GC-MS of EtSH and DEDS in wine.

### Oxidative faults
- Bueno, Marrufo-Curtido, Carrascón, Fernández-Zurbano, Escudero, Ferreira (2015) — *J. Agric. Food Chem.* — *Quantitative Analysis by GC-MS/MS of 18 Aroma Compounds Related to Oxidative Off-Flavor in Wines*. (Method-defining paper for the application chapter.)
- Bueno, Carrascón, Ferreira (2018) — *J. Agric. Food Chem.* — *Formation and Accumulation of Acetaldehyde and Strecker Aldehydes during Red Wine Oxidation* (PMC5817066 / Front. Chem. open version).
- Bueno et al. (2022) — *Foods* — Acetaldehyde Reactive Potential index (PMC8834303).
- Pereira et al. (Madeira; PMC6920768) — Is Sotolon Relevant to the Aroma of Madeira Wine Blends?
- Zea et al. (2015) — *Comprehensive Reviews in Food Science and Food Safety* — Acetaldehyde as Key Compound for the Authenticity of Sherry Wines.

### Application / multivariate
- Niu et al. (Foods 2022, PMC9801081) — HS-SPME/GC-MS + chemometrics across regions/varieties.
- Pan et al. (Foods 2022, PMC8750599) — Helan Mountains Cab Sauv 2013-2018 by HS-SPME-GC-MS + OPLS-DA + HCA.
- Souza Gonzaga et al. (2020) — *Food Res. Int.* / IVES — Regional Cab Sauv typicity (Coonawarra / Margaret River / Yarra Valley / Bordeaux).
- Casassa et al. (Mendoza Cab Sauv; ScienceDirect 2024 / 2026 in press) — geographical × vintage × plant material.
- Castro et al. — Brazilian Cab Sauv altitude/IBMP/norisoprenoid effects.
- Tao et al. — multivariate predictive modelling of oak-aged red wine aroma.

### Oak chemistry
- Pérez-Coello et al. — HS-SPME-GC-MS of oak wood volatiles (Pubmed 15537286).
- Spillman, Pollnitz, Liacopoulos, Skouroumounis, Sefton — vanillin barrel-aging accumulation (J. Agric. Food Chem. 1997).
- Pollnitz, Jones, Sefton — volatile composition of aged wine in used French/American oak barrels.
- Wilkinson, Cordingley, Boss, Capone — Cab Sauv aged in French oak with stave bending techniques.
- Spillman / Sefton oak lactone stereoisomer threshold paper.
- ETS Labs Oak Aroma Analysis library entry.
- ISC Barrels Trial 4 on bourbon whiskey lactone American vs European vs French oak.

### Other supporting
- Roujou de Boubée, Dubourdieu — Organoleptic impact of 2-methoxy-3-isobutylpyrazine on red Bordeaux and Loire wines (J. Agric. Food Chem. 2000).
- Sala et al. — IBMP location in Cab Sauv bunches and vinification extractability.
- Sacks, Gates, Ferry, Lavin, Kurtz, Acree (2012) — TDN sensory thresholds in Riesling vs non-Riesling.
- Pineau et al. — β-damascenone impact in red wines (J. Agric. Food Chem. 2007).
- Mendes-Pinto — C13-norisoprenoids review.

### Method-development references (cross-link to Chapters 3-5)
- A Simple GC-MS/MS Method for Determination of Smoke Taint-Related Volatile Phenols in Grapes (Metabolites 2020 — Tomasino group; PMC7407152).
- Agilent Application 5994-3161 — Analysis of Free Volatile Phenols in Smoke-Impacted Wines, SPME-8890 GC.
- Thermo Application 10752 — Rapid smoke-taint analysis of wine with SPME-GC-MS/MS.
- Cernicchiaro et al. (2024) — Fully Automated HS-SPME-GC-MS/MS Method for Odor-Active Carbonyls in Wines (J. Agric. Food Chem.).

---

**Hand-off note to writing team.** Sections 1-5 are concentration-fact tables ready to drop into the application chapter as Table 6.1-6.5. Section 6-8 frame the multivariate experimental design and the discriminant-compound shortlists Section 6 should aim to recover from real data. Section 7 oak table is the single most cited table in oak-related dissertation chapters and should be reproduced verbatim, with the dissertation's own measured values added as the rightmost columns.
