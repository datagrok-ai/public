# Box-Plot Analysis Patterns: Integrate as General Features, or Build a Specialized Toxicology Viewer?

## Executive Summary

**Bottom line: integrate all three patterns into the general-purpose box-plot viewer as a "structure-first" statistical core, and layer toxicology-specific behavior on top as an optional domain preset.** The three patterns drawn from the preclinical-toxicology cases are not toxicology inventions; they are mainstream, cross-disciplinary statistical and visualization practices that surface — under different names — in clinical trials, psychology, ecology, economics, A/B testing, genomics, agriculture, and more. Two of the three (stratified/matched comparison via faceting, and showing individual data points at small n) are essentially universal. The third (covariate adjustment via ratio/regress-out/ANCOVA) is very common but more methodologically contested and domain-flavored in its presentation. The only genuinely toxicology/pharmacology-leaning element is the specific *menu* of ordered dose-response trend tests (Jonckheere–Terpstra, Williams, Cochran–Armitage) and the control-vs-treatment Dunnett procedure — but even those have direct analogues in genetics, epidemiology, psychology, and education as "ordered-category trend tests." Existing tools implement the pieces in scattered, code-heavy ways; almost none offers an integrated point-and-click flow that combines faceting + matched baseline + covariate adjustment + the right test. That gap is the product opportunity, and it argues for a general core, not a niche viewer.

- **Universality:** Pattern 1 (stratify/facet) and Pattern 3's "show the data" display are near-universal; Pattern 2 (covariate adjustment) is widespread but contested; ordered trend tests are domain-flavored but cross-cutting.
- **Tooling gap:** Faceting and point display are widely available; matched-baseline + covariate-adjusted comparison + appropriate test as a single UI flow is essentially unserved.
- **Recommendation:** Build a general statistical/visual core (facet, point overlay, covariate adjust, control comparisons) with optional domain presets; do not silo this into a toxicology-only viewer.

---

## Pattern 1 — Stratified / Matched-Control Comparison (Faceting, Subgroup Analysis)

**Universality rating: 9.5/10 — essentially universal.**

Comparing treated groups to controls *within* a stratum rather than against a pooled control is one of the most deeply embedded ideas in applied statistics, appearing under the names stratification, subgroup analysis, blocking, faceting/small-multiples, and interaction/effect-modification analysis.

**Clinical trials / biostatistics.** Subgroup analysis is a first-class, regulated activity. The CONSORT reporting guideline explicitly addresses subgroup analyses; the NEJM "Statistics in Medicine" guidance defines subgroup analysis as "any evaluation of treatment effects for a specific end point in subgroups of patients defined by baseline characteristics" and notes that "most trial reports did include subgroup analyses." Methodological consensus stresses pre-specification, limiting the number of subgroups, and — crucially — using **interaction tests rather than separate within-subgroup p-values** (Pocock/Clayton: "an appropriate analysis typically calculates effect estimates and confidence intervals within such subgroups together with an overall interaction test"). Stratified randomization is standard, and "the validity of a subgroup finding is improved when stratified for at randomisation."

**The cross-cutting concept: Simpson's paradox.** The danger that pooling subgroups reverses or hides effects is a canonical statistical phenomenon. The Stanford Encyclopedia of Philosophy describes it as an association that "emerges, disappears or reverses when the population is divided into subpopulations," with implications "for a range of areas... including decision theory, causal inference, and evolutionary biology." Real-world examples span medicine (treatment efficacy reversing by group), social science (UC Berkeley gender-admissions case), economics (US median wage rising overall but falling within every education subgroup 2000–2012), epidemiology (COVID case-fatality rates between Italy and China by age), and A/B testing (inconsistent traffic splits across weeks). A dedicated paper in the Journal of Statistical Distributions and Applications is titled "The ubiquity of the Simpson's Paradox."

**Animal/preclinical research — the toxicology-adjacent mandate.** NIH notice NOT-OD-15-102 (issued June 2015, effective for applications submitted on or after January 25, 2016) states: "NIH expects that sex as a biological variable will be factored into research designs, analyses, and reporting in vertebrate animal and human studies. Strong justification ... must be provided for applications proposing to study only one sex." The ARRIVE guidelines and the SAGER (Sex and Gender Equity in Research) guidelines reinforce sex-disaggregated analysis. This is precisely the "compare males to male controls, females to female controls" stratification pattern, mandated across all of preclinical biology, not just toxicology.

**Visualization lineage.** Faceting is a foundational visualization principle: Tufte's "small multiples" and the Becker–Cleveland–Shyu "trellis displays" (1996), formalized as "faceting" by Wickham's ggplot2 and Wilkinson's Grammar of Graphics. Tufte's dictum — "At the heart of quantitative reasoning is a single question: Compared to what?... For a wide range of problems in data presentation, small multiples are the best design solution" — is the conceptual charter for this pattern.

**Verdict:** Faceting/stratified comparison is universal across every quantitative field. Building it into a general box-plot viewer is unambiguously justified.

---

## Pattern 2 — Covariate Adjustment / Normalization of One Column by Another (Ratio, Regress-Out, ANCOVA)

**Universality rating: 8/10 — very common, but methodologically contested and domain-flavored.**

Adjusting a response by a numeric covariate — by ratio (Y/covariate), by residualization (regress Y on covariate), or by ANCOVA — is a workhorse method across fields, but it is also one of the most debated, which has direct UI implications.

**Clinical trials.** ANCOVA is the standard tool for baseline adjustment. Trial statistical analysis plans routinely specify "change from baseline... analysed with an analysis of covariance (ANCOVA) model" with treatment and baseline covariates, reporting least-squares (adjusted) means. The biostatistics literature notes "the ANCOVA estimator is immune to such conditional bias, asymptotically," making covariate adjustment a recommended, even FDA-relevant, practice.

**Ecology / evolutionary biology — allometry.** Body-size correction is pervasive. A Journal of Experimental Biology review notes "commonly used body-size correction methods, such as ANCOVA and residual analysis, assume that the slopes of relationships between trait size and body size... are constant among treatment groups." García-Berthou (2001), "On the misuse of residuals in ecology: testing regression residuals vs. the analysis of covariance" (Journal of Animal Ecology 70:708–711), is a methods landmark: it documents that the "residual index" appeared in 8% of 1999 Journal of Animal Ecology papers and 2% of Ecology papers, argues the practice is "incorrect for at least four reasons," and concludes ANCOVA is the correct analysis.

**The cross-cutting cautionary concept: the ratio fallacy.** Richard A. Kronmal's "Spurious Correlation and the Fallacy of the Ratio Standard Revisited" (J. R. Statist. Soc. A, vol. 156, no. 3, pp. 379–392, 1993) is the canonical warning: "It is shown that the use of ratios in regression analyses can lead to incorrect or misleading inferences. A recommendation is made that the use of ratios in regression analyses be avoided." This traces to Pearson (1897) on spurious correlation of ratios with common components. The critique recurs across physiology (per-weight/per-surface-area standards), sociology, occupational health (effort-reward ratios), and morphometrics. This matters for the product: offering ratio mode without offering regress-out/ANCOVA would push users toward a method that statisticians broadly criticize — so the viewer should offer all three modes and ideally nudge toward regression adjustment.

**A/B testing / data science.** Covariate adjustment is now mainstream under the name CUPED (Controlled-experiment Using Pre-Experiment Data), introduced by Deng, Xu, Kohavi & Walker at Microsoft in their 2013 WSDM paper "Improving the Sensitivity of Online Controlled Experiments by Utilizing Pre-Experiment Data," which reports variance reduction of "more than 45%" when using a metric's pre-experiment value as the covariate; it has since been adopted at Netflix, Booking.com, Airbnb, the BBC, Optimizely and others, and practitioners report typical variance reductions of 20–50%. It is explicitly "a regression-based covariance adjustment method" — the same regress-out idea, repackaged for experimentation.

**Genomics / bioinformatics.** Regressing out technical covariates is standard. limma's `removeBatchEffect()` is documented as analogous to covariate adjustment: "The covariates argument allows correction for one or more continuous numeric effects, similar to the analysis of covariance method in statistics." (Caveat: limma's authors recommend including covariates in the design matrix for formal testing rather than pre-removing them — Aaron Lun, on the Bioconductor support forum: "we would recommend including additional covariates or factors in the design matrix rather than removing them with removeBatchEffect." The viewer should therefore treat regress-out as an exploratory/visualization transform.)

**The unifying concept: adjusted/estimated marginal means.** What R's emmeans calls estimated marginal means, SAS calls LSMEANS, SPSS reports as EMMEANS, and the ANCOVA literature calls "adjusted means" — all the same "model-based, equally weighted, covariate-adjusted predicted means." That this single concept has four names across four software ecosystems is itself evidence of cross-field universality.

**Verdict:** Covariate adjustment is broadly used and clearly belongs in a general viewer — but because ratios are statistically contested, the feature should offer ratio, regress-out, and ANCOVA modes with appropriate guidance, not ratio alone.

---

## Pattern 3 — Small-N Point Display + Ordered/Dose-Response Trend Tests

**Universality rating: display = 9/10 (near-universal); ordered trend tests = 6.5/10 (cross-cutting but domain-flavored).**

This pattern bundles two distinct ideas with different universality.

### 3a. "Show the data" at small n — near-universal

The movement against bar/box plots that hide raw data and toward dot/strip/beeswarm displays is one of the strongest cross-field methodological trends of the past decade. Weissgerber, Milic, Winham & Garovic, "Beyond Bar and Line Graphs: Time for a New Data Presentation Paradigm" (PLOS Biology 13(4):e1002128, 2015), conducted a systematic review of n = 703 articles in top physiology journals and found that "85.6% of papers included at least one bar graph," arguing that individual data points should be shown, especially for small samples. Major journals adopted policies requiring data-distribution displays — including **PLOS Biology, eLife, the Journal of Biological Chemistry, and Nature and affiliated journals.** Guidance is explicit and sample-size dependent: "Univariate scatterplot or dotplots showing the raw data points are the best option for very small samples (n ≤ 10 per group)... Dot plots should be used for small datasets, whereas dot plots may be combined with box plots or violin plots for medium-sized datasets." This applies everywhere small-n continuous data are compared — physiology, cell biology, pharmacology, psychology, materials testing — making it nearly universal.

### 3b. Ordered/trend tests and control-vs-treatment comparisons — cross-cutting but flavored

Tests for monotone trend across ordered groups, and comparisons of several treatments to a single control, are the most domain-specific element — but still far from toxicology-exclusive.

- **Dunnett's test** (many-to-one, comparing each treatment to a control) is a general-purpose multiple-comparison procedure available in essentially every stats package; the US National Toxicology Program recommends Dunnett and Williams, but Dunnett itself is used across biology, agriculture, and medicine.
- **Jonckheere–Terpstra** (nonparametric ordered-alternatives trend test) is "often used in toxicology" and recommended by the OECD for ecotoxicity data, but its documented application areas explicitly include "Medical Research," "Psychology," and "Education," and it is used in genetics as a nonparametric alternative to Cochran–Armitage.
- **Cochran–Armitage trend test** is, if anything, most prominent *outside* toxicology: it is a workhorse in epidemiology (dose-response of exposures like smoking) and a standard genotype-based test in **genome-wide association studies**.
- **Williams' test** (monotone dose-response) is more pharmacology/toxicology-specific; the NTP recommends the umbrella-protected Williams test combined with Dunnett.

So the *category* "ordered-category trend test + many-to-one comparison" is cross-disciplinary; the specific preferred test varies by field. This argues for a general "ordered/trend test" capability with field-appropriate defaults rather than a hardcoded toxicology test battery.

**Verdict:** Point display is universal and belongs in the core. Trend/Dunnett testing belongs in the core as a general capability, with the specific test menu surfaced via domain presets.

---

## Existing-Tool Support

Across leading tools, the individual pieces exist but are scattered, mostly code-driven, and rarely combined into one guided flow.

| Tool | Stratify / facet + matched control | Covariate adjust (ratio / regress-out / ANCOVA) | Small-N points + trend/Dunnett | UI model | Integrated "structure-first" flow? |
|---|---|---|---|---|---|
| **R: ggplot2 + ggpubr/rstatix + emmeans** | `facet_wrap`/`facet_grid`; `stat_compare_means`/`geom_pwc` with `ref.group` for control comparisons | emmeans/rstatix `emmeans_test` give ANCOVA adjusted means; ratios/residuals by manual transform | `geom_jitter`/`geom_dotplot`; trend tests via separate packages | Scripting | Partial — assemblable, but no single UI; faceted p-value placement is a known pain point |
| **Python: seaborn** | `catplot` with `col`/`row` facets, `hue` splits; `boxplot`/`stripplot`/`swarmplot` | None native (manual pandas transforms; statsmodels for ANCOVA) | `stripplot`/`swarmplot` overlay on `boxplot` | Scripting | No — visualization only, no built-in stats/tests |
| **GraphPad Prism** | Grouped tables; Dunnett, ANOVA "post-test for trend"; strong "plot the points" guidance | **No native ANCOVA** — official FAQ #607: "Neither Prism nor InStat do ANCOVA"; workaround via multiple regression (since Prism 8) or linear-regression slope/intercept comparison | One-click points-on-graph; Dunnett many-to-one; ANOVA test-for-trend | Point-and-click | Partial — best consumer-grade UI, but covariate adjustment is a gap and faceting is limited |
| **JMP / SAS** | Fit Model / GLM with interactions; LSMeans Dunnett; by-group analysis | Full ANCOVA (Fit Model, LSMEANS adjusted means); SAS PROC GLM/MIXED | Dunnett, Jonckheere, trend tests native | Point-and-click (JMP) / code (SAS) | Closest to integrated — but ANCOVA, faceting, and tests are separate platforms, not one box-plot flow |
| **SPSS / Stata** | GLM EMMEANS (SPSS); `margins`, subgroup, `nptrend` (Stata) | ANCOVA via GLM; adjusted means | Stata `nptrend` (Cochran–Armitage/Jonckheere-type); Dunnett | Menu + code | No single box-plot-centric flow |
| **Tableau / Power BI** | Small multiples / trellis (Power BI small multiples started as a preview, bar/column/line/area only; Tableau via calculated fields) | None (limited stats) | Limited; no built-in trend/Dunnett tests | Point-and-click | No — visualization-first, statistics minimal |
| **Spotfire / Origin / Minitab / Plotly** | Trellis/panels (Spotfire); facets vary | ANCOVA in Minitab/Origin; none in Plotly | Dot/jitter overlays vary; Dunnett in Minitab | Mixed | No integrated structure-first box-plot flow |

**Key tooling findings:**
- **Faceting/small-multiples is broadly available** (ggplot2, seaborn, Tableau, Power BI, Spotfire) — confirming Pattern 1's universality is matched by tool support.
- **Point overlay is everywhere** in plotting libraries (ggplot2, seaborn) and is one-click in Prism.
- **Covariate adjustment/ANCOVA is the biggest discoverability gap in consumer tools.** The most popular bench-science GUI, GraphPad Prism, still has no native ANCOVA command (confirmed against Prism 10/11 documentation and release notes); its official FAQ tells users to hand-build dummy-variable multiple regression or compare regression slopes/intercepts. Power BI/Tableau offer essentially none. ANCOVA lives mostly in code (emmeans, statsmodels) or heavyweight statistical platforms (JMP, SAS, SPSS, Stata, Minitab).
- **No tool offers the full "structure-first" box-plot flow** — facet by a stratifier + matched baseline + covariate adjust + appropriate control/trend test — as a single, guided UI. ggplot2+rstatix+emmeans can assemble it but only by scripting; JMP/SAS come closest but split it across separate platforms.

---

## Synthesis: General Feature vs. Specialized Viewer

**The evidence points decisively toward a general-purpose core, not a separate toxicology viewer.** Here is the weighing:

**Arguments for "general features":**
1. **Two of three patterns are universal.** Faceting/stratified comparison (Pattern 1) and "show the data points" (Pattern 3a) are recommended or required across virtually every quantitative discipline and are backed by reporting guidelines (CONSORT, ARRIVE, NIH SABV), foundational visualization theory (Tufte, Cleveland), and journal mandates (PLOS Biology, eLife, Nature, JBC).
2. **The cross-cutting *concepts* are field-neutral.** Simpson's paradox, ANCOVA/adjusted means, "plot the data," and ordered-category trend tests are taught as general statistics, not toxicology methods. The same construct appears under four names (emmeans/LSMEANS/EMMEANS/adjusted means) across four software ecosystems.
3. **Covariate adjustment (Pattern 2) is broad** — clinical baseline adjustment, allometry, CUPED in A/B testing, batch correction in genomics — even if its presentation is domain-flavored and statistically contested.
4. **There is a real, unfilled tooling gap.** No mainstream tool delivers the integrated structure-first flow through a UI. This is a horizontal opportunity across data science, not a vertical toxicology niche.

**Arguments for "specialized viewer" (and why they're weaker):**
1. The specific *test battery* (Jonckheere–Terpstra, Williams, umbrella-protected Williams, NTP-style Dunnett+Williams combinations) is genuinely pharmacology/toxicology-flavored. But the *category* (ordered trend + many-to-one) generalizes, and the field-specific tests can be encapsulated as a preset rather than as a separate app.
2. Some defaults (monotone dose ordering, regulatory conventions) are domain-specific. Again, these are configuration, not architecture.

**The hybrid resolution — general core + domain presets — is the strongest design.** Build:
- **General statistical/visual core (always on):** faceting by one or more stratifiers; choice of pooled vs. matched/within-stratum baseline; point overlay with sample-size-aware defaults (auto-show points at small n, e.g. n ≤ ~10–15); covariate adjustment offering ratio, regress-out, and ANCOVA modes with guidance steering users away from naive ratios; control-vs-treatment comparison (Dunnett) and a general ordered-trend test.
- **Optional domain presets (one click):** a "Toxicology/dose-response" preset that defaults to monotone dose ordering, Jonckheere–Terpstra/Williams trend tests, and NTP/OECD conventions; a "Clinical trial" preset (ANCOVA baseline adjustment, interaction-test framing for subgroups); an "A/B test" preset (CUPED-style covariate adjustment); a "Genomics" preset (regress-out as exploratory transform with the in-model caveat). Presets configure defaults and test menus on top of the same core engine.

This architecture maximizes addressable users (every field that compares groups), avoids fragmenting the product, and still serves toxicologists fully via a preset. A separate viewer would duplicate the universal 80% (faceting, points, control comparison, covariate adjustment) just to package the domain-specific 20% (the exact trend-test menu), which is poor engineering economics.

---

## Recommendations

**Stage 1 — Build the universal core first (highest ROI, lowest domain risk).**
- Ship faceting/stratified comparison with an explicit **pooled vs. matched-control toggle**, and surface an interaction/heterogeneity indicator (the CONSORT/Pocock best practice) plus a Simpson's-paradox warning when pooled and stratified conclusions diverge.
- Ship **sample-size-aware point display**: automatically overlay individual points (jitter/beeswarm) on top of the box at small n, citing the journal-policy rationale (n ≤ 10 per group → univariate scatter) in tooltips. This is low-risk, near-universally demanded, and a visible differentiator.
- *Benchmark that would change priorities:* if user telemetry shows faceting is rarely used, deprioritize the interaction-test layer — but evidence strongly predicts heavy use.

**Stage 2 — Add covariate adjustment as a guided feature.**
- Offer **three modes — ratio, regress-out (residualize), ANCOVA-adjusted means** — with inline guidance that flags the ratio fallacy (Kronmal 1993) and recommends regression adjustment by default. This directly fills the biggest gap left by Prism, Tableau, and Power BI.
- Treat regress-out as an exploratory/visualization transform and surface the limma-style caveat that formal inference should adjust in-model.
- *Threshold:* prioritize this for Stage 2 (not Stage 1) because it carries the most statistical-correctness risk; ship only with the guardrails.

**Stage 3 — Add control/trend testing with domain presets.**
- Implement **Dunnett (many-to-one)** and a **general ordered-trend test** in the core.
- Encapsulate the **toxicology/dose-response test battery** (Jonckheere–Terpstra, Williams, umbrella-protected Williams) as a one-click preset; add clinical, A/B, and genomics presets over time.
- *Threshold:* gate field-specific test menus behind presets so the core UI stays uncluttered for general users.

**Do not** build a standalone toxicology box-plot viewer. The domain-specific surface area is too small to justify a separate product, and doing so would forfeit the much larger horizontal market for an integrated structure-first box-plot experience that no incumbent currently offers.

---

## Caveats

- **"Universality ratings" are reasoned judgments,** not measured frequencies; they synthesize the strength and breadth of authoritative sources, but no single study quantifies cross-field prevalence of all three patterns side by side.
- **Tool capabilities evolve.** Power BI small multiples launched as a preview feature with documented limitations; Prism's lack of native ANCOVA is current as of Prism 10/11 (verified against GraphPad's FAQ and release notes) but could change. Verify against current release notes before finalizing competitive claims.
- **Covariate adjustment is statistically contested.** The ratio-vs-residual-vs-ANCOVA debate (Kronmal; García-Berthou; limma authors) means this feature must ship with guidance; a naive implementation could actively mislead users.
- **Trend tests carry assumptions** (monotonicity, ordered groups, distributional shape) that a general UI must not hide; the toxicology preset should expose these rather than auto-applying tests blindly.
- **Subgroup analysis is double-edged.** The same clinical-trial literature that mandates it also warns it is widely misused (multiplicity, post-hoc fishing). A general viewer that makes stratified comparison trivially easy should pair it with multiplicity warnings and pre-specification prompts to avoid amplifying a known abuse.