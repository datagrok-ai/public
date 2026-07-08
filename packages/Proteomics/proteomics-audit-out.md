# Numerical & Mathematical Correctness Audit — Datagrok Proteomics Plugin

Findings were verified against published reference implementations and primary literature. Three primary sources — Benjamini & Hochberg (1995), Huber et al. (2002, VSN), and Zhu et al. (2020, DEqMS) — were subsequently cross-checked **against the original PDFs**; the affected sections are marked "verified against source."

## Audit scope and method
Seven TypeScript files were audited **against published authoritative reference implementations and primary literature**, not merely for internal consistency. Each function received a line-level read of the described formulas, loop bounds, boundary/tie handling and numerical form. Severity scale: **Critical / High / Medium / Low / Informational**.

## Executive summary (severity-sorted)

| # | Severity | File / function | Finding |
|---|----------|-----------------|---------|
| 1 | **High** | spc.ts `rule4` | Nelson Test 4 is 14 points alternating in DIRECTION relative to the PREVIOUS point (a zigzag/sawtooth). The code derives "alternation" from whether points are above/below the centerline (mean). Wrong pattern → misses true sawtooth signals and can mis-fire on oscillation about the mean. |
| 2 | **Medium-High** | enrichment.ts `buildEnrichmentDf` | `r.intersections[j]` aligned to `queryGenes[j]` is only correct when every query gene maps 1:1 to a single Ensembl ID. The g:Profiler API aligns the `intersections` array positionally to the mapped `ensgs` list, not the raw query. Produces wrong member-gene names when unmapped/multi-mapped/background-filtered genes exist. |
| 3 | **Medium** | spc.ts `rule8` | Nelson Test 8 requires points on BOTH sides of the centerline (with none in Zone C). Code only checks `|z|>1` for all 8, not the both-sides condition. (Minitab uses an "either side" convention — defensible but non-literal vs Nelson 1984.) |
| 4 | **Medium** | imputation.ts Box-Muller | `u1` can be exactly 0 (`Math.random()` returns [0,1)), giving `log(0) = −Infinity`. Unguarded → occasional `Infinity`/`NaN` imputed values. |
| 5 | **Medium** | imputation.ts `imputeKnn` | Troyanskaya (2001) KNNimpute uses a similarity/inverse-distance-WEIGHTED average of neighbors; the code uses an UNWEIGHTED mean. Deviation from the cited reference. |
| 6 | **Low / Informational** | normalization.ts `quantileNormalize` | The fractional-rank interpolation for unequal non-missing counts is consistent with limma `normalizeQuantiles`' documented NA handling — a valid variant, not a bug (but differs from preprocessCore; document which reference is intended). |
| 7 | **Informational** | pca.ts | Covariance trick and score scaling are mathematically valid; Jacobi rotation matches the standard stable formula. |

All other functions audited (Pearson, median, Welford CV, MA-plot M/A, log2FC, BH FDR independence, MinProb parameters, control-chart constants, PCA centering, VSN scale) are **correct** and, where a primary source was supplied, verified against it.

---

## FILE 1 — spc.ts (Statistical Process Control — Nelson rules)

**Reference:** Lloyd S. Nelson, "The Shewhart Control Chart—Tests for Special Causes," *Journal of Quality Technology* 16(4):237–239 (1984); Western Electric *Statistical Quality Control Handbook* (1956/1958); AIAG/Montgomery/Wheeler control-chart constant tables.

**Control-chart constants (verified).** For n=2 moving-range charts, `WHEELER_D2 = 1.128` and `WHEELER_D4 = 3.267` are the correct standard values. The AIAG SPC constant table lists d2 = 1.128 and D4 = 3.267 for n=2, and Wheeler & Chambers (1992) *Understanding Statistical Process Control* is the correct citation basis for the moving-range unbiasing constants. **Verdict: Correct — Informational.**

- **rule1 — one point beyond 3σ (Test 1).** Substance correct. Uses strict `|z| > 3`, so a point exactly at 3σ is not flagged — a defensible boundary choice (Nelson's zones are open at 3σ). Evaluating only the tail point suits streaming detection but will not re-flag historical violations if the series is re-scanned. **Verdict: Correct — Low.**

- **rule2 — nine points same side (Test 2).** Correct. Nelson Test 2 = nine consecutive points on the same side of the centerline. Using strict `v > mean` / `v < mean` means a point landing exactly on the mean belongs to neither side and correctly terminates the run. **Verdict: Correct.**

- **rule3 — six steadily increasing/decreasing (Test 3).** Correct. "Six points in a row steadily increasing or decreasing" = six points / five intervals; the code's use of strict `>` enforces strict monotonicity so a tie correctly breaks the run. Verify the window is 6 points (5 comparisons), not 7. **Verdict: Correct.**

- **rule4 — fourteen points alternating up and down (Test 4). LIKELY BUG — High.** Nelson 1984 Test 4 is verbatim **"Fourteen points in a row alternating up and down"** — i.e. each point reverses direction *relative to the immediately preceding point* (a zigzag/sawtooth of 13 direction reversals). The code instead classifies each point by whether it is above or below the mean (`sides` derived from `v > mean` / `v < mean`) and requires that above/below label to alternate. That is a fundamentally different pattern: a genuine sawtooth that stays entirely above the mean would be MISSED, and a slow oscillation that crosses the mean every point but is not a local zigzag could FALSE-fire. **Fix:** compute `d[i] = sign(v[i] − v[i−1])` and require `d[i] = −d[i−1]` for 13 consecutive steps across 14 points; ignore/handle zero-differences explicitly. **Severity: High.**

- **rule5 — two of three beyond 2σ, same side (Test 5).** Substantially correct (Nelson Test 5; equivalent to Western Electric Rule 2). Strict `z > 2` / `z < −2` correctly encodes the Zone A boundary. The WECO/Nelson definition is "two of the **three most recent** points," and the code's current+prior counting with the most recent point participating matches the standard online semantics (the odd third point may be anywhere; the two Zone-A points must be on the same side). **Verdict: Correct — Low.**

- **rule6 — four of five beyond 1σ, same side (Test 6).** Matches Nelson Test 6 (Zone B, `z > 1` / `z < −1`, same side, most recent point participating). **Verdict: Correct.**

- **rule7 — fifteen points within 1σ (Test 7).** Correct. Nelson 1984 R7 is literally "Fifteen points in a row in Zone C (above and below center line)"; the code's `|v| < 1` for 15 points captures this. Strict `<` means a point sitting exactly on 1σ breaks the run — a minor boundary choice, acceptable. **Verdict: Correct — Low.**

- **rule8 — eight points beyond 1σ (Test 8). DEVIATION — Medium.** Nelson 1984 R8 is verbatim **"Eight points in a row on both sides of center line with none in Zone C"** (a mixture pattern). The code checks only `|v| > 1` for all 8 points and does **not** enforce that at least one point is above and at least one below the centerline. As written it will also fire on eight consecutive points beyond 1σ all on one side — which is a different (shift) signal. Note: Minitab's implementation deliberately states "Test 8: Eight points in a row more than 1σ from center line (**either side**)," so the code matches the Minitab convention but not Nelson's literal "both sides" definition. **Fix (to match Nelson literally):** additionally require the 8-point window to contain at least one point above and at least one below the mean. If the Minitab convention is intended, document it explicitly. **Severity: Medium.**

- **computeIterativeBaseline — iterative 3σ outlier removal (cap 2 iterations, sample SD n−1).** The double mean/SD recompute pattern (compute stats → drop |z|>3 → recompute) is sound and standard for a robust baseline; capping at 2 iterations bounds cost and avoids over-trimming. z is correctly recomputed against the updated mean/SD before the second drop. **Verdict: Correct — Low.**

- **meanSd** — sample SD (n−1) for n>1, sd=0 for n=1. **Correct.**
- **median** — standard even/odd handling. **Correct.**
- **pearson** — single-pass deviation-product form (Σ(x−x̄)(y−ȳ) over √Σ(x−x̄)²·Σ(y−ȳ)²), which avoids the catastrophic cancellation of the Σxy − nx̄ȳ form; returns NaN for n<3. **Verdict: Correct — good numerical form.**

---

## FILE 2 — differential-expression.ts

**References:** Welch's unequal-variances t-test; Benjamini & Hochberg (1995) *JRSS-B* 57(1):289–300 **(verified against source PDF)**; limma (Ritchie et al. 2015; Smyth 2004); DEqMS (Zhu et al. 2020, *Mol Cell Proteomics* 19:1047–1057) **(verified against source PDF)**.

- **log2FC = mean(vals2) − mean(vals1) on log2 data.** Correct: a difference of means in log2 space equals the log2 ratio of geometric-mean intensities. Sign convention = group2 − group1 (treatment − control). **Verdict: Correct.**

- **Minimum 2 per group for Welch.** This is the mathematical minimum for a within-group variance estimate; power at n=2 is very low and the Welch–Satterthwaite df is unstable, but it is not incorrect. **Verdict: Correct — Low (consider warning at n=2).**

- **BH FDR alpha-independence claim — VERIFIED CORRECT against Benjamini & Hochberg (1995).** The comment that "BH-adjusted p-values are alpha-independent, so alpha does not change corrected" is accurate. One subtlety is worth stating explicitly, because it can look like a contradiction: the *original BH procedure* in Benjamini & Hochberg (1995), Section 3.1, expression (1), is written as a **rejection rule** — "let k be the largest i for which P₍ᵢ₎ ≤ (i/m)·q*; then reject all H₍ᵢ₎, i = 1,…,k" — and there the threshold sequence (i/m)·q* **does** depend on the target rate q*. However, the *adjusted p-values* returned by `fdrcorrection` are a re-expression of that same procedure: p_adj₍ᵢ₎ = min over j≥i of (m/j)·p₍ⱼ₎ (cumulative-minimum of the step-up ratios), a quantity computed **purely from the sorted raw p-values and m**, with no α in it. `statsmodels.stats.multitest.fdrcorrection` with `method='i'`/`'indep'` maps to `fdr_bh`; the returned `pvalue-corrected` array does not depend on `alpha`, and `alpha` only determines the boolean `rejected` array (equivalent to comparing p_adj ≤ α, which reproduces the largest-k rule of expression (1)). Passing `pThreshold` as alpha therefore affects only which rows are flagged significant, not the corrected q-values. The independence caveat also matches the source: the 1995 proof (Theorem 1) establishes FDR control at q* **for independent test statistics**, so `method='i'` is the correct choice only under independence or positive dependence — the usual working assumption for per-protein DE tests. **Verdict: Correct.**

- **copyDEResultsToFrame — realign by 1-based `row` (target = row − 1).** Correct provided the R backend emits 1-based row indices (R is 1-indexed, so this is the right convention). The significance truthiness coverage (`true` / `1` / `'1'` / `'TRUE'`) is a sensible superset guarding against JSON/CSV type coercion. **Verdict: Correct — Low.**

- **DEqMS — peptide counts, default 1, Math.round — VERIFIED against Zhu et al. (2020).** DEqMS's `spectraCounteBayes` moderates the prior variance using the **per-protein PSM/peptide count** as the loess-regression covariate. The source is explicit: "we modified the IBMT function of Sartor et al. by changing the regression covariate from intensity to PSM or peptide count," and Equation (13) is `fitted(logVAR) = loess(logVAR ~ x)$fitted` with x = log2 of the peptide/PSM count. Integer counts are expected; defaulting a missing count to 1 (the minimum meaningful count → least variance shrinkage) and rounding to integer are safe and conservative. **One addition from the source:** for data sets spanning **multiple TMT experiments**, the paper's RSS analysis (data set D5, Fig. 4) shows the best prior-variance fit is obtained from the **minimum** PSM count across experiments ("The minimum number of quantified PSMs across multiple TMT sets was used in DEqMS"). The code passes the user-selected count column straight through to the server, which is fine for a single set but should apply or at least document the "minimum-count" rule for multi-set inputs to avoid a worse variance-count fit. **Verdict: Correct; add multi-TMT minimum-count note.**

---

## FILE 3 — enrichment.ts

**References:** g:Profiler — Raudvere et al. (2019) *NAR*; Kolberg et al. (2023) update; official g:GOSt API docs (biit.cs.ut.ee/gprofiler/page/apis); gprofiler2 CRAN package.

- **gGOSt params** (`significance_threshold_method='fdr'`, `domain_scope='custom'`, `background=backgroundGenes`, `ordered=false`, `measure_underrepresentation=false`). These produce FDR-corrected **over-representation** results against a custom background — the correct configuration for background-aware ORA. **Verdict: Correct.**

- **"Gene Ratio" → precision.** g:Profiler defines **`precision = intersection_size / query_size`** (verbatim from the gprofiler2 docs: "precision — the proportion of genes in the input list that are annotated to the function (defined as intersection_size/query_size)"). Note this is the ratio relative to the *query size*, whereas clusterProfiler's `GeneRatio = k/n` uses k = term hits and n = number of query genes with **any** annotation — the two coincide only when `query_size` equals the annotated-query count. **Verdict: Correct, with a terminology caveat — Low.** Recommend labeling the column "precision (intersection/query)" to avoid confusion with clusterProfiler GeneRatio.

- **FDR column from p_value.** With `significance_threshold_method='fdr'`, g:GOSt returns the **corrected** p-value in the `p_value` field, so populating FDR from `p_value` is correct. **Verdict: Correct.**

- **buildEnrichmentDf `memberGeneStrs` — intersections[j] aligned to queryGenes[j]. DEVIATION / likely bug — Medium-High.** The official g:Profiler API documentation states that the `intersections` array is **positionally aligned to the mapped Ensembl gene list** (`genes_metadata → query → <query_name> → ensgs`), NOT to the raw user query: *"The elements in this list correspond to query genes and are in the same order as the genes in genes_metadata → query → query_name → ensgs … Empty lists mean 'no intersection'; lists with elements mean 'this gene is part of the intersection'; null values mean the gene wasn't looked at (ordered queries)."* Indexing `r.intersections[j]` against `queryGenes[j]` (the raw input order) is therefore only correct when every input maps **1:1** to a single Ensembl ID, with no unmapped genes, no genes expanding to multiple ENSGs, and no background/annotation-domain filtering. When any of those occur, `ensgs` diverges from the raw query and the recovered member-gene names will be **misaligned**. The gprofiler2 R client does this correctly: it takes `ensgs`, selects positions `which(lengths(evcodes) > 0)`, then reverse-maps ENSG→symbol via the `mapping` table. **Fix:** align `intersections[j]` to `genes_metadata.query.<name>.ensgs[j]` and reverse-map to the original symbol via `mapping`, rather than to `queryGenes[j]`. **Severity: Medium-High.**

- **splitGenesByDirection** — up = fc>0, down = fc<0 among significant, background = all detected genes. This is standard practice for directional (up/down) enrichment with a detected-proteome background. **Verdict: Correct.**

- **applySmartPathwayFilter** — a heuristic ported verbatim from a Python tool; not a published statistical method. Internal logic only; results should not be presented as a statistically principled filter. **Verdict: Informational — label as a heuristic.**

---

## FILE 4 — normalization.ts

**References:** Bolstad et al. (2003) *Bioinformatics* 19(2):185–193 (quantile normalization); Huber et al. (2002) *Bioinformatics* 18 suppl 1:S96–S104 (VSN) **(verified against source PDF)**; limma `normalizeQuantiles`; preprocessCore `normalize.quantiles`.

- **medianNormalize** — subtracts the column median from each value (median centering in log space). Standard and correct for log-scale sample centering. **Verdict: Correct.**

- **quantileNormalize (scaled/interpolated rank alignment) — highest-priority check. VALID VARIANT — Informational.** The canonical Bolstad algorithm (sort each column, replace each rank with the cross-column mean at that rank, restore order) requires equal-length columns. The code handles columns with **different numbers of non-missing values** by mapping rank r∈[0..maxValid−1] to a fractional position `(r/(maxValid−1))·(len−1)` and interpolating. This is **mathematically legitimate and matches the limma approach**: limma's `normalizeQuantiles` explicitly "involves interpolating each column of non-missing values out to a full-length vector when computing the mean quantiles" (Smyth, Bioconductor), and limma was the first implementation to allow NAs. It is **not** identical to preprocessCore `normalize.quantiles` (Bolstad's own affy/preprocessCore variant handles NAs differently under a missing-at-random assumption), so results will differ slightly between the two references. **Assessment: valid, not a bug; no bias beyond the well-known limma-vs-preprocessCore NA-handling difference.** Recommendation: state in the docstring that the implementation follows the **limma interpolation** convention, and note that with heavy or non-random missingness (common in proteomics MNAR) neither variant is guaranteed unbiased. **Severity: Informational.**

- **vsnNormalize — VERIFIED against Huber et al. (2002).** Delegates to server-side R VSN and falls back to quantile normalization. The source confirms every claim in the flow: VSN operates on **raw (non-log) intensities** (`yₖᵢ`, Eqs (1),(6)) and outputs a generalized-log / arsinh scale via `h(y) = γ·arsinh(a + by)` (Eq (4)); the paper states that "for large intensities the transformation (4) becomes equivalent to the usual logarithmic transformation" while, "unlike the logarithm, it does not have a singularity at zero" (p. S98), and the difference statistic Δh (Eq (9)) coincides with the log-ratio at high intensity. So the described flow (operate on raw, then write log2/glog2) is exactly the reference behavior, and falling back to quantile is a reasonable degradation. **Verdict: Correct.** Verify the fallback still returns a log-scale output so downstream log2FC math stays valid.

- **VM/null-handling nuances (`raw[i] === DG.FLOAT_NULL` vs `col.isNone(i)`, `col.length` vs `raw.length`).** These are correctness-adjacent: if `raw` is a materialized copy, a mismatch between `col.isNone(i)` and the sentinel `DG.FLOAT_NULL`, or between `col.length` and `raw.length`, could cause a value to be treated as present/absent inconsistently. Not a math error but can corrupt which cells are normalized. **Verdict: Low — align null tests to a single source of truth.**

---

## FILE 5 — imputation.ts

**References:** Perseus/MaxQuant imputation (Tyanova et al. 2016, *Nature Protocols*; Cox & Mann); Troyanskaya et al. (2001) *Bioinformatics* 17(6):520–525 (KNNimpute).

- **imputeMinProb — VERIFIED CORRECT.** Draws from a normal with mean = colMean − downshift·colSD and sd = width·colSD, computed **per column** using that column's own stats, with Perseus defaults **downshift = 1.8, width = 0.3**. This matches Perseus exactly: "the missing values are replaced by random numbers drawn from a normal distribution of 1.8 standard deviation down-shift and with a width of 0.3 of each sample," deployed "in column-wise imputation mode (the default mode)." The per-column loop using `col.stats` is the correct (default) Perseus behavior. **Verdict: Correct.**

- **Box-Muller randomNormal — BUG (rare) — Medium.** `z = sqrt(−2·ln u1)·cos(2π·u2)` is the standard transform, but `Math.random()` returns a value in **[0, 1)**, so `u1` can be exactly 0, giving `ln(0) = −∞` and an `Infinity`/`NaN` draw. Over millions of imputed cells this will eventually occur. **Fix:** draw `u1` in (0,1] (e.g. `u1 = 1 − Math.random()`) or resample while `u1 === 0`. **Severity: Medium.**

- **imputeKnn — DEVIATION from Troyanskaya — Medium.**
  - *Orientation:* the code imputes per row (protein) using nearest rows, i.e. across genes/proteins — this **matches** Troyanskaya, whose KNNimpute finds neighbor **genes** (rows) and imputes gene-wise. Correct orientation.
  - *Distance:* neighbors are ranked by RMSE over shared non-missing columns, normalizing sumSq by `sharedCount` (mean-squared then sqrt). Normalizing by the overlap count to make distances comparable across differing overlaps is a reasonable missing-data adaptation of Euclidean distance (Troyanskaya used plain Euclidean on complete/log-transformed data). Acceptable deviation.
  - *Aggregation — the real issue:* Troyanskaya (2001) imputes each missing value with a **similarity-WEIGHTED average** of the k neighbors, where the contribution of each neighbor gene is weighted by the inverse of its Euclidean distance to the target ("the missing value is replaced with a weighted average of the equivalent value in those neighbors"). The code uses an **UNWEIGHTED** mean of neighbors. This is a genuine deviation from the cited reference and will slightly bias estimates toward more-distant neighbors. **Fix:** weight each neighbor by 1/distance (or a Gaussian kernel) as in the reference. **Severity: Medium.**

- **imputeMean / imputeMedian / imputeZero** — trivial column-statistic substitutions; correct use of column stats. **Verdict: Correct.**

---

## FILE 6 — pca.ts

**References:** Jolliffe, *Principal Component Analysis*; Golub & Van Loan, *Matrix Computations* (Jacobi); Turk & Pentland / Sirovich & Kirby (small-covariance trick).

- **Small-covariance trick — VALID.** With nSamples ≪ nProteins, forming `S = X Xᵀ / (nProteins − 1)` as an nSamples×nSamples matrix and eigendecomposing it is the standard Turk–Pentland/Sirovich–Kirby technique. The eigenvectors of the sample Gram/covariance matrix are exactly the sample-space principal directions; this yields correct PC structure for samples. **Verdict: Correct.**

- **Score scaling `scores = eigenvector[i][k]·sqrt(eigenvalue[k])` — VALID.** Scaling the (unit) sample-eigenvector components by √λ_k produces principal **coordinates** whose variance along PC k equals λ_k — the standard principal-coordinate/biplot scaling. This is a valid representation of sample positions in PC space (equivalent up to sign to projecting centered data onto the loading vectors). **Verdict: Correct.** (Caveat: eigenvector sign is arbitrary; if downstream code compares signs across runs, fix a sign convention.)

- **Centering.** Data is imputed with protein (row) means, then each protein is centered by subtracting its across-sample mean. Centering each feature (protein) is the correct PCA preprocessing. **Verdict: Correct.**

- **Jacobi rotation — VERIFIED CORRECT.** The formulas `θ = (A_qq − A_pp)/(2·A_pq)`, `t = sign(θ)/(|θ| + √(θ²+1))`, `c = 1/√(t²+1)`, `s = t·c` are exactly the standard numerically-stable Jacobi/Rutishauser rotation (identical to Golub & Van Loan §8.5: with μ = (a_qq − a_pp)/(2a_pq), t = 1/(μ + √(1+μ²)) for μ≥0 and the sign-mirrored branch otherwise, c = (1+t²)^(−1/2), s = tc). The `sign(θ)/(|θ|+…)` form is algebraically equivalent to the signed-branch form and picks the smaller rotation |θ|≤π/4, which is the correct, stable choice. Classical Jacobi (zero the largest off-diagonal each sweep), maxIter=100, tol 1e-10 are appropriate. **Verdict: Correct.**

- **Variance explained** — eigenvalues clamped ≥0 (guards tiny negative round-off), pct = λ/Σλ. **Correct.**

- **Divisor nProteins−1 vs nProteins.** Using (nProteins−1) yields the unbiased sample covariance; either choice only uniformly rescales all eigenvalues and leaves the variance-explained percentages and eigenvectors unchanged. **Verdict: Correct — Informational.**

---

## FILE 7 — qc-computations.ts

**References:** MA-plot (Dudoit/Yang; Bland–Altman mean-difference plot); Welford (1962) / Knuth *TAOCP* Vol 2 (online variance).

- **computeMA — CORRECT.** `M = g2Mean − g1Mean` on log2 data = log2 ratio (treatment − control); `A = (g1Mean + g2Mean)/2` = mean log2 intensity. This matches the standard MA-plot definition (M = log2(R/G) = log2 R − log2 G; A = ½(log2 R + log2 G), Dudoit/Yang). Sign convention: positive M = higher in treatment (g2) — clear and conventional. **Verdict: Correct.** (Just document that g1=control, g2=treatment so the M sign is interpreted correctly.)

- **computeCV — CORRECT, including the exponentiation.** The code correctly back-transforms log2 values to raw intensity via `Math.pow(2, v)` **before** computing CV. This is essential: coefficient of variation (sd/mean) is only meaningful on a **linear (raw) scale**, not on log-transformed data (where the mean can be near zero and CV is undefined/meaningless). The Welford update — `delta = val − mean; mean += delta/count; m2 += delta·(val − mean)` using the **new** mean — is exactly the Knuth/Welford recurrence (M2 += delta·delta2 with delta2 = val − updated mean), variance = m2/(count−1), CV = sd/mean. `count < 2 → 0` guards the variance. **Verdict: Correct.**

- **computeLoessTrend — NOT true LOESS — Informational.** The routine is a moving average over a window sorted by A, with `halfWindow = max(round(n·windowFraction/2), 5)`. The moving-average logic is internally sound, but this is **not** LOESS: there is no local weighted regression and no tricube weighting. The docstring's admission that it "approximates loess" understates the difference. **Fix (naming):** rename to `computeMovingAverageTrend` (or similar) to avoid implying local regression; if a true trend line is needed, implement local linear regression with tricube weights. **Severity: Informational.**

- **unpivotIntensities / createMissingnessMatrix / computeMissingBarData** — null/missing counting logic is consistent with treating `isNone`/null cells as missing. **Verdict: Correct.** (Same recommendation as File 4: use one consistent null test.)

---

## Consolidated fix list (priority order)
1. **spc.ts rule4 (High):** re-implement as direction-reversal test on consecutive first differences (13 reversals over 14 points).
2. **enrichment.ts buildEnrichmentDf (Med-High):** align `intersections[j]` to `genes_metadata.ensgs[j]` and reverse-map via `mapping`; do not index against the raw query.
3. **spc.ts rule8 (Med):** add the both-sides-of-centerline requirement to match Nelson (or document the Minitab "either side" convention).
4. **imputation.ts Box-Muller (Med):** guard `u1==0` (use `1 − Math.random()`).
5. **imputation.ts imputeKnn (Med):** switch to inverse-distance-weighted neighbor averaging per Troyanskaya.
6. **differential-expression.ts DEqMS (Low):** for multi-TMT inputs, apply/document the minimum-PSM-count-across-experiments rule (Zhu 2020, D5).
7. **normalization.ts quantile (Info):** document limma-interpolation convention and MNAR caveat.
8. **qc-computations.ts loess (Info):** rename to reflect moving-average, not LOESS.

## Source verification status
Verified **against the original PDFs supplied**: Nelson (1984) run-rule tables; Benjamini & Hochberg (1995) FDR procedure and adjusted-p-value independence; Huber et al. (2002) VSN arsinh/glog scale and raw-intensity input; Zhu et al. (2020) DEqMS PSM/peptide-count covariate and multi-TMT minimum-count recommendation.

Verified via open PDFs / PMC / official docs: Bolstad (2003), Troyanskaya (2001), limma documentation (Smyth/Bioconductor), gprofiler2 docs and g:Profiler API docs.

No primary source remains outstanding for the numerical claims in this audit.
