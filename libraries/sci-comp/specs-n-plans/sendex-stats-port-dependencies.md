# SENDEX Statistics Port — External Dependencies Analysis

Analysis of the Python dependencies of the five production files listed in
[sendex-stats-port-plan.md](sendex-stats-port-plan.md), which need to be
ported to TypeScript inside the `@datagrok-libraries/sci-comp` library.

**Decision:** **jstat 1.9.6** is used as the mathematical backend (already
declared under `dependencies` in [package.json](../package.json) and used by
[stats-utils.ts](../src/time-series/feature-extraction/stats-utils.ts) and
[moead.ts](../src/optimization/multi-objectives/moead/moead.ts)). jstat covers
all special functions, distributions, and linear algebra; the missing tests are
built on top.

---

## 1. Per-file dependency summary

| File | numpy | scipy.stats | stdlib |
|------|-------|-------------|--------|
| `statistics_fixed.py` | ✓ | `ttest_ind`, `mannwhitneyu`, `fisher_exact`, `norm.cdf`, `spearmanr`, `dunnett` | — |
| `trend_test_incidence_modified.py` | ✓ | `norm.cdf`, `norm.sf` | — |
| `williams_fixed.py` | ✓ | `t.ppf`, `t.cdf` | `dataclasses`, `typing` |
| `williams_tables.py` | — | — | `math.sqrt`, `bisect` |
| `ancova.py` | ✓ | `f.cdf`, `t.cdf` | — |

---

## 2. numpy — what is used

Trivially portable on top of `Float64Array` plus utility helpers.

- **Construction / shape**: `np.array`, `np.asarray`, `np.zeros`, `np.zeros_like`, `np.arange`, `np.concatenate`
- **Masks and NaN**: `np.isnan`, boolean indexing (`a[~np.isnan(a)]`)
- **Reductions**: `np.sum`, `np.mean`, `np.var(ddof=1)`, `np.all`, `np.sqrt`
- **Broadcasting**: outer subtraction `a[:, None] - b[None, :]` (used in `trend_test` for the JT statistic)
- **Linear algebra** (only in `ancova.py`): `np.linalg.lstsq`, `np.linalg.inv`, `np.linalg.pinv`, matrix multiplication `@`, transpose `.T` — covered through jstat (see §3.3).

---

## 3. jstat as a backend for scipy.stats

### 3.1. Distributions — fully covered

| scipy | jstat | Used by |
|-------|-------|---------|
| `norm.cdf`, `norm.sf` | `jStat.normal.cdf(x, 0, 1)`, `1 - jStat.normal.cdf(...)` | trend tests, JT, CA |
| `t.cdf` | `jStat.studentt.cdf(t, df)` | Williams, ANCOVA |
| `t.ppf` | `jStat.studentt.inv(p, df)` | Williams (critical values) |
| `f.cdf` | `jStat.centralF.cdf(f, df1, df2)` | ANCOVA (slope homogeneity) |
| `chi2.cdf` (reserve) | `jStat.chisquare.cdf(x, df)` | — |

### 3.2. Special functions — fully covered

`erf`, `erfc`, `erfcinv`, `gammafn`, `gammaln`, `loggam`, `gammap`, `gammapinv`,
`lowRegGamma`, `betafn`, `betaln`, `betacf`, **`ibeta`**, **`ibetainv`**,
`factorial`, `factorialln`, `combination`, **`combinationln`**, `permutation`.

→ **~250 lines of numerical-recipes code that we no longer need to write.**

### 3.3. Linear algebra — fully covered

jstat provides: `inv`, `det`, `lu`, **`cholesky`**, `gauss_jordan`,
`triaUpSolve` / `triaLowSolve`, `householder`, `jacobi` (eigen).

ANCOVA-OLS is implemented through the normal equations
`β = (XᵀX)⁻¹ Xᵀy`:

```ts
const XtX = matMul(transpose(X), X);
const Xty = matVec(transpose(X), y);
const beta = jStat.cholesky(XtX).then(L => triaLowSolve(L, Xty)...);  // or just jStat.inv(XtX) for p ≤ 6
```

For the dimensions encountered here (`p × p ≤ 6`, typical for ANCOVA), plain
inversion via `jStat.inv` is perfectly reliable.

### 3.4. Numerical methods — for Dunnett

jstat provides: `simpson`, `romberg`, `gauss_quadrature`, `richardson`. These
are used in the Dunnett-1955 implementation (see §4).

---

## 4. What needs to be implemented on top of jstat

| Test | Algorithm | Volume | Comment |
|------|-----------|--------|---------|
| `welchTTest` | Welch–Satterthwaite df + `jStat.studentt.cdf` | ~15 lines | `jStat.ttest` only handles paired/one-sample, so Welch is built ourselves |
| `rank()` (utility) | Sort + average rank for tied groups + `tieCorrection` | ~30 lines | Needed for both Mann-Whitney and Spearman |
| `mannWhitneyU` | Rank the pooled sample → U → normal approximation with tie correction → `jStat.normal.cdf` | ~30 lines | |
| `spearman` | `rank(x)`, `rank(y)` → Pearson → t-approximation → `jStat.studentt.cdf` | ~10 lines | |
| `fisherExact2x2` | `jStat.hypgeom` + `jStat.combinationln`, two-sided via "minlike" | ~30 lines | |
| **`dunnett`** | Dunnett (1955) formula: one-dimensional integration of the standard normal kernel, using `jStat.simpson` or `jStat.gauss_quadrature` | ~80 lines | The only test without direct jstat support |
| OLS wrapper (for ANCOVA) | `(XᵀX)⁻¹Xᵀy` via `jStat.inv` | ~10 lines | |
| `bonferroni` | Trivial iteration | ~5 lines | |

**Total custom code:** ~210 lines vs ~780 if we wrote everything from scratch.
**A 3.7× reduction.**

---

## 5. What changes compared with the pre-jstat analysis

❌ **No longer needed:**
- Lanczos approximation for `gamma` / `lgamma`.
- Regularised incomplete beta `I_x(a,b)` via continued fractions.
- Inverse `ibetainv` (one of the trickiest pieces).
- `erf`, `erfc`.
- All CDF/PPF for `normal`, `studentt`, `centralF`, `chisquare`.
- Cholesky / LU / matrix inversion.
- Numerical quadratures (needed for Dunnett).

✅ **Still to be written:**
- Aggregator tests (Welch, MWU, Spearman, Fisher, Dunnett, Bonferroni) on top of the ready-made distributions and special functions.
- The `rank()` utility with tie correction.
- A thin OLS wrapper for ANCOVA.

---

## 6. jstat caveats

1. **No TypeScript types** — `@ts-ignore: no types` is already in
   [stats-utils.ts:5](../src/time-series/feature-extraction/stats-utils.ts#L5).
   Solution: wrap every jstat call in typed functions inside `stats-utils.ts`
   (or a new `src/distributions.ts`) and expose only a clean TS API.
2. **Accuracy of `studentt.inv` for small df** — historically known rough
   spots. We verify with CSV fixtures against scipy in our tests.
3. **`jStat.ttest`** — does not cover Welch (only paired/one-sample), so we
   build Welch ourselves.
4. **Version 1.9.6** is locked in `package-lock.json` — sanity-check the
   stability of `ibetainv` for small α (it affects critical values).

---

## 7. Architectural decision

**jstat is used as an internal backend.** From the outside, `sci-comp`
exports only a clean TS API:

```
src/stats/
  distributions.ts   # typed wrappers: normalCdf, studentTCdf, fCdf, ibeta, ...
  rank.ts            # ranking utility with tie correction
  tests/
    welch-t.ts
    mann-whitney.ts
    spearman.ts
    fisher-exact.ts
    dunnett.ts
    jonckheere.ts
    cochran-armitage.ts
    williams.ts      # + williams-tables.ts (static tables)
    ancova.ts
  multiple-comparison/
    bonferroni.ts
```

Benefits:
- **Backend-independent** — jstat can be swapped for an in-house implementation
  without breaking the public API.
- **Consistent with the existing style** of [stats-utils.ts](../src/time-series/feature-extraction/stats-utils.ts) (`tdistCdf`, `twoTailPvalue`).
- **Validatable** — tests run against CSV fixtures generated by Python scripts
  derived from [validate_*.py](c:/Users/vmaka/Datagrok/SEND/send-data-browser/backend/services/analysis/validation/).
- **No new external dependencies** beyond the already-included jstat.

---

## 8. Internal Python dependencies (trivial)

- `math.sqrt`, `bisect.bisect_right` — available in JS out of the box.
- `dataclasses` → plain TS interfaces / classes.
- `typing.Optional` → `T | null` or `T | undefined`.
- `from williams_tables import lookup_1971, lookup_1972` — internal cross-file
  import.
