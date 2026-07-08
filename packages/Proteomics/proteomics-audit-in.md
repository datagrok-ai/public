# Proteomics — Numerical Methods Audit

Line-by-line correctness review of every numerical kernel in the package, including
the Nelson SPC rules in `src/analysis/spc.ts` and the library `tTest` / `fdrcorrection`
that differential expression relies on.

## Defects found

### 1. Nelson Rule 4 implemented incorrectly (real bug) — `src/analysis/spc.ts:291-301`

The canonical Nelson Rule 4 (1984) is *"14 points in a row **alternating in direction**
(up-down-up…)"* — a sawtooth / over-control pattern. It is detected from the **sign of
successive differences** `series[i] − series[i-1]`; the centerline is irrelevant to it.

The code instead checks alternation **relative to the centerline** (above/below CL):

```ts
const sides = t.map((v) => v > ctx.mean ? 1 : (v < ctx.mean ? -1 : 0));
```

**Divergence scenario:** a series `μ+5, μ+3, μ+5, μ+3, …` (14 points) is a true up-down
zigzag but is entirely **above** the CL. Canonical Rule 4 must trip; this implementation
does not (all `sides[i] === 1`, returns `false` on the first pair). Conversely, a slow
trend that repeatedly crosses the CL without direction oscillation is falsely flagged.

Mitigating factor: Rule 4 is **disabled by default** (`nelson_4: false`,
`src/analysis/spc.ts:244-245`), so it does not affect results until a user enables it.
As a named canonical algorithm it is still a defect. The comment
(`// 14 consecutive points alternating above/below CL`) shows the wrong semantics were
coded deliberately.

**Suggested fix:** detect the sign of consecutive differences and require it to alternate
across the last 14 points.

### 2. Quantile normalization yields NaN for a column with a single non-null value — `src/analysis/normalization.ts:90-97`

When `n === 1` (a sample column with exactly one non-null protein) and `maxValid > 1`:

```ts
const globalPos = maxValid === 1 ? 0 : (r / (n - 1)) * (maxValid - 1); // r=0, n=1 → 0/0 = NaN
```

`globalPos = NaN` → `NaN` is written to that value. Symmetrically, the `rankMeans`
computation at `src/analysis/normalization.ts:74` also produces `NaN` when `maxValid === 1`.
The case is degenerate (real proteomics has thousands of proteins), but it produces an
incorrect value instead of a no-op.

**Suggested fix:** guard with `if (n < 2) continue;` in the assignment loop (and skip the
`maxValid < 2` case up front).

## Notes (not bugs, but worth knowing)

- **MinProb is non-deterministic** — `src/analysis/imputation.ts:9-14` uses `Math.random()`
  with no seed, so imputation is not reproducible run-to-run, which means downstream DE
  results drift slightly. This matches Perseus behavior, but should be documented.
- **The MA-plot "LOESS" trend is a moving average**, not true loess
  (`src/viewers/qc-computations.ts:57-106`). Honestly labeled "approximates loess";
  adequate for a visual trend line, but it is not local regression.

## Verified correct

- **Median centering** (`src/analysis/normalization.ts:11-31`) — subtracts the column
  median; correct.
- **Box–Muller** (`src/analysis/imputation.ts:9-14`) — `√(−2 ln u₁)·cos(2π u₂)`, a correct
  standard normal; downshift/width 1.8/0.3 match the Perseus defaults.
- **kNN imputation** — distance `√(Σ diff²/sharedCount)` over shared observed columns only;
  zero-fills never enter the distance; neighbor averaging uses observed values only. Correct.
- **PCA (Jacobi)** (`src/analysis/pca.ts`) — derivation verified: sample-space covariance
  `S = XXᵀ/(n−1)`, eigenvectors = U from the SVD, scores `= U·√λ`. These are the true PC
  coordinates **up to a global factor** `1/√(nProteins−1)` that is identical for all axes,
  so the scatter-plot shape is unaffected. Variance-explained `λₖ/Σλ = σₖ²/Σσ²` is correct.
  Per-protein centering is consistent with the mean-imputation (imputed cells become exactly
  0 after centering).
- **DE / t-test** — log2FC = difference of means (data is in log space, correct); library
  `tTest` is Welch's (`devEqual=false`, separate variances, Welch–Satterthwaite df); the
  numerator/denominator direction is consistent with the dialog hint; BH-FDR is delegated to
  `@datagrok-libraries/statistics`.
- **CV** (`src/viewers/qc-computations.ts:110-161`) — Welford's algorithm on de-logged
  `2^log2` intensities; correctly avoids catastrophic cancellation at large intensities.
- **MA-plot** — `M = ḡ₂ − ḡ₁`, `A = (ḡ₁+ḡ₂)/2`, correct.
- **Nelson Rules 1, 2, 3, 5, 6, 7, 8** — checked against the canon (1 point beyond 3σ;
  9 on one side; 6 monotone; 2/3 beyond 2σ; 4/5 beyond 1σ; 15 within 1σ; 8 beyond 1σ) —
  implemented correctly. The iterative baseline (3σ pruning, ≤2 iterations, sample sd) is
  correct.
- **limma / DEqMS / VSN** (R) — correct use of `eBayes` / `spectraCounteBayes` / `justvsn`;
  results are copied back by the `row` key, not by position (load-bearing for DEqMS, whose
  output is sorted by significance).

## Summary

Of the numerical methods, one real defect — **Nelson Rule 4** (wrong semantics: alternation
of direction vs. alternation about the centerline), disabled by default — plus a degenerate
NaN edge in quantile normalization. Everything else is mathematically correct.
