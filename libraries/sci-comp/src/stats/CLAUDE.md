# CLAUDE.md — stats

Statistical tests, distributions, and multiple-comparison corrections. No dependency on `datagrok-api`.

## Architecture

```
src/stats/
  index.ts                        # Public entry — re-exports tests, types, distributions, multiple-comparison helpers
  README.md                       # Per-method tables, code snippets, run instructions for all 14 methods
  types.ts                        # NumericInput, Alternative, TestResult, SpearmanResult, FisherResult
  distributions.ts                # Typed wrappers around jstat (normal, t, F, chi², hypergeom, special functions)
  internal/
    normalize.ts                  # NumericInput → Float64Array, NaN stripping, mean / variance / std / sum
    rank.ts                       # Rank assignment with average-rank tie handling
    matrix.ts                     # Linear algebra helpers (matrix inverse via jstat)
    random.ts                     # mulberry32 seedable PRNG + shuffleInPlace (Fisher–Yates)
  tests/                          # One file per statistical test
    welch-t.ts, mann-whitney.ts, hedges-g.ts, spearman.ts, fisher-exact.ts
    welch-pairwise.ts, dunnett.ts, cochran-armitage.ts, ancova.ts
    williams.ts, williams-tables.ts
    jonckheere.ts                 # Approximate / permutation / exact, ±continuity, tie-corrected variance
    boschloo-exact.ts             # Unconditional exact: Fisher one-sided p as test stat, sup over nuisance π via grid + golden-section
  multiple-comparison/
    bonferroni.ts                 # bonferroniCorrect — multiplicity adjustment
  __tests__/
    helpers.ts                    # loadFixture, expectClose
    fixtures/                     # 18 JSON fixtures (179 cases), committed; produced out-of-band from scipy / numpy / R
  examples/                       # 9 runnable examples + _helpers.ts: npx tsx src/stats/examples/*.ts
```

## Key design patterns

- **Fixture-driven tests**: every method is validated against a JSON fixture in `__tests__/fixtures/`. Fixtures come from authoritative reference implementations, not hand-authored values:
  - scipy / numpy (most tests)
  - R `clinfun::jonckheere.test` → `jonckheere-clinfun.json` (144 cases)
  - R `PMCMRplus` permutation → `jonckheere-pmcmr.json` (144 cases)
  - Python `regressionpack` → `jonckheere-rp-{approximate,exact}.json`
  - `scipy.stats.boschloo_exact` → `boschloo-exact.json` (31 cases, hand-picked + 60 randomised)

  Fixtures are produced out-of-band against the reference tool and checked in. When you add a test, generate the fixture from the same authoritative source and commit the JSON.

- **`expectClose` gotcha** (`__tests__/helpers.ts`): the comparison uses `diff >= tol` (strict greater-or-equal), so `tol = 0` will fail on identical values. Pick `tol > 0` — typically `1e-12` for machine-precision matches.

- **Distribution wrappers** (`distributions.ts`): single source for normal / t / F / χ² / hypergeometric pdf-cdf-quantile and special functions, all typed wrappers around `jstat`. Tests and methods import from here, not from `jstat` directly — keeps the dependency swappable.

- **`NumericInput` normalization** (`internal/normalize.ts`): public APIs accept `number[]` or any TypedArray; the internal pipeline converts to `Float64Array` once and strips NaNs. Methods should not re-normalize.

## Adding a new test

1. Add `tests/<name>.ts` with the public function and TSDoc.
2. Re-export from `index.ts`.
3. Add a fixture generator (Python or R) and check the resulting `__tests__/fixtures/<name>.json` in.
4. Add `__tests__/<name>.test.ts` driven by the fixture via `loadFixture` + `expectClose`.
5. Update `README.md` (method table + code snippet) and the architecture tree above.
