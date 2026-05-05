# Patch: align `cochranArmitageBasic` degenerate semantics with new SEND `statistics.py`

**Target repo:** `sci-comp` (TS library, package `@datagrok-libraries/sci-comp`)

**Driver:** SENDEX commit `acecc1f8` ("promote modified trend_test_incidence") + `fed85601` ("promote tie-corrected statistics.py") — `statistics_fixed.py` was deleted; `cochranArmitageBasic` Python equivalent is now `statistics.trend_test_incidence(counts, totals)` with default args. The new Python returns `{statistic: 0.0, p_value: 1.0}` from `_degenerate_result(...)` for the cases that used to return `{None, None}`.

**Failing parity cases** (in `methods-validator`, suite `stats`):

- `cochranArmitageBasic.smoke.all_zero_counts` — `counts=[0,0,0,0], totals=[50,50,50,50]` → p̄ = 0
- `cochranArmitageBasic.smoke.all_full_counts` — `counts=[50,50,50,50], totals=[50,50,50,50]` → p̄ = 1

Python returns `{statistic: 0, pValue: 1}`; current TS returns `{statistic: null, pValue: null}`.

## Mapping per case

| Edge case | New Python | Recommended TS |
|---|---|---|
| `p̄ = 0` | `_degenerate_result` → `{0, 1.0}` | **Change** to `{0, 1}` |
| `p̄ = 1` | `_degenerate_result` → `{0, 1.0}` | **Change** to `{0, 1}` |
| `Sxx ≤ 0` (single non-zero group) | `_degenerate_result` → `{0, 1.0}` | **Change** to `{0, 1}` |
| `k < 2` | raises `ValueError` | leave `{null, null}` (full alignment would require throwing — separate API-breaking change) |
| `N = 0` | raises `ValueError` | leave `{null, null}` |

`cochranArmitage` (modified API) is already aligned — uses `degenerate(...)` helper for the same three cases.

## Patch — `src/stats/tests/cochran-armitage.ts`

```diff
 /**
  * Cochran-Armitage trend test for incidence data.
  *
  * Original references: Cochran (1954) Biometrics 10:417–451; Armitage (1955)
  * Biometrics 11:375–386. Modern textbook treatment: Agresti, Categorical
  * Data Analysis, 3rd ed., 2013, §6.4.
  *
- * - `cochranArmitageBasic` — original `statistics_fixed.py` API: returns
- *   `{statistic, p_value}` with `null` on degenerate input. Always uses
- *   default scores `0..k-1`, binomial variance, two-sided alternative.
+ * - `cochranArmitageBasic` — minimal API: returns `{statistic, p_value}`.
+ *   Always uses default scores `0..k-1`, binomial variance, two-sided
+ *   alternative. Returns `{statistic: 0, pValue: 1}` for degenerate inputs
+ *   (p̄ = 0, p̄ = 1, Sxx ≤ 0); returns `{null, null}` for invalid shape
+ *   (k < 2, N = 0).
  * - `cochranArmitage` — extended API with scores, alternative, choice of
```

```diff
 export function cochranArmitageBasic(
   counts: NumericInput,
   totals: NumericInput,
 ): TestResult {
   const c = toFloat64(counts);
   const t = toFloat64(totals);
   const k = c.length;
   if (k < 2 || sum(t) === 0) return {statistic: null, pValue: null};

   const n = sum(t);
   const pBar = sum(c) / n;
-  if (pBar === 0 || pBar === 1) return {statistic: null, pValue: null};
+  if (pBar === 0 || pBar === 1) return {statistic: 0, pValue: 1};

   const d = new Float64Array(k);
   for (let i = 0; i < k; i++) d[i] = i;

   let num = 0; let dT = 0;
   for (let i = 0; i < k; i++) {num += d[i] * c[i]; dT += d[i] * t[i];}
   num -= pBar * dT;

   const dBar = dT / n;
   let Sxx = 0;
   for (let i = 0; i < k; i++) {
     const dev = d[i] - dBar;
     Sxx += t[i] * dev * dev;
   }
   const denomSq = pBar * (1 - pBar) * Sxx;
-  if (denomSq <= 0) return {statistic: null, pValue: null};
+  if (denomSq <= 0) return {statistic: 0, pValue: 1};
```

## Patch — `src/stats/__tests__/fixtures/cochran-armitage.json`

Update the three "edge_none" cases that now have a defined degenerate result. Keep `k = 0`, `k = 1`, `N = 0` as `null/null`.

```diff
       "name": "p̄ = 0",
       "category": "edge_none",
       "inputs": {
         "counts": [0, 0, 0],
         "totals": [50, 50, 50]
       },
       "expected": {
-        "statistic": null,
-        "p_value": null
+        "statistic": 0,
+        "p_value": 1
       }
```

```diff
       "name": "p̄ = 1",
       "category": "edge_none",
       "inputs": {
         "counts": [50, 50, 50],
         "totals": [50, 50, 50]
       },
       "expected": {
-        "statistic": null,
-        "p_value": null
+        "statistic": 0,
+        "p_value": 1
       }
```

```diff
       "name": "single non-zero group (Sxx=0)",
       "category": "edge_none",
       "inputs": {
         "counts": [0, 10],
         "totals": [0, 50]
       },
       "expected": {
-        "statistic": null,
-        "p_value": null
+        "statistic": 0,
+        "p_value": 1
       }
```

## Verification

After applying both patches:

1. In `sci-comp`: `npm test -- cochran-armitage` should still pass (3 edge cases now expect `0, 1` instead of `null, null`).
2. In `methods-validator`: `npm run parity:stats` should clear the two `cochranArmitageBasic` FAILs:
   - `cochranArmitageBasic.smoke.all_zero_counts`
   - `cochranArmitageBasic.smoke.all_full_counts`
3. Suggested commit message in sci-comp:
   ```
   fix(stats): align cochranArmitageBasic degenerate semantics with SEND statistics.py

   Track SENDEX acecc1f8 / fed85601: statistics_fixed.py is gone; the new
   trend_test_incidence returns {statistic: 0, p_value: 1} (not null/null)
   for p̄ = 0, p̄ = 1, and Sxx ≤ 0. The modified CA API was already aligned
   via the degenerate() helper; this brings the basic API in line.
   Validation-error cases (k < 2, N = 0) keep the existing null/null
   return — full alignment to throw is a separate, API-breaking change.
   ```
