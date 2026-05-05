# Patch: align `fisherExact2x2.oddsRatio` degenerate semantics with SEND `incidence_exact_test`

**Target repo:** `sci-comp` (TS library, package `@datagrok-libraries/sci-comp`)

**Driver:** SENDEX `fed85601`/`acecc1f8`. SEND `incidence_exact_test` returns `null` for the odds ratio when either off-diagonal cell is zero (b = 0 or c = 0), citing "inf not JSON-serializable". TS currently returns `Infinity`/`0`. The validator surfaces ~9 sweep FAILs per run on degenerate tables.

## Failing parity cases (a sample)

```
fisherExact2x2.sweep.sweep_002_a1_b0_c3_d5  -- oddsRatio: py=None  ts=Inf
fisherExact2x2.sweep.sweep_004_a3_b0_c5_d2  -- oddsRatio: py=None  ts=Inf
fisherExact2x2.sweep.sweep_020_a5_b2_c0_d4  -- oddsRatio: py=None  ts=Inf
fisherExact2x2.sweep.sweep_036_a5_b5_c0_d1  -- oddsRatio: py=None  ts=Inf
…
```

(The smoke `diagonal` and `zero_cell` cases hit the same divergence.)

## Patch — `src/types.ts`

Allow `null` on `oddsRatio` (breaking type change for downstream consumers — they must check):

```diff
 export interface FisherResult {
-  oddsRatio: number;
+  /**
+   * `null` when either off-diagonal cell is zero (b = 0 or c = 0) — matches
+   * SEND `incidence_exact_test`, which avoids `±Infinity` because it is not
+   * JSON-serialisable. Callers should inspect the table's incidence rates
+   * directly in this case.
+   */
+  oddsRatio: number | null;
   pValue: number;
   pGreater: number;
   pLess: number;
 }
```

## Patch — `src/stats/tests/fisher-exact.ts`

```diff
- * `oddsRatio` is `(a·d)/(b·c)`. Returns `Infinity` when `b·c = 0` and `a·d > 0`,
- * `0` when `a·d = 0` and `b·c > 0`. Throws on negative inputs.
+ * `oddsRatio` is `(a·d)/(b·c)`, rounded to 6 decimals to match SEND's
+ * `round((a*d)/(b*c), 6)` in `incidence_exact_test`. Returns `null` when
+ * either off-diagonal cell is zero (b = 0 or c = 0); SEND avoids `±Infinity`
+ * because it is not JSON-serialisable. Throws on negative inputs.
  */
 export function fisherExact2x2(table: ReadonlyArray<ReadonlyArray<number>>): FisherResult {
```

```diff
-  let oddsRatio: number;
-  if (b * c === 0) oddsRatio = a * d > 0 ? Infinity : 0;
-  else oddsRatio = (a * d) / (b * c);
+  let oddsRatio: number | null;
+  if (b === 0 || c === 0) oddsRatio = null;
+  else oddsRatio = Math.round((a * d) / (b * c) * 1e6) / 1e6;
```

`incidenceExactBoth` (`boschloo-exact.ts`) already returns `oddsRatio` from the inner `fisherExact2x2` call, so the change propagates there with no further edit. The current `null` check on `boschloo.pValue` (the SEND policy fallback) keeps the `pValue` field a `number`.

## Patch — `src/stats/__tests__/fixtures/fisher-exact.json`

Any `b = 0` / `c = 0` fixtures that currently expect `Infinity` or `0` for `oddsRatio` need to be regenerated against `incidence_exact_test`. Regenerate via the existing `scripts/generate-fixtures.py` pipeline or update by hand:

```diff
       "expected": {
-        "oddsRatio": Infinity,
+        "oddsRatio": null,
         "pValue": …
       }
```

Non-degenerate fixtures should also be regenerated to pick up 6-decimal rounding (drift up to ~5e-7).

## Verification

After applying:

1. In `sci-comp`: `npm test -- fisher-exact` should pass against regenerated fixtures.
2. In `methods-validator`: `npm run parity:sweep` should clear the ~38 sweep `oddsRatio` FAILs (29 rounding + 9 None/Inf) and the 3 smoke ones (`diagonal`, `zero_cell`, `small_counts`). Note: the validator's `oddsRatio` tolerance was already loosened to `1e-6` to absorb SEND's rounding policy independently of this patch — both fixes together restore strict parity.
3. Suggested commit message in sci-comp:
   ```
   fix(stats): align fisherExact2x2 oddsRatio with SEND policy

   - Round (a*d)/(b*c) to 6 decimals (matches SEND round(..., 6)).
   - Return null instead of ±Infinity when b = 0 or c = 0 — SEND avoids
     non-JSON-serialisable infinities. Type changes from `number` to
     `number | null` (breaking).
   - Regenerate fixture expectations for degenerate and rounded cases.
   ```

## Open question — typing impact

`oddsRatio: number | null` is an API break. If preferable, `NaN` would keep the `number` type but would still need fixture regeneration and the validator-side tolerance would not catch it (NaN vs null is a no-match either way). Recommended path is `null` to match SEND exactly; if backwards compatibility is critical, consider a major version bump.
