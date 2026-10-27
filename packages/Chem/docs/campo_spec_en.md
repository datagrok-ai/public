# Specification: Correlation-aware (caMPO) weighting mode in MPO Profiles

Status: v8 (EN), implemented (§11). v8 allows `Correlation-aware` with every aggregation, with a warning balloon instead of the banner (§8). v5 rebased the spec on the existing MPO implementation (§1.1); v6 is a review against that code (`MpoProfileEditor`, `MpoProfileCreateView`, `MpoContextPanel`) that fixes how recomputation is triggered and adds edge cases; v7 closes the gaps between the controller code and the text (messages, tab activation, `d = 1`, no dataset, mode transitions, in-place edits). Correlation plot checked on a live instance (`campo_corrplot_check.js`, 2026-10-26).
Method reference: Garellick et al., *caMPO: Correlation-aware multiparameter optimization*, ChemRxiv 2026, doi:10.26434/chemrxiv.15008121/v1.

Decisions taken:

- the matrix `R` is computed by our own Pearson function on complete rows (§2.5); the built-in Correlation plot shows the same matrix and is not a source of data;
- matrix inversion uses Cholesky decomposition;
- the weight formula follows the paper, with negative values clipped and the rest renormalised;
- the Priority column is kept.

Decisions proposed in v5, to confirm:

- desirability arrays come from the existing `mapColumnDesirability` through the `MpoCalculator` cache, not from new code (§2.3);
- the property's existing `weight` always holds the **effective** weight, so the score, the MPO Score dialog, old clients and other packages need no change; caMPO inputs are stored next to it (§6);
- the mode viewers are one panel docked as a tab over the grid of a detached table view that fills the right half of the profile view (§3; replaced the tab control of v5–v7, §11);
- RESET keeps its current meaning, "revert to the last saved profile" (§4);
- the **Weighting** control exists only in the profile view; the MPO Score dialog and the profile preview use the stored weights (§1.2).

Decisions proposed in v6, to confirm:

- `editor.setProfile` fires no change event, so the view runs stage A itself after `attachLayout` and after RESET; `onChanged` covers only the user's edits (§4);
- a recomputation that does not change the weights does not mark the profile as modified; one that does (another dataset, edited data) does (§4);
- two properties mapped to the same column: the second one is excluded from `R` (§7);
- a property whose column is created by its compute function is picked up when the column appears (§4);
- on a stage A error the **Correlations** and **Weights** tabs are removed and the weights keep their last values (§4, §7).

Decisions proposed in v7, to confirm:

- stage A warnings (λ, few rows, rows less than half the dataset) are kept and shown by every successful stage B; an error replaces them (§4);
- the Weights tab becomes current only when the mode is switched on or RESET restores it, not when a saved profile is opened (§4);
- `hide()` guards `tabs.clear()` with a flag instead of relying on what `onTabRemoved` reports for it (§3.3);
- edits to the dataset's values do not recompute the weights in version 1 (§4);
- `Equal → Correlation-aware` keeps the stored priorities; a `weightColumn` property brings its static weight into the priority (§5.1).

### Reference artifacts (same folder)

| File | Purpose |
|---|---|
| `campo_demo.html` | Interactive demo of the feature: layout, weighting modes, priorities, pins, tooltips, viewers, tabs, aggregation warning (dropped, §8). Open in a browser. The reference for behaviour and look, with the differences listed below. |
| `campo_bench.html` | Verification bench for the computation: presets from the paper, generated data, intermediate values (s_i, raw, base, final, R⁻¹), the §9 test suite. |
| `campo_corrplot_check.js` | Datagrok script that checks what Correlation plot exposes: `getCorrelation` vs reference Pearson, missing values, hidden tab, filter, performance, close events. Its results are in §10. |

The math in both HTML files is the code from §2, ported 1:1 to JavaScript. Differences from the Datagrok implementation:

- data is generated: 400 compounds with drugs_props-like columns and fixed desirability curves (not editable in the demo);
- tabs and inputs are imitated with HTML; SAVE does not persist anything; the Data-driven method is not implemented;
- the demo's tab strip is at the bottom; `ui.tabControl` puts it on top (§3);
- the demo's RESET sets defaults; the real RESET reverts to the saved profile (§4);
- the demo offers three aggregations; the real list is `Average`, `Sum`, `Product`, `Geomean`, `Min`, `Max` (§8).

---

## 1. Goal and scope

The MPO profile view (**Manual** method) gets a choice of how weights are set, **Weighting**: `Equal` / `Manual` / `Correlation-aware`. `Manual` is today's behaviour and the default. In `Correlation-aware` mode, weights are computed from the correlations between desirability columns: redundant properties get less weight, independent ones get more. An expert can adjust the result in two ways: a priority (a multiplier on the weight) and pinning a specific weight value.

Out of scope: the Data-driven method (unchanged; **Weighting** is hidden there, as **Aggregation** already is).

### 1.1. Existing implementation

| File | What it does today | Change |
|---|---|---|
| `libraries/statistics/src/mpo/mpo-types.ts` | `DesirabilityProfile`, `PropertyDesirability` (`weight` 0–1, optional `weightColumn`, `missingValues`), `WEIGHTED_AGGREGATIONS` | `MpoWeighting`, `weighting`, `priority`, `pinnedWeight` (§6) |
| `libraries/statistics/src/mpo/mpo.ts` | `hoistColumns` + `mapColumnDesirability` (value → desirability per column, with missing values, log scale, inversion, qnum, categories), `reduceMpo` (aggregation), `MpoCalculator` (per-column cache) | `MpoCalculator.desirabilities()` (§2.3) |
| `libraries/statistics/src/mpo/campo.ts` | — | new: the math of §2 (pure functions, no UI) |
| `libraries/statistics/src/mpo/mpo-profile-editor.ts` | `MpoProfileEditor`: header Property / Weight / Desirability, one row per property, `aggregationInput`, column mapping (`columnMapping`), design mode (rename, add, delete, compute function), weight-from-column toggle, `onChanged` + `MPO_SCORE_CHANGED_EVENT` | weighting mode, Priority column, read-only weights, pins, name tooltip, `setWeights()` (§5) |
| `libraries/statistics/src/mpo/editors/mpo-property-row.ts` | `createMpoRow`: property cell, weight cell (`ui.input.float`, 0–1, `#0.000`), curve editor | priority cell, weight cell variants (§5) |
| `packages/Chem/src/mpo/mpo-create-profile.ts` | `MpoProfileCreateView`: title, description, Method / Dataset / Aggregation with hints, Save / Reset in the ribbon; content is `ui.splitH([editor, grid])`; Data-driven mode docks EDA viewers into a hidden `tableView` | Weighting input, tab control in the right half, caMPO controller (§3, §4) |
| `packages/Chem/src/mpo/campo-controller.ts`, `campo-viewers.ts` | — | new (§3.3, §4) |
| `packages/Chem/src/mpo/mpo-context-panel.ts` | Context panel: score histogram, best / worst compounds; owns the `MpoCalculator` | exposes its calculator to the controller |
| `packages/Chem/src/mpo/utils.ts` | `computeMpo` (maps property → column, tags columns, runs the calculator), `createProfileForDf` (up to `MAX_MPO_PROPERTIES` = 10, weight 1) | none |
| `packages/Chem/src/mpo/mpo-profile-store.ts`, `databases/mpo/` | Profiles in the `mpo.profile` table: `name`, `description`, `aggregation`, `format_version`, `properties` (jsonb) | `weighting` column, migration (§6) |
| `packages/Chem/src/analysis/mpo.ts` | `MpoProfileDialog`, "MPO Score" on a table: picks a saved profile, optional design mode, runs `computeMpo` | read-only weights for `Equal` / `Correlation-aware` profiles (§1.2) |
| `packages/Chem/src/mpo/mpo-profile-handler.ts` | Profile preview in the context panel (`MpoProfileEditor` in preview mode) | shows the stored weights, nothing else |
| `packages/Chem/src/ai-tools/mpo.ts` | AI functions of the profile view, including `setMpoPropertyWeight` | §1.2 |
| `packages/EDA/src/probabilistic-scoring/prob-scoring.ts`, `packages/PowerGrid/src/vlaaivis/editor.ts` | Use `MpoProfileEditor` / `createMpoRow` | none: the new behaviour is off unless a host turns it on |

### 1.2. Hosts other than the profile view

- **Score computation** does not change. `reduceMpo` with `Average` computes `Σ wᵢ·dᵢ / Σ wᵢ` over the properties a row has, so weights do not have to sum to 1, and the effective weights stored in `weight` (§6) give the same score everywhere.
- **MPO Score dialog** (`analysis/mpo.ts`) and **profile preview** (`mpo-profile-handler.ts`): no Weighting control. For a profile with `weighting` `equal` or `correlation`, the weight inputs are read-only with tooltip `Set by <Equal | Correlation-aware> weighting. Change it in MPO Profiles.` Weights are not recomputed on the dialog's table (§10, question 5).
- **AI functions** (`ai-tools/mpo.ts`): `getMpoProfile` reports `weighting`, priorities, pins and, in `Correlation-aware` mode, raw and final weights. `setMpoPropertyWeight` refuses in `Equal` / `Correlation-aware` mode with a message naming the mode. New functions: `setMpoWeighting(view, mode)`, `setMpoPriority(view, property, priority)`, `pinMpoWeight(view, property, weight | null)`.
- **EDA and PowerGrid** use the editor with its defaults: no Weighting, no Priority, editable weights, exactly as today.

---

## 2. Computation

### 2.1. Notation

| Symbol | Meaning |
|---|---|
| `n` | number of dataset rows |
| `d` | number of properties in the profile |
| `m` | number of properties included in `R` |
| `D_k` | desirability column of property `k` (values in [0, 1]) |
| `R` | `m×m` correlation matrix of the `D_k` columns (Pearson) |
| `R⁻¹` | inverse of `R` |
| `s_i` | sum of the elements of row `i` of `R⁻¹` |
| `S` | sum of all elements of `R⁻¹` (= `Σ s_i`) |
| `prio_i` | property priority (default 1) |
| `p_i` | pinned weight of the property (if pinned) |

### 2.2. Pipeline

The computation is split into two stages, so that editing a priority or a pin does not recompute the matrix.

```
Stage A (heavy)                                             Stage B (light)
MpoCalculator.desirabilities → D_k → complete rows,         base weights × priorities
active columns → Pearson → R → R⁻¹ → base weights        →  pins → final weights → editor.setWeights
(the same rows feed the Correlation plot, display only)

triggers: dataset change, property set or column            triggers: priority change,
mapping change, desirability change, mode switched on       pin change
```

All computation runs on typed arrays, as the existing MPO code does (`hoistColumns` reads `getRawData()` once per column).

### 2.3. Desirability columns

No new mapping code. `mapColumnDesirability` (`libraries/statistics/src/mpo/mpo.ts`) already turns a column into a `Float32Array` of desirabilities, with everything the score uses: missing-value strategy, log scale, inversion, qualifiers of qnum columns, categories. Its `D[i]` is `NaN` where the row does not contribute (strategy `exclude` or `skip`, unknown category) and the default score where the strategy is `default`. These are the values the score aggregates, so they are the values to correlate.

`MpoCalculator` caches these arrays per `df.id|column` and validates them with `desirabilityKey` (curve, missing values, scale; not the weight). A new method returns them from the same cache, so stage A does not map a column the score has already mapped:

```ts
// libraries/statistics/src/mpo/mpo.ts, class MpoCalculator
/// Desirability arrays of the given columns (tagged with 'desirabilityTemplate'), from the same cache compute() uses.
desirabilities(dataFrame: DG.DataFrame, columns: DG.Column[]): ColumnDesirability[] {
  const rowCount = dataFrame.rowCount;
  return hoistColumns(dataFrame, columns).map((h) => {
    const lineKey = desirabilityKey(h.template);
    const key = `${dataFrame.id}|${h.col.name}`;
    const hit = this.cache.get(key);
    if (hit && hit.lineKey === lineKey && hit.D.length === rowCount)
      return {D: hit.D, state: hit.state};
    const m = mapColumnDesirability(h, rowCount);
    this.cache.set(key, {lineKey, D: m.D, state: m.state});
    return m;
  });
}
```

`compute()` keeps its own loop or calls this method; either way one cache serves both. The controller uses the context panel's calculator (`MpoContextPanel` owns it today), and gets the columns the way `computeMpo` does: property → `columnMapping[property] ?? property` → `applyDesirabilityTags`. A property with no column in the dataset gets no `D_k` (§7).

The cache is an optimisation, not something to rely on. The view creates `MpoContextPanel` lazily in `setupGridAndContextPanel`, drops it in `closeContextPanel`, and the panel calls `calculator.release()` when the user leaves the view. So the controller reads `view.mpoContextPanel?.calculator` on every stage A run and keeps no reference; with no panel it uses a fresh `MpoCalculator`.

Known limitation: under the `default` strategy every missing value becomes the same score. With many missing values this lowers the column's variance and shifts its correlations towards that constant. These are still the values the score aggregates, so they are correlated as they are.

> Verified on a live instance: for `double` columns `getRawData()` returns a `Float32Array`, and `raw[i] === DG.FLOAT_NULL` finds every missing value (112 of 112). `hoistColumns` relies on this (`valNull`); keep it as a unit test.

### 2.4. Complete rows and constant columns

Two steps are performed before `R` is computed.

**Complete rows only.** `R` is computed only on rows where every `D_k` is a number (listwise deletion). Correlations on different row subsets (pairwise deletion) can produce a matrix that is not positive definite — Cholesky decomposition would then fail. The Correlation plot gets the same complete rows: it uses pairwise deletion (verified: with 5% missing values it matches the pairwise reference to 3e-9 and differs from the listwise one by 0.028), so on the full data it would show a matrix different from the one the weights come from.

The score is computed on a different set of rows: under the `skip` strategy `reduceMpo` averages a row over the properties it has. The weights come from the rows where every property is present. This is accepted for version 1. When complete rows are less than half the dataset, the form shows a warning (§5.4): `Weights are derived from <rows> of <n> rows that have every property.`

If there are not more complete rows than properties (`rows ≤ d`), the function returns at once with no active columns, so the controller reports `No rows with all properties present` (`rows = 0`) or `Not enough data…` instead of reading every column as constant (with `rows = 0` the mean is `NaN`).

**Constant columns.** If a desirability column has zero variance (for example, every row is 1.0), its correlation is undefined (the viewer returns `NaN`). Such a property is excluded from `R` and from the viewer, gets weight 0 and status `dropped`; the reason is shown in the name tooltip.

```ts
interface Prepared {
  cols: Float32Array[];  // complete rows of the active columns — input for pearsonMatrix
  frame: DG.DataFrame;   // the same values as columns — input for Correlation plot (display)
  active: Int32Array;    // indices (in the profile's property order) of columns included in R
  constant: Int32Array;  // indices excluded because variance == 0
  rows: number;          // number of complete rows
}

/** D[k] is null for a property with no column in the dataset. */
function prepareDesirabilityFrame(names: string[], D: (Float32Array | null)[], n: number): Prepared {
  const d = D.length;

  // complete-case mask
  const mask = new Uint8Array(n).fill(1);
  for (let k = 0; k < d; k++) {
    const a = D[k];
    if (!a) continue;
    for (let i = 0; i < n; i++) if (Number.isNaN(a[i])) mask[i] = 0;
  }
  let rows = 0;
  for (let i = 0; i < n; i++) rows += mask[i];
  if (rows <= d)
    return {cols: [], frame: DG.DataFrame.create(0), active: new Int32Array(0), constant: new Int32Array(0), rows};

  const activeList: number[] = [];
  const constantList: number[] = [];
  const arrays: Float32Array[] = [];
  const cols: DG.Column[] = [];

  for (let k = 0; k < d; k++) {
    const a = D[k];
    if (!a) continue;
    const out = new Float32Array(rows);
    let sum = 0;
    for (let i = 0, r = 0; i < n; i++) if (mask[i]) { out[r++] = a[i]; sum += a[i]; }
    const mean = sum / rows;
    let v = 0;
    for (let r = 0; r < rows; r++) { const t = out[r] - mean; v += t * t; }

    if (v > 1e-12) {
      activeList.push(k);
      arrays.push(out);
      cols.push(DG.Column.fromFloat32Array(names[k], out));
    } else {
      constantList.push(k);
    }
  }
  return {
    cols: arrays,
    frame: DG.DataFrame.fromColumns(cols),
    active: Int32Array.from(activeList),
    constant: Int32Array.from(constantList),
    rows,
  };
}
```

### 2.5. Matrix R: our own Pearson

`R` is computed from `prepared.cols` (complete rows, active columns) into a `Float64Array` (row-major, `m×m`). The Correlation plot is not used as a source of `R` (§10, questions 1 and 2):

- the viewer has no matrix to read: its `table` setting is the name of the bound table (`null` for a viewer built in code);
- `DG.CorrelationPlot.getCorrelation(c1, c2)` returns the right value (matches the reference Pearson to 4e-9, also in a hidden tab and right after the data changes), but it recomputes the pair on every call, with no cache: the full 20×20 matrix on 10⁵ rows takes ~1.5 s (measured). The code below does 190 dot products of length 10⁵ for the same matrix — tens of milliseconds by estimate. The §9 performance test covers all of stage A (desirability mapping on a cold cache, complete rows, `R`, inversion), so its 300 ms limit leaves room for more than this function.

On the same complete rows the viewer shows exactly this matrix, so what the user sees in **Correlations** is what the weights are computed from.

```ts
/** Pearson correlation matrix of columns with no missing values, all of length `rows`. */
function pearsonMatrix(cols: Float32Array[], rows: number): Float64Array {
  const m = cols.length;
  const centered = cols.map((a) => {
    let sum = 0;
    for (let r = 0; r < rows; r++) sum += a[r];
    const mean = sum / rows;
    const c = new Float64Array(rows);
    for (let r = 0; r < rows; r++) c[r] = a[r] - mean;
    return c;
  });
  const norm = centered.map((c) => {
    let s = 0;
    for (let r = 0; r < rows; r++) s += c[r] * c[r];
    return Math.sqrt(s);
  });

  const R = new Float64Array(m * m);
  for (let p = 0; p < m; p++) {
    R[p * m + p] = 1;
    for (let q = p + 1; q < m; q++) {
      const a = centered[p], b = centered[q];
      let s = 0;
      for (let r = 0; r < rows; r++) s += a[r] * b[r];
      const v = Math.max(-1, Math.min(1, s / (norm[p] * norm[q])));
      R[p * m + q] = R[q * m + p] = v;
    }
  }
  return R;
}
```

Constant columns are excluded in §2.4, so `norm[p] > 0`. The cost is `m(m−1)/2` dot products of length `rows`.

### 2.6. Inversion: Cholesky decomposition

`R` is symmetric and (on complete rows with no constant columns) positive definite, so we use Cholesky decomposition `R = L·Lᵀ`. It is about twice as fast as LU and more stable, and a failed decomposition immediately signals that the matrix is not positive definite or is nearly singular.

```ts
/** Lower-triangular L (row-major), or null if R is not positive definite. */
function cholesky(A: Float64Array, m: number): Float64Array | null {
  const L = new Float64Array(m * m);
  for (let i = 0; i < m; i++) {
    for (let j = 0; j <= i; j++) {
      let s = A[i * m + j];
      for (let k = 0; k < j; k++) s -= L[i * m + k] * L[j * m + k];
      if (i === j) {
        if (s <= 1e-12) return null;
        L[i * m + i] = Math.sqrt(s);
      } else {
        L[i * m + j] = s / L[j * m + j];
      }
    }
  }
  return L;
}

/** Solves L·Lᵀ·x = b. */
function choleskySolve(L: Float64Array, m: number, b: Float64Array): Float64Array {
  const y = new Float64Array(m);
  for (let i = 0; i < m; i++) {                     // L·y = b
    let s = b[i];
    for (let k = 0; k < i; k++) s -= L[i * m + k] * y[k];
    y[i] = s / L[i * m + i];
  }
  const x = new Float64Array(m);
  for (let i = m - 1; i >= 0; i--) {                // Lᵀ·x = y
    let s = y[i];
    for (let k = i + 1; k < m; k++) s -= L[k * m + i] * x[k];
    x[i] = s / L[i * m + i];
  }
  return x;
}

/** Full inverse, column by column. */
function choleskyInverse(L: Float64Array, m: number): Float64Array {
  const inv = new Float64Array(m * m);
  const e = new Float64Array(m);
  for (let c = 0; c < m; c++) {
    e.fill(0); e[c] = 1;
    const col = choleskySolve(L, m, e);
    for (let r = 0; r < m; r++) inv[r * m + c] = col[r];
  }
  return inv;
}
```

**Regularisation.** If `cholesky` returns `null` or the condition number exceeds `κ > 10⁶`, shrinkage `R_λ = (1−λ)·R + λ·I` is applied; `λ` is tried in order `[1e-3, 1e-2, 5e-2, 0.1]` until the decomposition succeeds. The fact of regularisation and the value of `λ` are shown in the form message (§5.4).

```ts
function regularise(R: Float64Array, m: number, lambda: number): Float64Array {
  const out = new Float64Array(m * m);
  for (let i = 0; i < m * m; i++) out[i] = (1 - lambda) * R[i];
  for (let k = 0; k < m; k++) out[k * m + k] += lambda;
  return out;
}

function norm1(A: Float64Array, m: number): number {
  let max = 0;
  for (let j = 0; j < m; j++) {
    let s = 0;
    for (let i = 0; i < m; i++) s += Math.abs(A[i * m + j]);
    if (s > max) max = s;
  }
  return max;
}

interface InverseResult { Rinv: Float64Array; lambda: number; kappa: number; }

function robustInverse(R: Float64Array, m: number): InverseResult | null {
  for (const lambda of [0, 1e-3, 1e-2, 5e-2, 0.1]) {
    const A = lambda === 0 ? R : regularise(R, m, lambda);
    const L = cholesky(A, m);
    if (!L) continue;
    const Rinv = choleskyInverse(L, m);
    const kappa = norm1(A, m) * norm1(Rinv, m);
    if (kappa <= 1e6) return { Rinv, lambda, kappa };
  }
  return null;
}
```

### 2.7. Base weights (element-wise)

The paper's formula `w = R⁻¹·1 / (1ᵀ·R⁻¹·1)` reduces to sums of elements of the inverse matrix:

- numerator for property `i`: `s_i = Σ_j R⁻¹[i][j]` — the sum of row `i`;
- denominator: `S = Σ_i Σ_j R⁻¹[i][j]` — the sum of all elements;
- raw weight: `raw_i = s_i / S`.

For a positive definite `R`, `S > 0` always holds. Raw weights can be negative; as in the paper, negative weights are set to zero and the rest are renormalised to sum to 1.

```ts
interface BaseWeights {
  raw: Float64Array;      // s_i / S, may be negative; NaN for properties outside R
  base: Float64Array;     // clipped and renormalised, sums to 1; 0 for properties outside R
  clipped: boolean;       // true if at least one raw weight was negative
}

function campoWeights(Rinv: Float64Array, m: number, active: Int32Array, d: number): BaseWeights {
  const s = new Float64Array(m);
  let S = 0;
  for (let i = 0; i < m; i++) {
    let row = 0;
    for (let j = 0; j < m; j++) row += Rinv[i * m + j];
    s[i] = row;
    S += row;
  }

  const raw = new Float64Array(d).fill(NaN);
  const base = new Float64Array(d);
  let clipped = false;
  let pos = 0;
  for (let i = 0; i < m; i++) {
    const r = s[i] / S;
    raw[active[i]] = r;
    if (r < 0) clipped = true; else pos += r;
  }
  for (let i = 0; i < m; i++) {
    const r = raw[active[i]];
    base[active[i]] = r > 0 ? r / pos : 0;
  }
  return { raw, base, clipped };
}
```

Diagnostics for the name tooltip: the share of a property explained by the others, `R²_i = 1 − 1/R⁻¹[i][i]` ("95% of MW is explained by other properties").

### 2.8. Priorities and pins (stage B)

Order of application:

1. **Priorities.** `b_i = base_i · prio_i`, then normalisation: `base'_i = b_i / Σ b`. If `Σ b = 0`, show the form error `All weighted priorities are 0`.
2. **Pins.** Pinned weights are fixed; unpinned weights are scaled so that the total is 1:

```
P = Σ of pinned values p_i
U = Σ base'_i over unpinned properties
w_i = p_i                         for pinned
w_i = base'_i · (1 − P) / U       for unpinned
```

The priority of a pinned property has no effect on the result.

```ts
type StageBResult = { final: Float64Array; error: null } | { final: null; error: string };

/**
 * prio[i] >= 0 (default 1).
 * pinned[i] is NaN when property i is not pinned (or its input is not valid yet).
 */
function finalWeights(base: Float64Array, prio: Float64Array, pinned: Float64Array): StageBResult {
  const d = base.length;

  // 1. priorities (pinned properties are skipped: their value is fixed anyway)
  const b = new Float64Array(d);
  let B = 0;
  for (let i = 0; i < d; i++) {
    if (!Number.isNaN(pinned[i])) continue;
    b[i] = base[i] * prio[i];
    B += b[i];
  }

  // 2. pins
  let P = 0;
  for (let i = 0; i < d; i++) if (!Number.isNaN(pinned[i])) P += pinned[i];
  if (P > 1 + 1e-9)
    return { final: null, error: `Pinned weights sum to ${P.toFixed(3)} (> 1)` };

  const rest = 1 - P;
  if (rest > 1e-9 && B <= 1e-12)
    return { final: null, error: 'Nothing to scale: all unpinned weights are 0' };

  const k = B > 0 ? rest / B : 0;
  const final = new Float64Array(d);
  for (let i = 0; i < d; i++)
    final[i] = Number.isNaN(pinned[i]) ? b[i] * k : pinned[i];
  return { final, error: null };
}
```

The normalisation after priorities is merged with the pin scaling: `b_i · (1 − P) / Σ b` is the same as `base'_i · (1 − P) / U`.

**Example.** Base weights after clipping: TPSA 0.258, TPSA_S 0.038, HBA 0, HBD 0.023, MW 0, nAtoms 0.103, cLogD 0.363, cLogP 0, mbpKa 0.215; all priorities 1. Pin HBA = 0.08: `P = 0.08`, `Σ b = 1`, `k = 0.92`; TPSA = 0.237, cLogD = 0.334, and so on; the sum is 1.

**Display statuses:** `pinned` — pinned; `dropped` — final weight 0; `kept` — all others. `pinned` is checked first, so a property pinned at 0 is `pinned`, not `dropped`: its weight is the user's choice, not the computation's.

The final weights go into the profile as the properties' `weight` (§6). They are in [0, 1] and sum to 1, which the weight input's 0–1 range and `reduceMpo` both accept.

> Scientific note. Simple scaling ignores the correlation between a pinned property and the others: if HBA is pinned, TPSA (r = 0.82 with HBA) keeps its share, and polarity is counted twice. A future alternative is conditional re-optimisation with fixed `w_P`: `w_U = μ·R_UU⁻¹·1 − R_UU⁻¹·R_UP·w_P`, where `μ` is chosen so that `Σ w_U = 1 − P`, followed by clipping. It also reduces to element-wise sums. Version 1 uses scaling.

---

## 3. Layout

Reference: `campo_demo.html`; the view is `MpoProfileCreateView` (`packages/Chem/src/mpo/mpo-create-profile.ts`), the same in create and edit mode.

- **Top** (`profileViewContainer`, unchanged except one input): editable title, editable description, then the toolbar with **Method**, **Dataset**, **Aggregation**, each with its hint below (`FIELD_DESCRIPTIONS`). **Weighting** is added after **Aggregation**, on the same row, with the hint `How property weights are set`; the Aggregation hint stays as it is (`How individual property scores combine into the final MPO score`). In Data-driven mode **Weighting** is hidden together with **Aggregation** (`chem-mpo-d-none` in `onMethodChanged`). **Save** and **Reset** stay in the ribbon.
- **Content**: `ui.splitH([editorPanel, gridPanel])`, as today.
  - **left** — `MpoProfileEditor` in design mode. Its header becomes **Property / Priority / Weight / Desirability**; the Priority column is shown only in `Correlation-aware` mode. The Weight column depends on the mode: `Equal` — read-only `1.000`; `Manual` — editable input and the weight-from-column toggle, as today; `Correlation-aware` — final weight (`0.000`, read-only) and the pin icon (§5.2). In `Equal` and `Correlation-aware` the weight-from-column toggle is hidden (§5.1).
  - **right** — today a grid (`df.plot.grid()`). It becomes a detached `DG.TableView` of the dataset (`DG.TableView.create(df, false)`, laid out as Flow embeds one): its grid shows the data as today, and in `Correlation-aware` mode one panel, **Weights**, is docked `FILL` on the grid (a tab over it). It is a vertical splitter:
    1. above — the Correlation plot on the complete rows (§2.4);
    2. below — bar charts "Raw" and "Final" side by side in a horizontal splitter, each filling its half, with the legend always shown at the bottom right. Each chart has a title with a tooltip: Raw — the weights straight from the correlations, before any adjustment (a negative weight means the property is redundant); Final — the weights used in the score, after clipping, priorities and pins, summing to 1.
    The MPO score column `MPO <profile name>` is maintained by the context panel's preview (`computeMpo`), which reruns on every weight change through `MPO_SCORE_CHANGED_EVENT`.
- If **no raw weight was negative** (`clipped === false`), the Raw chart is not shown and "Final" takes the whole lower half.
- The user can close the panel with its ×. It comes back the next time the mode is switched on, or on RESET when the restored mode is `Correlation-aware`.
- When Weighting is switched to `Equal` or `Manual`, the panel is undocked and only the grid stays.
- The context panel (score histogram, best and worst compounds) is unchanged.

**Bar chart titles.** "Raw weights (may be negative)" — `R⁻¹·1` normalised by `S`, before clipping; "Final weights (sum = 1)" — after clipping, priorities and pins, i.e. the weights used in the score. "Clipped & normalized" was rejected: the final weights also include priorities and pins.

> The Correlation plot is display only. Its tab being inactive or closed does not affect the computation.


### 3.1. Weights data frame for the bar charts

```ts
const EPS = 0.004; // display-only stub so zero bars keep their row (≈1% of axis 0–0.4)

function weightsFrame(names: string[], raw: Float64Array, final: Float64Array, status: string[]): DG.DataFrame {
  const df = DG.DataFrame.fromColumns([
    DG.Column.fromStrings('property', names),
    DG.Column.fromFloat32Array('raw', Float32Array.from(raw)),
    DG.Column.fromStrings('sign', Array.from(raw, (r) => (r < 0 ? 'negative' : 'positive'))),
    DG.Column.fromFloat32Array('final', Float32Array.from(final)),                               // real values
    DG.Column.fromFloat32Array('weight', Float32Array.from(final, (w) => (w === 0 ? EPS : w))),  // chart only
    DG.Column.fromStrings('status', status),
  ]);
  df.col('sign')!.meta.colors.setCategorical({ positive: 0xFF1F77B4, negative: 0xFFFF7F0E });
  df.col('status')!.meta.colors.setCategorical({ kept: 0xFF1F77B4, pinned: 0xFF2CA02C, dropped: 0xFF8A8F99 });
  return df;
}
```

> The `weight` column is for display only. Scoring and tooltips always use `final`. Properties outside `R` (`raw = NaN`) are not shown in the Raw chart: before building it, filter out rows with `NaN` in `raw` or write `DG.FLOAT_NULL`.

### 3.2. Creating the viewers

```ts
function createCorrPlot(df: DG.DataFrame, names: string[]): DG.Viewer {
  return DG.Viewer.fromType(DG.VIEWER.CORR_PLOT, df, {
    xColumnNames: names,
    yColumnNames: names,
    correlationType: 'Pearson',
    title: 'Correlations',
  });
}

function createBarChart(df: DG.DataFrame, value: 'raw' | 'weight', stack: 'sign' | 'status',
                        title: string, valueMin: number, valueMax: number): DG.Viewer {
  return DG.Viewer.fromType(DG.VIEWER.BAR_CHART, df, {
    splitColumnName: 'property',
    valueColumnName: value,
    valueAggrType: 'sum',
    stackColumnName: stack,
    orientation: 'Horizontal',
    barSortType: 'by value',
    barSortOrder: 'desc',
    valueMin, valueMax,
    showValueSelector: false,
    showCategorySelector: false,
    showStackSelector: false,
    title, showTitle: true,
  });
}

function addEqualLine(chart: DG.Viewer, d: number): void {
  chart.meta.formulaLines.addLine({
    formula: '${weight} = ' + (1 / d).toFixed(3),
    title: 'Equal', color: '#d62728', style: 'dashed',
  });
}
```

Axis ranges: Raw — from `min(−0.1, minRaw)` to `max(0.1, maxRaw)`, rounded outward to 0.1; Final — from 0 to `max(0.1, maxFinal)`, rounded to 0.1.

```ts
/** [valueMin, valueMax] of a chart over the finite values of `a`, rounded outward to 0.1. */
function axisRange(a: Float64Array, withNegative: boolean): [number, number] {
  let lo = 0, hi = 0.1;
  for (const v of a) if (Number.isFinite(v)) { lo = Math.min(lo, v); hi = Math.max(hi, v); }
  return [withNegative ? Math.floor(Math.min(lo, -0.1) * 10) / 10 : 0, Math.ceil(hi * 10) / 10];
}
```

> Verified on a live instance: horizontal orientation, negative values, stacking by category, legend, the Equal line, the `EPS` stub for zeros, and `correlationType` on Correlation plot. Not verified: hiding the selectors.

### 3.3. Tabs and viewers

The viewers are not docked: their roots are placed in tab panes, the way `MpoContextPanel` places its histogram. A viewer that was never docked throws on `close()` (§10, question 4), so `CampoViewers` never calls it; it drops the references and lets the panes go.

```ts
// packages/Chem/src/mpo/campo-viewers.ts
class CampoViewers {
  readonly tabs = ui.tabControl();
  private gridPane: DG.TabPane;
  private corr: DG.Viewer | null = null;
  private rawChart: DG.Viewer | null = null;
  private finalChart: DG.Viewer | null = null;
  private weightsHost = ui.div([], 'chem-mpo-campo-weights');
  private closedByUser = {corr: false, weights: false};
  private clearing = false;

  constructor(private gridRoot: HTMLElement, private gridTitle: string) {
    this.gridPane = this.tabs.addPane(gridTitle, () => gridRoot);
    this.tabs.onTabRemoved.subscribe((pane: DG.TabPane) => {   // the user pressed ×
      if (this.clearing) return;
      if (pane.name === 'Correlations') { this.closedByUser.corr = true; this.corr = null; }
      if (pane.name === 'Weights') { this.closedByUser.weights = true; this.rawChart = this.finalChart = null; }
    });
  }

  /** Mode switched on, or RESET restored it: both tabs come back. */
  open(): void {
    this.closedByUser = {corr: false, weights: false};
  }

  showCorr(frame: DG.DataFrame, names: string[]): void {
    if (this.closedByUser.corr) return;
    if (!this.corr) {
      this.corr = createCorrPlot(frame, names);
      this.tabs.addPane('Correlations', () => this.corr!.root, null, {allowClose: true});
    } else {
      this.corr.dataFrame = frame;
      this.corr.setOptions({xColumnNames: names, yColumnNames: names});
    }
  }

  showWeights(weightsDf: DG.DataFrame, raw: Float64Array, final: Float64Array, clipped: boolean,
              activate: boolean): void {
    if (this.closedByUser.weights) return;
    const first = !this.finalChart;
    const [fMin, fMax] = axisRange(final, false);
    if (!this.finalChart) {
      this.finalChart = createBarChart(weightsDf, 'weight', 'status', 'Final weights (sum = 1)', fMin, fMax);
      addEqualLine(this.finalChart, final.length);
    } else {
      this.finalChart.dataFrame = weightsDf;
      this.finalChart.setOptions({valueMin: fMin, valueMax: fMax});
    }

    const [rMin, rMax] = axisRange(raw, true);
    if (clipped && !this.rawChart)
      this.rawChart = createBarChart(weightsDf, 'raw', 'sign', 'Raw weights (may be negative)', rMin, rMax);
    else if (!clipped)
      this.rawChart = null;
    else {
      this.rawChart!.dataFrame = weightsDf;
      this.rawChart!.setOptions({valueMin: rMin, valueMax: rMax});
    }

    ui.empty(this.weightsHost);
    this.weightsHost.append(this.rawChart ?
      ui.splitH([this.rawChart.root, this.finalChart.root], {}, true) : this.finalChart.root);
    if (first) {
      const pane = this.tabs.addPane('Weights', () => this.weightsHost, null, {allowClose: true});
      if (activate) this.tabs.currentPane = pane;
    }
  }

  /** Mode switched to Equal or Manual, or stage A failed: only the grid stays. */
  hide(): void {
    this.clearing = true;
    try { this.tabs.clear(); }
    finally { this.clearing = false; }
    this.gridPane = this.tabs.addPane(this.gridTitle, () => this.gridRoot);
    this.corr = this.rawChart = this.finalChart = null;
  }
}
```

Axis ranges of the charts follow §3.2 (`axisRange`, applied with `setOptions({valueMin, valueMax})` on each update). `hide()` is also called when stage A fails (§4); the next successful stage A brings the tabs back, unless the user closed them. `MpoProfileCreateView.setupGridAndContextPanel` puts `campoViewers.tabs.root` where it puts `gridContainer` today; on dataset change the view is rebuilt as it is now (`attachLayout`), with a new `CampoViewers`.

> `hide()` runs `clear()` under the `clearing` flag, so whatever `onTabRemoved` reports for `clear()`, its removals are not taken for the user's. This matters because `hide()` also runs on every stage A error: without the flag, one error could mark both tabs as closed by the user, and they would not come back until the mode is switched on again. Whether `onTabRemoved` fires for `clear()` at all no longer needs checking.

---

## 4. Reactivity

All computation runs automatically; there is no Recompute button.

| Event | Where it comes from today | What is recomputed |
|---|---|---|
| Dataset change, view opened | `onDatasetChanged` → `attachLayout` → `setupGridAndContextPanel` | the view calls `campo.start(false)`: new controller, stage A + B, viewers rebuilt |
| Property added / removed / renamed, column remapped | `MpoProfileEditor.onChanged` | stage A + B |
| Desirability curve or missing-value handling edited | `MpoProfileEditor.onChanged` (debounced) | stage A + B |
| Column of a property created by its compute function | `df.onColumnsAdded` | stage A + B |
| `Correlation-aware` mode switched on | Weighting input | tabs opened, stage A + B |
| Priority changed (valid value) | Priority input | stage B only |
| Pin set / removed (valid value) | pin icon, pin input | stage B only |
| Switched to `Equal` / `Manual` | Weighting input | tabs removed, no computation (§5.1) |
| RESET | ribbon **Reset** → `resetProfile` | the saved profile is restored with its weighting, priorities and pins; if its mode is `Correlation-aware`, `resetProfile` reopens the tabs and calls `runStageA(true)` |

`editor.setProfile` (used by `attachLayout` and `resetProfile`) rebuilds the rows and fires **no** change event, and `resetProfile` does not go through `attachLayout`. So the two rows that start from it are wired in the view, not through `onChanged`. `runAllComputeFunctions` (called from `setProfile`) adds columns asynchronously and also fires nothing; hence the `onColumnsAdded` row.

`MpoProfileEditor.onChanged` fires for every edit, including the weights the controller itself writes. The controller tells the cases apart with a fingerprint of what stage A depends on, and runs stage A only when it changed. Whether the mapped column exists is part of the key, so a column that appears later is picked up:

```ts
/** Everything stage A depends on; weights, priorities and pins are not part of it. */
function stageAKey(df: DG.DataFrame, profile: DesirabilityProfile, mapping: Record<string, string | null>): string {
  return df.id + '|' + Object.entries(profile.properties).map(([name, p]) => {
    const col = mapping[name] ?? name;
    return `${name}>${col}>${df.col(col) ? 1 : 0}>${desirabilityKey(p)}`;
  }).join('|');
}
```

**Modified state.** `listenForProfileChanges` marks the profile modified on every `MPO_SCORE_CHANGED_EVENT` outside `updatingLayout`. A recomputation that gives the same weights must not do that, or a saved `Correlation-aware` profile would open as modified. So `editor.setWeights` compares the new weights with the stored ones at the display precision (`toFixed(3)`) and, if none changed, writes nothing and fires nothing. If they did change (the profile is opened on another dataset), it fires as usual and **Save** becomes active: the stored weights are no longer the ones the score would use.

**Edits to the data in place.** Changing values in the dataset (the same `df.id`) does not rerun stage A: nothing in §4 tracks cell edits, and `MpoCalculator` validates its cache by curve, not by data, so it would return the old desirabilities anyway. The MPO score preview has the same limitation today. Weights follow such edits when stage A next runs for another reason with a cold cache, or when the dataset is selected again. Tracking `df.onDataChanged` would need the cache invalidated per column and a guard against the score and desirability columns the preview itself writes; not in version 1.

```ts
// packages/Chem/src/mpo/campo-controller.ts
import {Subscription} from 'rxjs';
import {debounceTime} from 'rxjs/operators';

class CampoController {
  private baseW: BaseWeights | null = null;
  private keyA = '';
  private warningsA: string[] = [];   // stage A warnings; stage B shows them whenever it has no error
  private subs: Subscription[] = [];

  /** `activate`: make the Weights tab current when it is first shown (mode switched on, RESET into the mode). */
  /** `viewers` is null when there is no dataset (the view has no right half then). */
  constructor(private view: MpoProfileCreateView, private viewers: CampoViewers | null, private activate: boolean) {
    // curve dragging fires many events: debounce stage A
    this.subs.push(view.editor.onChanged.pipe(debounceTime(150)).subscribe(() => this.runStageA()));
    if (view.df)
      this.subs.push(view.df.onColumnsAdded.pipe(debounceTime(150)).subscribe(() => this.runStageA()));
  }

  private get names(): string[] { return Object.keys(this.view.profile.properties); }

  private fail(message: string): void {
    this.baseW = null;
    this.warningsA = [];
    this.viewers?.hide();
    this.view.showCampoMessage(message, 'error');
  }

  runStageA(force = false): void {
    const {df, profile, editor} = this.view;
    const names = this.names;
    if (names.length < 2) {
      this.keyA = '';
      editor.setWeights(Object.fromEntries(names.map((nm) => [nm, 1])), names.map(() => 'kept'), null);
      return this.fail('Needs at least two properties');
    }
    if (!df) {
      this.keyA = '';
      return this.fail('Select a dataset to compute weights');
    }
    const key = stageAKey(df, profile, editor.columnMapping);
    if (!force && key === this.keyA) return;
    this.keyA = key;

    // a column mapped to by an earlier property is not used again: one D_k per column
    const seen = new Set<DG.Column>();
    const columns = names.map((nm) => {
      const c = df.col(editor.columnMapping[nm] ?? nm);
      if (!c || seen.has(c)) return null;
      seen.add(c);
      return c;
    });
    const tagged = applyDesirabilityTags(df, Object.fromEntries(names
      .map((nm, k) => [columns[k]?.name, profile.properties[nm]]).filter(([c]) => c)), true);
    const calculator = this.view.mpoContextPanel?.calculator ?? new MpoCalculator();
    const maps = calculator.desirabilities(df, tagged);
    const D = columns.map((c) => (c ? maps[tagged.indexOf(c)].D : null));

    const prep = prepareDesirabilityFrame(names, D, df.rowCount);
    const m = prep.active.length;
    if (prep.rows === 0)
      return this.fail('No rows with all properties present');
    if (m < 2 || prep.rows <= m)
      return this.fail('Not enough data to compute correlation-aware weights');

    this.viewers!.showCorr(prep.frame, Array.from(prep.active, (k) => names[k]));
    const R = pearsonMatrix(prep.cols, prep.rows);
    const inv = robustInverse(R, m);
    if (!inv)
      return this.fail('Correlation matrix is singular');

    this.warningsA = [];
    if (prep.rows < 10 * m)
      this.warningsA.push(`Only ${prep.rows} complete rows for ${m} properties: weights may be unreliable.`);
    if (prep.rows < df.rowCount / 2)
      this.warningsA.push(`Weights are derived from ${prep.rows} of ${df.rowCount} rows that have every property.`);
    if (inv.lambda > 0)
      this.warningsA.push(`Correlation matrix is nearly singular: regularised with λ = ${inv.lambda}.`);
    this.baseW = campoWeights(inv.Rinv, m, prep.active, names.length);
    this.runStageB();
  }

  /** Cheap: priorities and pins only. */
  runStageB(): void {
    if (!this.baseW) return;
    const names = this.names;
    const props = names.map((nm) => this.view.profile.properties[nm]);
    const prio = Float64Array.from(props, (p) => p.priority ?? 1);
    const pinned = Float64Array.from(props, (p) => p.pinnedWeight ?? NaN);
    const {final, error} = finalWeights(this.baseW.base, prio, pinned);
    if (!final)
      return this.view.showCampoMessage(error, 'error');
    this.view.showCampoMessage(this.warningsA.length ? this.warningsA.join(' ') : null, 'warning');
    const status = names.map((_, i) =>
      !Number.isNaN(pinned[i]) ? 'pinned' : final[i] === 0 ? 'dropped' : 'kept');
    this.viewers!.showWeights(weightsFrame(names, this.baseW.raw, final, status), this.baseW.raw, final,
      this.baseW.clipped, this.activate);
    this.activate = false;
    this.view.editor.setWeights(Object.fromEntries(names.map((nm, i) => [nm, final[i]])), status,
      this.baseW.base);
  }

  dispose(): void { this.subs.forEach((s) => s.unsubscribe()); this.baseW = null; this.keyA = ''; }
}
```

- `editor.setWeights(weights, status, base)` writes `weight` on every property and refreshes the read-only cells without firing a change per property; `base` (the base weights before priorities, `null` when there are none) feeds the Priority hint (§5.2). If any weight changed (see **Modified state** above), it then fires `onChanged` and `MPO_SCORE_CHANGED_EVENT` once, so the score, the context panel and the view's "modified" state update as for any edit. The `onChanged` it fires reaches `runStageA`, which returns at once because the key is unchanged.
- Priority and pin inputs write `priority` / `pinnedWeight` on the property and call `runStageB()`.
- `MpoContextPanel.calculator` becomes a public read-only field; the controller reads it on every run (§2.3).
- **Messages.** One message line (§5.4). An error (stage A or stage B) replaces whatever is shown. Stage A warnings are kept in `warningsA` and shown by every successful stage B, so editing a priority or a pin does not wipe them; a stage B error hides them until it is fixed.
- The view owns the controller: `campo.start(activate)` disposes the old one, creates a new one and calls `runStageA(true)`. It is called with `activate = true` when the mode is switched on and in `resetProfile` when the restored mode is `Correlation-aware`, and with `false` at the end of `setupGridAndContextPanel` (opening the view, dataset change), so that opening a saved profile keeps the dataset tab current. `attachLayout` returns early when there is no dataset; in `Correlation-aware` mode that branch also calls `campo.start(false)`, with no `CampoViewers`, so the form says why there are no weights. `activate` is used once, by the first stage B that shows the Weights tab. The controller is disposed when the mode is switched off, before the dataset changes and in the view's `detach()`.

---

## 5. Form columns

All of this is in `MpoProfileEditor` and `createMpoRow`, behind a mode the host sets (`editor.setWeighting(mode)`); without it the editor behaves as today (EDA, PowerGrid, the MPO Score dialog).

### 5.1. Mode switching and weights

- **`Manual`** — today's editor: the weight input (0–1) and the weight-from-column toggle.
- **`Equal`** — every `weight` is set to 1; the cells show a read-only `1.000`.
- **`Correlation-aware`** — the cells show the final weights (§5.2).
- In `Equal` and `Correlation-aware` the weight-from-column toggle is hidden and `weightColumn` is removed from the properties: per-row weights would override the computed ones in `reduceMpo`.
- **Leaving `Manual`**: the manual weights and `weightColumn` values are stashed in the view for the session. **Back to `Manual`**: they are restored from the stash; without a stash (the profile was opened in another mode) Manual starts from the current effective weights.
- **`Manual → Correlation-aware`**: current manual weights are copied into the priorities (overwriting stored priorities). For a property with `weightColumn` the static `weight` is copied (the per-row column cannot be a priority); `weightColumn` goes to the stash with the manual weights.
- **`Equal → Correlation-aware`**: priorities are left as stored (absent = 1); in `Equal` every weight is 1, so copying them would only reset priorities to 1.
- **`Correlation-aware → Equal`** and **back**: `priority` and `pinnedWeight` stay on the properties untouched (§6), so the round trip restores them.
- `d = 1`: see §7. Aggregations other than `Average`: see §8.

### 5.2. Priority

- Shown only in `Correlation-aware` mode, as a column between Property and Weight (`createMpoRow` gets a `priorityCell`; the editor's header gets `Priority`).
- `ui.input.float`, min 0, no max, step 0.1, default 1; stored as the property's `priority` (§6).
- Validator: not empty → `Enter a priority`; ≥ 0 → `Must be 0 or greater`.
- For a pinned property the field stays editable: the priority is kept and applies once the weight is unpinned; while it is pinned and the priority ≠ 1 the field shows `Applies when the weight is unpinned.`
- If the base weight of a property is 0 and its priority ≠ 1, the field shows a hint: `No effect: computed weight is 0. Pin the weight instead.` — the priority multiplies the weight, and zero stays zero. The base weight comes from the third argument of `editor.setWeights` (§4); the final weight cannot tell this case apart from a priority of 0.

### 5.3. Weight and pinning

**Behaviour:**

1. Each row of the Weight column shows the final weight (format `0.000`, read-only) and a `ui.iconFA('thumbtack')` icon.
2. **Icon tooltip:** `Pin weight` (not pinned) / `Unpin weight` (pinned). A Datagrok tooltip (`ui.tooltip.bind`), not the browser `title`: shown on hover and on keyboard focus, hidden on click; after the click the next hover shows the new state.
3. **Clicking the icon on an unpinned property:** the number is replaced by a `ui.input.float` field holding the current final weight (rounded to 0.001), selected for editing. The property is pinned at once at that value, so the weights do not change until the user edits it.
4. **Entering a valid value:** `pinnedWeight` is set → stage B: unpinned weights are scaled (§2.8), the Weight column and the Final bar chart are updated.
5. **Clicking the icon on a pinned property:** `pinnedWeight` is removed, the field is replaced by the number → stage B.

**Input field:**

- min 0, max 1, step 0.01, format `0.000` (as the existing weight input: `#0.000`).
- Validator (the message is the first failed condition):
  1. the value is not empty → `Enter a weight`;
  2. the value is in [0, 1] → `Must be between 0 and 1`;
  3. the sum of all pinned values ≤ 1 → `Pinned weights sum to 1.07 (> 1)`.
- Condition 3 spans rows: when any pinned field changes, **all** pinned fields are revalidated.

```ts
interface PinState {
  inputs: Map<string, DG.InputBase<number | null>>;   // by property name
  final: (name: string) => number;
}

function pinnedSum(state: PinState): number {
  let s = 0;
  for (const inp of state.inputs.values()) {
    const v = inp.value;
    if (v != null && Number.isFinite(v)) s += v;
  }
  return s;
}

function revalidateAll(state: PinState): void {
  for (const inp of state.inputs.values()) inp.validate();
}

/** setPin(name, value | null) writes pinnedWeight and runs stage B. */
function weightCell(name: () => string, state: PinState, setPin: (name: string, v: number | null) => void): HTMLElement {
  const host = ui.divH([], 'statistics-mpo-campo-weight');

  const render = () => {
    ui.empty(host);
    const input = state.inputs.get(name());
    const valueEl = input ? input.root : ui.span([state.final(name()).toFixed(3)]);
    const icon = ui.iconFA('thumbtack', () => toggle());
    if (input) icon.classList.add('fas');          // filled when pinned, outline otherwise
    ui.tooltip.bind(icon, () => (state.inputs.has(name()) ? 'Unpin weight' : 'Pin weight'));
    host.append(valueEl, icon);
  };

  const toggle = () => {
    if (state.inputs.has(name())) {
      state.inputs.delete(name());
      revalidateAll(state);
      setPin(name(), null);
    } else {
      const input = ui.input.float('', {
        value: null, min: 0, max: 1, step: 0.01, format: '#0.000',
        onValueChanged: () => {
          revalidateAll(state);
          setPin(name(), input.validate() ? input.value : null);
        },
      });
      input.addValidator((s: string) => {
        if (s == null || s.trim() === '') return 'Enter a weight';
        const v = Number(s);
        if (!Number.isFinite(v) || v < 0 || v > 1) return 'Must be between 0 and 1';
        const sum = pinnedSum(state);
        if (sum > 1 + 1e-9) return `Pinned weights sum to ${sum.toFixed(2)} (> 1)`;
        return null;
      });
      state.inputs.set(name(), input);
      setTimeout(() => input.input.focus());
    }
    render();
  };

  render();
  return host;
}
```

`name` is a getter because the editor renames properties in place (`updatePropertyName`); the pin state follows the rename.

**Pins of a loaded profile.** `editor.setProfile` (opening the view, RESET, dataset change) creates a new `PinState`. For every property with `pinnedWeight`, it creates the pin input with that value before the rows are rendered, the same input `toggle()` creates, so the row opens pinned and the validator sees the sum of all pins. Nothing of the previous profile's pins is kept.

### 5.4. Form messages

The mode has one message at a time, not a line under the property table: a long table hid it below the fold. It is shown as an `exclamation-circle` icon (regular, not solid) right of the **Weighting** input (beside it, not in the input's options): red (`--red-3`) for an error, grey for a warning, with the message in its tooltip; the icon disappears when there is no message. When a message appears (or its kind changes), a warning balloon (`grok.shell.warning`) shows it too; a recomputation that keeps a message does not repeat the balloon, even if the text changes (e.g. the list of constant properties shrinks while the user shapes their curves) (v8):

- errors: `Desirability does not vary for …`, `No rows with all properties present`, `Not enough data…`, `Correlation matrix is singular`, the errors of `finalWeights` — the weights are not applied;
- warnings: regularisation used (`λ = 0.01`), few complete rows (§7), complete rows less than half the dataset (§2.4).

### 5.5. "Why this weight" tooltip

Hovering (or focusing) a property name in `Correlation-aware` mode shows a tooltip that explains the weight. In the profile view the editor is in design mode, so the name is an input (`buildPropertyNameInput`); the tooltip is bound to the property cell (`ctx.propertyCell`), which holds it.

- the property name and its final weight;
- the share of its desirability explained by the other properties, `R²_i = 1 − 1/R⁻¹[i][i]`;
- the reason, if any: `Pinned at 0.080.` or `Raw weight −0.08 < 0, so it was set to 0.`; for a constant column: `Constant in this dataset — excluded.`; for an unmapped property: `No column in this dataset — excluded.`; for a second property on the same column: `Column <name> is already used by <property> — excluded.` (§7; `WhyInput` gets a `sharedWith: Map<number, number>` next to `unmapped`);
- the three most correlated properties (by |r|) with their r values; the demo shows each with a small red (positive) or blue (negative) bar.

```ts
interface WhyInput {
  names: string[];
  R: Float64Array; m: number; active: Int32Array;   // R over active properties
  Rinv: Float64Array;
  raw: Float64Array; final: Float64Array; pinned: Float64Array;
  unmapped: Set<number>;
}

function whyTooltip(i: number, x: WhyInput): HTMLElement {
  const pos = Array.from(x.active).indexOf(i);
  const lines: HTMLElement[] = [ui.divText(`${x.names[i]} · final weight ${x.final[i].toFixed(3)}`, 'chem-mpo-why-title')];
  if (pos < 0) {
    lines.push(ui.divText(x.unmapped.has(i) ? 'No column in this dataset — excluded.' : 'Constant in this dataset — excluded.'));
    return ui.divV(lines);
  }
  const r2 = 1 - 1 / x.Rinv[pos * x.m + pos];
  lines.push(ui.divText(`${Math.round(r2 * 100)}% of its desirability is explained by the other properties.`));
  if (!Number.isNaN(x.pinned[i])) lines.push(ui.divText(`Pinned at ${x.pinned[i].toFixed(3)}.`));
  else if (x.raw[i] < 0) lines.push(ui.divText(`Raw weight ${x.raw[i].toFixed(2)} < 0, so it was set to 0.`));

  const partners = Array.from(x.active)
    .map((k, q) => ({ k, r: x.R[pos * x.m + q] }))
    .filter((o) => o.k !== i)
    .sort((a, b) => Math.abs(b.r) - Math.abs(a.r))
    .slice(0, 3);
  lines.push(ui.divText('Most correlated', 'chem-mpo-why-caption'));
  for (const o of partners) lines.push(ui.divText(`${x.names[o.k]}   ${o.r.toFixed(2)}`));
  return ui.divV(lines);
}
```

---

## 6. Persistence in the profile

The profile format (`mpo-types.ts`) gets three optional fields:

```ts
export type MpoWeighting = 'equal' | 'manual' | 'correlation';

type BasePropertyDesirability = {
  weight: number;          // the effective weight, in every mode
  weightColumn?: string;   // Manual mode only (§5.1)
  priority?: number;       // Correlation-aware: multiplier, absent = 1
  pinnedWeight?: number;   // Correlation-aware: pinned weight, absent = not pinned
  missingValues?: MissingValueConfig;
  computeFunction?: TemplateFunction;
};

export type DesirabilityProfile = {
  // ...existing fields
  weighting?: MpoWeighting;   // absent = 'manual'
};
```

- `weight` always holds the weight the score uses: 1 in `Equal`, the user's value in `Manual`, the final caMPO weight in `Correlation-aware`. That is the snapshot of §10 question 5: the MPO Score dialog, the profile preview, `mpoTransformFunction` / `mpoCalculate`, other packages and older clients compute the same score without knowing about caMPO.
- `priority` and `pinnedWeight` live on the property, not in name-keyed maps, so renaming a property in the editor carries them along. They are written only when they differ from the default and are kept, not applied, in other modes.
- `weighting` is written only when it is not `manual`.
- `CURRENT_MPO_VERSION` stays 2: every new field is optional, and a client that does not know them still reads a valid profile with the right weights.

**Storage.** `priority` and `pinnedWeight` are inside `properties` (jsonb) and need no schema change. `weighting` is profile-level and gets a column:

- `databases/mpo/schema.json`: `"weighting": {"type": "string", "choices": ["equal", "manual", "correlation"], "description": "How property weights are set; empty means manual."}`, schema version `1.2.0`;
- `databases/mpo/0002_weighting.sql` and `databases/mpo/down/0002_weighting.sql`, generated by the schema snapshot tool like `0001`;
- `src/generated/db.ts` regenerated;
- `mpo-profile-store.ts`: `toRow` writes `weighting`, `fromRow` reads it.

Export and import go through `toPlainProfile` / `parseMpoProfile` and carry the new fields with no change.

---

## 7. Edge cases

| Case | Behaviour |
|---|---|
| `d = 1` | `Correlation-aware` has nothing to compute: `runStageA` checks it first, shows `Needs at least two properties`, sets the weight to 1 and shows no tabs |
| Property with no column in the dataset | excluded from `R`, weight 0, status `dropped`, reason in the name tooltip |
| Two properties mapped to the same column | the first one in profile order uses the column; the others are excluded from `R`, weight 0, status `dropped`, tooltip `Column <name> is already used by <property> — excluded.` (otherwise `R` gets two identical columns, r = 1) |
| Column created by a compute function | until the column exists the property is treated as unmapped; stage A reruns when it is added (§4) |
| Constant desirability column | excluded from `R` and from the Correlation plot, weight 0, status `dropped`, reason in the name tooltip |
| Fewer than 2 properties in `R` | error `Not enough data…`, tabs are removed (`viewers.hide()`), weights keep their last values |
| No complete rows | error `No rows with all properties present`, tabs are removed, weights keep their last values |
| Few complete rows (`rows < 10·m`) | computation runs; the form message warns that weights may be unreliable |
| Complete rows less than half the dataset | computation runs; the form message names how many rows the weights come from (§2.4) |
| `rows ≤ d` | `prepareDesirabilityFrame` returns no active columns → error `Not enough data…`, weights keep their last values |
| Cholesky failed / `κ > 10⁶` | regularisation `R_λ` (§2.6); the fact and `λ` are shown in the form message |
| No negative raw weight | Raw bar chart is not shown |
| Sum of pinned = 1 | all unpinned weights are 0 (allowed) |
| Sum of pinned < 1, all unpinned weights × priorities = 0 | form error, result not applied |
| All priorities of unpinned properties = 0 and no pins | form error |
| No dataset selected | the profile view shows only the editor (as today); the view still starts the controller in `Correlation-aware` mode, and `runStageA` shows `Select a dataset to compute weights` and keeps the stored weights |
| Values edited in place | weights are not recomputed (§4, **Edits to the data in place**) |
| More than `MAX_MPO_PROPERTIES` (10) properties | no limit for caMPO; 10 is only the limit of auto-population |

---

## 8. Aggregations

caMPO weights are derived for a weighted average. The aggregations are `Average` (default), `Sum`, `Product`, `Geomean`, `Min`, `Max` (`WEIGHTED_AGGREGATIONS`); with any of them other than `Average` the weights mean less than caMPO computed (`Product`, `Min` and `Max` use them as exponents; `Geomean` averages `log dᵢ`, while the weights come from the correlations of `dᵢ`).

**Decision (v8):** `Correlation-aware` is allowed with every aggregation. The weights are computed the same way whatever the aggregation, and the score uses them with the selected one. The banner with **Switch to Average** / **Use Manual** proposed in v7 is dropped; instead, when the user (or the AI tools) switches Weighting to `Correlation-aware` while Aggregation is not `Average`, or picks another aggregation in that mode, a warning balloon is shown (`grok.shell.warning`): `Correlation-aware weights are derived for Average. With <aggregation>, they reduce redundancy only approximately.` Opening a saved profile shows none.

---

## 9. Tests

The math of §2 is in `libraries/statistics/src/mpo/campo.ts`; its tests and the integration tests go to the Chem package, next to the existing `src/tests/mpo-tests.ts`.

**Synthetic cases from the paper** (the matrix is passed directly to `robustInverse` + `campoWeights`):

| Matrix | Expected weights |
|---|---|
| 4 properties: A independent, B–D pairwise r = 0.98 | A ≈ 0.497, B = C = D ≈ 0.168 |
| 3 properties: A independent, r(B, C) = −0.35 | A ≈ 0.245, B = C ≈ 0.377 |
| `d×d` identity matrix | all weights = 1/d |

**Other:**

- `pearsonMatrix`: on complete rows, values match `DG.CorrelationPlot.getCorrelation` within `1e-6` (this also guards that **Correlations** shows the matrix the weights come from); the matrix is symmetric with 1 on the diagonal.
- `MpoCalculator.desirabilities`: returns the arrays `compute()` caches (same object on a second call with an unchanged curve and no `release()` in between); missing values of a `double` column come out as `NaN` under `exclude` and `skip`, as the default score under `default`.
- `prepareDesirabilityFrame` with no complete rows returns no active columns and no `NaN`.
- Two properties mapped to one column: the second is `dropped`, `R` has the column once.
- A compute-function property whose column is added after stage A gets a weight once the column appears.
- Opening a saved `Correlation-aware` profile on its dataset leaves it unmodified (Save inactive); on another dataset with different weights it is modified.
- RESET in `Correlation-aware` mode restores the saved pins (inputs filled) and recomputes; a stage A error removes the Correlations and Weights tabs, and the next successful stage A brings them back (a removal by `hide()` is not taken for the user's ×).
- Switching the mode on makes Weights current; opening a saved `Correlation-aware` profile leaves the dataset tab current.
- A regularisation warning survives a priority or pin edit; a stage B error replaces it, and fixing the error brings it back.
- `d = 1` → `Needs at least two properties`, weight 1; no dataset → `Select a dataset to compute weights`.
- The Priority hint appears for a property with base weight 0 and priority ≠ 1, and not for one with priority 0.
- `cholesky`: returns `null` on a non-positive-definite matrix; `robustInverse` finds a working `λ`.
- `finalWeights`: the example from §2.8; priority 2 on one property; sum of pinned > 1 → error; = 1 → unpinned are 0.
- A constant column and an unmapped property are excluded without exceptions and without `NaN` in final weights.
- The score with caMPO weights equals `Σ wᵢ·dᵢ` on complete rows (weights sum to 1) and matches `mpo()` with the same profile.
- Writing weights does not rerun stage A (`stageAKey` unchanged), and editing a curve does.
- Profile round trip through `mpoProfileStore.save` / `load`: `weighting`, `priority`, `pinnedWeight` and the effective `weight` survive; a profile without them loads as `manual`.
- Renaming a property keeps its priority and pin.
- Switching modes adds and removes tabs without leaking subscriptions; a tab closed with × stays closed until the mode is switched on again.
- The Raw chart appears / disappears as data changes.
- MPO Score dialog with a `correlation` profile: weights read-only, score equals the one in the profile view.
- Performance: `n = 10⁵`, `d = 20` — stage A (including `pearsonMatrix`) < 300 ms, excluding the viewers' rendering.

Changelogs: `## v.next` of `packages/Chem/CHANGELOG.md`; a version and changelog entry of `libraries/statistics` (now 1.12.16).

---

## 10. Open questions

Questions 1, 2 and 4 were checked with `campo_corrplot_check.js` on 2026-10-26.

1. ~~**Updating the `table` field of Correlation plot.**~~ Resolved. `DG.CorrelationPlot.getCorrelation(c1, c2)` computes the value synchronously on each call: correct in an inactive tab, right after an in-place edit and right after `dataFrame` is replaced, so there is nothing to wait for. It has no cache (1.5 s for 20×20 on 10⁵ rows), so `R` is computed by our own code (§2.5). Also found: the viewer uses pairwise deletion of missing values, returns `NaN` for a constant column, honours `correlationType`, and does not appear to apply `df.filter`. That last point is not confirmed: the script did not compare against the unfiltered data, and it does not matter here because the desirability frame has no filter.
2. ~~**Structure of `table`.**~~ Resolved: there is no matrix to read. `table` is the name of the bound table (`null` for a viewer built in code).
3. ~~**Docking `RIGHT` inside a `FILL` tab.**~~ No longer applies: the profile view has no dock manager in Manual mode, and the viewers live in a tab control (§3.3).
4. ~~**Closing tabs.**~~ Resolved for docked viewers: `grok.events.onViewerClosed` gives `args.args.viewer` (match by `.dart`), `dockManager.onClosed` gives `viewer.root`; `close()` on a viewer that was never docked throws. With the tab control (§3.3) the signal is `TabControl.onTabRemoved`; `hide()` guards `clear()` with a flag, so what it reports for `clear()` does not matter.
5. **Applying the profile to other tables.** Proposed (§1.2, §6): the stored effective weights are used as they are; weights are recomputed only in the profile view, on its dataset. To confirm.
6. **Pin semantics.** Version 1 uses simple scaling (§2.8). Consider conditional re-optimisation in version 2.
7. **Manual weights after a save in another mode.** With `weight` holding the effective weight (§6), the manual weights survive a mode switch only within the session (§5.1). If they must survive a save in `Equal` / `Correlation-aware` mode, add `manualWeight` to the property.
8. **AI function on a pinned property.** `setMpoPropertyWeight` refuses in `Correlation-aware` mode (§1.2). For that mode it could set `pinnedWeight` instead, which is what a user does by hand; or `pinMpoWeight` stays the only way. To decide.
9. **Rows of the weights vs rows of the score.** Weights come from complete rows, the score under `skip` uses every row (§2.4). Version 1 accepts this and warns below half the dataset; whether the threshold is right, and whether `skip` properties should be correlated on pairwise rows with a nearest positive-definite correction instead, is open.

---

## 11. Implementation notes

Implemented on 2026-10-27 (branch `vmakarichev/chem-ca-mpo`). The math is in `packages/Chem/src/mpo/campo-math.ts` for now, not in `libraries/statistics` as §1.1 says; it has no platform imports and moves to a library as is. Where the code differs from the text above:

- **Migration.** No `0002_weighting.sql`: a new optional column is an additive change that the domain-schema deploy applies itself. `schema.json` is at 1.2.0 and `src/generated/db.ts` is regenerated.
- **Flat curves.** When fewer than two desirabilities vary, the error names the properties with flat curves (`Desirability does not vary for …: shape their curves to compute weights`) instead of `Not enough data…`: a profile auto-populated from a dataset starts with flat curves (§7).
- **Computed weights in other hosts** are disabled, not read-only: a float input with a range has a slider that `readOnly` leaves draggable (§1.2).
- **"Why this weight" tooltip** is on the weight value and on the property cell around the name input; the name input keeps the hover for its own validation tooltip (§5.5).
- **`editor.setWeights(weights, base)`** takes no statuses; the editor reads pins from the properties (§4).
- **Stage B right after an edit.** A priority or pin edit, and every AI function, first reruns stage A when its key is stale (`CampoController.update()`), so a rename or deletion inside the stage A debounce cannot reach stage B with old names (§4).
- **Splitters.** The panel's splitters are flex-share splitters in `campo-viewers.ts` (`flexSplit`), not `ui.splitV` / `ui.splitH`: those fix pixel sizes when first inserted into the document, which happens before the dock has laid the panel out, so the panes stayed small.
- **Docked panel.** One panel (Correlation plot over both charts, in splitters) is docked `FILL` on the grid of a detached table view in the right half (§3), not shown in a tab control. `CampoViewers` owns the view; a detached view has no dock manager until it is added, so the view is added when its root enters the document, and shows requested before that are replayed then. Data updates change the viewers in place; Raw appearing or going only rebuilds the charts' half of the panel. The user's × is tracked through `dockManager.onClosed`. §3.3's code shows the earlier tab control.
- **Charts** are titled `Raw` and `Final`. The titles are drawn by `CampoViewers` above each chart, not by the viewers, so that they carry tooltips with the explanation (§3). When the chart is shown alone (no negative raw weights), it is titled `Weights`, with a tooltip that leaves out the clipping of negative weights.
- **Detached view layout.** A viewer root is `inline-flex`; in a block slot it sits on the text baseline and its bottom (the value axis) is clipped, so the slots are `display: flex`. The `div.ui-box > :not(...)` rule gives the view `overflow: auto !important`, and a pixel of rounding then added scrollbars under the tabs, so `overflow: hidden` is set inline with `important`. The JS `TableView` has no `_handleResize`; the dock follows the view's size through `DG.DockView.prototype._handleResize`.
- **Weighting input** is `MpoProfileEditor.weightingInput`, next to `aggregationInput`, and is set by `setProfile` / `setWeighting`; the view only places it and handles `onWeightingSelected` (moving the weights between modes).
- **Messages** of the mode are the icon next to Weighting and a balloon (§5.4, v8), and the aggregation warning is a balloon (§8, v8); there is no message line under the property table and no banner.
- **Weight column** is 96 px in Correlation-aware mode: the value sits in a box the size of the pin input, so it lines up with the row's inputs and the pin icon does not move on pinning. The header rules use `> div.` selectors because other packages bundle older copies of the editor's CSS whose `:nth-child` header rules would win otherwise.
- **Pins over 1** stay on the properties and stage B reports the error; only a pin whose own value is not a weight is treated as unpinned (§5.3). `setMpoPriority` and `pinMpoWeight` return `success: false` with the stage B error.
- **Statistics library** 1.12.17: `MpoWeighting`, `priority`, `pinnedWeight`, `weighting`, `MpoCalculator.desirabilities()`, the editor's weighting UI.
- **Tests** (`src/tests/campo-tests.ts`) cover §9 except the view-level cases (opening, RESET, tabs, messages), which were checked by hand on a local stand.
