# Specification: Correlation-aware (caMPO) weighting mode in MPO Profiles

Status: v5 (EN). Rebased on the existing MPO implementation (§1.1): the spec now describes changes to that code rather than a standalone feature. Correlation plot checked on a live instance (`campo_corrplot_check.js`, 2026-10-26).
Method reference: Garellick et al., *caMPO: Correlation-aware multiparameter optimization*, ChemRxiv 2026, doi:10.26434/chemrxiv.15008121/v1.

Decisions taken:

- the matrix `R` is computed by our own Pearson function on complete rows (§2.5); the built-in Correlation plot shows the same matrix and is not a source of data;
- matrix inversion uses Cholesky decomposition;
- the weight formula follows the paper, with negative values clipped and the rest renormalised;
- the Priority column is kept.

Decisions proposed in v5, to confirm:

- desirability arrays come from the existing `mapColumnDesirability` through the `MpoCalculator` cache, not from new code (§2.3);
- the property's existing `weight` always holds the **effective** weight, so the score, the MPO Score dialog, old clients and other packages need no change; caMPO inputs are stored next to it (§6);
- the mode viewers live in a `ui.tabControl` in the right half of the profile view, which today is a plain `ui.splitH`, not a table view with a dock manager (§3);
- RESET keeps its current meaning, "revert to the last saved profile" (§4);
- the **Weighting** control exists only in the profile view; the MPO Score dialog and the profile preview use the stored weights (§1.2).

### Reference artifacts (same folder)

| File | Purpose |
|---|---|
| `campo_demo.html` | Interactive demo of the feature: layout, weighting modes, priorities, pins, tooltips, viewers, tabs, aggregation warning. Open in a browser. The reference for behaviour and look, with the differences listed below. |
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

Out of scope: the Data-driven method (unchanged; **Weighting** is hidden there, as **Aggregation** already is) and aggregations other than `Average` (§8).

### 1.1. Existing implementation

| File | What it does today | Change |
|---|---|---|
| `libraries/statistics/src/mpo/mpo-types.ts` | `DesirabilityProfile`, `PropertyDesirability` (`weight` 0–1, optional `weightColumn`, `missingValues`), `WEIGHTED_AGGREGATIONS` | `MpoWeighting`, `weighting`, `priority`, `pinnedWeight` (§6) |
| `libraries/statistics/src/mpo/mpo.ts` | `hoistColumns` + `mapColumnDesirability` (value → desirability per column, with missing values, log scale, inversion, qnum, categories), `reduceMpo` (aggregation), `MpoCalculator` (per-column cache) | `MpoCalculator.desirabilities()` (§2.3) |
| `libraries/statistics/src/mpo/campo.ts` | — | new: the math of §2 (pure functions, no UI) |
| `libraries/statistics/src/mpo/mpo-profile-editor.ts` | `MpoProfileEditor`: header Property / Weight / Desirability, one row per property, `aggregationInput`, column mapping (`columnMapping`), design mode (rename, add, delete, compute function), weight-from-column toggle, `onChanged` + `MPO_SCORE_CHANGED_EVENT` | weighting mode, Priority column, read-only weights, pins, name tooltip, `setWeights()` (§5) |
| `libraries/statistics/src/mpo/editors/mpo-property-row.ts` | `createMpoRow`: property cell, weight cell (`ui.input.float`, 0–1, `#0.000`), curve editor | priority cell, weight cell variants (§5) |
| `packages/Chem/src/mpo/mpo-create-profile.ts` | `MpoProfileCreateView`: title, description, Method / Dataset / Aggregation with hints, Save / Reset in the ribbon; content is `ui.splitH([editor, grid])`; Data-driven mode docks EDA viewers into a hidden `tableView` | Weighting input, tab control in the right half, caMPO controller, aggregation banner (§3, §4, §8) |
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

> Verified on a live instance: for `double` columns `getRawData()` returns a `Float32Array`, and `raw[i] === DG.FLOAT_NULL` finds every missing value (112 of 112). `hoistColumns` relies on this (`valNull`); keep it as a unit test.

### 2.4. Complete rows and constant columns

Two steps are performed before `R` is computed.

**Complete rows only.** `R` is computed only on rows where every `D_k` is a number (listwise deletion). Correlations on different row subsets (pairwise deletion) can produce a matrix that is not positive definite — Cholesky decomposition would then fail. The Correlation plot gets the same complete rows: it uses pairwise deletion (verified: with 5% missing values it matches the pairwise reference to 3e-9 and differs from the listwise one by 0.028), so on the full data it would show a matrix different from the one the weights come from.

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
- `DG.CorrelationPlot.getCorrelation(c1, c2)` returns the right value (matches the reference Pearson to 4e-9, also in a hidden tab and right after the data changes), but it recomputes the pair on every call, with no cache: the full 20×20 matrix on 10⁵ rows takes ~1.5 s (measured). The code below does 190 dot products of length 10⁵ for the same matrix — tens of milliseconds by estimate; the §9 performance test confirms it.

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
type StageBResult = { final: Float64Array; error: string | null };

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
    return { final: base, error: `Pinned weights sum to ${P.toFixed(3)} (> 1)` };

  const rest = 1 - P;
  if (rest > 1e-9 && B <= 1e-12)
    return { final: base, error: 'Nothing to scale: all unpinned weights are 0' };

  const k = B > 0 ? rest / B : 0;
  const final = new Float64Array(d);
  for (let i = 0; i < d; i++)
    final[i] = Number.isNaN(pinned[i]) ? b[i] * k : pinned[i];
  return { final, error: null };
}
```

The normalisation after priorities is merged with the pin scaling: `b_i · (1 − P) / Σ b` is the same as `base'_i · (1 − P) / U`.

**Example.** Base weights after clipping: TPSA 0.258, TPSA_S 0.038, HBA 0, HBD 0.023, MW 0, nAtoms 0.103, cLogD 0.363, cLogP 0, mbpKa 0.215; all priorities 1. Pin HBA = 0.08: `P = 0.08`, `Σ b = 1`, `k = 0.92`; TPSA = 0.237, cLogD = 0.334, and so on; the sum is 1.

**Display statuses:** `pinned` — pinned; `dropped` — final weight 0; `kept` — all others.

The final weights go into the profile as the properties' `weight` (§6). They are in [0, 1] and sum to 1, which the weight input's 0–1 range and `reduceMpo` both accept.

> Scientific note. Simple scaling ignores the correlation between a pinned property and the others: if HBA is pinned, TPSA (r = 0.82 with HBA) keeps its share, and polarity is counted twice. A future alternative is conditional re-optimisation with fixed `w_P`: `w_U = μ·R_UU⁻¹·1 − R_UU⁻¹·R_UP·w_P`, where `μ` is chosen so that `Σ w_U = 1 − P`, followed by clipping. It also reduces to element-wise sums. Version 1 uses scaling.

---

## 3. Layout

Reference: `campo_demo.html`; the view is `MpoProfileCreateView` (`packages/Chem/src/mpo/mpo-create-profile.ts`), the same in create and edit mode.

- **Top** (`profileViewContainer`, unchanged except one input): editable title, editable description, then the toolbar with **Method**, **Dataset**, **Aggregation**, each with its hint below (`FIELD_DESCRIPTIONS`). **Weighting** is added after **Aggregation**, on the same row, with the hint `How property weights are set`; the Aggregation hint stays as it is (`How individual property scores combine into the final MPO score`). In Data-driven mode **Weighting** is hidden together with **Aggregation** (`chem-mpo-d-none` in `onMethodChanged`). The aggregation banner (§8) goes at the bottom of the toolbar. **Save** and **Reset** stay in the ribbon.
- **Content**: `ui.splitH([editorPanel, gridPanel])`, as today.
  - **left** — `MpoProfileEditor` in design mode. Its header becomes **Property / Priority / Weight / Desirability**; the Priority column is shown only in `Correlation-aware` mode. The Weight column depends on the mode: `Equal` — read-only `1.000`; `Manual` — editable input and the weight-from-column toggle, as today; `Correlation-aware` — final weight (`0.000`, read-only) and the pin icon (§5.2). In `Equal` and `Correlation-aware` the weight-from-column toggle is hidden (§5.1).
  - **right** — today a grid (`df.plot.grid()`). It becomes a `DG.TabControl` (`ui.tabControl`):
    1. `<dataset name>` — the grid, as today. The MPO score column `MPO <profile name>` is maintained by the context panel's preview (`computeMpo`), which reruns on every weight change through `MPO_SCORE_CHANGED_EVENT`.
    2. **Correlations** — Correlation plot on the complete rows (§2.4), `allowClose: true`.
    3. **Weights** — bar chart "Raw weights (may be negative)" and bar chart "Final weights (sum = 1)" side by side in a `ui.splitH`, `allowClose: true`.
- If **no raw weight was negative** (`clipped === false`), the Raw chart is not shown and **Weights** holds only "Final weights (sum = 1)".
- When `Correlation-aware` is switched on, the **Weights** tab becomes current.
- The user can close **Correlations** or **Weights** with the tab's ×. Both tabs reopen the next time the mode is switched on, or on RESET when the restored mode is `Correlation-aware`.
- When Weighting is switched to `Equal` or `Manual`, both tabs are removed and the grid tab becomes current.
- The context panel (score histogram, best and worst compounds) is unchanged.

**Bar chart titles.** "Raw weights (may be negative)" — `R⁻¹·1` normalised by `S`, before clipping; "Final weights (sum = 1)" — after clipping, priorities and pins, i.e. the weights used in the score. "Clipped & normalized" was rejected: the final weights also include priorities and pins.

> The Correlation plot is display only. Its tab being inactive or closed does not affect the computation.

> `TabControl` puts its tab strip on top; the demo draws it at the bottom. Keep the platform's look.

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

  constructor(private gridRoot: HTMLElement, private gridTitle: string) {
    this.gridPane = this.tabs.addPane(gridTitle, () => gridRoot);
    this.tabs.onTabRemoved.subscribe((pane: DG.TabPane) => {   // the user pressed ×
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

  showWeights(weightsDf: DG.DataFrame, d: number, clipped: boolean, activate: boolean): void {
    if (this.closedByUser.weights) return;
    const first = !this.finalChart;
    if (!this.finalChart) {
      this.finalChart = createBarChart(weightsDf, 'weight', 'status', 'Final weights (sum = 1)', 0, 0.4);
      addEqualLine(this.finalChart, d);
    } else
      this.finalChart.dataFrame = weightsDf;

    if (clipped && !this.rawChart)
      this.rawChart = createBarChart(weightsDf, 'raw', 'sign', 'Raw weights (may be negative)', -0.4, 0.6);
    else if (!clipped)
      this.rawChart = null;
    if (this.rawChart) this.rawChart.dataFrame = weightsDf;

    ui.empty(this.weightsHost);
    this.weightsHost.append(this.rawChart ?
      ui.splitH([this.rawChart.root, this.finalChart.root], {}, true) : this.finalChart.root);
    if (first) {
      const pane = this.tabs.addPane('Weights', () => this.weightsHost, null, {allowClose: true});
      if (activate) this.tabs.currentPane = pane;
    }
  }

  /** Mode switched to Equal or Manual: only the grid stays. */
  hide(): void {
    this.tabs.clear();
    this.gridPane = this.tabs.addPane(this.gridTitle, () => this.gridRoot);
    this.corr = this.rawChart = this.finalChart = null;
  }
}
```

Axis ranges of the charts follow §3.2 (set them with `setOptions({valueMin, valueMax})` on each update). `MpoProfileCreateView.setupGridAndContextPanel` puts `campoViewers.tabs.root` where it puts `gridContainer` today; on dataset change the view is rebuilt as it is now (`attachLayout`), with a new `CampoViewers`.

> To check on a live instance: `onTabRemoved` fires for a pane closed with × and not for `clear()`; if it fires for `clear()` too, `hide()` must set a flag around the call so that the closures are not taken for the user's.

---

## 4. Reactivity

All computation runs automatically; there is no Recompute button.

| Event | Where it comes from today | What is recomputed |
|---|---|---|
| Dataset change | `onDatasetChanged` → `attachLayout` → `editor.setProfile` | stage A + B, viewers rebuilt |
| Property added / removed / renamed, column remapped | `MpoProfileEditor.onChanged` | stage A + B |
| Desirability curve or missing-value handling edited | `MpoProfileEditor.onChanged` (debounced) | stage A + B |
| `Correlation-aware` mode switched on | Weighting input | tabs opened, stage A + B |
| Priority changed (valid value) | Priority input | stage B only |
| Pin set / removed (valid value) | pin icon, pin input | stage B only |
| Switched to `Equal` / `Manual` | Weighting input | tabs removed, no computation (§5.1) |
| RESET | ribbon **Reset** → `resetProfile` | the saved profile is restored with its weighting, priorities and pins; if its mode is `Correlation-aware`, tabs reopen and stage A + B run |

`MpoProfileEditor.onChanged` fires for every edit, including the weights the controller itself writes. The controller tells the cases apart with a fingerprint of what stage A depends on, and runs stage A only when it changed:

```ts
/** Everything stage A depends on; weights, priorities and pins are not part of it. */
function stageAKey(df: DG.DataFrame, profile: DesirabilityProfile, mapping: Record<string, string | null>): string {
  return df.id + '|' + Object.entries(profile.properties)
    .map(([name, p]) => `${name}>${mapping[name] ?? name}>${desirabilityKey(p)}`).join('|');
}
```

```ts
// packages/Chem/src/mpo/campo-controller.ts
import {Subscription} from 'rxjs';
import {debounceTime} from 'rxjs/operators';

class CampoController {
  private baseW: BaseWeights | null = null;
  private keyA = '';
  private sub: Subscription;

  constructor(private view: MpoProfileCreateView, private viewers: CampoViewers,
              private calculator: MpoCalculator) {
    // curve dragging fires many events: debounce stage A
    this.sub = view.editor.onChanged.pipe(debounceTime(150)).subscribe(() => this.runStageA());
  }

  private get names(): string[] { return Object.keys(this.view.profile.properties); }

  runStageA(force = false): void {
    const {df, profile, editor} = this.view;
    if (!df) return;
    const key = stageAKey(df, profile, editor.columnMapping);
    if (!force && key === this.keyA) return;
    this.keyA = key;

    const names = this.names;
    const columns = names.map((nm) => df.col(editor.columnMapping[nm] ?? nm));
    const tagged = applyDesirabilityTags(df, Object.fromEntries(names
      .map((nm, k) => [columns[k]?.name, profile.properties[nm]]).filter(([c]) => c)), true);
    const maps = this.calculator.desirabilities(df, tagged);
    const D = columns.map((c) => (c ? maps[tagged.indexOf(c)].D : null));

    const prep = prepareDesirabilityFrame(names, D, df.rowCount);
    const m = prep.active.length;
    if (m < 2 || prep.rows <= m) {
      this.view.showCampoMessage('Not enough data to compute correlation-aware weights', 'error');
      return;
    }

    this.viewers.showCorr(prep.frame, Array.from(prep.active, (k) => names[k]));
    const R = pearsonMatrix(prep.cols, prep.rows);
    const inv = robustInverse(R, m);
    if (!inv) {
      this.view.showCampoMessage('Correlation matrix is singular', 'error');
      return;
    }
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
    this.view.showCampoMessage(error, 'error');
    if (error) return;
    const status = names.map((_, i) =>
      !Number.isNaN(pinned[i]) ? 'pinned' : final[i] === 0 ? 'dropped' : 'kept');
    this.viewers.showWeights(weightsFrame(names, this.baseW.raw, final, status), names.length,
      this.baseW.clipped, false);
    this.view.editor.setWeights(Object.fromEntries(names.map((nm, i) => [nm, final[i]])), status);
  }

  dispose(): void { this.sub.unsubscribe(); this.baseW = null; this.keyA = ''; }
}
```

- `editor.setWeights` writes `weight` on every property and refreshes the read-only cells without firing a change per property; it then fires `onChanged` and `MPO_SCORE_CHANGED_EVENT` once, so the score, the context panel and the view's "modified" state update as for any edit. The `onChanged` it fires reaches `runStageA`, which returns at once because the key is unchanged.
- Priority and pin inputs write `priority` / `pinnedWeight` on the property and call `runStageB()`.
- The controller is created when the mode is switched on and disposed when it is switched off, on dataset change and when the view is detached (`detach()`).

---

## 5. Form columns

All of this is in `MpoProfileEditor` and `createMpoRow`, behind a mode the host sets (`editor.setWeighting(mode)`); without it the editor behaves as today (EDA, PowerGrid, the MPO Score dialog).

### 5.1. Mode switching and weights

- **`Manual`** — today's editor: the weight input (0–1) and the weight-from-column toggle.
- **`Equal`** — every `weight` is set to 1; the cells show a read-only `1.000`.
- **`Correlation-aware`** — the cells show the final weights (§5.2).
- In `Equal` and `Correlation-aware` the weight-from-column toggle is hidden and `weightColumn` is removed from the properties: per-row weights would override the computed ones in `reduceMpo`.
- **Leaving `Manual`**: the manual weights and `weightColumn` values are stashed in the view for the session. **Back to `Manual`**: they are restored from the stash; without a stash (the profile was opened in another mode) Manual starts from the current effective weights.
- **`Manual → Correlation-aware`**: current manual weights are copied into the priorities (overwriting stored priorities).
- `d = 1` and aggregations other than `Average`: see §7 and §8.

### 5.2. Priority

- Shown only in `Correlation-aware` mode, as a column between Property and Weight (`createMpoRow` gets a `priorityCell`; the editor's header gets `Priority`).
- `ui.input.float`, min 0, no max, step 0.1, default 1; stored as the property's `priority` (§6).
- Validator: not empty → `Enter a priority`; ≥ 0 → `Must be 0 or greater`.
- For a pinned property the field is disabled, with tooltip `Weight is pinned`.
- If the base weight of a property is 0 and its priority ≠ 1, the field shows a hint: `No effect: computed weight is 0. Pin the weight instead.` — the priority multiplies the weight, and zero stays zero.

### 5.3. Weight and pinning

**Behaviour:**

1. Each row of the Weight column shows the final weight (format `0.000`, read-only) and a `ui.iconFA('thumbtack')` icon.
2. **Icon tooltip:** `Pin weight` (not pinned) / `Unpin weight` (pinned). A Datagrok tooltip (`ui.tooltip.bind`), not the browser `title`: shown on hover and on keyboard focus, hidden on click; after the click the next hover shows the new state.
3. **Clicking the icon on an unpinned property:** the number is replaced by a `ui.input.float` field **with no value**, and the field gets focus. Until a value is entered and passes validation, the property is treated as unpinned (the computation does not change).
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

### 5.4. Form messages

Below the property table (`profileEditorContainer`) there is one message line for the mode, built like the view's other messages (`chem-mpo-error-message`):

- errors: `Not enough data…`, `Correlation matrix is singular`, the errors of `finalWeights` — the weights are not applied;
- warnings: regularisation used (`λ = 0.01`), few complete rows (§7).

### 5.5. "Why this weight" tooltip

Hovering (or focusing) a property name in `Correlation-aware` mode shows a tooltip that explains the weight. In the profile view the editor is in design mode, so the name is an input (`buildPropertyNameInput`); the tooltip is bound to the property cell (`ctx.propertyCell`), which holds it.

- the property name and its final weight;
- the share of its desirability explained by the other properties, `R²_i = 1 − 1/R⁻¹[i][i]`;
- the reason, if any: `Pinned at 0.080.` or `Raw weight −0.08 < 0, so it was set to 0.`; for a constant column: `Constant in this dataset — excluded.`; for an unmapped property: `No column in this dataset — excluded.`;
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
| `d = 1` | `Correlation-aware` has nothing to compute: selecting it shows `Needs at least two properties` and keeps the weight at 1 |
| Property with no column in the dataset | excluded from `R`, weight 0, status `dropped`, reason in the name tooltip |
| Constant desirability column | excluded from `R` and from the Correlation plot, weight 0, status `dropped`, reason in the name tooltip |
| Fewer than 2 properties in `R` | error `Not enough data…`, tabs are not shown, weights keep their last values |
| Few complete rows (`rows < 10·m`) | computation runs; the form message warns that weights may be unreliable |
| `rows ≤ m` | `R` is singular → error, weights keep their last values |
| Cholesky failed / `κ > 10⁶` | regularisation `R_λ` (§2.6); the fact and `λ` are shown in the form message |
| No negative raw weight | Raw bar chart is not shown |
| Sum of pinned = 1 | all unpinned weights are 0 (allowed) |
| Sum of pinned < 1, all unpinned weights × priorities = 0 | form error, result not applied |
| All priorities of unpinned properties = 0 and no pins | form error |
| No dataset selected | the profile view shows only the editor (as today); `Correlation-aware` shows `Select a dataset to compute weights` and keeps the stored weights |
| More than `MAX_MPO_PROPERTIES` (10) properties | no limit for caMPO; 10 is only the limit of auto-population |

---

## 8. Aggregation constraint

caMPO weights are derived for a weighted average. The aggregations are `Average` (default), `Sum`, `Product`, `Geomean`, `Min`, `Max` (`WEIGHTED_AGGREGATIONS`); with any of them other than `Average` the weights no longer mean what caMPO computed (`Product`, `Min` and `Max` use them as exponents).

`ui.input.choice` cannot disable one item, so the mode is not blocked. Whenever Weighting is `Correlation-aware` and Aggregation is not `Average` — whichever of the two was changed last — a warning is shown at the bottom of the toolbar: `Correlation-aware weights are derived for a weighted average. With this aggregation they no longer remove redundancy.` with two buttons: **Switch to Average** (sets `editor.aggregationInput`) and **Use Manual** (switches Weighting to `Manual`). Until the user chooses, the score is computed with the selected aggregation and the current weights.

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
- `MpoCalculator.desirabilities`: returns the arrays `compute()` caches (same object on a second call with an unchanged curve); missing values of a `double` column come out as `NaN` under `exclude` and `skip`, as the default score under `default`.
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
4. ~~**Closing tabs.**~~ Resolved for docked viewers: `grok.events.onViewerClosed` gives `args.args.viewer` (match by `.dart`), `dockManager.onClosed` gives `viewer.root`; `close()` on a viewer that was never docked throws. With the tab control (§3.3) the signal is `TabControl.onTabRemoved`; what it reports for `clear()` is still to check.
5. **Applying the profile to other tables.** Proposed (§1.2, §6): the stored effective weights are used as they are; weights are recomputed only in the profile view, on its dataset. To confirm.
6. **Pin semantics.** Version 1 uses simple scaling (§2.8). Consider conditional re-optimisation in version 2.
7. **Manual weights after a save in another mode.** With `weight` holding the effective weight (§6), the manual weights survive a mode switch only within the session (§5.1). If they must survive a save in `Equal` / `Correlation-aware` mode, add `manualWeight` to the property.
