# Specification: Correlation-aware (caMPO) weighting mode in MPO Profiles

Status: v3 (EN). Updated to match the reference demo `campo_demo.html`.
Method reference: Garellick et al., *caMPO: Correlation-aware multiparameter optimization*, ChemRxiv 2026, doi:10.26434/chemrxiv.15008121/v1.

Decisions taken:

- the matrix `R` is read from the `table` field of the built-in Correlation plot viewer;
- matrix inversion uses Cholesky decomposition;
- the weight formula follows the paper, with negative values clipped and the rest renormalised;
- viewers are docked in `DG.DOCK_TYPE.FILL` mode;
- the Priority column is kept.

### Reference artifacts (same folder)

| File | Purpose |
|---|---|
| `campo_demo.html` | Interactive demo of the feature with the approved UI: layout, weighting modes, priorities, pins, tooltips, viewers, bottom tabs, aggregation warning. Open in a browser. This is the reference for behaviour and look. |
| `campo_bench.html` | Verification bench for the computation: presets from the paper, generated data, intermediate values (s_i, raw, base, final, R⁻¹), the §9 test suite. |

The math in both files is the code from §2, ported 1:1 to JavaScript. Differences from the Datagrok implementation:

- `R` is computed by our own Pearson function on complete rows — a stand-in for reading the Correlation plot `table` field (§2.5);
- data is generated: 400 compounds with drugs_props-like columns and fixed desirability curves (not editable in the demo);
- docking, tabs and inputs are imitated with HTML; SAVE does not persist anything; the Data-driven method is not implemented.

---

## 1. Goal and scope

The MPO profile editor (**Manual** method) gets a choice of how weights are set, **Weighting**: `Equal` / `Manual` / `Correlation-aware`. In `Correlation-aware` mode, weights are computed from the correlations between desirability columns: redundant properties get less weight, independent ones get more. An expert can adjust the result in two ways: a priority (a multiplier on the weight) and pinning a specific weight value.

Out of scope: the Data-driven method (unchanged) and aggregations other than `Average` (see §8).

---

## 2. Computation

### 2.1. Notation

| Symbol | Meaning |
|---|---|
| `n` | number of dataset rows |
| `d` | number of properties in the profile |
| `m` | number of properties included in `R` (constant columns excluded) |
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
dataset + curves → D_k → desirability data frame →          base weights × priorities
Correlation plot → R (from table) → R⁻¹ → base weights   →  pins → final weights

triggers: dataset change, property set change,              triggers: priority change,
desirability curve edit, mode switched on                   pin change
```

All computation runs on typed arrays. From Datagrok we read `col.getRawData()` once per column; no `col.get(i)` and no iteration over `df.rows`. Results go back to Datagrok with one `DG.Column.fromFloat32Array` call per column.

### 2.3. Desirability columns

```ts
interface PropSpec {
  name: string;
  column: DG.Column;                     // source property column
  desirability: (x: number) => number;   // desirability function from the profile
}

/** One Float32Array per property; missing values become NaN. */
function desirabilityArrays(props: PropSpec[], n: number): Float32Array[] {
  return props.map((p) => {
    const raw = p.column.getRawData();   // Float32Array | Float64Array | Int32Array
    const nullValue = p.column.type === DG.TYPE.INT ? DG.INT_NULL : DG.FLOAT_NULL;
    const out = new Float32Array(n);
    for (let i = 0; i < n; i++)
      out[i] = raw[i] === nullValue ? NaN : p.desirability(raw[i]);
    return out;
  });
}
```

> Cover with a test: for `double` columns `getRawData()` returns a `Float32Array`, and comparison with `DG.FLOAT_NULL` matches exactly.

### 2.4. Complete rows and constant columns

Two steps are performed before the data is passed to the viewer.

**Complete rows only.** The desirability data frame contains only rows with no missing values in any desirability column. Otherwise the viewer would compute correlations on different row subsets (pairwise deletion), and the matrix could end up not positive definite — Cholesky decomposition would then fail.

**Constant columns.** If a desirability column has zero variance (for example, every row is 1.0), its correlation is undefined. Such a property is not passed to the viewer, gets weight 0 and status `dropped`; the reason is shown in the row tooltip.

```ts
interface Prepared {
  frame: DG.DataFrame;   // complete rows, active columns only — input for Correlation plot
  active: Int32Array;    // indices (in props) of columns included in R
  constant: Int32Array;  // indices excluded because variance == 0
  rows: number;          // number of complete rows
}

function prepareDesirabilityFrame(names: string[], D: Float32Array[], n: number): Prepared {
  const d = D.length;

  // complete-case mask
  const mask = new Uint8Array(n).fill(1);
  for (let k = 0; k < d; k++) {
    const a = D[k];
    for (let i = 0; i < n; i++) if (Number.isNaN(a[i])) mask[i] = 0;
  }
  let rows = 0;
  for (let i = 0; i < n; i++) rows += mask[i];

  const activeList: number[] = [];
  const constantList: number[] = [];
  const cols: DG.Column[] = [];

  for (let k = 0; k < d; k++) {
    const a = D[k];
    const out = new Float32Array(rows);
    let sum = 0;
    for (let i = 0, r = 0; i < n; i++) if (mask[i]) { out[r++] = a[i]; sum += a[i]; }
    const mean = sum / rows;
    let v = 0;
    for (let r = 0; r < rows; r++) { const t = out[r] - mean; v += t * t; }

    if (v > 1e-12) {
      activeList.push(k);
      cols.push(DG.Column.fromFloat32Array(names[k], out));
    } else {
      constantList.push(k);
    }
  }
  return {
    frame: DG.DataFrame.fromColumns(cols),
    active: Int32Array.from(activeList),
    constant: Int32Array.from(constantList),
    rows,
  };
}
```

### 2.5. Matrix R from Correlation plot

The Correlation plot is built on `prepared.frame`. The matrix computed by the viewer is available in its `table` field; we read `R` from it into a `Float64Array` (row-major, `m×m`), matching rows and columns by property names.

```ts
/**
 * Reads R from Correlation plot's `table` field.
 * Assumed layout (to verify): one string column with row labels
 * plus one numeric column per property.
 */
function matrixFromCorrPlot(corr: DG.Viewer, names: string[]): Float64Array {
  const t = (corr as unknown as { table: DG.DataFrame }).table;
  const m = names.length;

  const labelCol = t.columns.toList().find((c) => c.type === DG.TYPE.STRING)!;
  const labels = labelCol.toList() as string[];            // m values, cheap
  const rowOf = new Map(labels.map((l, i) => [l, i]));

  const R = new Float64Array(m * m);
  for (let q = 0; q < m; q++) {
    const raw = t.col(names[q])!.getRawData();
    for (let p = 0; p < m; p++) {
      const v = raw[rowOf.get(names[p])!];
      R[p * m + q] = v === DG.FLOAT_NULL ? NaN : v;
    }
  }
  for (let k = 0; k < m; k++) R[k * m + k] = 1;             // diagonal may hold no value
  return R;
}
```

**Checks after reading:** the matrix is symmetric (`|R[p][q] − R[q][p]| < 1e-6`), all off-diagonal values are in [−1, 1], and there is no `NaN`. If a check fails, log an error to the console and show a message in the form tooltip; weights are not recomputed.

**When to read.** `R` is read after the viewer has recomputed the matrix for the new data frame. The notification mechanism needs to be confirmed (see §10, question 1). For now we wait until the table matches the current set of columns:

```ts
async function waitCorrTable(corr: DG.Viewer, names: string[], timeoutMs = 2000): Promise<boolean> {
  const start = performance.now();
  while (performance.now() - start < timeoutMs) {
    const t = (corr as unknown as { table?: DG.DataFrame }).table;
    if (t && names.every((nm) => t.col(nm) != null) && t.rowCount === names.length) return true;
    await new Promise((r) => requestAnimationFrame(r));
  }
  return false;
}
```

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

**Regularisation.** If `cholesky` returns `null` or the condition number exceeds `κ > 10⁶`, shrinkage `R_λ = (1−λ)·R + λ·I` is applied; `λ` is tried in order `[1e-3, 1e-2, 5e-2, 0.1]` until the decomposition succeeds. The fact of regularisation and the value of `λ` are shown in the form tooltip.

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
  raw: Float64Array;      // s_i / S, may be negative; NaN for constant columns
  base: Float64Array;     // clipped and renormalised, sums to 1; 0 for constant columns
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

Diagnostics for the row tooltip: the share of a property explained by the others, `R²_i = 1 − 1/R⁻¹[i][i]` ("95% of MW is explained by other properties").

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

> Scientific note. Simple scaling ignores the correlation between a pinned property and the others: if HBA is pinned, TPSA (r = 0.82 with HBA) keeps its share, and polarity is counted twice. A future alternative is conditional re-optimisation with fixed `w_P`: `w_U = μ·R_UU⁻¹·1 − R_UU⁻¹·R_UP·w_P`, where `μ` is chosen so that `Σ w_U = 1 − P`, followed by clipping. It also reduces to element-wise sums. Version 1 uses scaling.

---

## 3. Layout

Reference: `campo_demo.html`.

- **Top panel** (full width): profile title and description; **Method**; below it **Dataset**; below it **Aggregation** and **Weighting** on one row. Each input keeps its italic hint (Aggregation: `How property scores combine into the MPO score`; Weighting: `How property weights are set`). The aggregation warning (§8) is shown at the bottom of this panel.
- **Bottom panel**, split in two:
  - **left** — the property table **Property / Priority / Weight / Desirability**, with its own scroll. The Priority column is shown only in `Correlation-aware` mode. The Weight column depends on the mode: `Equal` — read-only `1.000`; `Manual` — editable input, as today; `Correlation-aware` — final weight (`0.000`, read-only) and the pin icon (§5.2).
  - **right** — the table view: the grid and the mode-specific viewers, docked with `DG.DOCK_TYPE.FILL`, so they appear as tabs on the standard Datagrok tab strip at the bottom:
    1. `<dataset name>` (for example `drugs_props`) — the grid. The MPO score is the last column, `MPO <profile name>`, recomputed on every weight change.
    2. **Correlations** — Correlation plot, `FILL` relative to the grid.
    3. **Weights** — bar chart "Raw weights (may be negative)" docked `FILL` relative to the grid, and bar chart "Final weights (sum = 1)" docked `RIGHT` of it, in the same tab.
- If **no raw weight was negative** (`clipped === false`), the Raw chart is not shown, and "Final weights (sum = 1)" is docked `FILL` relative to the grid as the **Weights** tab.
- When `Correlation-aware` is switched on, the **Weights** tab becomes active.
- The user can close **Correlations** or **Weights** with the tab's ×. The viewer is closed; both tabs reopen the next time the mode is switched on, or on RESET.
- When Weighting is switched to `Equal` or `Manual`, both tabs are closed and the grid tab becomes active.

**Bar chart titles.** "Raw weights (may be negative)" — `R⁻¹·1` normalised by `S`, before clipping; "Final weights (sum = 1)" — after clipping, priorities and pins, i.e. the weights used in the score. "Clipped & normalized" was rejected: the final weights also include priorities and pins.

> `campo_demo.html` shows the expected result of docking `RIGHT` inside a `FILL` tab: Raw and Final side by side in the **Weights** tab. Confirm on a live instance (§10, question 3).

> The Correlation plot may be in an inactive (not rendered) tab. Confirm that its `table` field is recomputed for a hidden viewer as well (§10, question 1).

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

> The `weight` column is for display only. Scoring and tooltips always use `final`. Constant columns (`raw = NaN`) are not shown in the Raw chart: before building it, filter out rows with `NaN` in `raw` or write `DG.FLOAT_NULL`.

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

> Verified on a live instance: horizontal orientation, negative values, stacking by category, legend, the Equal line, and the `EPS` stub for zeros. Not verified: `correlationType` on Correlation plot, hiding the selectors.

### 3.3. Docking and closing

```ts
class CampoViewers {
  corr: DG.Viewer | null = null;
  private rawChart: DG.Viewer | null = null;
  private finalChart: DG.Viewer | null = null;
  private rawNode: DG.DockNode | null = null;
  private finalNode: DG.DockNode | null = null;

  constructor(private tv: DG.TableView) {}

  private get gridNode(): DG.DockNode {
    return this.tv.dockManager.findNode(this.tv.grid.root);
  }

  openCorr(desirDf: DG.DataFrame, names: string[]): DG.Viewer {
    if (!this.corr) {
      this.corr = createCorrPlot(desirDf, names);
      this.tv.dockManager.dock(this.corr, DG.DOCK_TYPE.FILL, this.gridNode, 'Correlations');
    } else {
      this.corr.dataFrame = desirDf;
      this.corr.setOptions({ xColumnNames: names, yColumnNames: names });
    }
    return this.corr;
  }

  updateWeights(weightsDf: DG.DataFrame, d: number, clipped: boolean): void {
    const dm = this.tv.dockManager;

    // Raw chart appears / disappears depending on clipping
    if (clipped && !this.rawChart) {
      this.finalChart?.close();                       // re-dock Final to the right of Raw
      this.finalChart = null;
      this.rawChart = createBarChart(weightsDf, 'raw', 'sign', 'Raw weights (may be negative)', -0.4, 0.6);
      this.rawNode = dm.dock(this.rawChart, DG.DOCK_TYPE.FILL, this.gridNode, 'Weights');
    } else if (!clipped && this.rawChart) {
      this.rawChart.close();
      this.rawChart = null;
      this.rawNode = null;
      this.finalChart?.close();                       // re-dock Final on its own
      this.finalChart = null;
    }

    if (!this.finalChart) {
      this.finalChart = createBarChart(weightsDf, 'weight', 'status', 'Final weights (sum = 1)', 0, 0.4);
      addEqualLine(this.finalChart, d);
      this.finalNode = this.rawNode
        ? dm.dock(this.finalChart, DG.DOCK_TYPE.RIGHT, this.rawNode, 'Final weights (sum = 1)', 0.5)
        : dm.dock(this.finalChart, DG.DOCK_TYPE.FILL, this.gridNode, 'Weights');
    } else {
      this.finalChart.dataFrame = weightsDf;
    }
    if (this.rawChart) this.rawChart.dataFrame = weightsDf;
  }

  /** Called when the user closed a viewer's tab with ×. */
  forget(v: DG.Viewer): void {
    if (v === this.corr) this.corr = null;
    if (v === this.rawChart) { this.rawChart = null; this.rawNode = null; }
    if (v === this.finalChart) { this.finalChart = null; this.finalNode = null; }
  }

  close(): void {
    for (const v of [this.rawChart, this.finalChart, this.corr]) v?.close();
    this.corr = this.rawChart = this.finalChart = null;
    this.rawNode = this.finalNode = null;
  }
}
```

When the user closes a tab with ×, the controller must forget the closed viewer, so that the next update does not touch it and the next switch to `Correlation-aware` reopens it:

```ts
// The user closed a tab with ×: forget the viewer; it reopens on the next switch to Correlation-aware.
// Event name to confirm against the current datagrok-api.
grok.events.onViewerClosed.subscribe((args: DG.EventData) => {
  const v = args.args.viewer as DG.Viewer;
  campoViewers.forget(v);   // sets the matching field (corr / rawChart / finalChart) and its node to null
});
```

After `forget`, `updateWeights` must not recreate a closed chart until the mode is switched on again (keep a `closedByUser` flag per tab, as `campo_demo.html` does).

> Closing **Correlations** removes the source of `R` (§2.5). Options: (a) closing the tab only hides it — the viewer stays alive and keeps computing; (b) after the tab is closed, stage A falls back to our own Pearson computation on the same complete rows. See §10, question 4.

---

## 4. Reactivity

All computation runs automatically; there is no Recompute button.

| Event | What is recomputed |
|---|---|
| Dataset change | stage A + B, all viewers updated |
| Property added / removed | stage A + B |
| Desirability curve edited | stage A + B (debounced) |
| `Correlation-aware` mode switched on | viewers opened, stage A + B |
| Priority changed (valid value) | stage B only |
| Pin set / removed (valid value) | stage B only |
| Switched to `Equal` / `Manual` | viewers closed, no computation |
| RESET | priorities = 1, pins cleared, aggregation `Average`, mode `Correlation-aware`, closed tabs reopened; stage A + B |

```ts
import { Subject } from 'rxjs';
import { debounceTime } from 'rxjs/operators';

class CampoController {
  private baseW: BaseWeights | null = null;
  private names: string[] = [];
  prio = new Float64Array(0);       // default 1
  pinned = new Float64Array(0);     // NaN = not pinned
  private recomputeA$ = new Subject<void>();
  private viewers: CampoViewers;
  onWeights?: (final: Float64Array, status: string[]) => void;
  onError?: (message: string | null) => void;

  constructor(private tv: DG.TableView, private props: () => PropSpec[]) {
    this.viewers = new CampoViewers(tv);
    // curve dragging fires many events: debounce stage A
    this.recomputeA$.pipe(debounceTime(150)).subscribe(() => this.runStageA());
  }

  requestStageA(): void { this.recomputeA$.next(); }

  private async runStageA(): Promise<void> {
    const props = this.props();
    const d = props.length;
    this.names = props.map((p) => p.name);
    if (this.prio.length !== d) this.prio = new Float64Array(d).fill(1);
    if (this.pinned.length !== d) this.pinned = new Float64Array(d).fill(NaN);

    const n = this.tv.dataFrame.rowCount;
    const D = desirabilityArrays(props, n);
    const prep = prepareDesirabilityFrame(this.names, D, n);
    const activeNames = Array.from(prep.active, (k) => this.names[k]);
    const m = activeNames.length;

    if (m < 2 || prep.rows <= m) {
      this.onError?.('Not enough data to compute correlation-aware weights');
      return;
    }

    const corr = this.viewers.openCorr(prep.frame, activeNames);
    if (!(await waitCorrTable(corr, activeNames))) {
      this.onError?.('Correlation matrix is not available');
      return;
    }
    const R = matrixFromCorrPlot(corr, activeNames);
    const inv = robustInverse(R, m);
    if (!inv) {
      this.onError?.('Correlation matrix is singular');
      return;
    }
    this.baseW = campoWeights(inv.Rinv, m, prep.active, d);
    this.runStageB();
  }

  /** Cheap: priorities and pins only. */
  runStageB(): void {
    if (!this.baseW) return;
    const { final, error } = finalWeights(this.baseW.base, this.prio, this.pinned);
    this.onError?.(error);
    if (error) return;
    const status = this.names.map((_, i) =>
      !Number.isNaN(this.pinned[i]) ? 'pinned' : final[i] === 0 ? 'dropped' : 'kept');
    const wDf = weightsFrame(this.names, this.baseW.raw, final, status);
    this.viewers.updateWeights(wDf, this.names.length, this.baseW.clipped);
    this.onWeights?.(final, status);
  }

  setPriority(i: number, value: number): void { this.prio[i] = value; this.runStageB(); }
  setPin(i: number, value: number | null): void { this.pinned[i] = value ?? NaN; this.runStageB(); }
  close(): void { this.viewers.close(); this.baseW = null; }
}
```

> Connect the editor events (curve change, dataset change, property set change) to `requestStageA()` in the existing profile editor code.

---

## 5. Form columns

### 5.1. Priority

- Shown only in `Correlation-aware` mode.
- `ui.input.float`, min 0, no max, step 0.1, default 1.
- Validator: not empty → `Enter a priority`; ≥ 0 → `Must be 0 or greater`.
- For a pinned property the field is disabled, with tooltip `Weight is pinned`.
- If the computed weight of a property is 0 and its priority ≠ 1, the field shows a hint tooltip: `No effect: computed weight is 0. Pin the weight instead.` — the priority multiplies the weight, and zero stays zero.
- **Switching `Manual → Correlation-aware`:** current manual weights are copied into priorities. **Switching back:** manual weights are restored unchanged.

### 5.2. Weight and pinning

**Behaviour:**

1. Each row of the Weight column shows the final weight (format `0.000`, read-only) and a `ui.iconFA('thumbtack')` icon.
2. **Icon tooltip:** `Pin weight` (not pinned) / `Unpin weight` (pinned). A Datagrok tooltip (`ui.tooltip.bind`), not the browser `title`: shown on hover and on keyboard focus, hidden on click; after the click the next hover shows the new state.
3. **Clicking the icon on an unpinned property:** the number is replaced by a `ui.input.float` field **with no value**, and the field gets focus. Until a value is entered and passes validation, the property is treated as unpinned (the computation does not change).
4. **Entering a valid value:** the property is pinned → stage B: unpinned weights are scaled (§2.8), the Weight column and the Final bar chart are updated.
5. **Clicking the icon on a pinned property:** the pin is removed, the field is replaced by the number → stage B.

**Input field:**

- min 0, max 1, step 0.01, format `0.000`.
- Validator (the message is the first failed condition):
  1. the value is not empty → `Enter a weight`;
  2. the value is in [0, 1] → `Must be between 0 and 1`;
  3. the sum of all pinned values ≤ 1 → `Pinned weights sum to 1.07 (> 1)`.
- Condition 3 spans rows: when any pinned field changes, **all** pinned fields are revalidated.
- The form-level error from `finalWeights` ("all unpinned weights are 0, nothing to distribute the rest to") is shown below the property table.

```ts
interface PinState {
  inputs: Map<number, DG.InputBase<number | null>>;
  final: () => Float64Array;
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

function weightCell(i: number, ctl: CampoController, state: PinState): HTMLElement {
  const host = ui.divH([], { style: { alignItems: 'center', gap: '8px' } });

  const render = () => {
    ui.empty(host);
    const input = state.inputs.get(i);
    const valueEl = input ? input.root : ui.span([state.final()[i].toFixed(3)]);
    const icon = ui.iconFA('thumbtack', () => toggle());
    if (input) icon.classList.add('fas');          // filled when pinned, outline otherwise
    ui.tooltip.bind(icon, () => (state.inputs.has(i) ? 'Unpin weight' : 'Pin weight'));
    host.append(valueEl, icon);
  };

  const toggle = () => {
    if (state.inputs.has(i)) {
      state.inputs.delete(i);
      revalidateAll(state);
      ctl.setPin(i, null);
    } else {
      const input = ui.input.float('', {
        value: null, min: 0, max: 1, step: 0.01, nullable: false,
        onValueChanged: () => {
          revalidateAll(state);
          ctl.setPin(i, input.validate() ? input.value : null);
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
      state.inputs.set(i, input);
      setTimeout(() => input.input.focus());
    }
    render();
  };

  render();
  return host;
}
```

### 5.3. "Why this weight" tooltip

Hovering (or focusing) a property name in `Correlation-aware` mode shows a tooltip that explains the weight:

- the property name and its final weight;
- the share of its desirability explained by the other properties, `R²_i = 1 − 1/R⁻¹[i][i]`;
- the reason, if any: `Pinned at 0.080.` or `Raw weight −0.08 < 0, so it was set to 0.`; for a constant column: `Constant in this dataset — excluded.`;
- the three most correlated properties (by |r|) with their r values; the demo shows each with a small red (positive) or blue (negative) bar.

```ts
interface WhyInput {
  names: string[];
  R: Float64Array; m: number; active: Int32Array;   // R over active properties
  Rinv: Float64Array;
  raw: Float64Array; final: Float64Array; pinned: Float64Array;
}

function whyTooltip(i: number, x: WhyInput): HTMLElement {
  const pos = Array.from(x.active).indexOf(i);
  const lines: HTMLElement[] = [ui.divText(`${x.names[i]} · final weight ${x.final[i].toFixed(3)}`, { style: { fontWeight: '500' } })];
  if (pos < 0) {
    lines.push(ui.divText('Constant in this dataset — excluded.'));
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
  lines.push(ui.divText('Most correlated', { style: { color: 'var(--grey-4)', marginTop: '6px' } }));
  for (const o of partners) lines.push(ui.divText(`${x.names[o.k]}   ${o.r.toFixed(2)}`));
  return ui.divV(lines);
}

// in the property table: ui.tooltip.bind(nameElement, () => whyTooltip(i, currentWhyInput()));
```

---

## 6. Persistence in the profile

The following fields are added to the profile JSON:

```json
{
  "weighting": "correlation",
  "priority": { "cLogD_ACD_v15": 1.5 },
  "pinned": { "HBA": 0.08 }
}
```

- `weighting`: `"equal"` | `"manual"` | `"correlation"`.
- `priority`: only values different from 1.
- `pinned`: only valid pins.
- `priority` and `pinned` are used only in `correlation` mode; when switching to other modes they are kept but not applied.
- Manual weights of the `Manual` mode are stored as today and are not overwritten when switching modes.

---

## 7. Edge cases

| Case | Behaviour |
|---|---|
| `d = 1` | `Correlation-aware` mode is unavailable (weight = 1): the option is disabled with a tooltip |
| Constant desirability column | not passed to the viewer, weight 0, status `dropped`, reason in the tooltip |
| Fewer than 2 non-constant columns | error `Not enough data…`, viewers are not shown |
| Few complete rows (`rows < 10·m`) | computation runs; the form tooltip warns that weights may be unreliable |
| `rows ≤ m` | `R` is singular → mode unavailable for this dataset, message shown |
| Cholesky failed / `κ > 10⁶` | regularisation `R_λ` (§2.6); the fact and `λ` are shown in the tooltip |
| Matrix from the viewer failed the checks (§2.5) | weights are not recomputed, message in the tooltip |
| No negative raw weight | Raw bar chart is not shown |
| Sum of pinned = 1 | all unpinned weights are 0 (allowed) |
| Sum of pinned < 1, all unpinned weights × priorities = 0 | form error, result not applied |
| All priorities of unpinned properties = 0 and no pins | form error |

---

## 8. Aggregation constraint

caMPO weights are derived for a weighted average. With any aggregation other than `Average`, the `Correlation-aware` mode is unavailable: the option is disabled with tooltip `Available for Average aggregation only`. If the user changes the aggregation while the mode is active, a warning is shown at the bottom of the top panel: `Correlation-aware weights are derived for a weighted average. With this aggregation they no longer remove redundancy.` with two buttons: **Switch to Average** (sets the aggregation back) and **Use Manual** (switches Weighting to `Manual`). Until the user chooses, the score is computed with the selected aggregation and the current weights.

---

## 9. Tests

**Synthetic cases from the paper** (the matrix is passed directly to `robustInverse` + `campoWeights`):

| Matrix | Expected weights |
|---|---|
| 4 properties: A independent, B–D pairwise r = 0.98 | A ≈ 0.497, B = C = D ≈ 0.168 |
| 3 properties: A independent, r(B, C) = −0.35 | A ≈ 0.245, B = C ≈ 0.377 |
| `d×d` identity matrix | all weights = 1/d |

**Other:**

- `matrixFromCorrPlot`: on a known data frame, values match a reference Pearson computation within `1e-6`; the matrix is symmetric.
- `cholesky`: returns `null` on a non-positive-definite matrix; `robustInverse` finds a working `λ`.
- `finalWeights`: the example from §2.8; priority 2 on one property; sum of pinned > 1 → error; = 1 → unpinned are 0.
- A constant column is excluded without exceptions and without `NaN` in final weights.
- Switching modes closes and reopens viewers without leaking subscriptions.
- The Raw chart appears / disappears as data changes: Final stays to the right of Raw or takes its place.
- Performance: `n = 10⁵`, `d = 20` — stage A < 300 ms, excluding the viewer's own computation.

---

## 10. Open questions

1. **Updating the `table` field of Correlation plot.** Is there an event signalling that the matrix was recomputed (for now we wait for a matching table, §2.5)? Is `table` recomputed when the viewer sits in an inactive `FILL` tab? If not, fallback: our own Pearson computation on the same complete rows (code from spec v1).
2. **Structure of `table`.** Confirm the layout: one string column with row labels + one numeric column per property; what is stored on the diagonal.
3. **Docking `RIGHT` inside a `FILL` tab** (§3): the expected result is shown in `campo_demo.html`; confirm on a live instance.
4. **Closing tabs** (§3.3): confirm the event name (`grok.events.onViewerClosed`) and its payload; choose what closing **Correlations** means for the source of `R` — hide the viewer or fall back to our own Pearson computation.
5. **Applying the profile to other tables.** Recompute weights on the new dataset or use the stored ones? Recomputing on small tables and virtual molecules is unreliable; the proposal is to store a snapshot of final weights in the profile.
6. **Pin semantics.** Version 1 uses simple scaling (§2.8). Consider conditional re-optimisation in version 2.
