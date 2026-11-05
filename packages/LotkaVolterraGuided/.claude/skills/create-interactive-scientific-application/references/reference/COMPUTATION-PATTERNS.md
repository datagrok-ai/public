# Computation Patterns

Reference for implementing computational methods in Datagrok packages.

> **End-to-end example:** core `../examples/lotka-volterra-spec/code/src/lotka-volterra/core.ts` + `../examples/lotka-volterra-spec/code/src/lotka-volterra/model.ts`.

For worker-based methods, see `WORKER-GUIDE.md`.
For array allocation and reuse patterns, see `ARRAY-OPERATIONS.md`.
For architectural context (core tasks, pipelines, specification structure), see `../guide.md` and `../spec-template.md`.

## Raw Typed Arrays

Access column data via `col.getRawData()` instead of per-element `col.get(i)`. This returns the underlying
typed array (`Float32Array`, `Float64Array`, `Int32Array`, `Uint32Array`) and avoids boxing/unboxing overhead on every iteration.

**IMPORTANT:** The raw array's `.length` may be larger than the column's element count (due to internal buffer allocation). Always use `col.length` for iteration bounds, never `rawData.length`.

```typescript
const vals: Float32Array = features.getRawData();
const cats: Int32Array = categories.getRawData();
const len = features.length; // use column length, NOT vals.length

for (let i = 0; i < len; i++) {
  // direct access — no per-element API calls
  const value = vals[i];
  const category = cats[i];
}
```

Reference: `anova/anova-tools.ts` (`FactorizedData` constructor), `missing-values-imputation/knn-imputer.ts` (`KnnImputer.transform`).

---

## Missing Values Strategy

**Before implementing any new method, define the missing values strategy.** Different methods require different approaches:

| Strategy | When to use | Example |
|----------|-------------|---------|
| **Skip** | Aggregation, statistics | ANOVA skips rows where factor or value is null |
| **Impute before computation** | Methods that require complete data (e.g., matrix operations) | KNN imputation, mean/median fill |
| **Propagate** | Result column should reflect original nulls | Copy null sentinel to output at the same index |
| **Reject** | Method cannot handle nulls at all | Throw error if `missingValueCount > 0` |

Document the chosen strategy in the method's JSDoc or function header. When multiple input columns are involved, specify per-column behavior (e.g., ANOVA: skip if factor OR value is null; KNN: skip feature columns with nulls at the target row but impute the target).

---

## Null Handling in Loops

Before processing, check `col.stats.missingValueCount`. If there are no missing values, skip all null checks
in the loop — this eliminates a branch per iteration and significantly speeds up computation.

```typescript
const hasMissing = col.stats.missingValueCount > 0;

if (hasMissing) {
  const nullValue = getNullValue(col);
  for (let i = 0; i < len; i++) {
    if (raw[i] === nullValue) continue;
    // process raw[i]
  }
} else {
  for (let i = 0; i < len; i++) {
    // process raw[i] — no null checks needed
  }
}
```

When nulls are present, use `getNullValue(col)` from `utils.ts` to obtain the sentinel value and compare
against it directly in loops. Do not use platform null-checking APIs in hot paths.

| Column type | Sentinel | Notes |
|-------------|----------|-------|
| `int`, `string`, `bool` | `-2147483648` | Min 32-bit int |
| `float`, `datetime`, `qnum` | `2.6789344063684636e-34` | Special float constant |

```typescript
import {getNullValue} from '../utils';

const nullValue = getNullValue(col);
const raw = col.getRawData();

for (let i = 0; i < col.length; i++) {
  if (raw[i] === nullValue) continue; // skip missing
  // process raw[i]
}
```

For categorical (string) columns, raw data stores integer category indices. Check for null categories separately:

```typescript
const categoriesNull = categories.stats.missingValueCount > 0 ? getNullValue(categories) : -1;

for (let i = 0; i < size; i++) {
  if ((cats[i] === categoriesNull) || (vals[i] === featuresNull)) continue;
  // process non-null pair
}
```

Reference: `anova/anova-tools.ts` (`FactorizedData.setStats`), `missing-values-imputation/knn-imputer.ts` (`KnnImputer.transform`).

---

## Single-Pass Aggregation

Compute all required statistics in one loop over the data. Pre-allocate output buffers as typed arrays.

```typescript
const K = uniqueCategoryCount;
const sums = new Float64Array(K).fill(0);
const sumsOfSquares = new Float64Array(K).fill(0);
const subSampleSizes = new Int32Array(K).fill(0);

for (let i = 0; i < size; i++) {
  const cat = cats[i];
  if (vals[i] === nullValue) continue;

  sums[cat] += vals[i];
  sumsOfSquares[cat] += vals[i] ** 2;
  ++subSampleSizes[cat];
}
```

Reference: `anova/anova-tools.ts`, `FactorizedData.setStats()`.

---

## Bool Column Handling

Bool columns are stored as packed bit arrays. Extract individual bits via bitwise operations:

```typescript
const raw = boolCol.getRawData(); // Uint32Array with packed bits
let catIdx = 0;
let shift = 0;
let packed = raw[0];
const MAX_SHIFT = 8 * raw.BYTES_PER_ELEMENT - 1;

for (let i = 0; i < size; i++) {
  const bit = 1 & (packed >> shift);
  // use `bit` as 0 or 1

  ++shift;
  if (shift > MAX_SHIFT) {
    shift = 0;
    ++catIdx;
    packed = raw[catIdx];
  }
}
```

Reference: `anova/anova-tools.ts` (`FactorizedData.setStats`, bool branch).

---

## Data Locality

Typed arrays store elements contiguously in memory. Sequential access maximizes CPU cache utilization
and enables hardware prefetching. This is a key reason to prefer typed arrays over `number[]` or per-element API calls.

### Why it matters

- **Cache lines**: CPU loads data in 64-byte blocks. One cache line holds 16 `float32` or 8 `float64` values.
  Sequential access means every loaded cache line is fully utilized.
- **Prefetching**: CPU detects sequential access patterns and preloads next cache lines automatically.
  This hides memory latency almost entirely for linear traversals.
- **No boxing**: `number[]` stores boxed values as heap-allocated objects (pointer → header → value).
  Typed arrays store raw values inline — more useful data per cache line.

### Access pattern guidelines

| Pattern | Cache behavior | Use when |
|---------|---------------|----------|
| Sequential typed array traversal | Optimal — prefetcher active, full cache line utilization | Aggregation, statistics, transforms |
| Multiple typed arrays in parallel (`vals[i]`, `cats[i]`) | Good — each array has its own prefetch stream | Multi-column single-pass (ANOVA, KNN distances) |
| Random access to typed array | Cache miss per access — up to 100x slower than sequential | Avoid; restructure if possible |
| `col.get(i)` in a loop | Method call + potential unboxing per element | Avoid in hot loops |

### Column-major vs row-major

When building matrices from multiple columns, the layout determines which access patterns are cache-friendly:

- **Column-major** (`data[i + j * nRows]`): optimal when processing columns independently
  (e.g., centering, scaling, per-feature statistics)
- **Row-major** (`data[i * nCols + j]`): optimal when accessing all features of one row
  (e.g., distance computation, KNN, nearest neighbor search)

Choose the layout that matches the method's primary access pattern. See `WORKER-GUIDE.md`
for `toFlatColumnMajor` and `toFlatRowMajor` helper functions.

### Pre-allocate output buffers

Allocate result arrays once before the loop to avoid repeated allocations and garbage collection:

```typescript
// Good: single allocation
const result = new Float64Array(len);
for (let i = 0; i < len; i++)
  result[i] = vals[i] * scale;

// Bad: growing array triggers re-allocation and copying
const result: number[] = [];
for (let i = 0; i < len; i++)
  result.push(vals[i] * scale);
```

---

## Module Structure

Separate computation from UI into distinct files:

| File pattern | Purpose | Dependencies |
|-------------|---------|-------------|
| `*-tools.ts` | Pure computation on raw data | `datagrok-api` types, `utils.ts`, math libraries |
| `*-ui.ts` or `ui.ts` | Dialog, inputs, validation, visualization | `datagrok-api` UI, computation module |
| `*-constants.ts` or `ui-constants.ts` | Enums, error messages, UI labels | None |

The computation module must not import UI components. This keeps it testable and potentially reusable in workers.
