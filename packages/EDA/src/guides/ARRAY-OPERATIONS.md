# Array Operations Guide

Reference for implementing efficient array operations in the EDA package.
For raw data access and null handling, see `COMPUTATION-PATTERNS.md`.
For worker-specific patterns, see `WORKER-GUIDE.md`.

## Pre-allocate and Reuse

The core principle: allocate buffers once before the loop, reuse them across iterations.
Every `new Float32Array(n)` inside a loop is a hidden cost — allocation + eventual GC pause.

```typescript
// Bad: allocation per iteration
for (let iter = 0; iter < maxIter; iter++) {
  const temp = new Float32Array(n); // GC pressure grows with maxIter
  // ... use temp ...
}

// Good: single allocation, reused across iterations
const temp = new Float32Array(n);
for (let iter = 0; iter < maxIter; iter++) {
  // ... use temp — same memory, zero allocations ...
}
```

Reference: `probabilistic-scoring/nelder-mead.ts` (lines 103-106) — `centroid`, `reflectionPoint`,
`expansionPoint`, `contractionPoint` are allocated once and reused across all Nelder-Mead iterations.

---

## Out-Parameter Pattern

Write results into a caller-provided array instead of allocating and returning a new one.
This gives the caller control over allocation and enables buffer reuse.

```typescript
function add(a: Float32Array, b: Float32Array, out: Float32Array, len: number): void {
  for (let i = 0; i < len; i++) out[i] = a[i] + b[i];
}

const buf = new Float32Array(n);
add(x, y, buf, n);      // buf = x + y
scale(buf, 2.0, buf, n); // buf = 2 * (x + y), in-place
```

Reference: `probabilistic-scoring/nelder-mead.ts` — `fillPoint` and `fillCentroid` write into
pre-allocated arrays, called repeatedly inside the optimization loop.

---

## Scratch Buffers for Iterative Algorithms

When an algorithm runs many iterations, declare all temporary arrays before the loop.

```typescript
// softmax-worker.ts: Z, dZ, dW allocated once before training loop
const Z = new Array<Float32Array>(m);
for (let i = 0; i < m; i++) Z[i] = new Float32Array(c);
const dZ = new Array<Float32Array>(c);
for (let i = 0; i < c; i++) dZ[i] = new Float32Array(m);

for (let iter = 0; iter < iterations; iter++) {
  // Forward/backward pass writes into Z, dZ — zero allocations per iteration
}
```

Reference: `workers/softmax-worker.ts` (lines 48-56).

---

## Local Aliases for Inner Loops

Store a reference to a sub-array in a local variable before the inner loop.
The primary benefit is **readability and reduced index errors**: `wBuf[k] * xBuf[k]` is
clearer than `params[i][k] * X[j][k]`, and there is less chance of mixing up `i`/`j` indices.

> **Note on performance:** Modern V8 often hoists loop-invariant array lookups automatically
> (loop-invariant code motion), so the performance gain may be minimal. Use this pattern
> primarily for clarity in multi-level loops.

```typescript
// Before: dense indexing, easy to confuse i/j
for (let j = 0; j < m; j++)
  for (let k = 0; k < n; k++)
    sum += params[i][k] * X[j][k];

// After: meaningful names, less index juggling
for (let j = 0; j < m; j++) {
  const xBuf = X[j];       // alias, not copy
  const wBuf = params[i];
  for (let k = 0; k < n; k++)
    sum += wBuf[k] * xBuf[k];
}
```

Reference: `workers/softmax-worker.ts` (lines 62-76) — `xBuf`, `wBuf`, `zBuf` aliases in forward propagation.

---

## Accumulation into Pre-allocated Output

Allocate the output array once and accumulate contributions in-place.

```typescript
// regression.ts: prediction = bias + sum(weight_j * feature_j)
const prediction = new Float32Array(samplesCount);
let rawData = features.byIndex(0).getRawData();
const bias = params[featuresCount];

for (let i = 0; i < samplesCount; i++)
  prediction[i] = bias + params[0] * rawData[i];

for (let j = 1; j < featuresCount; j++) {
  rawData = features.byIndex(j).getRawData();
  for (let i = 0; i < samplesCount; i++)
    prediction[i] += params[j] * rawData[i];
}
```

Reference: `regression.ts` (`getPredictionByLinearRegression`, lines 105-120).

---

## Logical Length vs Physical Length

Pre-allocated buffers may have a fixed physical size but a variable logical length.
Track the logical length separately and use it for all iteration bounds.

```typescript
const properIndices = new Uint32Array(featuresCount);
let properIndicesCount = 0;

const getProperIndices = (idx: number) => {
  properIndicesCount = 0; // reset logical length
  for (let i = 0; i < featuresCount; i++) {
    if (featureSource[i][idx] !== featureNullVal[i])
      properIndices[properIndicesCount++] = i;
  }
};

// Later: iterate only over valid elements
for (let i = 0; i < properIndicesCount; i++)
  sum += bufferVector[properIndices[i]];
```

Reference: `missing-values-imputation/knn-imputer.ts` (lines 139-162).

---

## In-Place Transforms

When the input array is no longer needed after the transform, write results directly into it.

```typescript
function normalizeInPlace(arr: Float32Array, len: number, avg: number, stdev: number): void {
  for (let i = 0; i < len; i++) arr[i] = (arr[i] - avg) / stdev;
}
```

Reference: `regression.ts` (`getTestDatasetForLinearRegression`, lines 148-154).

**Caution:** Only use in-place transforms when you own the array. Never modify arrays obtained
via `col.getRawData()` on user data — this mutates the underlying DataFrame column.

---

## Bulk Copy with TypedArray.set()

Use the built-in `set()` instead of a manual loop — engines optimize it to a memcpy-like path.

```typescript
dst.set(src);                        // full copy
dst.set(src, offset);                // copy into dst starting at offset
dst.set(src.subarray(start, end));   // copy a slice (subarray is a zero-copy view)

// Clone raw column data for safe mutation
const clone = new Float32Array(col.getRawData().length);
clone.set(col.getRawData());
```

**Tip:** `subarray(start, end)` returns a zero-copy view — use it to pass a logical slice
to `set()` or to functions that accept a typed array, without allocating.

---

## Array Pool for Variable-Size Buffers

When buffer sizes vary between calls, a pool recycles previously created arrays.
Interface: `acquire(minLen)` returns a buffer of at least `minLen` (contents uninitialized),
`release(arr)` returns it to the pool, `clear()` drops all pooled arrays.

```typescript
const pool = new Float32Pool();

function processChunk(chunkSize: number): void {
  const tmp = pool.acquire(chunkSize);
  // ... compute into tmp ...
  pool.release(tmp);
}

pool.clear(); // after all work is done
```

Guidelines:
- **Always release** — otherwise it degrades to plain allocation.
- **Never read stale contents** — treat as uninitialized, `arr.fill(0)` if needed.
- **Scope the lifetime** — create per invocation and `clear()` when done.
- **Keep it simple** — for fixed-size buffers, plain pre-allocation is better.

---

## Ring Buffer for Fixed-Length History

When an algorithm needs a sliding window of the last N values, pre-allocate N arrays
and use a modular head index — O(1) per step, zero allocations.

```typescript
const HIST_LEN = 5;
const history: Float64Array[] = [];
for (let i = 0; i < HIST_LEN; i++)
  history[i] = new Float64Array(dim);
let head = 0;

computeValues(history[head]);

for (let step = 0; step < totalSteps; step++) {
  const newest = history[head];
  const oldest = history[(head - (HIST_LEN - 1) + HIST_LEN) % HIST_LEN];

  // Advance: overwrite oldest slot — O(1)
  head = (head + 1) % HIST_LEN;
  computeValues(history[head]);
}
```

**Alternative — reference shift** (O(N) per step): when consumers expect `[0]` = newest,
`[N-1]` = oldest, shift references instead. Acceptable for small N.

```typescript
const recycled = history[HIST_LEN - 1];
for (let j = HIST_LEN - 1; j > 0; --j) history[j] = history[j - 1];
history[0] = recycled;
computeValues(history[0]);
```

Reference: `diff-grok` library, `solver-tools/ab5-method.ts` (lines 343-347) — reference
shift with N=5.

---

## Multi-Purpose Scratch Buffers

The same buffer can serve different purposes at different stages within one iteration.
Each stage must fully overwrite the buffer before reading it.

```typescript
const scratch0 = new Float64Array(dim);
const scratch1 = new Float64Array(dim);

while (solving) {
  // Stage 1: Jacobian — fills scratch0, scratch1 entirely
  jacobian(t, y, f, eps, scratch0, scratch1, W);

  // Stage 2: time derivative — overwrites all elements
  tDerivative(t, y, f, eps, scratch0, scratch1, hdT);

  // Stage 3: scratch0 reused as RHS for linear solve
  for (let i = 0; i < dim; i++) scratch0[i] = f0[i] + hdT[i];
  luSolve(L, U, scratch0, luBuf, k1, dim);
}
```

For many stages, use **stage-scoped aliases**: `const rhs = scratch0;` gives semantic
context without misleading names. Both point to the same memory — zero overhead.

> **Aliasing hazard:** Never pass the same buffer as both `src` and `dst` of a single call.
> If the function reads `src` while writing `dst`, aliasing corrupts the result.
> When unsure, use separate buffers — the cost is negligible vs a silent data corruption bug.

Guidelines:
- **Document the reuse** with comments at each stage.
- **Never read previous-stage contents** — each stage must fully overwrite before reading.
- **Watch for aliasing** — never pass the same buffer as both source and destination.

Reference: `diff-grok` library, `solver-tools/mrt-method.ts` (lines 218-219, 257-283).

---

## Summary

| Pattern | When to use | Saves |
|---------|-------------|-------|
| **Pre-allocate and reuse** | Iterative algorithms | N allocations per loop |
| **Out-parameter** | Utility functions called repeatedly | 1 allocation per call |
| **Scratch buffers** | Multi-step computations in a loop | All intermediate arrays per iteration |
| **Local aliases** | Nested loops with array-of-arrays | Index errors and readability |
| **Accumulation into output** | Aggregation from multiple sources | Intermediate result arrays |
| **Logical length** | Variable-size subsets of a fixed buffer | Re-allocation on size change |
| **In-place transforms** | Input no longer needed after transform | 1 output array |
| **Bulk copy (set())** | Copying blocks between typed arrays | Loop overhead; engine-optimized |
| **Array pool** | Variable-size temporary buffers | Repeated allocation of similar arrays |
| **Ring buffer** | Sliding window / fixed-length history | O(1) advance with modular index |
| **Multi-purpose scratch** | Multi-stage algorithms | Extra buffer per stage |
