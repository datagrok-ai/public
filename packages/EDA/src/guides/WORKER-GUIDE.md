# Worker Implementation Guide for EDA ML Methods

Reference for implementing in-worker ML and data analysis methods using the `worker-utils` infrastructure.

## Worker-Utils Infrastructure

### Definitions (`worker-defs.ts`) — no dependencies, safe to import in workers

```typescript
type RawData = Int32Array | Float32Array | Float64Array | Uint32Array;
type ColumnType = 'int' | 'float32' | 'float64' | 'string' | 'bool' | 'datetime' | 'qnum' | 'bigint';

interface WorkerColumnStats {
  totalCount: number;
  missingValueCount: number;
  uniqueCount: number;
  valueCount: number;
  min: number;       max: number;
  sum: number;       avg: number;
  stdev: number;     variance: number;
  skew: number;      kurt: number;
  med: number;
  q1: number;        q2: number;        q3: number;
  nullValue: number; // INT_NULL (-2147483648) or FLOAT_NULL (2.6789344063684636e-34)
}

interface WorkerColumn {
  name: string;
  type: ColumnType;
  length: number;
  rawData: RawData;
  stats: WorkerColumnStats;
  categories?: string[]; // only for type === 'string'
}

interface WorkerDataFrame {
  name: string;
  rowCount: number;
  columns: WorkerColumn[];
}
```

### Transforms (`worker-transforms.ts`) — requires `datagrok-api`, main-thread only

| Function | Signature | Direction |
|----------|-----------|-----------|
| `toWorkerColumn` | `(col: DG.Column) => WorkerColumn` | DG -> Worker |
| `toWorkerColumns` | `(columns: DG.ColumnList) => WorkerColumn[]` | DG -> Worker |
| `toWorkerDataFrame` | `(df: DG.DataFrame) => WorkerDataFrame` | DG -> Worker |
| `fromWorkerColumn` | `(wc: WorkerColumn) => DG.Column` | Worker -> DG |
| `fromWorkerDataFrame` | `(wdf: WorkerDataFrame) => DG.DataFrame` | Worker -> DG |

### Null Sentinel Values

| ColumnType | Sentinel | Constant |
|------------|----------|----------|
| `int`, `string`, `bool` | -2147483648 | `INT_NULL` |
| `float32`, `float64`, `datetime`, `qnum` | 2.6789344063684636e-34 | `FLOAT_NULL` |

---

## Missing Values Strategy

Before implementing any new worker-based method, define the missing values strategy (skip, impute, propagate, or reject). See `COMPUTATION-PATTERNS.md` (Missing Values Strategy section) for the full decision table and per-column behavior guidelines.

---

## Working with WorkerColumn Inside a Worker

**IMPORTANT:** The `rawData` array's `.length` may be larger than `col.length` (due to internal buffer allocation). Always use `col.length` for iteration bounds, never `rawData.length`.

### Reading numerical data

Check `col.stats.missingValueCount` before processing. If there are no missing values, skip all null checks
in the loop — this eliminates a branch per iteration and significantly speeds up computation.

```typescript
// worker.ts
import {WorkerColumn} from './worker-defs';

onmessage = (e: MessageEvent) => {
  const col: WorkerColumn = e.data;
  const raw = col.rawData as Float32Array;
  const n = col.length;

  if (col.stats.missingValueCount > 0) {
    const nullVal = col.stats.nullValue;
    for (let i = 0; i < n; i++) {
      if (raw[i] === nullVal) continue; // skip missing
      // process raw[i]
    }
  } else {
    for (let i = 0; i < n; i++) {
      // process raw[i] — no null checks needed
    }
  }
};
```

### Centering / scaling using stats

```typescript
function centerAndScale(col: WorkerColumn): Float32Array {
  const raw = col.rawData as Float32Array;
  const result = new Float32Array(col.length);
  const nullVal = col.stats.nullValue;
  const avg = col.stats.avg;
  const stdev = col.stats.stdev;

  for (let i = 0; i < col.length; i++) {
    if (raw[i] === nullVal)
      result[i] = nullVal;
    else
      result[i] = (raw[i] - avg) / stdev;
  }
  return result;
}
```

### Building a feature matrix from WorkerColumn[]

Choose the matrix layout based on the method's primary access pattern — this directly affects
CPU cache utilization:

- **Column-major**: sequential access within each column. Optimal when columns are processed
  independently (centering, scaling, per-feature statistics, WASM interop).
- **Row-major flat**: sequential access across features of each row. Optimal for distance
  computation, KNN, nearest neighbor search.
- **Row-major typed**: same access pattern as row-major flat, but each row is a separate
  `Float32Array`. Use as a drop-in replacement for `number[][]`.

Avoid random access patterns — a cache miss per access can be up to 100x slower than sequential traversal.
For more details on data locality, see `COMPUTATION-PATTERNS.md` (Data Locality section).

```typescript
// Row-major flat Float32Array: data[i * nCols + j]
// Single allocation, contiguous memory, no boxing overhead.
function toFlatRowMajor(cols: WorkerColumn[]): Float32Array {
  const nRows = cols[0].length;
  const nCols = cols.length;
  const data = new Float32Array(nRows * nCols);
  for (let j = 0; j < nCols; j++) {
    const raw = cols[j].rawData;
    for (let i = 0; i < nRows; i++)
      data[i * nCols + j] = raw[i];
  }
  return data;
}

// Row-major Float32Array[]: data[i][j]
// Drop-in replacement for number[][] with unboxed typed rows.
function toTypedRowMajor(cols: WorkerColumn[]): Float32Array[] {
  const nRows = cols[0].length;
  const nCols = cols.length;
  const data: Float32Array[] = new Array(nRows);
  for (let i = 0; i < nRows; i++)
    data[i] = new Float32Array(nCols);
  for (let j = 0; j < nCols; j++) {
    const raw = cols[j].rawData;
    for (let i = 0; i < nRows; i++)
      data[i][j] = raw[i];
  }
  return data;
}

// Column-major flat Float32Array: data[i + j * nRows]
// Optimal when columns are processed independently (WASM, matrix ops).
function toFlatColumnMajor(cols: WorkerColumn[]): Float32Array {
  const nRows = cols[0].length;
  const nCols = cols.length;
  const data = new Float32Array(nRows * nCols);
  for (let j = 0; j < nCols; j++) {
    const raw = cols[j].rawData;
    const offset = j * nRows;
    for (let i = 0; i < nRows; i++)
      data[i + offset] = raw[i];
  }
  return data;
}
```

### Creating a result WorkerColumn

```typescript
function makeResultColumn(name: string, data: Float32Array): WorkerColumn {
  return {
    name: name,
    type: 'float32',
    length: data.length,
    rawData: data,
    stats: computeStats(data), // compute in worker or leave zeros if not needed
  };
}
```

---

## Worker Lifecycle Pattern

### Main thread (caller)

```typescript
import {toWorkerColumns, fromWorkerColumn} from './worker-utils/worker-transforms';
import {WorkerColumn} from './worker-utils/worker-defs';

async function runInWorker(
  features: DG.ColumnList, components: number
): Promise<DG.Column[]> {
  const workerFeatures = toWorkerColumns(features);

  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL('./workers/my-worker.ts', import.meta.url));

    worker.postMessage({features: workerFeatures, components});

    worker.onmessage = (e) => {
      worker.terminate();
      if (e.data.success) {
        const result = e.data.data.columns as WorkerColumn[];
        resolve(result.map(fromWorkerColumn));
      } else {
        reject(new Error(e.data.error));
      }
    };

    worker.onerror = (e) => {
      worker.terminate();
      reject(new Error(e.message));
    };
  });
}
```

### Web worker

```typescript
import {WorkerColumn} from '../worker-utils/worker-defs';

interface MyWorkerInput {
  features: WorkerColumn[];
  components: number;
}

interface MyWorkerOutput {
  success: true;
  data: {columns: WorkerColumn[]};
} | {
  success: false;
  error: string;
}

onmessage = (e: MessageEvent<MyWorkerInput>) => {
  try {
    const {features, components} = e.data;

    // Access raw data directly:
    const nRows = features[0].length;
    const nCols = features.length;

    // Use stats:
    for (const f of features) {
      const avg = f.stats.avg;
      const stdev = f.stats.stdev;
      const nullVal = f.stats.nullValue;
      // ...
    }

    // Build result columns:
    const resultCols: WorkerColumn[] = [];
    for (let c = 0; c < components; c++) {
      const data = new Float32Array(nRows);
      // ... fill data ...
      resultCols.push({
        name: `Component ${c + 1}`,
        type: 'float32',
        length: nRows,
        rawData: data,
        stats: { /* fill or leave defaults */ } as any,
      });
    }

    postMessage({success: true, data: {columns: resultCols}} satisfies MyWorkerOutput);
  } catch (err) {
    postMessage({success: false, error: String(err)});
  }
};
```

---

## Parallel Execution

For distributing independent computations across multiple workers, see `PARALLEL-EXECUTION.md`.
