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

## Working with WorkerColumn Inside a Worker

### Reading numerical data

```typescript
// worker.ts
import {WorkerColumn} from './worker-defs';

onmessage = (e: MessageEvent) => {
  const col: WorkerColumn = e.data;
  const raw = col.rawData as Float32Array;
  const nullVal = col.stats.nullValue;
  const n = col.length;

  for (let i = 0; i < n; i++) {
    if (raw[i] === nullVal) continue; // skip missing
    // process raw[i]
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

For distributing independent computations across multiple workers.

### Worker count

```typescript
import {MIN_WORKERS_COUNT, WORKERS_COUNT_DOWNSHIFT} from './worker-utils/worker-defs';

const workerCount = Math.max(MIN_WORKERS_COUNT, navigator.hardwareConcurrency - WORKERS_COUNT_DOWNSHIFT);
```

### Fan-out / fan-in pattern

```typescript
async function runParallel<TInput, TOutput>(
  inputs: TInput[],
  workerUrl: URL,
): Promise<TOutput[]> {
  const nWorkers = Math.min(
    Math.max(MIN_WORKERS_COUNT, navigator.hardwareConcurrency - WORKERS_COUNT_DOWNSHIFT),
    inputs.length,
  );

  // Distribute inputs round-robin
  const chunks: TInput[][] = Array.from({length: nWorkers}, () => []);
  for (let i = 0; i < inputs.length; i++)
    chunks[i % nWorkers].push(inputs[i]);

  const promises = chunks.map((chunk) =>
    new Promise<TOutput[]>((resolve, reject) => {
      const worker = new Worker(workerUrl);
      worker.postMessage(chunk);
      worker.onmessage = (e) => {
        worker.terminate();
        if (e.data.success)
          resolve(e.data.data);
        else
          reject(new Error(e.data.error));
      };
      worker.onerror = (e) => {
        worker.terminate();
        reject(new Error(e.message));
      };
    }),
  );

  const results = await Promise.all(promises);
  return results.flat();
}
```

---
