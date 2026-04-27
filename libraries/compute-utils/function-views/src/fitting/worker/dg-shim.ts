/**
 * Worker-side `DG` shim.
 *
 * Compiled JS bodies of fittable scripts construct outputs by calling
 * `DG.DataFrame.fromColumns([DG.Column.from*Array(...)])`. Inside a web
 * worker there is no Dart-bridge `DG`; this module supplies the minimum
 * factory surface those bodies use, returning LiteDataFrame / LiteColumn
 * objects backed directly by typed arrays.
 *
 * Pure TypeScript, no `datagrok-api` imports — safe to use inside a worker.
 */

import {LiteColumn, LiteColumnList, LiteColumnStats, LiteColumnType, LiteDataFrame, LiteRawData}
  from './types';

// Mirror DG.FLOAT_NULL / DG.INT_NULL from js-api/src/const.ts. Hard-coded so
// this module has no DG dependency.
export const FLOAT_NULL = 2.6789344063684636e-34;
export const INT_NULL = -2147483648;

const COLUMN_TYPE = Object.freeze({
  INT: 'int',
  FLOAT: 'float',
  DATE_TIME: 'datetime',
  STRING: 'string',
  BOOL: 'bool',
  BIG_INT: 'bigint',
} as const);

function makeNumericColumn(name: string, type: LiteColumnType, raw: Float64Array | Float32Array | Int32Array,
  nullSentinel: number): LiteColumn {
  let cached: LiteColumnStats | null = null;
  return {
    name,
    type,
    length: raw.length,
    getRawData: () => raw,
    toList: () => Array.from(raw as ArrayLike<number>),
    get: (i) => {
      const v = raw[i];
      return v === nullSentinel ? null : v;
    },
    get stats() {
      if (cached) return cached;
      let min = Number.POSITIVE_INFINITY;
      let max = Number.NEGATIVE_INFINITY;
      for (let i = 0; i < raw.length; ++i) {
        const v = raw[i];
        if (v === nullSentinel) continue;
        if (v < min) min = v;
        if (v > max) max = v;
      }
      // Match DG.Column.stats behaviour for empty/all-null columns: min/max
      // collapse to 0. Avoids ±Infinity in scaling logic downstream.
      if (min === Number.POSITIVE_INFINITY) min = 0;
      if (max === Number.NEGATIVE_INFINITY) max = 0;
      cached = {min, max};
      return cached;
    },
  };
}

function makeBigIntColumn(name: string, raw: BigInt64Array): LiteColumn {
  let cached: LiteColumnStats | null = null;
  return {
    name,
    type: COLUMN_TYPE.BIG_INT,
    length: raw.length,
    getRawData: () => raw,
    toList: () => Array.from(raw),
    get: (i) => raw[i],
    get stats() {
      if (cached) return cached;
      // Stats are reported as numbers; lossy for values > Number.MAX_SAFE_INTEGER
      // but matches DG.Column.stats which also returns numeric min/max.
      let min = Number.POSITIVE_INFINITY;
      let max = Number.NEGATIVE_INFINITY;
      for (let i = 0; i < raw.length; ++i) {
        const v = Number(raw[i]);
        if (v < min) min = v;
        if (v > max) max = v;
      }
      if (min === Number.POSITIVE_INFINITY) min = 0;
      if (max === Number.NEGATIVE_INFINITY) max = 0;
      cached = {min, max};
      return cached;
    },
  };
}

function makeStringColumn(name: string, values: string[]): LiteColumn {
  // Stats are not meaningful for strings; expose 0/0 to match DG behaviour.
  const cached: LiteColumnStats = {min: 0, max: 0};
  return {
    name,
    type: COLUMN_TYPE.STRING,
    length: values.length,
    getRawData: () => values,
    toList: () => values.slice(),
    get: (i) => values[i],
    get stats() {
      return cached;
    },
  };
}

function makeBoolColumn(name: string, values: boolean[]): LiteColumn {
  const cached: LiteColumnStats = {min: 0, max: 0};
  return {
    name,
    type: COLUMN_TYPE.BOOL,
    length: values.length,
    getRawData: () => values,
    toList: () => values.slice(),
    get: (i) => values[i],
    get stats() {
      return cached;
    },
  };
}

function makeDataFrame(cols: LiteColumn[]): LiteDataFrame {
  const byName = new Map<string, LiteColumn>();
  for (const c of cols) byName.set(c.name, c);
  const rowCount = cols.length === 0 ? 0 : cols[0].length;
  const columns: LiteColumnList = {
    length: cols.length,
    byName: (n) => {
      const c = byName.get(n);
      if (!c) throw new Error(`column not found: ${n}`);
      return c;
    },
    names: () => cols.map((c) => c.name),
    [Symbol.iterator]: function* () {
      yield* cols;
    },
  };
  return {
    columns,
    rowCount,
    col: (n) => byName.get(n) ?? null,
  };
}

export interface WorkerDG {
  DataFrame: {
    fromColumns(cols: LiteColumn[]): LiteDataFrame;
  };
  Column: {
    fromFloat64Array(name: string, arr: Float64Array): LiteColumn;
    fromFloat32Array(name: string, arr: Float32Array): LiteColumn;
    fromInt32Array(name: string, arr: Int32Array): LiteColumn;
    fromBigInt64Array(name: string, arr: BigInt64Array): LiteColumn;
    fromList(type: LiteColumnType, name: string, values: unknown[]): LiteColumn;
    fromStrings(name: string, values: string[]): LiteColumn;
  };
  COLUMN_TYPE: typeof COLUMN_TYPE;
  FLOAT_NULL: number;
  INT_NULL: number;
}

export function createWorkerDG(): WorkerDG {
  return {
    DataFrame: {
      fromColumns: makeDataFrame,
    },
    Column: {
      // DG.Column.fromFloat64Array reports `.type === 'double'` — match that.
      // fromFloat32Array reports `.type === 'float'`.
      fromFloat64Array: (name, arr) => makeNumericColumn(name, 'double', arr, FLOAT_NULL),
      fromFloat32Array: (name, arr) => makeNumericColumn(name, 'float', arr, FLOAT_NULL),
      fromInt32Array: (name, arr) => makeNumericColumn(name, COLUMN_TYPE.INT, arr, INT_NULL),
      fromBigInt64Array: (name, arr) => makeBigIntColumn(name, arr),
      fromStrings: (name, values) => makeStringColumn(name, values),
      fromList: (type, name, values) => {
        switch (type) {
        case COLUMN_TYPE.FLOAT: {
          const buf = new Float64Array(values.length);
          for (let i = 0; i < values.length; ++i)
            buf[i] = values[i] == null ? FLOAT_NULL : Number(values[i]);
          return makeNumericColumn(name, COLUMN_TYPE.FLOAT, buf, FLOAT_NULL);
        }
        case COLUMN_TYPE.INT: {
          const buf = new Int32Array(values.length);
          for (let i = 0; i < values.length; ++i)
            buf[i] = values[i] == null ? INT_NULL : Number(values[i]) | 0;
          return makeNumericColumn(name, COLUMN_TYPE.INT, buf, INT_NULL);
        }
        case COLUMN_TYPE.STRING:
          return makeStringColumn(name, values.map((v) => v == null ? '' : String(v)));
        case COLUMN_TYPE.BOOL:
          return makeBoolColumn(name, values.map((v) => Boolean(v)));
        case COLUMN_TYPE.DATE_TIME: {
          const buf = new Float64Array(values.length);
          for (let i = 0; i < values.length; ++i) {
            const v: any = values[i];
            if (v == null) {
              buf[i] = FLOAT_NULL;
              continue;
            }
            // Accept Date | dayjs | number (epoch ms). DG stores datetime as
            // microseconds since epoch in raw data.
            const ms = typeof v === 'number' ? v :
              (typeof v.valueOf === 'function' ? Number(v.valueOf()) : Number(v));
            buf[i] = ms * 1000;
          }
          return makeNumericColumn(name, COLUMN_TYPE.DATE_TIME, buf, FLOAT_NULL);
        }
        case COLUMN_TYPE.BIG_INT: {
          const buf = new BigInt64Array(values.length);
          for (let i = 0; i < values.length; ++i)
            buf[i] = typeof values[i] === 'bigint' ? values[i] as bigint : BigInt(Number(values[i] ?? 0));
          return makeBigIntColumn(name, buf);
        }
        default:
          throw new Error(`unsupported column type for fromList: ${type}`);
        }
      },
    },
    COLUMN_TYPE,
    FLOAT_NULL,
    INT_NULL,
  };
}

// Re-export the raw-data type so callers can narrow without re-importing.
export type {LiteRawData} from './types';
