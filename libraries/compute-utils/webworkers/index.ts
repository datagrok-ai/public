// Generic, worker-safe building blocks shared across compute-utils workers.
//
//   - dg-lite: a minimal DG.DataFrame / DG.Column shim (LiteDataFrame /
//     LiteColumn) and an Arrow IPC → LiteDataFrame deserializer. Lets worker
//     code construct and consume DG-shaped values without `datagrok-api`.
//   - script-runner: compiler for Datagrok JS-script bodies (`//input:` /
//     `//output:` header convention) that wraps a body into a `(DG, ...inputs)`
//     function whose return is `{ <output>: ... }`.

export type {
  LiteColumn, LiteColumnList, LiteColumnStats, LiteColumnType, LiteDataFrame, LiteRawData,
} from './dg-lite/types';
export {createWorkerDG, FLOAT_NULL, INT_NULL} from './dg-lite/dg-shim';
export type {WorkerDG} from './dg-lite/dg-shim';
export {arrowIpcToLite} from './dg-lite/arrow-to-lite';

export {
  compileBody, clearCompileCache, _getCompileStats, _setCompileCacheCap,
  createWorkerFuncCall,
} from './script-runner/func-call-shim';
export type {WorkerFuncCall} from './script-runner/func-call-shim';
