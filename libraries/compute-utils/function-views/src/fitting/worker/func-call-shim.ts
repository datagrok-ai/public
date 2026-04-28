// In-worker FuncCall shim.
//
// Compiled JS bodies of fittable scripts read `funcCall.inputs.X` and assign
// to `funcCall.outputs.Y` (or to bare locals declared in the //output block —
// `grok api` rewrites those to `funcCall.outputs.Y` at registration time, but
// here we hand the body the *names* directly via `new Function(...names, body)`
// and recover assignments via the returned context). The shim covers only the
// surface cost-functions.ts:42,45,75,79 actually reads.
//
// Compile cache is module-level so it hits across NM iterations within a task
// AND across tasks within a long-lived worker — pays the body's parse cost
// exactly once per source. Bounded LRU so accidental script churn over a
// tab-lifetime worker can't grow memory without bound. Active sessions hold
// strong references to their compiled function via `createWorkerFuncCall`'s
// closure, so cache eviction never invalidates an in-flight fit.

import {LRUCache} from 'lru-cache';
import {createWorkerDG, WorkerDG} from './dg-shim';

// Mirror of serialize.MAX_FN_SOURCE_BYTES — duplicated here so the worker
// entry doesn't pull serialize.ts (which transitively imports
// @datagrok-libraries/arrow → datagrok-api/dg) into the worker bundle.
// 1 MB covers Diff Studio bodies and multi-equation models; raise again
// only with measured workloads.
const MAX_FN_SOURCE_BYTES = 1024 * 1024;

function checkSourceSize(source: string): void {
  if (source.length > MAX_FN_SOURCE_BYTES)
    throw new Error(`function source exceeds ${MAX_FN_SOURCE_BYTES} bytes (got ${source.length})`);
}

type CompiledFn = (...args: any[]) => any;

const DEFAULT_COMPILE_CACHE_CAP = 32;
let compileCache = new LRUCache<string, CompiledFn>({max: DEFAULT_COMPILE_CACHE_CAP});

let compileHits = 0;
let compileMisses = 0;

export function clearCompileCache(): void {
  compileCache.clear();
  compileHits = 0;
  compileMisses = 0;
}

// Test-only: read hit/miss counters to assert cache behavior.
export function _getCompileStats(): {hits: number; misses: number} {
  return {hits: compileHits, misses: compileMisses};
}

// Test-only: rebuild cache with a smaller cap so eviction can be observed
// without compiling 33+ bodies. Resets counters too.
export function _setCompileCacheCap(cap: number): void {
  compileCache = new LRUCache<string, CompiledFn>({max: cap});
  compileHits = 0;
  compileMisses = 0;
}

// Compile a JS body that reads named inputs and writes named outputs. The body
// runs as the script body of a function with signature `(_DG, ...inputs)` and
// returns an object whose keys are the //output param names.
//
// We rewrite the body to:
//   1. Hoist `let <output> = undefined;` per output param so the body's bare
//      assignments (`y = ...`, `simulation = DG.DataFrame...`) succeed without
//      ReferenceError.
//   2. Append `return { <outputs> };` to capture them.
//
// `DG` is bound to the worker shim (`createWorkerDG()`) — fittable bodies
// touch DG.DataFrame.fromColumns and DG.Column.from*Array, all covered.
export function compileBody(source: string, paramList: string[], outputNames: string[]): CompiledFn {
  checkSourceSize(source);
  const key = `${paramList.join(',')}|${outputNames.join(',')}|${source}`;
  let fn = compileCache.get(key);
  if (fn) {
    ++compileHits;
    return fn;
  }
  ++compileMisses;
  const declOutputs = outputNames.map((n) => `let ${n} = undefined;`).join('');
  const returnLine = `return { ${outputNames.join(', ')} };`;
  const wrapped = `${declOutputs}\n${source}\n${returnLine}`;
  fn = new Function('DG', ...paramList, wrapped) as CompiledFn;
  compileCache.set(key, fn);
  return fn;
}

export interface WorkerFuncCall {
  inputs: Record<string, any>;
  outputs: Record<string, any>;
  getParamValue(name: string): any;
  call(varied: Record<string, number>): void;
}

export function createWorkerFuncCall(args: {
  source: string;
  paramList: string[]; // ordered list: all inputs (fixed + varied) by declaration order
  outputNames: string[];
  fixedInputs: Record<string, any>;
  variedInputNames: string[];
  dg?: WorkerDG;
}): WorkerFuncCall {
  const dg = args.dg ?? createWorkerDG();
  const compiled = compileBody(args.source, args.paramList, args.outputNames);
  const inputs: Record<string, any> = {...args.fixedInputs};
  let outputs: Record<string, any> = {};

  return {
    get inputs() {return inputs;},
    get outputs() {return outputs;},
    getParamValue(name: string): any {
      if (name in outputs) return outputs[name];
      if (name in inputs) return inputs[name];
      return undefined;
    },
    call(varied: Record<string, number>): void {
      for (const n of args.variedInputNames)
        inputs[n] = varied[n];
      const orderedArgs = args.paramList.map((n) => inputs[n]);
      outputs = compiled(dg, ...orderedArgs) ?? {};
    },
  };
}
