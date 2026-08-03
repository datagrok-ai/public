// In-worker Datagrok JS-script body runner.
//
// Compiles JS bodies that read named inputs and write named outputs (the
// `//input:` / `//output:` script header convention) into a `(DG, ...inputs)`
// function whose return is `{ <output>: ... }`. The body is handed the names
// directly via `new Function(...names, body)` and outputs are recovered via a
// typeof-guarded `return { … }` trailer that mirrors the platform's own JS
// script wrap (see grok_shared.dart.js:121177 — `runScript(...) { <body>;
// return [outputs] }`).
//
// Compile cache is module-level so it hits across repeated invocations within
// a long-lived worker — pays the body's parse cost exactly once per source.
// Bounded LRU so accidental script churn over a tab-lifetime worker can't
// grow memory without bound. Live consumers (e.g. an in-flight fit holding a
// `WorkerFuncCall`) keep their compiled function alive via closure, so cache
// eviction never invalidates work in progress.

import {LRUCache} from 'lru-cache';
import {createWorkerDG, WorkerDG} from '../dg-lite/dg-shim';

// 1 MB covers Diff Studio bodies and multi-equation models; raise only with
// measured workloads.
const MAX_FN_SOURCE_BYTES = 1024 * 1024;

function checkSourceSize(source: string): void {
  if (source.length > MAX_FN_SOURCE_BYTES)
    throw new Error(`function source exceeds ${MAX_FN_SOURCE_BYTES} bytes (got ${source.length})`);
}

type CompiledFn = ((...args: any[]) => any) & {__isAsync?: boolean};

// AsyncFunction isn't a global; pull the constructor off an async lambda
// parsed at runtime so TS down-levelling (ES6 target → __awaiter helpers)
// can't substitute a regular Function for the async one.
const AsyncFunction = new Function(
  'return Object.getPrototypeOf(async function(){}).constructor;')() as
  new (...args: string[]) => CompiledFn;

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
// runs as the script body of a function with signature `(DG, ...inputs)` and
// returns an object whose keys are the //output param names.
//
// We append `return { <out>: typeof <out> !== 'undefined' ? <out> : undefined,
// ... };` and rely on JS engine scope resolution at the return site:
//   - `const X = ...` / `let X = ...` in the body resolve as same-scope
//     lexical bindings.
//   - bare `X = ...` becomes an implicit global (non-strict `new Function`),
//     resolved via global lookup.
//   - outputs the body never sets resolve as `undefined` via the typeof guard,
//     instead of throwing ReferenceError.
// No predeclaration → no `let X` vs `const X` SyntaxError when a body declares
// an output as `const` (e.g. `const tempDiff` in object-cooling.js).
//
// `DG` is bound to the worker shim (`createWorkerDG()`) — fittable bodies
// touch DG.DataFrame.fromColumns and DG.Column.from*Array, all covered.
//
// Cache key is a 128-bit FNV-1a fingerprint over (paramList, outputNames,
// source) — four parallel 32-bit passes with different seeds, length-prefixed
// per part so ['foo','bar'] and ['fooba','r'] hash differently. Avoids
// per-call ~1 MB ConsString allocation + flatten + memcmp inside lru-cache's
// underlying Map; collision probability for cap-32 cache is ~5×10⁻³⁵ per
// pair — operationally zero.
// Try `new Function` first; on SyntaxError (top-level `await` in the body
// is the only common case here) fall back to `new AsyncFunction`. The JS
// parser is the source of truth — no source-scanning heuristic that could
// false-positive on `// await` comments or `"await"` string literals. The
// `__isAsync` flag is stamped on the compiled function so cache hits also
// recover it.
export function compileBody(source: string, paramList: string[], outputNames: string[]): CompiledFn {
  checkSourceSize(source);
  const key = fnv1a128([paramList.join(','), outputNames.join(','), source]);
  let fn = compileCache.get(key);
  if (fn) {
    ++compileHits;
    return fn;
  }
  ++compileMisses;
  const captures = outputNames
    .map((n) => `${n}: typeof ${n} !== 'undefined' ? ${n} : undefined`)
    .join(', ');
  const wrapped = `${source}\nreturn { ${captures} };`;
  try {
    fn = new Function('DG', ...paramList, wrapped) as CompiledFn;
    fn.__isAsync = false;
  } catch (e) {
    if (!(e instanceof SyntaxError)) throw e;
    fn = new AsyncFunction('DG', ...paramList, wrapped);
    fn.__isAsync = true;
  }
  compileCache.set(key, fn);
  return fn;
}

// Four parallel 32-bit FNV-1a passes with different seeds, each fed every
// part length-prefixed (4 little-endian bytes) so adjacent parts can't be
// re-sliced into the same byte stream. Packed as four base36 digits joined
// by `-` (~26 chars total). Math.imul is the JIT-friendly 32-bit signed
// multiply primitive; `>>> 0` re-interprets the result as unsigned.
function fnv1a128(parts: string[]): string {
  const PRIME = 0x01000193;
  let a = 0x811c9dc5, b = 0xdeadbeef, c = 0xcafebabe, d = 0x12345678;
  for (const s of parts) {
    const len = s.length;
    for (let lb = 0; lb < 4; ++lb) {
      const ch = (len >>> (lb * 8)) & 0xff;
      a = Math.imul(a ^ ch, PRIME);
      b = Math.imul(b ^ ch, PRIME);
      c = Math.imul(c ^ ch, PRIME);
      d = Math.imul(d ^ ch, PRIME);
    }
    for (let i = 0; i < len; ++i) {
      const ch = s.charCodeAt(i);
      a = Math.imul(a ^ ch, PRIME);
      b = Math.imul(b ^ ch, PRIME);
      c = Math.imul(c ^ ch, PRIME);
      d = Math.imul(d ^ ch, PRIME);
    }
  }
  return `${(a >>> 0).toString(36)}-${(b >>> 0).toString(36)}-` +
    `${(c >>> 0).toString(36)}-${(d >>> 0).toString(36)}`;
}

// Split sync/async into a discriminated union: callers branch once on
// `kind` and then run a monomorphic `call`, which keeps the sync hot path
// byte-identical to the pre-async-support implementation and lets the JIT
// specialise each variant.
interface WorkerFuncCallBase {
  inputs: Record<string, any>;
  outputs: Record<string, any>;
  getParamValue(name: string): any;
}

export interface SyncWorkerFuncCall extends WorkerFuncCallBase {
  kind: 'sync';
  call(varied: Record<string, number>): void;
}

export interface AsyncWorkerFuncCall extends WorkerFuncCallBase {
  kind: 'async';
  call(varied: Record<string, number>): Promise<void>;
}

export type WorkerFuncCall = SyncWorkerFuncCall | AsyncWorkerFuncCall;

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

  const base: WorkerFuncCallBase = {
    get inputs() {return inputs;},
    get outputs() {return outputs;},
    getParamValue(name: string): any {
      if (name in outputs) return outputs[name];
      if (name in inputs) return inputs[name];
      return undefined;
    },
  };

  if (compiled.__isAsync) {
    return {
      ...base,
      kind: 'async',
      async call(varied: Record<string, number>): Promise<void> {
        for (const n of args.variedInputNames)
          inputs[n] = varied[n];
        const orderedArgs = args.paramList.map((n) => inputs[n]);
        outputs = (await compiled(dg, ...orderedArgs)) ?? {};
      },
    };
  }

  return {
    ...base,
    kind: 'sync',
    call(varied: Record<string, number>): void {
      for (const n of args.variedInputNames)
        inputs[n] = varied[n];
      const orderedArgs = args.paramList.map((n) => inputs[n]);
      outputs = compiled(dg, ...orderedArgs) ?? {};
    },
  };
}
