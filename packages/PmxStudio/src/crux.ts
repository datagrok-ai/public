/*
 * Crux (https://github.com/datagrok-ai/crux) in the browser: the wasm modules this package ships in
 * crux-modules/ (nonclinical and stats, the ones its functions are routed to; README "Crux modules"), loaded once
 * by Crux's own module loader
 * (crux-modules.js, renamed from .mjs: a Web Worker per module, functions routed by manifest.json) and called
 * through its leaf ABI: plain columns in, the result's columns out.
 */
import {_package} from './package';

/** The part of crux-modules.js's call result this package reads. */
export interface CruxResult {
  main: {rows: number, columns: string[], col: {[name: string]: any}};
  named?: {[name: string]: any};
}

interface CruxModules {
  call(name: string, table: {[col: string]: any}, params?: object): Promise<CruxResult>;
}

/** Sets of workers, opened as calls overlap: a set's module runs one call at a time, and every call is seeded, so
 *  which set runs it changes nothing. Each set costs its modules' memory, hence the cap. */
const pool: {modules: Promise<CruxModules>, busy: number}[] = [];
export const MAX_SETS = Math.max(1, Math.min(15, (navigator.hardwareConcurrency ?? 2) - 1));

let first: Promise<CruxModules> | null = null;

/** A new set of workers: the loader and its manifest are fetched for the first only (a slow server pays once). */
function load(): Promise<CruxModules> {
  if (first)
    return first.then((m: any) => new m.constructor(m.manifest, m.options));
  const base = new URL('crux-modules/', _package.webRoot.endsWith('/') ? _package.webRoot : `${_package.webRoot}/`);
  first = import(/* webpackIgnore: true */ new URL('crux-modules.js', base).href)
    .then((m) => m.CruxModules.fromBase(base) as Promise<CruxModules>);
  first.catch(() => first = null);
  return first;
}

/** A float column of `n` nulls (a function's "solve for this" input): the leaf ABI's explicit layout, since
 *  a typed array's NaN is a value and an all-null plain array has no type to infer. */
export function nullFloats(n: number): {format: string, values: Float64Array, validity: Uint8Array} {
  return {format: 'g', values: new Float64Array(n), validity: new Uint8Array(Math.ceil(n / 8))};
}

const utf8 = new TextDecoder();

/** Row i of a decoded leaf column, as crux-leaf.js's toArray reads it: null where the validity bitmap says so,
 *  a boolean from its bit, a 64-bit integer as a number, a UTF-8 string decoded, a list of numbers as a view. */
function at(c: any, i: number): any {
  const k = (c.offset ?? 0) + i;
  if (c.validity && !((c.validity[k >> 3] >> (k & 7)) & 1))
    return null;
  if (c.format === '+l' || c.format === '+L') {
    const ch = c.children[0];
    return ch.values.subarray((ch.offset ?? 0) + Number(c.offsets[k]), (ch.offset ?? 0) + Number(c.offsets[k + 1]));
  }
  if (c.format === 'b')
    return !!((c.values[k >> 3] >> (k & 7)) & 1);
  if (c.format === 'u')
    return utf8.decode(c.values.subarray(c.offsets[k], c.offsets[k + 1]));
  const v = c.values[k];
  return typeof v === 'bigint' ? Number(v) : v;
}

/** A decoded leaf table ({rows, columns}), e.g. a result's main table, as row objects. */
export const rows = (t: {rows: number, columns: any[]}) => Array.from({length: t.rows}, (_, i) =>
  Object.fromEntries(t.columns.map((c) => [c.name, at(c, i)])) as {[k: string]: any});

/** A result's named output table as row objects. */
export const named = (r: CruxResult, name: string) => rows(r.named![name].table);

/** Calls the crux-fn `name` (e.g. 'stats.map_prior') on `table` (column name -> values) with `params`. */
export async function crux(name: string, table: {[col: string]: any}, params: object = {}): Promise<CruxResult> {
  let set = pool.find((s) => s.busy === 0);
  if (!set && pool.length < MAX_SETS) {
    const opened = {modules: load(), busy: 0};
    pool.push(opened);
    opened.modules.catch(() => pool.splice(pool.indexOf(opened), 1)); // the next call retries
    set = opened;
  }
  set ??= pool.reduce((a, b) => b.busy < a.busy ? b : a);
  set.busy++;
  try {
    return await (await set.modules).call(name, table, params);
  } finally {
    set.busy--;
  }
}
