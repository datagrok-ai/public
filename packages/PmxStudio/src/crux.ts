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
  /** `columns`: each column as decoded (its name and values); `col`: the values by name */
  main: {rows: number, columns: {name: string, [k: string]: any}[], col: {[name: string]: any}};
  named?: {[name: string]: any};
}

export interface CruxModules {
  call(name: string, table: {[col: string]: any}, params?: object): Promise<CruxResult>;
  /** the modules loaded so far, by name; `load` loads one without calling it; `terminate` stops their workers */
  loaded(): string[];
  load(module: string): Promise<void>;
  terminate(): Promise<void>;
}

export const MAX_SETS = Math.max(1, Math.min(15, (navigator.hardwareConcurrency ?? 2) - 1));
/** How long the engine stays idle before letting its sets of workers go: longer than a pause to read, so only a user
 *  who has gone (another tab, a meeting) pays the reopening, and a gesture that wakes it (`cruxWake`) hides that. */
export const IDLE_MS = 180000;

/** Sets of workers, opened as calls overlap: a set's module runs one call at a time, and every call is seeded, so
 *  which set runs it changes nothing. A module's wasm memory never shrinks (it keeps what its largest call took),
 *  so after `idleMs` with no call every set is stopped and one new one is opened, its modules loaded again: the
 *  next call starts warm, the memory of the calls before is given back. `wake` (a gesture: someone is there) reopens,
 *  at rest, as many sets as there were, so the calls the gesture is about to make find them open; open, it starts the
 *  countdown again. */
export function cruxPool(open: () => Promise<CruxModules>, maxSets: number, idleMs: number) {
  type Set = {modules: Promise<CruxModules>, busy: number};
  const pool: Set[] = [];
  let timer: ReturnType<typeof setTimeout> | undefined;
  let generation = 0; // a call since the pause began: the pause is not over
  let rested: {sets: number, modules: string[]} | null = null; // what the last pause closed, until a call or a wake
  const countdown = () => {
    clearTimeout(timer);
    timer = setTimeout(() => void rest(), idleMs);
  };
  const fresh = (modules: string[]): Set => {
    const set: Set = {busy: 0, modules: open().then(async (m) => {
      await Promise.all(modules.map((x) => m.load(x).catch(() => {}))); // (a module that fails loads on its call)
      return m;
    })};
    set.modules.catch(() => { // a set that fails to open leaves the pool: the next call opens another
      const i = pool.indexOf(set);
      if (i >= 0)
        pool.splice(i, 1);
    });
    return set;
  };
  const rest = async () => {
    const at = generation;
    if (!pool.length || pool.some((s) => s.busy > 0)) // (none open, as after a set that failed: nothing to keep)
      return;
    // the modules the sets had loaded, read before the swap: the new set loads them, so the next call finds them
    const modules = (await Promise.all(pool.map((s) => s.modules.then((m) => m.loaded(), () => []))));
    if (generation !== at || pool.some((s) => s.busy > 0)) // (a call came in meanwhile: the next pause)
      return;
    rested = {sets: pool.length, modules: [...new Set(modules.flat())]};
    const old = pool.splice(0, pool.length, fresh(rested.modules));
    old.forEach((s) => s.modules.then((m) => m.terminate(), () => {}));
  };
  return {
    sets: () => pool.length, // (for the package's tests)
    wake(): void {
      // (never opened, or calculating: the calculation's end counts down; at rest with nothing to reopen, a countdown
      // would only swap the one set for the same)
      if (!pool.length || pool.some((s) => s.busy > 0) || (rested && pool.length >= rested.sets))
        return;
      generation++;
      while (rested && pool.length < rested.sets)
        pool.push(fresh(rested.modules));
      rested = null;
      countdown();
    },
    async call(name: string, table: {[col: string]: any}, params: object = {}): Promise<CruxResult> {
      clearTimeout(timer);
      generation++;
      rested = null;
      let set = pool.find((s) => s.busy === 0);
      if (!set && pool.length < maxSets) {
        set = fresh([]);
        pool.push(set);
      }
      set ??= pool.reduce((a, b) => b.busy < a.busy ? b : a);
      set.busy++;
      try {
        return await (await set.modules).call(name, table, params);
      } finally {
        set.busy--;
        if (pool.every((s) => s.busy === 0))
          countdown();
      }
    },
  };
}

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

const engine = cruxPool(load, MAX_SETS, IDLE_MS);
let calls = 0;
/** The calls made so far (for the package's tests: what opens the engine). */
export const cruxCalls = () => calls;

/** A failure of the engine itself, its text as Crux gives it: technical, so not shown to the user (ui.ts's
 *  plainError). */
export class EngineError extends Error {}

/** Not a failure: what is still missing before a stage can show anything (shown plainly, not as an error). */
export class EmptyState extends Error {}

/** Calls the crux-fn `name` (e.g. 'stats.map_prior') on `table` (column name -> values) with `params`. */
export const crux = (name: string, table: {[col: string]: any}, params: object = {}): Promise<CruxResult> => {
  calls++;
  return engine.call(name, table, params).catch((e) => {
    throw new EngineError(String(e?.message ?? e));
  });
};

/** A gesture that will compute (a field changed, a stage opened, the tab back): the engine at rest reopens its sets
 *  now, not on the first call. Nothing when it is open, or never was. */
export const cruxWake = () => engine.wake();
