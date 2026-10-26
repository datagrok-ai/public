/**
 * The lazy loader of Chem's one Crux file, `vendor/chem-wasm` (RDKit.js-shaped calls, the search Collection,
 * descriptors, the scaffold tree): fetched and compiled once, on the first use; instantiated on the main thread on the
 * first use there (rendering, the search's query check); the compiled Module is what the Crux workers instantiate.
 * In RDKit mode only what runs on Crux in both modes loads it (molzip; FCFP and chirality fingerprints). A load not
 * done within a minute of its fetch's start (a stalled fetch) has failed; a failed load is tried once more a minute
 * later (a network hiccup), then it is final for the session; a wasm trap on the main thread
 * (`WebAssembly.RuntimeError`) is final at once. The callers take their RDKit paths meanwhile.
 */
import * as DG from 'datagrok-api/dg';
import {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import {getRdKitWebRoot} from '../utils/chem-common-rdkit';

/** The vendored wasm, emitted under its own name into `dist/` (rspack.config.js `emit`). */
export const CRUX_WASM_FILE = 'chem_wasm_bg.wasm';
/** A first failed load is tried again once, this long after. */
const RETRY_AFTER_MS = 60000;
/** The fetch and compile of the file (and the glue's instance) within this long of the fetch's start: a 2 MB file on a
 * slow link takes some 20 s; past it the load has failed. */
const LOAD_LIMIT_MS = 60000;

/** The glue's module namespace: RDKit.js's module functions (`get_mol`, `get_qmol`, ...) plus Crux's own exports. */
export type ChemWasm = typeof import('../../vendor/chem-wasm/chem_wasm.js');

let compiled: Promise<WebAssembly.Module> | null = null;
let instantiated: Promise<ChemWasm> | null = null;
let ready: ChemWasm | null = null;
let memory: WebAssembly.Memory | null = null;
let failure: Error | null = null;
/** The failure is for the session: the second failed load, or a trap. */
let failureFinal = false;
let failedAt = 0;
let retried = false;
/** Tests: while set, the module counts as still loading (see {@link holdCruxForTests}). */
let held: Promise<void> | null = null;
/** Tests: the file's URL in place of the real one (see {@link resetCruxLoaderForTests}). */
let urlOverride: string | null = null;
/** Tests: the fetch in place of the real one, and the load's limit (see {@link setCruxFetchForTests}). */
let fetchOverride: ((url: string, init?: RequestInit) => Promise<Response>) | null = null;
let loadLimitMs = LOAD_LIMIT_MS;
/** Who asked for the module, in order ('engine', 'search', 'descriptors'): tests read that RDKit mode asks nothing. */
const requests: string[] = [];

/** Chem's web root as RDKit's loader has it, else (the package-test bundle before any init) the published Chem's. */
function chemWebRoot(): string {
  return getRdKitWebRoot() ?? DG.Func.find({package: 'Chem', name: 'getRdKitModule'})[0]?.package?.webRoot ?? '';
}

export function cruxWasmUrl(): string {
  return urlOverride ?? `${chemWebRoot()}/dist/${CRUX_WASM_FILE}`;
}

/** The failure in effect: a first load failure is cleared, once, a minute after it happened. */
function currentFailure(): Error | null {
  if (failure && !failureFinal && !retried && performance.now() - failedAt >= RETRY_AFTER_MS) {
    retried = true;
    failure = null;
    compiled = null;
    instantiated = null;
  }
  return failure;
}

/** The file, fetched and compiled once (streaming where the server sends application/wasm). */
export function compileCrux(): Promise<WebAssembly.Module> {
  const failed = currentFailure();
  if (failed)
    return Promise.reject(failed);
  const abort = new AbortController(); // a load past its limit stops its download: a retry starts the only one
  compiled ??= withinLimit((async () => {
    const url = cruxWasmUrl();
    const get = (u: string) => fetchOverride ? fetchOverride(u, {signal: abort.signal}) :
      fetch(u, {signal: abort.signal});
    return WebAssembly.compileStreaming(get(url)).catch(async () => {
      const response = await get(url);
      if (!response.ok)
        throw new Error(`Crux wasm is not available: ${response.status} ${response.statusText}`);
      return WebAssembly.compile(await response.arrayBuffer());
    });
  })(), () => abort.abort()).catch((e) => fail(e));
  return compiled;
}

/** The main thread's instance, created once from the compiled file, set up as Chem sets up RDKit. */
export function loadCrux(reason: string): Promise<ChemWasm> {
  requests.push(reason);
  const failed = currentFailure();
  if (failed)
    return Promise.reject(failed);
  if (held)
    return held.then(() => loadCrux(reason));
  instantiated ??= withinLimit((async () => {
    const module = await compileCrux();
    const glue: ChemWasm = await import('../../vendor/chem-wasm/chem_wasm.js');
    memory = (await glue.default({module_or_path: module})).memory;
    glue.prefer_coordgen(false);
    glue.use_legacy_stereo_perception(false);
    ready = glue;
    console.log('Chem | Crux module initialized');
    return glue;
  })()).catch((e) => fail(e));
  return instantiated;
}

/** `work`, or a failure once the load's limit has passed since it started (a fetch that never settles); `onTimeout`
 * then stops what it started. */
function withinLimit<T>(work: Promise<T>, onTimeout?: () => void): Promise<T> {
  const limitMs = loadLimitMs;
  let timer: ReturnType<typeof setTimeout> | undefined;
  const limit = new Promise<never>((_, reject) => timer = setTimeout(() => {
    onTimeout?.();
    reject(new Error(`the Crux wasm did not load within ${limitMs / 1000} s`));
  }, limitMs));
  return Promise.race([work, limit]).finally(() => clearTimeout(timer));
}

/** The main thread's instance as an RDKit.js module, or null until {@link loadCrux} has finished. */
export function cruxModuleSync(): RDModule | null {
  if (held || !ready)
    return null;
  if (!parseTrap)
    return ready as unknown as RDModule;
  const trap = () => {
    throw new WebAssembly.RuntimeError('simulated: unreachable');
  };
  return {...ready, get_mol: trap, get_qmol: trap} as unknown as RDModule;
}

/** Tests: the main thread's Crux module traps when it parses a molecule (get_mol, get_qmol), until called with
 * false. */
export function setCruxParseTrapForTests(on: boolean): void {
  parseTrap = on;
}
let parseTrap = false;

/** The main thread's Crux heap in bytes (0 before the module is loaded): tests watch it for growth. */
export function cruxHeapBytes(): number {
  return memory?.buffer.byteLength ?? 0;
}

export function cruxLoadFailed(): boolean {
  return currentFailure() !== null;
}

export function cruxLoadRequests(): readonly string[] {
  return requests;
}

function fail(e: any): never {
  if (!failure) {
    failure = e instanceof Error ? e : new Error(String(e));
    failedAt = performance.now();
    failureFinal = retried;
    console.warn(`Chem | the Crux module failed to load, RDKit is used ${failureFinal ? 'for this session' :
      'until it is tried again in a minute'}: ${failure.message}`);
  }
  throw failure;
}

/**
 * A main-thread Crux call threw `e`: a wasm trap (`WebAssembly.RuntimeError`: a Rust panic or an out-of-bounds access)
 * may have left the instance broken, so Crux is off for the session and RDKit does its work; logged once.
 */
export function noteCruxError(e: unknown): void {
  if (!(e instanceof WebAssembly.RuntimeError) || failureFinal)
    return;
  failure = e;
  failureFinal = true;
  console.warn(`Chem | the Crux module trapped, RDKit is used for this session: ${e.message}`);
}

/**
 * Tests: makes the loader fail as a missing file would (`error`), or clears a failure set so (null). A real failure
 * stays: the session has no Crux then.
 */
export function setCruxLoadFailureForTests(error: Error | null): void {
  if (error) {
    failure = error;
    simulated.add(error);
  } else if (failure && simulated.has(failure))
    failure = null;
}
const simulated = new WeakSet<Error>();

/**
 * Tests: forgets the loaded module and any failure, and loads from `url` next (null: the real file). A load from a
 * missing file then fails as it would in a session; the real file loads again after the next reset.
 */
export function resetCruxLoaderForTests(url: string | null): void {
  urlOverride = url;
  compiled = null;
  instantiated = null;
  ready = null;
  failure = null;
  failureFinal = false;
  retried = false;
}

/**
 * Tests: the loader fetches with `fetchFile` (a stalled fetch: a promise that never settles) and gives up after
 * `limitMs` (default: the real limit); null restores the real fetch and limit. It applies to the next load: reset the
 * loader first ({@link resetCruxLoaderForTests}).
 */
export function setCruxFetchForTests(fetchFile: ((url: string) => Promise<Response>) | null, limitMs?: number):
  void {
  fetchOverride = fetchFile;
  loadLimitMs = fetchFile ? limitMs ?? LOAD_LIMIT_MS : LOAD_LIMIT_MS;
}

/**
 * Tests: the module counts as still loading (a switch to Crux then draws with RDKit) until the returned function is
 * called; the load then completes as a first load does.
 */
export function holdCruxForTests(): () => void {
  let release!: () => void;
  const hold = new Promise<void>((resolve) => release = resolve);
  held = hold;
  return () => {
    if (held === hold)
      held = null;
    release();
  };
}
