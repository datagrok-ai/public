// Every function Chem has moved from RDKit.js (or OpenChemLib) to Crux, measured on a large dataset against the engine
// Chem uses for it in RDKit mode, with Chem's own code (its TypeScript transpiled on the fly) on each engine:
//
//   node compare100k.mjs [rows=100000] [only=<name,...>] [--compare-only]
//
// (--compare-only: the comparison and the report again from the outputs in OUT/raw, nothing run)
// names: parse, drawing, fingerprints, alerts, pharmacophore, search, searchtypes, searchtypesfiles, descriptors,
//        scaffoldtree, properties, reactions, mmpa, rgroups, fmcs, mces, molzip
// Engines: `rdkit` Chem's RDKit.js 2024.09 (RDKit mode), `crux` Chem's Crux file (vendor/chem-wasm) as Crux mode runs it
// (with the RDKit fallback Chem has behind it, counted), `rdkit2026` RDKit.js 2026.03.6 (the RDKit Crux reproduces),
// `py` Python RDKit 2026.3.6, `chempy` the chem-chem container's chem.py in that Python, `ocl` OpenChemLib (Chem's
// worker code), `sg` scripts/scaffold_tree.py on ScaffoldGraphDG 1.1.8. Each engine run is its own process (this file
// with --child), one thread; times are warm: each child first runs the call on a few items (excluded) after loading and
// compiling its module (excluded).
//
// OUT        output folder (default <tmpdir>/chem-compare100k): results.json, results.md, raw/ (each child's outputs)
// DATASET    SMILES csv (default <reddata>/data/demo/chem/chembl/chembl-100k.csv)
// PYTHON     a Python with RDKit 2026.3.6 (default <crux-core>/.venv/bin/python)
// SG_PYTHON  a Python with ScaffoldGraphDG 1.1.8, networkx, pandas and RDKit (the scaffold tree's oracle; without it the
//            row times Crux alone)
// SCAFFOLD_SIZES  the scaffold tree's sizes (default 1000,2000,5000,10000; Crux also at `rows`): the script's time
//            grows ~n^1.7 (6.7 s at 1000, 21 s at 2000 on an M5 Pro)
// SCAFFOLD_BUDGET_MIN  minutes the script may take at one size (default 45): past it it is killed and larger sizes skipped
// ENGINES    comma-separated engines to run (e.g. `crux` after a new Crux file); the others' outputs in OUT/raw are reused
// RDKITJS    a folder with node_modules/@rdkit/rdkit 2026.3.6 (default <crux-core>/tools/rdkitjs; `npm ci` there once)
// CRUX_BENCH crux-bench checkout (the search's query sets; default <grok>/crux/crux-bench)
// CONTAINER  the chem-chem container's URL (e.g. http://localhost:55005): its descriptors timed on 2000 molecules
// Outputs: OUT/results.json, OUT/results.md, OUT/processes.log (each process: exit, time, peak memory), OUT/inputs/
// (the data subsets: compare100k.py subsets), OUT/raw/<function>/<engine>.json
// SLICE3D    molecules for the 3D descriptors (default 1000)
import {registerHooks} from 'module';
import {spawnSync} from 'child_process';
import {Worker} from 'worker_threads';
import {existsSync, mkdirSync, readFileSync, statSync, writeFileSync} from 'fs';
import {createRequire} from 'module';
import {tmpdir} from 'os';
import {dirname, join, resolve as resolvePath} from 'path';
import {performance} from 'perf_hooks';
import {fileURLToPath, pathToFileURL} from 'url';
import {CHEM_DIR, DATA_DIR, loadRdkit, readSmilesCsv, requireFromChem} from './node-env.mjs';

const HERE = dirname(fileURLToPath(import.meta.url));
const SRC = join(CHEM_DIR, 'src');
const GROK = resolvePath(CHEM_DIR, '..', '..', '..', '..', '..');
const CRUX_CORE = join(GROK, 'crux', 'crux-core');
const OUT = process.env.OUT ?? join(tmpdir(), 'chem-compare100k');
const RAW = join(OUT, 'raw');
const DATASET = process.env.DATASET ?? join(DATA_DIR, 'demo', 'chem', 'chembl', 'chembl-100k.csv');
const PYTHON = process.env.PYTHON ?? join(CRUX_CORE, '.venv', 'bin', 'python');
const RDKITJS = process.env.RDKITJS ?? join(CRUX_CORE, 'tools', 'rdkitjs');
const CRUX_BENCH = process.env.CRUX_BENCH ?? join(GROK, 'crux', 'crux-bench');
const PY_HELPER = join(HERE, 'compare100k.py');
const VENDOR_GLUE = join(CHEM_DIR, 'vendor', 'chem-wasm', 'chem_wasm.js');
const WARMUP = 50;

// ---- Chem's TypeScript, transpiled on import; datagrok-api and the main thread's engine modules stubbed -------------

const ts = requireFromChem('typescript');
/** Module files replaced by stubs (their callers' exports only), by path under src/. */
const STUB_FILES = {
  'engine/chem-engine.ts': `
    export var ChemicalEngine; (function (e) { e.RDKit = 'RDKit'; e.Crux = 'Crux'; })(ChemicalEngine || (ChemicalEngine = {}));
    export const engineInEffect = () => globalThis.__engine.mode;
    export const isCruxRequested = () => globalThis.__engine.mode === 'Crux';
    export const engineGeneration = () => 0;
    const warn = (call, e) => globalThis.__engine.fallbacks.push([call, String(e?.message ?? e)]);
    export async function byEngineAsync(call, crux, rdkit) {
      if (globalThis.__engine.mode !== 'Crux') return rdkit();
      try { return await crux(); } catch (e) { warn(call, e); return rdkit(); } }
    export const renderingModule = (rdkit) => globalThis.__engine.mode === 'Crux' ?
      {module: globalThis.__engine.crux, engine: 'Crux'} : {module: rdkit, engine: 'RDKit'};
    export function byEngineModule(call, rdkit, job, fallback) {
      if (globalThis.__engine.mode !== 'Crux') return job(rdkit, false);
      let threw = false, crux = null;
      try { crux = job(globalThis.__engine.crux, true); } catch (e) { warn(call, e); threw = true; }
      if (crux !== null) return crux;
      const r = job(rdkit, false);
      if (!threw && r !== null) warn(call, fallback ?? 'Crux read no molecule, RDKit was asked');
      return r; }
    export function withMol(module, crux, molecule, read, asQuery = false) {
      let mol = null;
      try { mol = asQuery ? module.get_qmol(molecule) : module.get_mol(molecule); return mol ? read(mol) : null; }
      catch (e) { if (crux) throw e; return null; } finally { mol?.delete(); } }`,
  'engine/crux-module.ts': `
    export const cruxModuleSync = () => globalThis.__engine.mode === 'Crux' ? globalThis.__engine.crux : null;
    export const cruxLoadFailed = () => false;
    export const noteCruxError = () => {};
    export const compileCrux = async () => null;
    export const loadCrux = async () => globalThis.__engine.crux;`,
  'engine/fallback.ts': `
    export function warnCruxFallback(call, e) { globalThis.__engine.fallbacks.push([call, String(e?.message ?? e)]); }
    export function cruxOrRdKit(call, crux, rdkit) {
      try { return crux(); } catch (e) { warnCruxFallback(call, e); return rdkit(); } }
    export async function cruxOrRdKitAsync(call, crux, rdkit) {
      try { return await crux(); } catch (e) { warnCruxFallback(call, e); return rdkit(); } }`,
  'utils/chem-common-rdkit.ts': `
    export const getRdKitService = async () => globalThis.__engine.rdkitService;
    export const getRdKitModule = () => globalThis.__engine.rdkit;
    export const RDKIT_COMMON_RENDER_OPTS = ${readFileSync(join(SRC, 'utils', 'chem-common-rdkit.ts'), 'utf8')
    .match(/export const RDKIT_COMMON_RENDER_OPTS[^=]*= (\{[\s\S]*?\n\});/)[1]};`,
  // reactions.ts canonicalizes its inputs by engine (the Crux workers in Crux mode); the harness gives every engine the
  // canonical SMILES it computed beforehand, outside the timings, as the RDKit pool's answer
  'engine/worker-jobs.ts': `
    export const convertMolNotationByEngine = (molecules, notation, kekulize) =>
      globalThis.__engine.rdkitService.convertMolNotation(molecules, notation, kekulize);`,
  'rendering/rdkit-reaction-renderer.ts': `export const BRANCH_DELIMITER = '--**--';`,
  'rendering/rdkit-cell-renderer.ts': `export const _addColorsToBondsAndAtoms = () => {};`,
};
/** Exports added to modules that keep the functions measured here to themselves. */
const EXTRA_EXPORTS = {
  'crux/crux.worker.ts': 'export {rdkitShapedWorker};',
  'analysis/deprotect.ts': 'export {deprotectOne};',
};
const STUB_ANY = `
  const any = (name) => new Proxy(function () {}, {
    get: (t, p) => p === 'prototype' ? t.prototype : p === Symbol.toPrimitive ? () => name : any(name + '.' + String(p)),
    apply: () => any(name + '()'), construct: () => ({})});`;
const STUB_MODULES = {
  'datagrok-api/dg': STUB_ANY + `
    export class LruCache { constructor() { this.m = new Map(); }
      getOrCreate(k, f) { if (!this.m.has(k)) this.m.set(k, f(k)); return this.m.get(k); } }
    export const SEMTYPE = {MOLECULE: 'Molecule'}, TYPE = any('TYPE'), COLUMN_TYPE = any('COLUMN_TYPE'),
      FuncCallEditor = any('FuncCallEditor'), GridCellRenderer = any('GridCellRenderer'), Widget = any('Widget'),
      Column = any('Column'), DataFrame = any('DataFrame'), TaskBarProgressIndicator = any('Progress'),
      Color = any('Color'), chem = {isMolBlock: (s) => s.includes('M  END'), Notation: {Smiles: 'smiles',
      MolBlock: 'molblock', Smarts: 'smarts', Unknown: 'unknown', V3KMolBlock: 'v3Kmolblock'}},
      FLOAT_NULL = -3.4028234663852886e+38, UNITS = any('UNITS'), Func = any('Func'), InputBase = any('InputBase'),
      FILTER_TYPE = any('FILTER_TYPE'), BitSet = any('BitSet'), Viewer = any('Viewer'), Utils = any('Utils');`,
  'datagrok-api/grok': STUB_ANY + `
    export const shell = any('shell'), chem = {Notation: {Smiles: 'smiles', MolBlock: 'molblock'}}, data = any('data'),
      functions = any('functions'), events = any('events'), dapi = any('dapi'), userSettings = any('userSettings'),
      log = any('log'), decorators = any('decorators');`,
  'datagrok-api/ui': STUB_ANY + `
    export const div = any('div'), divText = any('divText'), divV = any('divV'), divH = any('divH'), input = any('input'),
      tooltip = any('tooltip'), iconFA = any('iconFA'), tools = any('tools'), dialog = any('dialog'), label = any('label'),
      button = any('button'), canvas = any('canvas'), tableFromMap = any('tableFromMap'), setUpdateIndicator = any('u');`,
};

const fileExists = (p) => existsSync(p) && statSync(p).isFile();
function candidates(base) {
  return [base, `${base}.ts`, `${base}.js`, `${base}.mjs`, join(base, 'index.ts'), join(base, 'index.js'),
    base.endsWith('.js') ? base.replace(/\.js$/, '.ts') : null].filter(Boolean);
}
function nodeModulesDirs(from) {
  const dirs = [];
  for (let d = from; ; d = dirname(d)) {
    dirs.push(join(d, 'node_modules'));
    if (dirname(d) === d)
      break;
  }
  return dirs;
}

registerHooks({
  resolve(specifier, context, nextResolve) {
    if (STUB_MODULES[specifier])
      return {url: `stub:${specifier}`, shortCircuit: true};
    const parent = context.parentURL?.startsWith('file:') ? fileURLToPath(context.parentURL) : null;
    let found = null;
    if (parent && (specifier.startsWith('./') || specifier.startsWith('../')))
      found = candidates(resolvePath(dirname(parent), specifier)).find(fileExists);
    else if (parent && specifier.startsWith('@datagrok-libraries/') || specifier.startsWith('datagrok-api/src/')) {
      for (const dir of nodeModulesDirs(dirname(parent))) {
        found = candidates(join(dir, specifier)).find(fileExists);
        if (found)
          break;
      }
    }
    if (found)
      return {url: pathToFileURL(found).href, shortCircuit: true, format: found.endsWith('.ts') ? 'module' : undefined};
    return nextResolve(specifier, context);
  },
  load(url, context, nextLoad) {
    if (url.startsWith('stub:'))
      return {format: 'module', source: STUB_MODULES[url.slice(5)], shortCircuit: true};
    if (!url.startsWith('file:'))
      return nextLoad(url, context);
    const path = fileURLToPath(url);
    const underSrc = path.startsWith(SRC + '/') ? path.slice(SRC.length + 1) : null;
    if (underSrc && STUB_FILES[underSrc] !== undefined)
      return {format: 'module', source: STUB_FILES[underSrc], shortCircuit: true};
    if (path.endsWith('.ts')) {
      let code = readFileSync(path, 'utf8');
      if (underSrc && EXTRA_EXPORTS[underSrc])
        code += '\n' + EXTRA_EXPORTS[underSrc];
      const js = ts.transpileModule(code, {fileName: path, compilerOptions: {module: ts.ModuleKind.ESNext,
        target: ts.ScriptTarget.ES2022, preserveConstEnums: true, experimentalDecorators: true,
        useDefineForClassFields: false}}).outputText;
      return {format: 'module', source: js, shortCircuit: true};
    }
    if (path === VENDOR_GLUE || path.includes('/datagrok-api/src/') || path.includes('/js-api/src/'))
      return {format: 'module', source: readFileSync(path, 'utf8'), shortCircuit: true};
    return nextLoad(url, context);
  },
});

// a worker's globals: crux.worker.ts registers its message handler at import, OpenChemLib reads `self`
globalThis.self ??= globalThis;
globalThis.addEventListener ??= () => {};

/** Imports a module of Chem's src/ (its TypeScript transpiled). */
const chem = (path) => import(pathToFileURL(join(SRC, path)).href);

// ---- the engines ---------------------------------------------------------------------------------------------------

/** The loaded engines' wasm memories (heapBytes: a child's restart when its engine's heap has grown too far). */
let chemRdkitModule = null;
let rdkit2026Module = null;
let cruxMemory = null;

/** Chem's Crux file, the instance the Crux worker's module (crux/crux.worker.ts) uses. */
async function loadCruxModule() {
  const glue = await import(pathToFileURL(VENDOR_GLUE).href);
  const wasm = await glue.default({module_or_path: readFileSync(join(CHEM_DIR, 'vendor', 'chem-wasm',
    'chem_wasm_bg.wasm'))});
  cruxMemory ??= wasm?.memory ?? null;
  glue.prefer_coordgen(false); // as engine/crux-module.ts sets up the main thread's instance
  glue.use_legacy_stereo_perception(false); // as crux.worker.ts and crux-module.ts
  return glue;
}

let rdkit2026 = null;
async function loadRdkit2026() {
  return rdkit2026 ??= loadRdkit2026Once();
}
async function loadRdkit2026Once() {
  const require = createRequire(join(RDKITJS, 'package.json'));
  const init = require(join(RDKITJS, 'node_modules', '@rdkit', 'rdkit', 'dist', 'RDKit_minimal.js'));
  const rdkit = await init();
  rdkit.prefer_coordgen(false);
  rdkit.use_legacy_stereo_perception(false);
  rdkit2026Module = rdkit;
  return rdkit;
}

let chemRdkit = null;
async function loadChemRdkit() {
  return chemRdkit ??= loadChemRdkitOnce();
}
async function loadChemRdkitOnce() {
  const rdkit = await loadRdkit(); // Chem's RDKit.js 2024.09, use_legacy_stereo_perception(false) as Chem's workers
  rdkit.prefer_coordgen(false); // as Chem's main thread (utils/chem-common-rdkit.ts)
  chemRdkitModule = rdkit;
  return rdkit;
}

/** The modules a child needs: `module` is the engine's, `rdkit` Chem's RDKit (Crux mode's fallback). */
async function engineModules(engine) {
  const rdkit = await loadChemRdkit();
  const crux = engine === 'crux' ? await loadCruxModule() : null;
  const module = engine === 'crux' ? crux : engine === 'rdkit2026' ? await loadRdkit2026() : rdkit;
  globalThis.__engine = {mode: engine === 'crux' ? 'Crux' : 'RDKit', crux, rdkit, fallbacks: [], rdkitService: null};
  return {module, rdkit, crux};
}

// ---- utilities -----------------------------------------------------------------------------------------------------

// A yield to the event loop (`setTimeout(r, 0)`, which Chem's worker and reaction code take every 50 molecules) costs
// no wait here: Node's 0 ms timer is 1 ms, which would add ~20 us per molecule that is no computation.
const realSetTimeout = globalThis.setTimeout;
globalThis.setTimeout = (fn, ms, ...rest) => !ms ? setImmediate(fn, ...rest) : realSetTimeout(fn, ms, ...rest);

/** A 64-bit hash (two 32-bit lanes) and the length of a string or bytes: outputs too large to keep whole. */
function hash(x) {
  let h1 = 0x811c9dc5;
  let h2 = 0x2b992ddf;
  const str = typeof x === 'string';
  for (let i = 0; i < x.length; i++) {
    const c = str ? x.charCodeAt(i) : x[i];
    h1 = Math.imul(h1 ^ c, 0x01000193);
    h2 = Math.imul(h2 ^ c, 0x5bd1e995);
    h2 ^= h2 >>> 13;
  }
  return (h1 >>> 0).toString(16).padStart(8, '0') + (h2 >>> 0).toString(16).padStart(8, '0') + ':' + x.length;
}

/** RFC 4180 CSV (quoted fields, doubled quotes) as objects by header. */
function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = '';
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"') {
        field += '"';
        i++;
      } else if (c === '"')
        quoted = false;
      else
        field += c;
    } else if (c === '"')
      quoted = true;
    else if (c === ',') {
      row.push(field);
      field = '';
    } else if (c === '\n' || c === '\r') {
      if (c === '\r' && text[i + 1] === '\n')
        i++;
      row.push(field);
      rows.push(row);
      row = [];
      field = '';
    } else
      field += c;
  }
  if (field || row.length) {
    row.push(field);
    rows.push(row);
  }
  const [header, ...body] = rows.filter((r) => r.length > 1 || r[0] !== '');
  return body.map((r) => Object.fromEntries(header.map((h, k) => [h, r[k] ?? ''])));
}
const chemFile = (name) => parseCsv(readFileSync(join(CHEM_DIR, 'files', name), 'utf8'));

let molCache = null;
/** The dataset's first `rows` SMILES. */
function molecules(rows) {
  molCache ??= readSmilesCsv(DATASET, rows);
  return molCache;
}

const tally = (values) => values.reduce((t, v) => (t[v] = (t[v] ?? 0) + 1, t), {});

// Crux panics (a Rust panic: console_error_panic_hook writes "panicked at ...", then the call traps with
// RuntimeError: unreachable), counted by the child's console.error hook; `traps`: the items whose call panicked
const PANICS = [];
/** Per panic, what the run was doing (globalThis.__panicContext(), where a runner sets one). */
const PANIC_CTX = [];
const TRAPS = [];
const panicCount = () => PANICS.length;
/** `one` on each item, the items whose call panicked recorded (TRAPS). */
function trackItem(i, run) {
  const before = PANICS.length;
  try {
    return run();
  } finally {
    if (PANICS.length > before)
      TRAPS.push(i);
  }
}
/** The items of a chunk whose call panics, each run alone (untimed; their panics not counted again). */
async function trapsIn(chunk, start, runOne) {
  const before = PANICS.length;
  const found = [];
  for (let j = 0; j < chunk.length; j++) {
    const k = PANICS.length;
    try {
      await runOne([chunk[j]]);
    } catch (e) {}
    if (PANICS.length > k)
      found.push(start + j);
  }
  PANICS.length = before;
  PANIC_CTX.length = before;
  return found;
}

/** `one` over the items, timed after `warmup` items (excluded; their outputs dropped). */
function timeSync(items, one, warmup = WARMUP) {
  for (let i = 0; i < Math.min(warmup, items.length); i++)
    one(items[i], i);
  TRAPS.length = 0;
  const out = new Array(items.length);
  const t0 = performance.now();
  for (let i = 0; i < items.length; i++)
    out[i] = trackItem(i, () => one(items[i], i));
  return {ms: performance.now() - t0, out, traps: [...TRAPS]};
}

const inputs = (name) => join(OUT, 'inputs', name);
const rawFile = (fn, engine, ext = 'json') => join(RAW, fn, `${engine}.${ext}`);
/** A Python output's path without its extensions (its folder made). */
function rawPrefix(fn, name) {
  mkdirSync(join(RAW, fn), {recursive: true});
  return join(RAW, fn, name);
}
function saveRaw(fn, engine, obj) {
  mkdirSync(join(RAW, fn), {recursive: true});
  writeFileSync(rawFile(fn, engine), JSON.stringify(obj));
}
function loadRaw(fn, engine) {
  const path = rawFile(fn, engine);
  return existsSync(path) ? JSON.parse(readFileSync(path, 'utf8')) : null;
}
function writeJson(path, obj) {
  mkdirSync(dirname(path), {recursive: true});
  writeFileSync(path, JSON.stringify(obj));
}
const readJson = (path) => JSON.parse(readFileSync(path, 'utf8'));
/** A file of float64s (its own copy: a small file's Buffer is a slice of a shared pool). */
function readF64(path) {
  const b = readFileSync(path);
  return new Float64Array(b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength));
}

/** Crux mode's fallbacks so far (engine/fallback.ts's warnCruxFallback, stubbed to count). */
const fallbacks = () => globalThis.__engine?.fallbacks?.length ?? 0;

// ---- an engine that aborts, traps, raises or dies: counted outcomes, not the end of a function --------------------
// RDKit.js aborts (Emscripten's abort(): "Aborted(...)" on stderr, then a trap: RuntimeError: unreachable) are counted
// by the child's console.error hook with the run's context (ABORTS), like Crux's panics. A child whose engine heap has
// grown past HEAP_RESTART_BYTES between two units of work asks for a fresh process (exit RESTART_EXIT): RDKit.js
// 2024.09 leaks in the reaction runs (its wasm heap reaches its 4 GB limit within one 100K run, after which every call
// aborts), and a run that only measured that would compare out-of-memory aborts. A child that dies (a trap or an
// exception nobody catches, out of memory) is started again by the parent (child()), which records the unit it died
// in (the child's `current` marker); the next process skips it: such a unit is run again in chunks, a crashed chunk item
// by item, a crashed item is a counted outcome ("crashed"). The units done are kept in raw/<fn>/<engine>.state-<run>/.
const RESTART_EXIT = 75;
const wanted = (engine) => !process.env.ENGINES || process.env.ENGINES.split(',').includes(engine);
const HEAP_RESTART_BYTES = parseFloat(process.env.HEAP_RESTART_GB ?? '1.5') * 2 ** 30;
const MAX_RESTARTS = 400;
const ABORTS = [];
/** Exceptions an RDKit.js reaction's run_reactants threw (a C++ exception: Chem's loop catches it as an error), with the
 * run's context: RDKit.js 2026.03.6 throws where Crux panics, and its instance can fail later in the process. */
const ENGINE_EXCEPTIONS = [];
/** The largest wasm heap of the engines loaded in this process (bytes). */
function heapBytes() {
  return Math.max(chemRdkitModule?.HEAP8?.length ?? 0, rdkit2026Module?.HEAP8?.length ?? 0,
    cruxMemory?.buffer?.byteLength ?? 0);
}
const gb = (b) => +(b / 2 ** 30).toFixed(2);
/** A child's units of work done so far, across its processes (the parent passes RESUME_DIR). */
class Resumable {
  constructor(fn, engine) {
    this.dir = process.env.RESUME_DIR ?? join(RAW, fn, `${engine}.state-${process.pid}`);
    mkdirSync(this.dir, {recursive: true});
    this.crashed = existsSync(join(this.dir, 'crashed.json')) ? readJson(join(this.dir, 'crashed.json')) : [];
  }
  file(key) {
    return join(this.dir, encodeURIComponent(key).replace(/%/g, '_') + '.json');
  }
  has(key) {
    return existsSync(this.file(key));
  }
  get(key) {
    return readJson(this.file(key));
  }
  /** The parent's record of a process that died in `key` (null: none). */
  crash(key) {
    return this.crashed.find((c) => c.key === key) ?? null;
  }
  start(key) {
    writeFileSync(join(this.dir, 'current.json'), JSON.stringify({key, heap: heapBytes(), aborts: ABORTS.length}));
  }
  finish(key, value) {
    this.doneHere = (this.doneHere ?? 0) + 1;
    writeFileSync(this.file(key), JSON.stringify(value));
    writeFileSync(join(this.dir, 'current.json'), JSON.stringify({key: null}));
    return value;
  }
  /** A fresh process from here when the engine's heap has grown past HEAP_RESTART_BYTES (exits). */
  restartIfHeapHigh(next) {
    const heap = heapBytes();
    // a Crux panic leaves its instance in an unknown state (Chem drops the module; here the rest of the unit ran on,
    // the next unit gets a fresh process)
    const panicked = PANICS.length > 0 || ENGINE_EXCEPTIONS.length > 0;
    if ((heap <= HEAP_RESTART_BYTES && !panicked) || !this.doneHere) // a fresh process does one unit at least
      return;
    const log = join(this.dir, 'restarts.json');
    const restarts = existsSync(log) ? readJson(log) : [];
    restarts.push({before: next, heapGb: gb(heap), ...(panicked ? {afterPanics: PANICS.length + ENGINE_EXCEPTIONS.length}
      : {})});
    writeFileSync(log, JSON.stringify(restarts));
    console.log(`  restart before ${next}: heap ${gb(heap)} GB${panicked ? `, ${PANICS.length} Crux panics` : ''}`);
    process.exit(RESTART_EXIT);
  }
  restarts() {
    const log = join(this.dir, 'restarts.json');
    return existsSync(log) ? readJson(log) : [];
  }
}
/** A unit that hangs (an engine looping in wasm, which no JS timer can interrupt) is ended by a watchdog thread: past
 * the unit's limit in CPU time (this machine is shared: wall time says little) it writes hung.json and kills the
 * process; the parent records the unit as crashed ("hung"). */
const WATCH = new Int32Array(new SharedArrayBuffer(12)); // [watching, unit start (process CPU ms), its limit (CPU s)]
let watchdog = null;
const cpuMs = () => {
  const u = process.cpuUsage();
  return Math.floor((u.user + u.system) / 1000);
};
function watch(state, key) {
  watchdog ??= new Worker(`
    const {workerData} = require('worker_threads');
    const fs = require('fs');
    const w = new Int32Array(workerData.sab);
    const cpuMs = () => {
      const u = process.cpuUsage();
      return Math.floor((u.user + u.system) / 1000);
    };
    setInterval(() => {
      const limit = Atomics.load(w, 2);
      if (Atomics.load(w, 0) && (cpuMs() - Atomics.load(w, 1)) / 1000 > limit) {
        fs.writeFileSync(workerData.file, JSON.stringify({seconds: limit}));
        process.kill(process.pid, 'SIGKILL');
      }
    }, 2000);`, {eval: true, workerData: {sab: WATCH.buffer, file: join(state.dir, 'hung.json')}});
  watchdog.unref();
  // limits in CPU s: ~4x the slowest normal unit here (a 100K reaction unit 71 s on RDKit.js 2024.09; an FMCS pair is
  // bounded by its 60 s timeout)
  const limit = /^pair /.test(key) ? 200 : /#\d+$/.test(key) ? 30 : /#\d+-\d+$/.test(key) ? 60 : 300;
  Atomics.store(WATCH, 2, limit);
  Atomics.store(WATCH, 1, cpuMs());
  Atomics.store(WATCH, 0, 1);
}
const unwatch = () => Atomics.store(WATCH, 0, 0);
/** `run(key)` for `key` unless done in an earlier process (its saved result) or crashed there (null). */
async function unitOnce(state, key, run) {
  if (state.has(key))
    return state.get(key);
  if (state.crash(key))
    return null;
  state.restartIfHeapHigh(key);
  state.start(key);
  if (process.env.CRASH_TEST && new RegExp(process.env.CRASH_TEST).test(key)) // the harness's own test: a death here
    process.kill(process.pid, 'SIGKILL');
  let value;
  try {
    watch(state, key);
    value = await run();
    unwatch();
  } catch (e) { // the engine left in an unknown state: the parent records the unit as crashed and starts a process
    console.error(`Error: ${key} threw: ${String(e?.message ?? e).slice(0, 300)}`);
    process.exit(70);
  }
  return state.finish(key, value);
}
/** Aborts (RDKit.js) since `before`: count, the heap then, and up to 3 contexts. */
function abortsSince(before) {
  const list = ABORTS.slice(before);
  return list.length ? {count: list.length, heapGb: gb(Math.max(...list.map((a) => a.heap))),
    contexts: list.slice(0, 3).map((a) => a.ctx)} : null;
}

// ---- the runs, one per function and engine (each in its own process) -----------------------------------------------

const RUN = {};
/** Per function and engine, the call on one item with its whole output: the parent's examples of differences. */
const ONE = {};
export {ONE, RUN, chem, engineModules, loadCruxModule, loadRdkit2026, loadChemRdkit, molecules, svgGeometry, roundSvg};

// parse: what each engine reads, the way Chem reads a column (getMolSafe: sanitized, then kekulize: false, then a query)
ONE.parse = async (engine) => {
  const {module} = await engineModules(engine);
  const {getMolSafe} = await chem('utils/mol-creation_rdkit.ts');
  return (s) => {
    if (!s)
      return {status: 'empty', out: ''};
    const ctx = getMolSafe(s, {}, module);
    if (!ctx.mol)
      return {status: 'failed', out: ''};
    try {
      return {status: ctx.isQMol ? 'query' : ctx.kekulize ? 'molecule' : 'molecule (kekulize: false)',
        out: ctx.isQMol ? ctx.mol.get_smarts() : ctx.mol.get_smiles()};
    } catch (e) {
      return {status: 'threw', out: 'THROWS'};
    } finally {
      ctx.mol.delete();
    }
  };
};
RUN.parse = async (engine, rows) => {
  const one = await ONE.parse(engine);
  const mols = molecules(rows);
  const {ms, out, traps} = timeSync(mols, one);
  return {n: mols.length, ms, counts: tally(out.map((o) => o.status)), traps, items: out.map((o) => o.out)};
};

// drawing: the cell renderer's SVG (rdkit-cell-renderer.ts _molToSvg: _fetchMolGetOrCreate's depiction, _infoToSvg's
// details), plain and highlighted by a scaffold (benzene, as "highlight by scaffold" matches it); Crux mode draws a
// molecule Crux cannot with RDKit
const SVG_W = 200;
const SVG_H = 100;
const HIGHLIGHT = 'c1ccccc1';
const roundSvg = (svg) => svg.replace(/-?\d+\.\d+/g, (x) => (+x).toFixed(1));
/** The picture's geometry alone: each path's and label's coordinates as written, a straight line's two ends in a fixed
 * order, sorted (drawing order, atom and bond numbering, a line's direction and style attributes aside). */
const svgGeometry = (svg) => [...svg.matchAll(/ d='([^']*)'|<text[^>]*x='([^']*)'[^>]*y='([^']*)'[^>]*>([^<]*)</g)]
  .map((m) => {
    if (m[1] === undefined)
      return `${m[2]},${m[3]},${m[4]}`;
    const line = m[1].match(/^M ([-\d.]+,[-\d.]+) L ([-\d.]+,[-\d.]+)$/);
    return line ? [line[1], line[2]].sort().join(' ') : m[1];
  }).sort().join('|');
ONE.drawing = async (engine) => {
  const {module, rdkit} = await engineModules(engine);
  const {getMolSafe} = await chem('utils/mol-creation_rdkit.ts');
  const {RDKIT_COMMON_RENDER_OPTS} = await chem('utils/chem-common-rdkit.ts');
  /** The atoms' places after the depiction (element and x, y from the MOL block, sorted): the layout. */
  const layoutOf = (mol) => {
    const lines = mol.get_molblock().split('\n');
    const na = Number(lines[3].slice(0, 3));
    return lines.slice(4, 4 + na).map((l) => `${l.slice(31, 34).trim()}${l.slice(0, 20).replace(/\s+/g, ',')}`).sort()
      .join(';');
  };
  const svg = (mod, s, highlight, layout = false) => {
    const ctx = getMolSafe(s, {}, mod, true);
    const mol = ctx.mol;
    if (!mol)
      return null;
    try {
      const substruct = {atoms: [], bonds: [], highlightAtomColors: {}, highlightBondColors: {}};
      if (highlight) { // the scaffold read as the renderer reads one (isSubstructure: get_qmol, aromatic form)
        const query = mod.get_qmol(HIGHLIGHT);
        try {
          query.convert_to_aromatic_form();
          const matches = JSON.parse(mol.get_substruct_matches(query)); // '{}' for none
          for (const match of Array.isArray(matches) ? matches : []) {
            for (const a of match.atoms ?? []) {
              if (!substruct.atoms.includes(a))
                substruct.atoms.push(a);
              substruct.highlightAtomColors[a] = [1.0, 0.7, 0.7, 1.0];
            }
            for (const b of match.bonds ?? []) {
              if (!substruct.bonds.includes(b))
                substruct.bonds.push(b);
              substruct.highlightBondColors[b] = [1.0, 0.7, 0.7, 1.0];
            }
          }
        } finally {
          query.delete();
        }
      }
      let own = mol.has_coords() > 0;
      const useMolBlockWedging = own;
      if (mol.has_coords() === 0) {
        mol.set_new_coords(false);
        mol.normalize_depiction(1);
        mol.straighten_depiction();
        own = false;
      }
      mol.normalize_depiction(own ? 0 : 1);
      mol.straighten_depiction(own);
      if (!mol.is_valid())
        return null;
      if (layout)
        return layoutOf(mol);
      const details = {...RDKIT_COMMON_RENDER_OPTS, width: SVG_W, height: SVG_H, ...substruct};
      if (!ctx.kekulize)
        details.kekulize = false;
      if (useMolBlockWedging) {
        details.useMolBlockWedging = true;
        details.wedgeBonds = false;
        details.addChiralHs = false;
      }
      return mol.get_svg_with_highlights(JSON.stringify(details));
    } finally {
      mol.delete();
    }
  };
  const draw = (s) => {
    const out = {fallback: false};
    for (const [key, hl] of [['plain', false], ['highlighted', true]]) {
      let picture = null;
      try {
        picture = svg(module, s, hl);
      } catch (e) {
        picture = engine === 'crux' ? null : 'THROWS';
      }
      if (picture === null && engine === 'crux') { // _fetchMol / _molToSvg: RDKit draws what Crux cannot
        out.fallback = true;
        try {
          picture = svg(rdkit, s, hl);
        } catch (e) {
          picture = 'THROWS';
        }
      }
      out[key] = picture;
    }
    return out;
  };
  draw.layout = (s) => { // the layout the plain picture is drawn from (Crux mode: RDKit's where RDKit drew it)
    let l = null;
    try {
      l = svg(module, s, false, true);
    } catch (e) {}
    if (l === null && engine === 'crux') {
      try {
        l = svg(rdkit, s, false, true);
      } catch (e) {}
    }
    return l ?? '';
  };
  return draw;
};
RUN.drawing = async (engine, rows) => {
  const one = await ONE.drawing(engine);
  const mols = molecules(rows);
  const {ms, out, traps} = timeSync(mols, one);
  const status = out.map((o) => o.plain === null ? 'failed' : o.plain === 'THROWS' ? 'threw' : 'drawn');
  const layouts = mols.map((s) => one.layout(s)); // untimed: the atoms' places, to tell layout from drawing
  return {n: mols.length, ms, counts: tally(status), traps, fallbacks: out.flatMap((o, i) => o.fallback ? [i] : []), layouts,
    items: out.map((o) => [o.plain, o.highlighted].map((svg) => svg === null ? '' :
      `${hash(svg)} ${hash(roundSvg(svg))} ${hash(svgGeometry(svg))}`))};
};

/** Whether two layouts (layoutOf) place the same elements at the same places, within 0.01 (MOL block units, Å). */
function sameLayout(a, b) {
  if (a === b)
    return true;
  const parse = (l) => l.split(';').map((x) => x.split(',')).map(([el, x, y]) => [el, +x, +y]);
  const [pa, pb] = [parse(a), parse(b)];
  if (pa.length !== pb.length)
    return false;
  const used = new Uint8Array(pb.length);
  return pa.every(([el, x, y]) => {
    const k = pb.findIndex(([e2, x2, y2], j) => !used[j] && e2 === el && Math.abs(x - x2) <= 0.01 && Math.abs(y - y2) <= 0.01);
    if (k < 0)
      return false;
    used[k] = 1;
    return true;
  });
}

// fingerprints: the RDKit worker's getFingerprints (rdkit-service-worker-similarity.ts) on each module; Crux mode as
// crux.worker.ts (crux/jobs/batch.ts: the whole chunk in one call) + engine/worker-jobs.ts run it (the RDKit pool for
// the molecules Crux cannot read, but for FCFP and chirality, which are Crux's in both modes)
const FP_VARIANTS = [
  ['Morgan', 'Morgan', null], ['Pattern', 'Pattern', null], ['RDKit', 'RDKit', null], ['MACCS', 'MACCS', null],
  ['AtomPair', 'AtomPair', null], ['TopologicalTorsion', 'TopologicalTorsion', null],
  ['Morgan radius 3', 'Morgan', {radius: 3}], ['Morgan FCFP', 'Morgan', {useFeatures: true}],
  ['Morgan chirality', 'Morgan', {useChirality: true}], ['AtomPair chirality', 'AtomPair', {useChirality: true}],
  ['TopologicalTorsion chirality', 'TopologicalTorsion', {useChirality: true}],
];
const cruxOnlyFp = (options) => !!options?.useFeatures || !!options?.useChirality;
ONE.fingerprints = async (engine) => {
  const {module, rdkit, crux} = await engineModules(engine);
  const {RdKitServiceWorker} = await chem('rdkit-service/rdkit-service-worker.ts');
  const fo = await chem('utils/fingerprint-options.ts');
  const batch = engine === 'crux' ? await chem('crux/jobs/batch.ts') : null;
  return async (type, options, mols) => {
    const opts = fo.fingerprintOptions(type, options);
    if (engine !== 'crux')
      return {fps: (await new RdKitServiceWorker(module, '').getFingerprints(type, mols, false, opts)).fps, fallback: []};
    fo.assertFingerprintOptions(crux, type, opts);
    const {fps, unread} = batch.fingerprints(crux, type, mols, false, opts);
    const fallback = [];
    if (unread.length && !fo.needsCrux(opts)) {
      const patch = await new RdKitServiceWorker(rdkit, '').getFingerprints(type, unread.map((i) => mols[i]), false, opts);
      unread.forEach((row, k) => {
        fps[row] = patch.fps[k];
        if (patch.fps[k])
          fallback.push(row);
      });
    }
    return {fps, fallback};
  };
};
RUN.fingerprints = async (engine, rows) => {
  const compute = await ONE.fingerprints(engine);
  const mols = molecules(rows);
  const variants = {};
  for (const [name, type, options] of FP_VARIANTS) {
    if (engine === 'rdkit' && cruxOnlyFp(options))
      continue; // Crux's in both modes
    await compute(type, options, mols.slice(0, WARMUP));
    const t0 = performance.now();
    const c0 = process.cpuUsage();
    const before = panicCount();
    const {fps, fallback} = await compute(type, options, mols);
    const ms = performance.now() - t0;
    const cpu = process.cpuUsage(c0);
    const cpuMs = (cpu.user + cpu.system) / 1000; // the process's CPU time: less moved by a busy machine than ms
    const traps = panicCount() > before ? await trapsIn(mols, 0, (one) => compute(type, options, one)) : [];
    variants[name] = {ms, cpuMs, counts: tally(fps.map((f) => f ? 'computed' : 'none')), fallbacks: fallback, traps,
      items: fps.map((f) => f ? hash(f) : '')};
  }
  return {n: mols.length, variants};
};

// structural alerts and pharmacophore features: the RDKit worker's getStructuralAlerts, Crux mode crux/jobs/batch.ts's
// structuralAlerts (the whole chunk in one call; + the RDKit pool's for what Crux cannot read, engine/worker-jobs.ts)
function alertRules() {
  const rules = {};
  for (const r of chemFile('alert-collection.csv'))
    (rules[r.rule_set_name] ??= []).push(r.smarts);
  return rules;
}
function pharmacophoreRules() {
  const families = {D: 'Donor', A: 'Acceptor', H: 'Hydrophobic', a: 'Aromatic', P: 'Positive', N: 'Negative',
    X: 'Halogen Bond'};
  const rules = {};
  for (const r of chemFile('pharmacophore-features.csv')) {
    if (families[r.family] && r.smarts)
      (rules[families[r.family]] ??= []).push(r.smarts);
  }
  return rules;
}
async function alertsRunner(engine) {
  const {module, rdkit, crux} = await engineModules(engine);
  const {RdKitServiceWorker} = await chem('rdkit-service/rdkit-service-worker.ts');
  const batch = engine === 'crux' ? await chem('crux/jobs/batch.ts') : null;
  return (rules, mols) => {
    if (engine !== 'crux')
      return {found: new RdKitServiceWorker(module, '').getStructuralAlerts(rules, mols), fallback: []};
    const {alerts: found, unread} = batch.structuralAlerts(crux, rules, mols);
    const fallback = [];
    if (unread.length) {
      const patch = new RdKitServiceWorker(rdkit, '').getStructuralAlerts(rules, unread.map((i) => mols[i]));
      for (const rule of Object.keys(found)) {
        unread.forEach((row, k) => {
          found[rule][row] = patch[rule]?.[k] ?? found[rule][row];
        });
      }
      fallback.push(...unread);
    }
    return {found, fallback};
  };
}
const flagsOf = (found, i) => Object.keys(found).filter((rule) => found[rule][i]).join(',');
for (const [fn, rulesOf] of [['alerts', alertRules], ['pharmacophore', pharmacophoreRules]]) {
  ONE[fn] = async (engine) => {
    const run = await alertsRunner(engine);
    const rules = rulesOf();
    return (s) => flagsOf(run(rules, [s]).found, 0);
  };
  RUN[fn] = async (engine, rows) => {
    const run = await alertsRunner(engine);
    const rules = rulesOf();
    const mols = molecules(rows);
    run(rules, mols.slice(0, WARMUP));
    const t0 = performance.now();
    const c0 = process.cpuUsage();
    const before = panicCount();
    const {found, fallback} = run(rules, mols);
    const ms = performance.now() - t0;
    const cpu = process.cpuUsage(c0);
    const cpuMs = (cpu.user + cpu.system) / 1000; // the process's CPU time: less moved by a busy machine than ms
    const traps = panicCount() > before ? await trapsIn(mols, 0, (one) => run(rules, one)) : [];
    const items = mols.map((_, i) => flagsOf(found, i));
    return {n: mols.length, ms, cpuMs, rules: Object.fromEntries(Object.entries(rules).map(([r, s]) => [r, s.length])),
      perRule: Object.fromEntries(Object.keys(found).map((r) => [r, found[r].filter((x) => x).length])),
      fallbacks: fallback, traps, counts: {'molecules with a match': items.filter((x) => x).length}, items};
  };
}

// descriptors: crux's calcDescriptors (crux.worker.ts, chunks of 1000 as crux-service.ts sends them) vs chem.py
async function descriptorNames(dims) {
  const crux = await loadCruxModule();
  const tree = JSON.parse(crux.descriptorTree());
  return Object.values(tree).flatMap((g) => g.descriptors.filter((d) => d.tags.includes(dims)).map((d) => d.name));
}
async function cruxDescriptors(rows, dims, count) {
  const {crux} = await engineModules('crux');
  const names = await descriptorNames(dims);
  const mols = molecules(rows).slice(0, count);
  const calc = (chunk) => JSON.parse(crux.calcDescriptors(chunk, names));
  calc(mols.slice(0, WARMUP));
  const values = new Float64Array(names.length * mols.length).fill(NaN);
  const types = {};
  const chunksPanicked = [];
  const t0 = performance.now();
  for (let start = 0; start < mols.length; start += 1000) {
    const before = panicCount();
    const res = calc(mols.slice(start, start + 1000));
    if (panicCount() > before)
      chunksPanicked.push(start);
    names.forEach((name, k) => {
      types[name] = res[name].type;
      res[name].value.forEach((v, j) => {
        if (v !== null)
          values[k * mols.length + start + j] = v;
      });
    });
  }
  const ms = performance.now() - t0;
  const traps = [];
  for (const start of chunksPanicked)
    traps.push(...await trapsIn(mols.slice(start, start + 1000), start, (one) => calc(one)));
  writeFileSync(rawFile(dims === '2D' ? 'descriptors' : 'descriptors3d', 'crux', 'f64'), Buffer.from(values.buffer));
  // a row of nulls only: Chem asks the container for it (descriptors-engine.ts fromContainerWhereUnread)
  const unread = mols.flatMap((m, i) => m && names.every((_, k) => Number.isNaN(values[k * mols.length + i])) ? [i] : []);
  return {n: mols.length, ms, names, types, fallbacks: unread, traps};
}
RUN.descriptors = async (engine, rows) => {
  mkdirSync(join(RAW, 'descriptors'), {recursive: true});
  return cruxDescriptors(rows, '2D', Infinity);
};
RUN.descriptors3d = async (engine, rows) => {
  mkdirSync(join(RAW, 'descriptors3d'), {recursive: true});
  return cruxDescriptors(rows, '3D', parseInt(process.env.SLICE3D ?? '1000'));
};

// the scaffold tree: crux's scaffoldTree (the Crux worker's call) at each size, ringCutoff 10, no discharge
// (getScaffoldTree's defaults)
// (OUT/scaffold-sizes.txt, when present, overrides SCAFFOLD_SIZES: a run already started can be given fewer sizes)
const scaffoldSizes = (rows) => (existsSync(join(OUT, 'scaffold-sizes.txt')) ? readFileSync(join(OUT,
  'scaffold-sizes.txt'), 'utf8').trim() : process.env.SCAFFOLD_SIZES ?? '1000,2000,5000,10000').split(',').map(Number)
  .filter((n) => n < rows).concat([rows]);
RUN.scaffoldtree = async (engine, rows) => {
  mkdirSync(join(RAW, 'scaffoldtree'), {recursive: true});
  const {crux} = await engineModules('crux');
  const mols = molecules(rows);
  crux.scaffoldTree(mols.slice(0, WARMUP), 10, false);
  const sizes = {};
  for (const size of scaffoldSizes(mols.length)) {
    const t0 = performance.now();
    try {
      const tree = crux.scaffoldTree(mols.slice(0, size), 10, false);
      sizes[size] = {ms: performance.now() - t0};
      writeFileSync(rawFile('scaffoldtree', `crux-${size}`), tree);
    } catch (e) { // Crux reproduces the script's own exceptions (e.g. recurs_append_nodes' IndexError): an outcome
      sizes[size] = {ms: performance.now() - t0, raised: String(e?.message ?? e).slice(0, 300)};
    }
  }
  return {n: mols.length, sizes};
};

// the property table: crux's chemProperties (crux-service.ts, chunks; molecules over MAX_SMILES_LENGTH as null, as
// properties.ts sends them) vs OpenChemLib through Chem's OCL worker (ocl-worker.ts, the RDKit-mode path)
const PROPERTY_NAMES = ['MW', 'HBA', 'HBD', 'LogP', 'LogS', 'PSA', 'Rotatable bonds', 'Stereo centers',
  'Molecule charge', 'Molecular formula'];
ONE.properties = async (engine) => {
  if (engine === 'ocl') {
    globalThis.onmessage = null;
    let posted = null;
    globalThis.postMessage = (m) => posted = m;
    await chem('open-chem/ocl-service/workers/ocl-worker.ts');
    const handler = globalThis.onmessage;
    return (chunk) => {
      handler({data: {op: 'chem_properties', data: chunk, argList: PROPERTY_NAMES, notationType: 'smiles'}});
      return posted.res;
    };
  }
  const {crux} = await engineModules('crux');
  const {MAX_SMILES_LENGTH} = await chem('utils/chem-constants.ts');
  return (chunk) => JSON.parse(crux.chemProperties(chunk.map((m) => typeof m === 'string' && m !== '' &&
    (m.includes('M  END') || m.length <= MAX_SMILES_LENGTH) ? m : null), PROPERTY_NAMES));
};
RUN.properties = async (engine, rows) => {
  const calc = await ONE.properties(engine);
  const mols = molecules(rows);
  calc(mols.slice(0, WARMUP));
  const values = Object.fromEntries(PROPERTY_NAMES.map((p) => [p, []]));
  const chunksPanicked = [];
  const t0 = performance.now();
  for (let start = 0; start < mols.length; start += 1000) {
    const before = panicCount();
    const res = calc(mols.slice(start, start + 1000));
    if (panicCount() > before)
      chunksPanicked.push(start);
    for (const p of PROPERTY_NAMES)
      values[p].push(...res[p]);
  }
  const ms = performance.now() - t0;
  const traps = [];
  for (const start of chunksPanicked)
    traps.push(...await trapsIn(mols.slice(start, start + 1000), start, (one) => calc(one)));
  return {n: mols.length, ms, traps, counts: tally(values.MW.map((v) => v === null || v === undefined ? 'none' : 'computed')),
    values};
};

// reactions (utils/reactions/reactions.ts, analysis/deprotect.ts, utils/reaction-enumeration/enumerate.ts) on the
// engine's module through Chem's EngineReaction: Crux mode on the Crux module with RDKit behind it per molecule, RDKit
// mode on Chem's RDKit, the reference with RDKit.js 2026.03.6 as the module. The canonical SMILES Chem reacts
// (RdKitService.convertMolNotation, RDKit in both modes) are computed once beforehand by Chem's RDKit worker
// (inputs/reaction-inputs.json), outside the engines' times.
const DEFAULT_REACTIONS = () => readJson(join(CHEM_DIR, 'files', 'reactions', 'default-reactions.json'));
const reactantCount = (smarts) => smarts.split('>>')[0].split('.').length;
const variablesOf = (r) => Object.fromEntries(Object.entries(r.variables ?? {}).map(([k, v]) => [k, v.defaultValue]));
const DEPROTECT_FRAGMENTS = {Fmoc: 'O=C([*:1])OCC1c2ccccc2-c2ccccc21', Boc: 'CC(C)(C)OC(=O)[*:1]'};

async function reactionSetup(engine) {
  const {module, rdkit} = await engineModules(engine);
  const rx = await chem('utils/reactions/reactions.ts');
  const dp = await chem('analysis/deprotect.ts');
  const {EngineReaction} = await chem('utils/reactions/reaction-engine.ts');
  const en = await chem('utils/reaction-enumeration/enumerate.ts');
  const cfg = await chem('utils/reaction-enumeration/config.ts');
  const {hasNewLines} = await chem('utils/chem-common.ts');
  const {MAX_SMILES_LENGTH} = await chem('utils/chem-constants.ts');
  const lastRead = []; // the modules' last two get_mol inputs: the molecule(s) a panic or an abort came from
  const noted = (s) => {
    lastRead.push(s);
    if (lastRead.length > 2)
      lastRead.shift();
  };
  // the `rdkit` argument (RDKit mode's module), its get_mol inputs noted (an RDKit.js abort's context)
  const rdkitModule = engine === 'rdkit2026' ? module : rdkit;
  const wrapRxn = (rxn) => new Proxy(rxn, {get: (t, p) => p === 'run_reactants' ? (list, max) => {
    try {
      return t.run_reactants(list, max);
    } catch (e) {
      ENGINE_EXCEPTIONS.push({ctx: globalThis.__panicContext?.() ?? null, what: String(e?.constructor?.name ?? e)});
      throw e;
    }
  } : typeof t[p] === 'function' ? t[p].bind(t) : t[p]});
  const rdkitParam = new Proxy(rdkitModule, {get: (t, p) => p === 'get_mol' ? (s, d) => {
    noted(s);
    return d === undefined ? t.get_mol(s) : t.get_mol(s, d);
  } : p === 'get_rxn' ? (s, d) => {
    const rxn = d === undefined ? t.get_rxn(s) : t.get_rxn(s, d);
    return rxn ? wrapRxn(rxn) : rxn;
  } : t[p]});
  if (engine === 'crux') { // as Chem's own test hook wraps it (reaction-engine.ts setCruxUnreadableReactantsForTests)
    const cruxModule = globalThis.__engine.crux;
    globalThis.__engine.crux = {...cruxModule, get_mol: (smiles, details) => {
      noted(smiles);
      return cruxModule.get_mol(smiles, details);
    }, get_rxn: (smarts, details) => {
      // a panic inside run_reactants leaves the reaction object borrowed (wasm-bindgen: "attempted to take ownership
      // of Rust value while it was borrowed"), so every later molecule of the run would fall back to RDKit: the
      // reaction is made again after a panic, and the molecules after it are still Crux's
      let rxn = cruxModule.get_rxn(smarts, details);
      if (!rxn)
        return rxn;
      return new Proxy({}, {get: (_, prop) => {
        if (prop === 'run_reactants') {
          return (list, max) => {
            const before = panicCount();
            try {
              return rxn.run_reactants(list, max);
            } finally {
              if (panicCount() > before)
                rxn = cruxModule.get_rxn(smarts, details);
            }
          };
        }
        if (prop === 'delete')
          return () => rxn.delete();
        const v = rxn[prop];
        return typeof v === 'function' ? v.bind(rxn) : v;
      }});
    }};
  }
  const canonical = new Map(); // RdKitService.convertMolNotation's answers, by the list asked for
  globalThis.__engine.rdkitService = {convertMolNotation: async (list) => {
    if (!canonical.has(list))
      throw new Error('a list without its canonical SMILES');
    return canonical.get(list);
  }};
  const deprotect = (fragment, mols) => { // addDeprotectedColumn's work, its column aside
    const fragmentMol = rdkitParam.get_mol(fragment);
    fragmentMol.add_hs_in_place();
    const hFragment = fragmentMol.get_smiles();
    fragmentMol.delete();
    const reaction = EngineReaction.create('deprotection', rdkitParam, `${hFragment}>>[*:1]`);
    try {
      return mols.map((s) => {
        try {
          if (s && !hasNewLines(s) && s.length > MAX_SMILES_LENGTH)
            return s;
          return reaction.run((mod, rxn) => dp.deprotectOne(mod, rxn, s));
        } catch (_) {
          return s;
        }
      });
    } finally {
      reaction.delete();
    }
  };
  const enumeration = async (variant) => {
    const dir = join(CHEM_DIR, 'files', 'enumerations');
    const config = variant === 'config.yaml' ? cfg.configFromYaml(readFileSync(join(dir, 'config.yaml'), 'utf8')) :
      cfg.cloneConfig(cfg.DEFAULT_CONFIG);
    const [tFile, bFile] = variant === 'amino acids' ? ['aa_reactions.csv', 'aa_bb.csv'] : ['reactions.csv', 'bb.csv'];
    const templates = parseCsv(readFileSync(join(dir, tFile), 'utf8')).filter((r) => r.reaction_smarts?.trim())
      .map((r) => ({smarts: r.reaction_smarts.trim(), reactionName: r.reaction_name ?? '', propagated: '',
        blockingSmartsList: (r.blocking_fg ?? '').split(/[;|]/).map((s) => s.trim()).filter((s) => s.length > 0)}));
    const buildingBlocks = parseCsv(readFileSync(join(dir, bFile), 'utf8')).map((r) => r.SMILES)
      .filter((s) => s && s.trim());
    const exclusionSmarts = parseCsv(readFileSync(join(dir, 'ex_smarts.csv'), 'utf8')).map((r) => r.SMARTS)
      .filter((s) => s && s.trim());
    return (await en.enumerate({rdkit: rdkitParam, config, templates, buildingBlocks, exclusionSmarts, reagents: []}))
      .rows;
  };
  return {rx, rdkitParam, canonical, deprotect, enumeration, lastRead};
}

/** A unit over items [0, n): whole; if a process died in it, again in chunks of 1000; a chunk it died in, item by item;
 * an item it died on: 'CRASHED' (past 30 deaths in one unit the rest is 'NOT RUN'). `run(a, b)` gives the outputs of
 * items [a, b) ({ms, items, numbers...}); the parts merged (items concatenated, numbers summed). */
async function rangedUnit(state, key, n, run) {
  const whole = await unitOnce(state, key, () => run(0, n));
  if (whole)
    return whole;
  const parts = [];
  const deaths = () => state.crashed.filter((c) => c.key === key || c.key.startsWith(`${key}#`)).length;
  const crashedItems = [];
  for (let a = 0; a < n; a += 1000) {
    const b = Math.min(n, a + 1000);
    if (deaths() > 30) {
      parts.push({items: new Array(b - a).fill('NOT RUN'), notRun: b - a});
      continue;
    }
    const chunk = await unitOnce(state, `${key}#${a}-${b}`, () => run(a, b));
    if (chunk) {
      parts.push(chunk);
      continue;
    }
    for (let i = a; i < b; i++) {
      const one = deaths() > 30 ? {items: ['NOT RUN'], notRun: 1} : await unitOnce(state, `${key}#${i}`, () => run(i, i + 1));
      if (!one)
        crashedItems.push(i);
      parts.push(one ?? {items: ['CRASHED']});
    }
  }
  const out = {items: [], crashedItems, crashes: state.crashed.filter((c) => c.key === key || c.key.startsWith(`${key}#`))};
  for (const p of parts) {
    for (const [k, v] of Object.entries(p)) {
      if (k === 'items')
        out.items.push(...v);
      else if (k === 'abortContexts' || k === 'trapContexts' || k === 'exceptionContexts')
        (out[k] ??= []).push(...v);
      else if (typeof v === 'number')
        out[k] = (out[k] ?? 0) + v;
    }
  }
  return out;
}

RUN.reactions = async (engine, rows) => {
  const state = new Resumable('reactions', engine);
  const {rx, rdkitParam, canonical, deprotect, enumeration, lastRead} = await reactionSetup(engine);
  let current = null;
  globalThis.__panicContext = () => ({part: current?.[0], key: current?.[1], inputs: [...lastRead]});
  const prepared = readJson(inputs(`reaction-inputs-${rows}.json`));
  const mols = prepared.molecules;
  const warm = mols.slice(0, WARMUP);
  canonical.set(warm, prepared.canonical.slice(0, WARMUP));
  const out = {n: mols.length, transformation: {}, twoComponent: {}, deprotect: {}, enumeration: {}};
  const reactions = DEFAULT_REACTIONS();
  const warmed = new Set(); // the units warmed up in this process
  /** `call()` timed, with the fallbacks, the aborts and their contexts (up to 200) it gave. */
  const measured = async (part, key, call) => {
    current = [part, key];
    const before = fallbacks();
    const ab = ABORTS.length;
    const pc = PANIC_CTX.length;
    const ex = ENGINE_EXCEPTIONS.length;
    const t0 = performance.now();
    const res = await call();
    const ms = performance.now() - t0;
    const aborts = ABORTS.slice(ab);
    const traps = PANIC_CTX.slice(pc).filter((c) => c);
    const exceptions = ENGINE_EXCEPTIONS.slice(ex);
    return {ms, res, fallbacks: fallbacks() - before, aborts: aborts.length, ...(traps.length ? {traps: traps.length,
      trapContexts: traps.slice(0, 500)} : {}), ...(exceptions.length ? {exceptions: exceptions.length,
      exceptionContexts: exceptions.slice(0, 500).map((x) => x.ctx)} : {}),
      ...(aborts.length ? {abortContexts: aborts.slice(0, 200).map((a) => a.ctx), abortHeapGb: gb(Math.max(...aborts.map(
        (a) => a.heap)))} : {})};
  };
  const strip = ({res, ...m}) => m;
  // transformation: every one-reactant reaction over the dataset (removeSaltsAndWater, the variables' defaults)
  for (const r of reactions.filter((r) => reactantCount(r.reactionSmarts) === 1)) {
    const opts = {removeSaltsAndWater: true, resolvedVariables: variablesOf(r)};
    out.transformation[r.id] = await rangedUnit(state, `transformation/${r.id}`, mols.length, async (a, b) => {
      if (!warmed.has(r.id)) {
        current = ['warm-up', r.id];
        await rx.runTransformationReaction(rdkitParam, r.reactionSmarts, warm, opts);
        warmed.add(r.id);
      }
      const list = a === 0 && b === mols.length ? mols : mols.slice(a, b);
      canonical.set(list, a === 0 && b === mols.length ? prepared.canonical : prepared.canonical.slice(a, b));
      const m = await measured('transformation', r.id, () => rx.runTransformationReaction(rdkitParam, r.reactionSmarts,
        list, opts));
      canonical.delete(list);
      // products: the product where it reacted, '=' where it gave back its input (no product, an error)
      return {...strip(m), errors: m.res.errors, noProducts: m.res.noProducts, empty: m.res.emptyMolecules,
        items: m.res.products.map((p, j) => p === (prepared.cleaned[a + j] ?? list[j]) ? '=' : p)};
    });
  }
  // two-component, pairwise: reactant lists from the dataset matching each reactant template
  for (const [id, pair] of Object.entries(prepared.twoComponent)) {
    const r = reactions.find((x) => x.id === id);
    out.twoComponent[id] = await rangedUnit(state, `twoComponent/${id}`, pair.r1.length, async (a, b) => {
      if (!warmed.has(id)) {
        const w1 = pair.r1.slice(0, 5);
        const w2 = pair.r2.slice(0, 5);
        canonical.set(w1, pair.c1.slice(0, 5));
        canonical.set(w2, pair.c2.slice(0, 5));
        current = ['warm-up', id];
        await rx.runTwoComponentReaction(rdkitParam, r.reactionSmarts, w1, w2, {mode: 'pairwise'});
        warmed.add(id);
      }
      const all = a === 0 && b === pair.r1.length;
      const [l1, l2] = all ? [pair.r1, pair.r2] : [pair.r1.slice(a, b), pair.r2.slice(a, b)];
      canonical.set(l1, all ? pair.c1 : pair.c1.slice(a, b));
      canonical.set(l2, all ? pair.c2 : pair.c2.slice(a, b));
      const m = await measured('twoComponent', id, () => rx.runTwoComponentReaction(rdkitParam, r.reactionSmarts, l1, l2,
        {mode: 'pairwise'}));
      return {...strip(m), errors: m.res.errors, noProducts: m.res.noProducts, items: m.res.products};
    });
  }
  // the demo's matrix: amide coupling of files/carb_acids.csv x files/amines.csv
  {
    const key = 'amide-coupling (demo matrix)';
    out.twoComponent[key] = await unitOnce(state, `twoComponent/${key}`, async () => {
      const acids = chemFile('carb_acids.csv').map((r) => r.canonical_smiles);
      const amines = chemFile('amines.csv').map((r) => r.canonical_smiles);
      canonical.set(acids, prepared.demo.acids);
      canonical.set(amines, prepared.demo.amines);
      const r = reactions.find((x) => x.id === 'amide-coupling');
      const m = await measured('twoComponent', key, () => rx.runTwoComponentReaction(rdkitParam, r.reactionSmarts, acids,
        amines, {mode: 'matrix'}));
      return {...strip(m), errors: m.res.errors, noProducts: m.res.noProducts, items: m.res.products};
    }) ?? {crashed: state.crash(`twoComponent/${key}`), items: []};
  }
  // Deprotect over the dataset: Chem's default fragment (Fmoc) and Boc
  for (const [name, fragment] of Object.entries(DEPROTECT_FRAGMENTS)) {
    out.deprotect[name] = await rangedUnit(state, `deprotect/${name}`, mols.length, async (a, b) => {
      if (!warmed.has(name)) {
        current = ['warm-up', name];
        deprotect(fragment, warm);
        warmed.add(name);
      }
      const list = a === 0 && b === mols.length ? mols : mols.slice(a, b);
      const m = await measured('deprotect', name, async () => deprotect(fragment, list));
      return {...strip(m), items: m.res.map((p, j) => p === list[j] ? '=' : p)};
    });
  }
  // Reaction Enumeration: the bundled demo (reactions.csv x bb.csv, the default config), the bundled config.yaml, and
  // the amino-acid files
  for (const variant of ['default config', 'config.yaml', 'amino acids']) {
    out.enumeration[variant] = await unitOnce(state, `enumeration/${variant}`, async () => {
      current = ['warm-up', variant];
      await enumeration(variant);
      const m = await measured('enumeration', variant, () => enumeration(variant));
      return {...strip(m), items: m.res.map((row) => JSON.stringify(row))};
    }) ?? {crashed: state.crash(`enumeration/${variant}`), items: []};
  }
  // across the processes: the fallbacks and Crux's panics in this process only (as before), the aborts of every unit
  out.fallbackCalls = tally(globalThis.__engine.fallbacks.map(([call]) => call));
  const units = Object.values(out).flatMap((p) => p && typeof p === 'object' && !Array.isArray(p) ? Object.values(p) :
    []).filter((u) => u && typeof u === 'object' && Array.isArray(u.items));
  // Crux's panics: every unit's (kept across processes), this process's where a unit has none recorded
  const unitTraps = units.flatMap((u) => u.trapContexts ?? []);
  out.trapContexts = unitTraps.length ? unitTraps : PANIC_CTX.filter((c) => c);
  out.trapCount = units.reduce((n, u) => n + (u.traps ?? 0), 0);
  out.aborts = {count: units.reduce((n, u) => n + (u.aborts ?? 0), 0),
    heapGb: Math.max(0, ...units.map((u) => u.abortHeapGb ?? 0)) || null};
  out.restarts = state.restarts();
  out.crashed = state.crashed;
  return out;
};

// the search: parity.mjs runs both of Chem's engines; for the queries where they differ, RDKit.js 2026.03.6 searches the
// same way as Chem's RDKit search (pattern fingerprint screen, then the query getQueryMolSafe builds matched against the
// canonical SMILES read unsanitized: parity.mjs rdkitSearch)
RUN.search = async (engine, rows) => {
  const module = await loadRdkit2026();
  const chemRdkit = await loadChemRdkit();
  const {getMolSafe, getQueryMolSafe} = await chem('utils/mol-creation_rdkit.ts');
  const mols = molecules(rows);
  const diffs = readJson(inputs(`search-diffs-${rows}.json`));
  const parityText = readFileSync(join(HERE, 'parity.mjs'), 'utf8');
  const smilesQueries = eval(parityText.match(/const SMILES_QUERIES = (\[[\s\S]*?\]);/)[1]);
  const setFile = (set) => ['queries.txt', 'queries.smi', `${set}.smarts`].map((f) => join(CRUX_BENCH, 'queries', set, f))
    .find((f) => existsSync(f));
  const resolve = ({set, label}) => {
    if (set === 'smiles')
      return smilesQueries.find((q) => q.slice(0, 36).trim() === label);
    if (set === 'molblock') {
      const mol = chemRdkit.get_mol(smilesQueries.find((q) => q.slice(0, 36).trim() === label));
      try {
        return mol.get_molblock();
      } finally {
        mol.delete();
      }
    }
    if (set === 'self')
      return mols[+label.replace('row ', '')];
    return readFileSync(setFile(set), 'utf8').split(/\r?\n/).map((l) => l.split(/[\t ]/)[0])
      .find((q) => q && q.slice(0, 36).trim() === label);
  };
  const unsanitized = JSON.stringify({sanitize: false, removeHs: false, assignStereo: false});
  const t0 = performance.now();
  const fps = new Array(mols.length).fill(null);
  const targets = new Array(mols.length).fill(null);
  for (let i = 0; i < mols.length; i++) {
    const mol = mols[i] ? getMolSafe(mols[i], {}, module).mol : null;
    if (!mol)
      continue;
    try {
      fps[i] = mol.get_pattern_fp_as_uint8array();
      targets[i] = module.get_mol(mol.get_smiles(), unsanitized);
    } catch (e) {
    } finally {
      mol.delete();
    }
  }
  const prepareMs = performance.now() - t0;
  const items = diffs.map((d) => {
    const q = resolve(d);
    const qmol = q ? getQueryMolSafe(q, '', module) : null;
    if (!qmol)
      return {...d, query: q ?? null, rdkit2026: null};
    const qfp = qmol.get_pattern_fp_as_uint8array();
    let hits = 0;
    outer:
    for (let i = 0; i < targets.length; i++) {
      if (!fps[i] || !targets[i])
        continue;
      for (let j = 0; j < qfp.length; j++) {
        if ((fps[i][j] & qfp[j]) !== qfp[j])
          continue outer;
      }
      if (targets[i].get_substruct_match(qmol) !== '{}')
        hits++;
    }
    qmol.delete();
    return {...d, query: q, rdkit2026: hits};
  });
  targets.forEach((t) => t?.delete());
  return {n: diffs.length, prepareMs, items};
};

// the search types beside Contains (Exact, Stereo agnostic, Included in, Not included in, Similar), one hit set per
// query and type. RDKit mode: Chem's RDKit search as its pool runs it on a column whose rows each worker caches (a
// worker caches its first 10,000 molecules, as getFingerprints reads them: a 100K column over 16 workers, an 18-thread machine's):
// RdKitServiceWorker.getFingerprints (pattern fingerprints, canonical SMILES) per 10,000 rows, RdKitService's pattern
// screen (filterByPatternFps), RdKitServiceWorker.searchSubstructure on the rows it passes; Similar: getFingerprints of
// the type and Chem's Tanimoto. Crux mode: crux/jobs/search-types.ts over the column (one segment), RDKit's answer for
// the rows Crux does not read (the same code on Chem's RDKit), a query Crux does not read on RDKit (counted); Similar:
// the RDKit worker's getFingerprints on the Crux module (engine/worker-jobs.ts), the same Tanimoto.
// The queries: SMILES_QUERIES (parity.mjs) as SMILES and as Chem's MOL blocks (the sketcher's), 200 of the column's
// rows as given and as MOL blocks (with their wedges), 20 SMARTS of crux-bench's bdb-1k; Not included in on a fifth of
// them, Stereo agnostic on the rows and a fifth of the others (RDKit's search scans (nearly) every row for them); Similar for 100 rows (Morgan at 0.8, 0.5, 0.3; the other types at 0.5
// for 40).
const SEARCH_FILES = ['smiles.csv', 'sar-small.csv', 'mol1K.csv', 'molblocks.csv', 'stereo-agnostic-test.csv',
  'v3000_sample.csv', 'wedging.csv', 'chem_standards.csv', 'test.csv', 'activity_cliffs.csv', 'r-groups.csv'];
/** The molecules of Chem's test datasets, each file's molecule column (searchtypesfiles). */
function chemFileMolecules() {
  return SEARCH_FILES.flatMap((f) => {
    const table = chemFile(f);
    const col = ['smiles', 'canonical_smiles', 'molecule'].find((c) => c in (table[0] ?? {}));
    return table.map((r) => r[col] ?? '');
  });
}
async function searchTypeQueries(rdkit, mols) {
  const parityText = readFileSync(join(HERE, 'parity.mjs'), 'utf8');
  const smilesQueries = eval(parityText.match(/const SMILES_QUERIES = (\[[\s\S]*?\]);/)[1]);
  const molblock = (s) => {
    const mol = s ? rdkit.get_mol(s) : null;
    try {
      return mol ? mol.get_molblock() : null;
    } finally {
      mol?.delete();
    }
  };
  const queries = smilesQueries.map((q) => ({set: 'smiles', label: q, q}));
  queries.push(...smilesQueries.map((q) => ({set: 'molblock', label: q, q: molblock(q)})).filter((q) => q.q));
  const self = [];
  for (let k = 0; k < 200 && mols.length; k++) {
    const i = Math.floor((k + 0.5) * mols.length / Math.min(200, mols.length));
    if (mols[i] && !self.includes(i))
      self.push(i);
  }
  queries.push(...self.map((i) => ({set: 'self', label: `row ${i}`, q: mols[i]})));
  queries.push(...self.map((i) => ({set: 'selfmol', label: `row ${i}`, q: mols[i].includes('M  END') ? mols[i] :
    molblock(mols[i])})).filter((q) => q.q));
  const smartsFile = join(CRUX_BENCH, 'queries', 'bdb-1k', 'queries.txt');
  if (existsSync(smartsFile)) {
    queries.push(...readFileSync(smartsFile, 'utf8').split(/\r?\n/).map((l) => l.split(/[\t ]/)[0]).filter((q) => q)
      .slice(0, 20).map((q) => ({set: 'smarts', label: q, q})));
  }
  return queries;
}
/** A hit set as compared: its size and hash; the rows themselves (or the rows of a Not type outside it) when few. */
function hitSetItem(hits, n, not) {
  const sorted = Int32Array.from(hits).sort();
  const key = `${sorted.length}:${hash(Array.from(sorted).join(','))}`;
  const outside = not ? (() => {
    const inSet = new Uint8Array(n);
    for (const h of sorted)
      inSet[h] = 1;
    const rest = [];
    for (let i = 0; i < n && rest.length <= 3000; i++) {
      if (!inSet[i])
        rest.push(i);
    }
    return rest;
  })() : null;
  return {key, rows: !not && sorted.length <= 3000 ? Array.from(sorted) : null, outside};
}
RUN.searchtypes = async (engine, rows, source = 'dataset') => {
  const {module, rdkit} = await engineModules(engine);
  const {getMolSafe, getQueryMolSafe} = await chem('utils/mol-creation_rdkit.ts');
  const {SubstructureSearchType: T, AVAILABLE_FPS} = await chem('constants.ts');
  const {RdKitServiceWorker} = await chem('rdkit-service/rdkit-service-worker.ts');
  const {rdKitFingerprintToBitArray} = await chem('utils/chem-common.ts');
  const {tanimotoSimilarity} = await import('@datagrok-libraries/ml/src/distance-metrics-methods');
  const mols = source === 'files' ? chemFileMolecules() : molecules(rows);
  const n = mols.length;
  const queries = await searchTypeQueries(rdkit, mols);
  const types = [T.EXACT_MATCH, T.STEREO_AGNOSTIC, T.INCLUDED_IN, T.NOT_INCLUDED_IN];
  const t0 = performance.now();
  let search;
  const counts = {queries: queries.length, rowsCruxDoesNotRead: 0, queriesOnRdkit: 0};
  if (engine === 'crux') {
    const {SearchTypeKeys, searchTypeHits} = await chem('crux/jobs/search-types.ts');
    const {cruxReadsCxExtension} = await chem('crux/crux-service.ts');
    const sent = mols.map((m) => cruxReadsCxExtension(m) ? m : '');
    const keys = new SearchTypeKeys(sent);
    keys.patternFps(module);
    const unread = [...keys.unread(module), ...sent.flatMap((m, i) => !m && mols[i] ? [i] : [])].sort((a, b) => a - b);
    counts.rowsCruxDoesNotRead = unread.length;
    const rdkitKeys = new SearchTypeKeys(unread.map((i) => mols[i]));
    search = (q, type) => {
      let cruxQuery = null;
      try {
        cruxQuery = cruxReadsCxExtension(q) ? getQueryMolSafe(q, '', module) : null;
      } catch (e) {}
      if (!cruxQuery)
        return null; // RDKit's search (RDKit mode's hits)
      cruxQuery.delete();
      return [...searchTypeHits(module, keys, q, '', type),
        ...searchTypeHits(rdkit, rdkitKeys, q, '', type).map((k) => unread[k])];
    };
  } else {
    // the pool's workers, each with the rows it caches (at most MAX_MOL_CACHE_SIZE)
    const per = 10000;
    const workers = [];
    for (let start = 0; start < n; start += per) {
      const worker = new RdKitServiceWorker(module, '');
      const fps = await worker.getFingerprints('Pattern', mols.slice(start, start + per), true);
      workers.push({worker, start, fps});
    }
    search = async (q, type) => {
      const qmol = getQueryMolSafe(q, '', module);
      if (!qmol)
        return null;
      const qfp = qmol.get_pattern_fp_as_uint8array();
      qmol.delete();
      const hits = [];
      for (const {worker, start, fps} of workers) {
        // RdKitService.searchSubstructureWithFps: the screen (filterByPatternFps), then the worker's search
        const not = type === T.NOT_INCLUDED_IN;
        const passed = [];
        for (let i = 0; i < fps.fps.length; i++) {
          if (not) {
            passed.push(i);
            continue;
          }
          const fp = fps.fps[i];
          if (!fp)
            continue;
          let ok = true;
          for (let j = 0; j < 256 && ok; j++)
            ok = (fp[j] & qfp[j]) === (type === T.INCLUDED_IN ? fp[j] : qfp[j]);
          if (ok)
            passed.push(i);
        }
        const buf = await worker.searchSubstructure(q, '', passed.map((i) => fps.smiles[i]), type);
        const bits = new Uint32Array(buf);
        passed.forEach((i, k) => {
          if (bits[k >>> 5] & (1 << (k & 31)))
            hits.push(start + i);
        });
      }
      return hits;
    };
  }
  const prepareMs = performance.now() - t0;
  const items = [];
  const searchStart = performance.now();
  for (const type of types) {
    for (const [k, {set, label, q}] of queries.entries()) {
      // Not included in on a fifth of the queries, Stereo agnostic on the rows and a fifth of the others: RDKit's search
      // scans (nearly) every row for them
      if ((type === T.NOT_INCLUDED_IN || (type === T.STEREO_AGNOSTIC && !set.startsWith('self'))) && k % 5 !== 0)
        continue;
      const t = performance.now();
      const hits = await search(q, type);
      if (hits === null)
        counts.queriesOnRdkit++;
      items.push({type, set, label, ms: performance.now() - t, onRdkit: hits === null,
        ...(hits === null ? {key: null} : hitSetItem(hits, n, type === T.NOT_INCLUDED_IN))});
    }
  }
  const searchMs = performance.now() - searchStart;
  // Similar: the column's fingerprints (the RDKit pool's getFingerprints; Crux mode the same code on Crux, RDKit's for
  // the rows Crux does not read), the query's as chemGetFingerprint computes it, Chem's Tanimoto
  const simStart = performance.now();
  const fpWorker = new RdKitServiceWorker(module, '');
  fpWorker._cacheCounter = 10000; // no cache (crux.worker.ts rdkitShapedWorker; the RDKit pool's cache adds nothing)
  const rdkitFpWorker = new RdKitServiceWorker(rdkit, '');
  rdkitFpWorker._cacheCounter = 10000;
  const selfRows = queries.filter((q) => q.set === 'self');
  // chem-searches.ts getRDKitFpAsUint8Array (that module imports the package): the query's fingerprint
  const getRDKitFpAsUint8Array = (mol, type) => {
    if (type === 'Morgan') {
      mol.remove_hs_in_place();
      return mol.get_morgan_fp_as_uint8array(JSON.stringify({radius: 2, nBits: 2048}));
    }
    return {Pattern: () => mol.get_pattern_fp_as_uint8array(), AtomPair: () => mol.get_atom_pair_fp_as_uint8array(),
      MACCS: () => mol.get_maccs_fp_as_uint8array(), RDKit: () => mol.get_rdkit_fp_as_uint8array(),
      TopologicalTorsion: () => mol.get_topological_torsion_fp_as_uint8array()}[type]();
  };
  for (const fp of AVAILABLE_FPS) {
    const cutoffs = fp === 'Morgan' ? [0.8, 0.5, 0.3] : [0.5];
    const qs = fp === 'Morgan' ? selfRows.slice(0, 100) : selfRows.slice(0, 40);
    const res = await fpWorker.getFingerprints(fp, mols, false);
    if (engine === 'crux') { // the molecules Crux cannot read: the RDKit pool's fingerprints (engine/worker-jobs.ts)
      const unreadRows = mols.flatMap((m, i) => {
        if (!m || res.fps[i])
          return [];
        const mol = getMolSafe(m, {}, module).mol;
        mol?.delete();
        return mol ? [] : [i];
      });
      const patch = await rdkitFpWorker.getFingerprints(fp, unreadRows.map((i) => mols[i]), false);
      unreadRows.forEach((row, k) => res.fps[row] = patch.fps[k]);
    }
    const bits = res.fps.map((f) => f ? rdKitFingerprintToBitArray(f) : null);
    for (const {label, q} of qs) {
      // chemGetFingerprint: the engine's, RDKit's where it reads none
      const queryFp = (m) => {
        const mol = getMolSafe(q, {}, m).mol;
        try {
          return mol ? getRDKitFpAsUint8Array(mol, fp) : null;
        } finally {
          mol?.delete();
        }
      };
      const qfp = queryFp(module) ?? queryFp(rdkit);
      const qbits = rdKitFingerprintToBitArray(qfp);
      const sims = bits.map((b) => b ? tanimotoSimilarity(b, qbits) : -1);
      for (const cutoff of cutoffs) {
        const hits = sims.flatMap((x, i) => x >= cutoff ? [i] : []);
        items.push({type: `${T.IS_SIMILAR} ${fp} ${cutoff}`, set: 'self', label, ...hitSetItem(hits, n, false)});
      }
    }
  }
  return {n: items.length, rows: n, prepareMs, searchMs, similarMs: performance.now() - simStart, counts,
    fallbacks: [], items};
};
RUN.searchtypesfiles = (engine, rows) => RUN.searchtypes(engine, rows, 'files');

// MMPA fragmentation: RdKitServiceWorker.mmpGetFragments; Crux mode crux-analysis-jobs.ts's mmpGetFragments (+ the RDKit
// pool for its `unread` rows, engine/analysis-jobs.ts), in chunks of 1000
async function analysisSetup(engine) {
  const {module, rdkit, crux} = await engineModules(engine);
  const {RdKitServiceWorker} = await chem('rdkit-service/rdkit-service-worker.ts');
  const {ANALYSIS_JOBS} = await chem('crux/crux-analysis-jobs.ts');
  const cw = engine === 'crux' ? await chem('crux/crux.worker.ts') : null;
  const rdkitWorker = new RdKitServiceWorker(rdkit, '');
  const worker = engine === 'crux' ? cw.rdkitShapedWorker() : new RdKitServiceWorker(module, '');
  // a Crux worker job's context (crux-api.ts WorkerJobContext)
  const ctx = {crux: module, module, rdkitWorker: () => worker, transfer: []};
  return {module, rdkit, crux, ANALYSIS_JOBS, cw, rdkitWorker, worker, ctx};
}
ONE.mmpa = async (engine) => {
  const {ANALYSIS_JOBS, rdkitWorker, worker, ctx} = await analysisSetup(engine);
  return (chunk) => {
    if (engine !== 'crux')
      return {...worker.mmpGetFragments(chunk), fallback: []};
    const {frags, smiles, unread} = ANALYSIS_JOBS.mmpGetFragments(ctx, chunk);
    if (unread.length) {
      const patch = rdkitWorker.mmpGetFragments(unread.map((i) => chunk[i]));
      unread.forEach((row, k) => {
        frags[row] = patch.frags[k];
        smiles[row] = patch.smiles[k];
      });
    }
    return {frags, smiles, fallback: unread};
  };
};
RUN.mmpa = async (engine, rows) => {
  const frag = await ONE.mmpa(engine);
  const mols = molecules(rows);
  frag(mols.slice(0, WARMUP));
  const items = [];
  const fallback = [];
  const chunksPanicked = [];
  const t0 = performance.now();
  for (let start = 0; start < mols.length; start += 1000) {
    const before = panicCount();
    const {frags, smiles, fallback: fb} = frag(mols.slice(start, start + 1000));
    if (panicCount() > before)
      chunksPanicked.push(start);
    for (let j = 0; j < smiles.length; j++)
      items.push(JSON.stringify([smiles[j] ?? null, frags[j] ?? null]));
    fallback.push(...fb.map((j) => start + j));
  }
  const ms = performance.now() - t0;
  const traps = [];
  for (const start of chunksPanicked)
    traps.push(...await trapsIn(mols.slice(start, start + 1000), start, (one) => frag(one)));
  // '["",[]]': mmpGetFragments' catch (get_mmpa_frags gave no lists: no single bond to cut, or a refusal)
  const status = items.map((s) => s.startsWith('[null') ? 'not read' : s === '["",[]]' ? 'no fragments (no cut)' :
    s.endsWith(',[]]') || s.endsWith(',null]') ? 'read, no pairs' : 'fragmented');
  return {n: mols.length, ms, counts: tally(status), fallbacks: fallback, traps, items};
};

// R-group decomposition: RdKitServiceWorker.rGroupAnalysis with the R-group dialog's options and core (a MOL block of
// the core as RDKit writes it from its SMARTS reading); Crux mode crux-analysis-jobs.ts's rGroupAnalysis, the RDKit pool
// for a series it cannot do (engine/analysis-jobs.ts rGroupsByEngine)
const RGD_OPTIONS = JSON.stringify({matchingStrategy: 'Greedy', includeTargetMolInResults: true,
  onlyMatchAtRGroups: false});
const rgdText = (r) => r === null ? 'null' : JSON.stringify({colNames: r.colNames, smiles: r.smiles,
  atoms: r.atomsToHighLight.map((c) => Array.from(c, (a) => a ? Array.from(a) : null)),
  bonds: r.bondsToHighLight.map((c) => Array.from(c, (b) => b ? Array.from(b) : null))});
ONE.rgroups = async (engine) => {
  const {ANALYSIS_JOBS, rdkitWorker, worker, ctx} = await analysisSetup(engine);
  return async (mols, core, coreIsQMol) => {
    if (engine !== 'crux')
      return {result: await worker.rGroupAnalysis(mols, core, coreIsQMol, RGD_OPTIONS), fallback: false};
    try {
      const result = await ANALYSIS_JOBS.rGroupAnalysis(ctx, mols, core, coreIsQMol, RGD_OPTIONS);
      if (result)
        return {result, fallback: false};
    } catch (e) {}
    return {result: await rdkitWorker.rGroupAnalysis(mols, core, coreIsQMol, RGD_OPTIONS), fallback: true};
  };
};
RUN.rgroups = async (engine, rows) => {
  const decompose = await ONE.rgroups(engine);
  const mols = molecules(rows);
  const series = readJson(inputs(`rgroup-series-${rows}.json`));
  await decompose(series[0].rows.slice(0, 20).map((i) => mols[i]), series[0].coreMolBlock, series[0].coreIsQMol);
  const items = [];
  const fallback = [];
  const status = [];
  const times = [];
  const t0 = performance.now();
  for (let k = 0; k < series.length; k++) {
    const s = series[k];
    const t1 = performance.now();
    const before = panicCount();
    try {
      const {result, fallback: fb} = await decompose(s.rows.map((i) => mols[i]), s.coreMolBlock, s.coreIsQMol);
      items.push(rgdText(result));
      status.push(result.colNames.length ? 'decomposed' : 'no columns');
      if (fb)
        fallback.push(k);
    } catch (e) {
      items.push(`THROWS ${e?.message ?? e}`);
      status.push('threw');
    }
    times.push(performance.now() - t1);
    if (panicCount() > before)
      TRAPS.push(k);
  }
  return {n: series.length, molecules: series.reduce((n, s) => n + s.rows.length, 0), ms: performance.now() - t0,
    times, counts: tally(status), fallbacks: fallback, traps: [...TRAPS], items};
};

// FMCS: Cluster MCS (performClusterMCS: exact atoms and bonds; the RDKit pool gives '' for a cluster past 45 s, the
// Crux pool FMCS's Timeout of 45 s) and its SMILES (utils/mol-creation_rdkit.ts clusterMcsSmiles)
ONE.fmcs = async (engine) => {
  const {module, rdkit, ANALYSIS_JOBS, rdkitWorker, worker, ctx} = await analysisSetup(engine);
  const {clusterMcsSmiles} = await chem('utils/mol-creation_rdkit.ts');
  const {CLUSTER_MCS_TIMEOUT_S} = await chem('crux/crux-api.ts');
  return async (mols) => {
    let smarts;
    let fallback = false;
    if (engine !== 'crux')
      smarts = worker.mostCommonStructure(mols, true, true, CLUSTER_MCS_TIMEOUT_S);
    else {
      [smarts] = await ANALYSIS_JOBS.clusterMCS(ctx, [mols], true, true);
      if (smarts === null) {
        smarts = rdkitWorker.mostCommonStructure(mols, true, true, CLUSTER_MCS_TIMEOUT_S);
        fallback = true;
      }
    }
    let [smiles] = clusterMcsSmiles([smarts], [mols], module);
    if (engine === 'crux' && smarts && !smiles) // analysis-jobs.ts: an MCS the Crux module does not write: RDKit's
      [smiles] = clusterMcsSmiles([smarts], [mols], rdkit);
    return {smarts, smiles, fallback};
  };
};
RUN.fmcs = async (engine, rows) => {
  const mcs = await ONE.fmcs(engine);
  const mols = molecules(rows);
  const clusters = readJson(inputs(`fmcs-clusters-${rows}.json`));
  await mcs(clusters[clusters.length - 1].rows.map((i) => mols[i]));
  const items = [];
  const fallback = [];
  const times = [];
  const t0 = performance.now();
  for (let k = 0; k < clusters.length; k++) {
    const t1 = performance.now();
    const before = panicCount();
    const {smarts, smiles, fallback: fb} = await mcs(clusters[k].rows.map((i) => mols[i]));
    times.push(performance.now() - t1);
    if (panicCount() > before)
      TRAPS.push(k);
    items.push(JSON.stringify([smarts, smiles]));
    if (fb)
      fallback.push(k);
  }
  return {n: clusters.length, molecules: clusters.reduce((n, c) => n + c.rows.length, 0), ms: performance.now() - t0,
    times, counts: tally(items.map((s) => s.startsWith('[""') ? 'no MCS (timed out or none)' : 'MCS')),
    fallbacks: fallback, traps: [...TRAPS], items};
};

// Rascal MCES (an activity cliff's MCS in Crux mode): get_mces with crux-service.ts's options on molecules read by
// getMolSafe (the parity), and the whole Crux job (crux.worker.ts pairUncommonAtomsAndBonds: the read, get_mces and the
// uncommon atoms and bonds; the time Chem spends)
function mcesOptions() {
  const src = readFileSync(join(SRC, 'crux', 'crux-service.ts'), 'utf8');
  const num = (name) => Number(src.match(new RegExp(`const ${name} = ([\\d.]+);`))[1]);
  return {similarityThreshold: num('PAIR_MCS_SIMILARITY_THRESHOLD'),
    maxBondMatchPairs: num('PAIR_MCS_MAX_BOND_MATCH_PAIRS'), timeout: num('PAIR_MCS_TIMEOUT_S')};
}
const f64 = new Float64Array(1);
const u8 = new Uint8Array(f64.buffer);
const floatBits = (x) => {
  f64[0] = x;
  return Array.from(u8).reverse().map((b) => b.toString(16).padStart(2, '0')).join('');
};
RUN.mces = async (engine, rows) => {
  const {crux} = await engineModules('crux');
  const {getMolSafe} = await chem('utils/mol-creation_rdkit.ts');
  const {ANALYSIS_JOBS} = await chem('crux/crux-analysis-jobs.ts');
  const pairs = readJson(inputs(`mces-pairs-${rows}.json`));
  const details = JSON.stringify(mcesOptions());
  const parity = (p) => {
    const m1 = getMolSafe(p.smiles1, {}, crux).mol;
    const m2 = getMolSafe(p.smiles2, {}, crux).mol;
    try {
      if (!m1 || !m2)
        return {error: 'parse'};
      const t = performance.now();
      const res = JSON.parse(crux.get_mces(m1, m2, details));
      const s = (performance.now() - t) / 1000;
      return {s, results: res.map((r) => ({smarts: r.smartsString, atomMatches: r.atomMatches,
        bondMatches: r.bondMatches, similarity: floatBits(r.similarity), tier1Sim: floatBits(r.tier1Sim),
        tier2Sim: floatBits(r.tier2Sim), timedOut: r.timedOut, numFragments: r.numFragments,
        largestFragmentSize: r.largestFragmentSize}))};
    } catch (e) {
      return {error: String(e?.message ?? e)};
    } finally {
      m1?.delete();
      m2?.delete();
    }
  };
  const {ms, out, traps} = timeSync(pairs, parity, 5);
  const job = timeSync(pairs, (p) => {
    try {
      return ANALYSIS_JOBS.pairUncommonAtomsAndBonds({module: crux}, p.smiles1, p.smiles2, details) === null ? 'below threshold' : 'MCES';
    } catch (e) {
      return 'threw';
    }
  }, 5);
  return {n: pairs.length, ms, options: mcesOptions(), jobMs: job.ms, counts: tally(job.out), traps, items: out};
};
// FMCS of the same pairs: what RDKit mode gives an activity cliff (getMCS(pair, exact atoms, exact bonds) on the RDKit
// pool, no time limit there: 60 s here as a guard), and Crux's FMCS (Crux mode's fallback for a cliff)
// Pair by pair, each pair's outcome kept across processes (Resumable): a pair that throws (RDKit.js 2024.09 raised a
// C++ exception here, WebAssembly.Exception, with its heap at its 4 GB limit) is 'THREW', one the process dies on
// 'CRASHED', and a process whose heap has grown past HEAP_RESTART_BYTES gives way to a fresh one between two pairs.
RUN['mces-fmcs'] = async (engine, rows) => {
  const state = new Resumable('mces-fmcs', engine);
  const {worker} = await analysisSetup(engine);
  const pairs = readJson(inputs(`mces-pairs-${rows}.json`));
  let current = null;
  globalThis.__panicContext = () => ({pair: current});
  const mcs = (p) => worker.mostCommonStructure([p.smiles1, p.smiles2], true, true, 60);
  let warmed = false;
  const items = [];
  const times = [];
  const outcomes = [];
  // MCES_FMCS_PAIRS: only the first n pairs (the rest 'NOT RUN'), for a time budget
  // (OUT/mces-fmcs-pairs.json {engine: n}, when present, overrides it: a run already started can be given fewer)
  const limits = existsSync(join(OUT, 'mces-fmcs-pairs.json')) ? readJson(join(OUT, 'mces-fmcs-pairs.json')) : {};
  const limit = parseInt(String(limits[engine] ?? process.env.MCES_FMCS_PAIRS ?? pairs.length));
  for (let i = 0; i < pairs.length; i++) {
    if (i >= limit) {
      items.push('NOT RUN');
      times.push(null);
      outcomes.push(undefined);
      continue;
    }
    const r = await unitOnce(state, `pair ${i}`, () => {
      if (!warmed) { // the first 5 pairs once per process, untimed
        for (const p of pairs.slice(0, 5)) {
          try {
            mcs(p);
          } catch (_) {}
        }
        warmed = true;
      }
      current = i;
      const before = panicCount();
      const ab = ABORTS.length;
      const t = performance.now();
      let smarts;
      let threw = null;
      try {
        smarts = mcs(pairs[i]);
      } catch (e) {
        threw = `${e?.constructor?.name ?? 'Error'}: ${String(e?.message ?? e).slice(0, 120)}`;
      }
      return {ms: performance.now() - t, smarts: threw ? null : smarts, threw, trapped: panicCount() > before,
        aborts: ABORTS.length - ab, heapGb: gb(heapBytes())};
    });
    items.push(r ? r.threw ? 'THREW' : r.smarts : 'CRASHED');
    times.push(r?.ms ?? null);
    outcomes.push(r);
  }
  const done = outcomes.filter((r) => r && !r.threw);
  const threw = outcomes.flatMap((r, i) => r?.threw ? [{pair: i, error: r.threw, heapGb: r.heapGb}] : []);
  return {n: pairs.length, ms: done.reduce((t, r) => t + r.ms, 0), times: done.map((r) => r.ms), pairsRun: Math.min(limit,
    pairs.length), counts: tally(items.map((s) => s === 'THREW' ? 'threw' : s === 'CRASHED' ? 'the process died' :
      s === 'NOT RUN' ? 'not run (time budget)' : s ? 'MCS' : 'none or past 60 s')), traps: outcomes.flatMap((r, i) => r?.trapped ? [i] : []), items,
    threw, crashed: state.crashed, restarts: state.restarts(),
    aborts: {count: outcomes.reduce((n, r) => n + (r?.aborts ?? 0), 0)}};
};

// molzip (Chem:molzip, in both modes on Crux: crux-analysis-jobs.ts molzip, details {} = RDKit's defaults) vs RDKit.js
// 2026.03.6's molzip(a, b) of the same pieces and Python's Chem.molzip
ONE.molzip = async (engine) => {
  const {module, ANALYSIS_JOBS, ctx} = await analysisSetup(engine);
  const {zippedSmiles} = await chem('rdkit-service/rdkit-service-worker-substructure.ts');
  const {getMolSafe} = await chem('utils/mol-creation_rdkit.ts');
  if (engine === 'crux')
    return (chunk) => ANALYSIS_JOBS.molzip(ctx, chunk, '{}').smiles;
  return (chunk) => chunk.map((s) => {
    const [first, ...rest] = s.split('.');
    const a = getMolSafe(first, {}, module).mol;
    const b = getMolSafe(rest.join('.'), {}, module).mol;
    try {
      return a && b ? zippedSmiles(module.molzip(a, b, '{}'), module) : '';
    } catch (e) {
      return '';
    } finally {
      a?.delete();
      b?.delete();
    }
  });
};
RUN.molzip = async (engine, rows) => {
  const zip = await ONE.molzip(engine);
  const list = readJson(inputs(`molzip-inputs-${rows}.json`)).map((x) => x.smiles);
  zip(list.slice(0, WARMUP));
  const items = [];
  const chunksPanicked = [];
  const t0 = performance.now();
  for (let start = 0; start < list.length; start += 1000) {
    const before = panicCount();
    items.push(...zip(list.slice(start, start + 1000)));
    if (panicCount() > before)
      chunksPanicked.push(start);
  }
  const ms = performance.now() - t0;
  const traps = [];
  for (const start of chunksPanicked)
    traps.push(...await trapsIn(list.slice(start, start + 1000), start, (one) => zip(one)));
  return {n: list.length, ms, counts: tally(items.map((s) => s ? 'joined' : 'refused or unread')), traps, items};
};

// ---- the parent: inputs, the runs, the comparisons, the report -----------------------------------------------------

const LABEL = {
  rdkit: 'RDKit.js 2024.09.1pre (Chem, RDKit mode)', crux: 'Crux (Chem\'s vendored file, crux-core a8a136474)',
  rdkit2026: 'RDKit.js 2026.03.6', py: 'Python RDKit 2026.3.6', chempy: 'chem.py on Python RDKit 2026.3.6 (one process)',
  container: 'chem-chem container (RDKit 2026.03.5; 0.5 CPU and 1 GB on this stand)', ocl: 'OpenChemLib 8.21.0 (Chem\'s OCL worker)',
  sg: 'scripts/scaffold_tree.py on ScaffoldGraphDG 1.1.8 (Python RDKit 2026.3.6)',
};

/** A process, its peak memory measured (macOS /usr/bin/time -l), written to OUT/processes.log. */
function spawnMeasured(label, cmd, args, opts = {}) {
  const t0 = performance.now();
  const timed = process.platform === 'darwin' && existsSync('/usr/bin/time');
  const r = spawnSync(timed ? '/usr/bin/time' : cmd, timed ? ['-l', cmd, ...args] : args,
    {stdio: ['inherit', 'inherit', 'pipe'], maxBuffer: 1 << 26, ...opts});
  const err = r.stderr?.toString() ?? '';
  const rss = err.match(/(\d+)\s+maximum resident set size/);
  // time's own lines dropped (the counters after the timing line)
  const timeLines = /^\s*\d+(\.\d+)?\s+(real|user|sys|maximum resident|average|page|swaps|block|messages|signals|voluntary|involuntary|instructions retired|cycles elapsed|peak memory)/;
  process.stderr.write(err.split('\n').filter((l) => !timeLines.test(l)).join('\n').slice(-4000));
  const line = `${new Date().toISOString()} ${label}: exit ${r.status}${r.signal ? ` (${r.signal})` : ''} in ` +
    `${((performance.now() - t0) / 1000).toFixed(0)} s, peak ${rss ? (rss[1] / 1e9).toFixed(2) + ' GB' : '?'}`;
  console.log(`  ${line}`);
  mkdirSync(OUT, {recursive: true});
  writeFileSync(join(OUT, 'processes.log'), line + '\n', {flag: 'a'});
  return r;
}
/** One engine's run in a child process: started again where it asks for a fresh process (RESTART_EXIT) or dies in a
 * unit of work it marked (the unit recorded as crashed, the next process skips it); the processes in `raw/<fn>/<engine>
 * .processes.json`. */
function child(fn, engine, rows) {
  if (!wanted(engine))
    return;
  const dir = join(RAW, fn, `${engine}.state-${Date.now()}`);
  mkdirSync(dir, {recursive: true});
  const processes = [];
  for (let k = 0; k <= MAX_RESTARTS; k++) {
    const t0 = performance.now();
    const r = spawnMeasured(`${fn}/${engine}${k ? ` (process ${k + 1})` : ''}`, process.execPath,
      ['--max-old-space-size=5500', fileURLToPath(import.meta.url), '--child', fn, engine, String(rows)],
      {env: {...process.env, OUT, RESUME_DIR: dir}});
    const current = existsSync(join(dir, 'current.json')) ? readJson(join(dir, 'current.json')) : {key: null};
    processes.push({exit: r.status, signal: r.signal, s: +((performance.now() - t0) / 1000).toFixed(1),
      in: current.key ?? null});
    writeJson(join(RAW, fn, `${engine}.processes.json`), processes);
    if (r.status === 0)
      return;
    if (r.status === RESTART_EXIT)
      continue;
    if (!current.key)
      throw new Error(`${fn}/${engine} failed (exit ${r.status}${r.signal ? `, ${r.signal}` : ''})`);
    const crashed = existsSync(join(dir, 'crashed.json')) ? readJson(join(dir, 'crashed.json')) : [];
    const err = r.stderr?.toString() ?? '';
    const hung = existsSync(join(dir, 'hung.json')) ? readJson(join(dir, 'hung.json')) : null;
    if (hung)
      writeJson(join(dir, 'hung.json'), null);
    crashed.push({key: current.key, exit: r.status, signal: r.signal, heapGb: current.heap ? gb(current.heap) : null,
      abortsBefore: current.aborts ?? null, error: hung ? `hung: not done within ${hung.seconds} s of CPU (killed by ` +
        'the watchdog)' : err.split('\n').filter((l) => /Error|Exception|Abort|unreachable|panick/.test(l)).slice(0, 4)
        .join(' | ').slice(0, 400) || 'killed (no message)'});
    writeJson(join(dir, 'crashed.json'), crashed);
    console.log(`  ${fn}/${engine}: the process died in ${current.key} (exit ${r.status}): recorded, run again past it`);
  }
  throw new Error(`${fn}/${engine}: more than ${MAX_RESTARTS} restarts`);
}
function python(args, {py = PYTHON, timeoutMs} = {}) {
  if (!wanted({descriptors: 'chempy', scaffoldtree: 'sg'}[args[0]] ?? 'py'))
    return true;
  const r = spawnMeasured(`python ${args[0]}`, py, ['-I', PY_HELPER, ...args], {timeout: timeoutMs,
    killSignal: 'SIGKILL'});
  return r.status === 0;
}

/** The data subsets (compare100k.py subsets) and each function's prepared inputs, made once per row count. */
async function prepareInputs(rows) {
  mkdirSync(inputs(''), {recursive: true});
  const mols = molecules(rows);
  const molFile = inputs(`molecules-${rows}.json`);
  if (!existsSync(molFile))
    writeJson(molFile, mols);
  const subsetDir = inputs(`subsets-${rows}`);
  if (!existsSync(join(subsetDir, 'molzip-inputs.json')))
    python(['subsets', molFile, subsetDir]);
  const rdkit = await loadChemRdkit();
  const {getMolSafe} = await chem('utils/mol-creation_rdkit.ts');
  const {RdKitServiceWorker} = await chem('rdkit-service/rdkit-service-worker.ts');
  if (!existsSync(inputs(`rgroup-series-${rows}.json`))) {
    // the dialog's core: convertMolNotation(core, Smarts, MolBlock) on Chem's RDKit (getMolSafe, get_molblock)
    writeJson(inputs(`rgroup-series-${rows}.json`), readJson(join(subsetDir, 'rgroup-series.json')).map((s) => {
      const mol = getMolSafe(s.core, {}, rdkit).mol;
      const block = mol.get_molblock();
      mol.delete();
      const coreIsQMol = block.includes('M  ALS') || block.includes('M  RAD');
      return {...s, coreMolBlock: coreIsQMol ? s.core : block, coreIsQMol};
    }));
  }
  if (!existsSync(inputs(`fmcs-clusters-${rows}.json`)))
    writeJson(inputs(`fmcs-clusters-${rows}.json`), readJson(join(subsetDir, 'fmcs-clusters.json')));
  if (!existsSync(inputs(`mces-pairs-${rows}.json`))) {
    writeJson(inputs(`mces-pairs-${rows}.json`), readJson(join(subsetDir, 'mces-pairs.json'))
      .map((p) => ({...p, smiles1: mols[p.a], smiles2: mols[p.b]})));
  }
  if (!existsSync(inputs(`molzip-inputs-${rows}.json`)))
    writeJson(inputs(`molzip-inputs-${rows}.json`), readJson(join(subsetDir, 'molzip-inputs.json')));
  if (!existsSync(inputs(`reaction-inputs-${rows}.json`))) {
    // RdKitService.convertMolNotation(mols, Smiles) on Chem's RDKit worker (both modes), then the salts removed
    const worker = new RdKitServiceWorker(rdkit, '');
    const t0 = performance.now();
    const canonical = await worker.convertMolNotation(mols, 'smiles');
    const canonicalMs = performance.now() - t0;
    const {removeWaterAndSaltsSingle} = await chem('utils/reactions/reactions.ts');
    const cleaned = canonical.map((s) => removeWaterAndSaltsSingle(s));
    // two-component reactant lists: the first 1000 molecules matching each reactant template (Chem's RDKit)
    const {splitSmartsByReactants} = await chem('utils/reaction-enumeration/enumerate.ts');
    const twoComponent = {};
    const read = cleaned.map((s) => s ? rdkit.get_mol(s) : null);
    for (const r of DEFAULT_REACTIONS().filter((r) => reactantCount(r.reactionSmarts) === 2)) {
      const lists = splitSmartsByReactants(r.reactionSmarts.split('>>')[0]).map((smarts) => {
        const q = rdkit.get_qmol(smarts);
        const hits = [];
        for (let i = 0; i < read.length && hits.length < 1000; i++) {
          if (read[i] && read[i].get_substruct_match(q) !== '{}')
            hits.push(i);
        }
        q.delete();
        return hits;
      });
      if (lists.every((l) => l.length)) {
        const len = Math.min(...lists.map((l) => l.length));
        const [r1, r2] = lists.map((l) => l.slice(0, len).map((i) => mols[i]));
        twoComponent[r.id] = {r1, r2, c1: await worker.convertMolNotation(r1, 'smiles'),
          c2: await worker.convertMolNotation(r2, 'smiles'), rows1: lists[0].slice(0, len), rows2: lists[1].slice(0, len)};
      }
    }
    read.forEach((m) => m?.delete());
    const demo = {acids: await worker.convertMolNotation(chemFile('carb_acids.csv').map((r) => r.canonical_smiles), 'smiles'),
      amines: await worker.convertMolNotation(chemFile('amines.csv').map((r) => r.canonical_smiles), 'smiles')};
    writeJson(inputs(`reaction-inputs-${rows}.json`), {molecules: mols, canonical, cleaned, canonicalMs, twoComponent,
      demo});
  }
}

// -- comparison helpers

const pct = (k, n) => n ? `${(100 * k / n).toFixed(3)} %` : '-';
const us = (ms, n) => n ? 1000 * ms / n : null;
// stereo marks out: @/@@, / and \\, and the brackets they leave ([C@@H] -> C)
const stripStereo = (s) => s.replace(/@+/g, '').replace(/[/\\]/g, '').replace(/\[([CNSP])H?\]/g, '$1');
const clip = (s, k = 160) => s === null || s === undefined ? String(s) : s.length > k ? s.slice(0, k) + '…' : s;
function firstDifference(a, b) {
  let k = 0;
  while (k < a.length && k < b.length && a[k] === b[k])
    k++;
  return `at char ${k}: …${a.slice(Math.max(0, k - 30), k + 50)}… vs …${b.slice(Math.max(0, k - 30), k + 50)}…`;
}

/** Crux's outputs against RDKit mode's and the references', the differences grouped by cause (rows for examples). */
const TRAP_CAUSE = 'Crux panicked (a trap)';
function compareItems({crux, base, refs = {}, fallback = new Set(), traps = new Set(), n, classify,
  refName = 'rdkit2026', override = null, forced = null}) {
  const res = {n, identical: base ? 0 : null, identicalRefs: Object.fromEntries(Object.keys(refs).map((r) => [r, 0])),
    causes: new Map(), refCauses: new Map()};
  const add = (map, cause, i) => {
    const c = map.get(cause) ?? {count: 0, rows: []};
    c.count++;
    if (c.rows.length < 3)
      c.rows.push(i);
    map.set(cause, c);
  };
  const ref = refs[refName] ?? null;
  for (let i = 0; i < n; i++) {
    const c = crux[i];
    if (forced?.has(i)) { // an engine gave no output for it (its process died there): never identical
      add(base ? res.causes : res.refCauses, override(i), i);
      continue;
    }
    for (const [r, arr] of Object.entries(refs)) {
      if (c === arr[i])
        res.identicalRefs[r]++;
    }
    if (base && c === base[i]) {
      res.identical++;
      if (ref && c !== ref[i]) {
        const sub = traps.has(i) ? classify?.(i) : null;
        add(res.refCauses, traps.has(i) ? (sub ? `${TRAP_CAUSE}: ${sub}` : TRAP_CAUSE) :
          'Crux = Chem\'s RDKit 2024.09, not RDKit.js 2026.03.6', i);
      }
      continue;
    }
    let cause;
    const special = override?.(i) ?? null; // an engine that aborted or died on the item
    if (special) {
      add(base ? res.causes : res.refCauses, special, i);
      continue;
    }
    if (traps.has(i))
      cause = TRAP_CAUSE;
    else if (fallback.has(i))
      cause = 'Crux mode gave RDKit\'s (Crux did not read it)';
    else if (!base)
      cause = ref && c !== ref[i] ? 'differs from the reference' : null;
    else if (ref)
      cause = c === ref[i] ? 'RDKit 2024.09 vs 2026.03 (Crux = RDKit.js 2026.03.6)' : base[i] === ref[i] ?
        'Crux differs from both RDKit releases' : 'Crux, RDKit 2024.09 and 2026.03.6 all differ';
    else
      cause = 'differs';
    if (!cause)
      continue;
    const sub = classify?.(i);
    add(base ? res.causes : res.refCauses, sub ? `${cause}: ${sub}` : cause, i);
  }
  return res;
}
async function withExamples(causes, show) {
  const out = [];
  for (const [cause, {count, rows}] of [...causes].sort((a, b) => b[1].count - a[1].count)) {
    const examples = [];
    for (const i of rows)
      examples.push(await show(i));
    out.push({cause, count, examples});
  }
  return out;
}
function engineStats(engine, raw, n, extra = {}) {
  return {engine, label: LABEL[engine], ms: raw?.ms ?? null, usPerItem: raw ? us(raw.ms, n) : null,
    counts: raw?.counts ?? null, fallbacks: Array.isArray(raw?.fallbacks) ? raw.fallbacks.length : raw?.fallbacks ?? null,
    ...(raw?.traps?.length ? {traps: raw.traps.length} : {}), ...extra};
}
/** What Chem does on a Crux trap: on the main thread (drawing, reactions) crux-module.ts noteCruxError leaves RDKit for
 * the session; in a Crux worker the pool restarts and the call falls back to RDKit (molzip: the call fails). */
function trapNote(raw, where) {
  const traps = raw?.traps?.length ?? 0;
  if (!traps && !raw?.panics?.count)
    return [];
  const consequence = where === 'main' ? 'Chem marks the Crux module failed (crux-module.ts noteCruxError): RDKit ' +
    'for the rest of the session' : where === 'molzip' ? 'the Crux worker traps, the pool restarts and the call ' +
    'fails (molzip has no RDKit behind it)' : 'the Crux worker traps, the pool restarts and the whole call runs on RDKit';
  return [`Crux panicked on ${traps} items (${Object.entries(raw.panics?.messages ?? {}).map(([m, k]) => `${k}× ${m}`)
    .join('; ')}). In Chem: ${consequence}. Here every other item still ran on Crux (the outputs of the items after a ` +
    'panic come from the same instance).'];
}

// -- the functions: engines, and how they compare

const ROWS_OUT = [];
const pushRow = (row) => {
  ROWS_OUT.push(row);
  console.log(`  = ${row.function}${row.variant ? ` / ${row.variant}` : ''}: n ${row.n}; ` +
    `${row.identical ? `= RDKit mode ${row.identical.count}/${row.identical.of}` : ''} ` +
    `${(row.refs ?? []).map((r) => `= ${r.engine} ${r.count}/${r.of}`).join(' ')}`);
};
function baseRow(fn, variant, n, unit, engines, cmp, extra = {}) {
  return {function: fn, variant, n, unit, engines,
    identical: cmp && cmp.identical !== null ? {vs: 'rdkit', count: cmp.identical, of: n, pct: pct(cmp.identical, n)} :
      null,
    refs: cmp ? Object.entries(cmp.identicalRefs).map(([engine, count]) => ({engine, count, of: n, pct: pct(count, n)})) :
      [], ...extra};
}

const FUNCTIONS = {};

FUNCTIONS.parse = {engines: ['rdkit', 'crux', 'rdkit2026'], async compare(rows) {
  const [rd, cx, r6] = ['rdkit', 'crux', 'rdkit2026'].map((e) => loadRaw('parse', e));
  const mols = molecules(rows);
  const cmp = compareItems({crux: cx.items, base: rd.items, refs: {rdkit2026: r6.items}, n: cx.n,
    traps: new Set(cx.traps), classify: (i) => !cx.items[i] || !rd.items[i] ? 'read by one engine only' :
      /\[N@/.test(cx.items[i]) && !/\[N@/.test(r6.items[i]) ? 'a stereo N Crux keeps (as Python RDKit 2026.3.6; ' +
        'RDKit.js drops N stereocenters)' : stripStereo(cx.items[i]) === stripStereo(rd.items[i]) ? 'stereo marks only' :
        'other'});
  const show = (i) => `${clip(mols[i], 90)} → Crux ${clip(cx.items[i], 90)} | RDKit 2024.09 ${clip(rd.items[i], 90)} | ` +
    `2026.03.6 ${clip(r6.items[i], 90)}`;
  pushRow(baseRow('parse', 'getMolSafe + canonical SMILES', cx.n, 'molecule',
    ['rdkit', 'crux', 'rdkit2026'].map((e) => engineStats(e, loadRaw('parse', e), cx.n)), cmp,
    {causes: await withExamples(cmp.causes, show), refCauses: await withExamples(cmp.refCauses, show),
      notes: ['Output: the canonical SMILES of what getMolSafe reads (a query: its SMARTS). Crux mode has no fallback ' +
        'here: this is what Crux itself reads.', ...trapNote(cx, 'main')]}));
}};

FUNCTIONS.drawing = {engines: ['rdkit', 'crux', 'rdkit2026'], async compare(rows) {
  const raw = Object.fromEntries(['rdkit', 'crux', 'rdkit2026'].map((e) => [e, loadRaw('drawing', e)]));
  const mols = molecules(rows);
  const n = raw.crux.n;
  let draw = null;
  for (const [k, variant] of [[0, `plain, ${SVG_W}x${SVG_H}`], [1, `highlighted by ${HIGHLIGHT}, ${SVG_W}x${SVG_H}`]]) {
    const key = (e) => raw[e].items.map((x) => x[k].split(' ')[0] ?? '');
    const rounded = (e) => raw[e].items.map((x) => x[k].split(' ')[1] ?? '');
    const geometry = (e) => raw[e].items.map((x) => x[k].split(' ')[2] ?? '');
    const [rdR, cxR, r6R] = ['rdkit', 'crux', 'rdkit2026'].map(rounded);
    const [rdG, cxG, r6G] = ['rdkit', 'crux', 'rdkit2026'].map(geometry);
    const [rdL, cxL, r6L] = ['rdkit', 'crux', 'rdkit2026'].map((e) => raw[e].layouts);
    const cmp = compareItems({crux: key('crux'), base: key('rdkit'), refs: {rdkit2026: key('rdkit2026')}, n,
      fallback: new Set(raw.crux.fallbacks), traps: new Set(raw.crux.traps),
      classify: (i) => cxR[i] === rdR[i] ? 'equal at 0.1 px (sub-pixel numbers)' :
        cxG[i] === rdG[i] ? 'the same lines and labels (drawing order, numbering, line direction, styles differ)' :
          sameLayout(cxL[i], rdL[i]) ? 'the same atom places (within 0.01 Å), other lines or labels at 0.1 px ' +
            '(Kekulé double bonds, wedges, label positions, or 0.1 px rounding)' : 'other atom places (another layout)'});
    const layoutSame = cxL.filter((x, i) => sameLayout(x, rdL[i])).length;
    const layoutSame6 = cxL.filter((x, i) => sameLayout(x, r6L[i])).length;
    const roundedSame = cxR.filter((x, i) => x === rdR[i]).length;
    const roundedSame6 = cxR.filter((x, i) => x === r6R[i]).length;
    const geometrySame = cxG.filter((x, i) => x === rdG[i]).length;
    const geometrySame6 = cxG.filter((x, i) => x === r6G[i]).length;
    const show = async (i) => {
      draw ??= Object.fromEntries(await Promise.all(['rdkit', 'crux', 'rdkit2026'].map(async (e) => [e, await ONE.drawing(e)])));
      const pics = Object.fromEntries(Object.entries(draw).map(([e, f]) => [e, f(mols[i])[k ? 'highlighted' : 'plain'] ?? '']));
      return `${clip(mols[i], 80)}: Crux vs RDKit 2024.09 ${firstDifference(pics.crux, pics.rdkit)}; vs 2026.03.6 ` +
        (pics.crux === pics.rdkit2026 ? 'identical' : firstDifference(pics.crux, pics.rdkit2026));
    };
    pushRow(baseRow('drawing', variant, n, 'molecule (both SVGs timed together)',
      ['rdkit', 'crux', 'rdkit2026'].map((e) => engineStats(e, raw[e], n)), cmp,
      {causes: await withExamples(cmp.causes, show), refCauses: await withExamples(cmp.refCauses, show),
        extra: {equalAtTenthPx: {rdkit: roundedSame, rdkit2026: roundedSame6},
          sameLinesAndLabels: {rdkit: geometrySame, rdkit2026: geometrySame6},
          sameAtomPlaces: {rdkit: layoutSame, rdkit2026: layoutSame6}},
        notes: ['Output: the exact SVG text of get_svg_with_highlights with the renderer\'s options ' +
          '(RDKIT_COMMON_RENDER_OPTS, its depiction calls); "equal at 0.1 px": the SVGs with every decimal number ' +
          'rounded to 0.1; "the same lines and labels": every path\'s and label\'s coordinates as written, a line\'s ' +
          'two ends in a fixed order, sorted ' +
          `(${geometrySame} of ${n} vs RDKit mode, ${geometrySame6} vs RDKit.js 2026.03.6); "the same atom places": ` +
          'the depicted molecule\'s atoms (element, x, y of its MOL block) matched within 0.01, computed apart from the ' +
          'timed run ' +
          `(${layoutSame} of ${n} vs RDKit mode, ${layoutSame6} vs RDKit.js 2026.03.6). One time covers both ` +
          'pictures of a molecule.', ...trapNote(raw.crux, 'main')]}));
  }
}};

FUNCTIONS.fingerprints = {engines: ['rdkit', 'crux', 'rdkit2026', 'py'], async run(rows) {
  for (const e of ['rdkit', 'crux', 'rdkit2026'])
    child('fingerprints', e, rows);
  python(['fcfp', inputs(`molecules-${rows}.json`), rawPrefix('fingerprints', 'py-fcfp')]);
}, async compare(rows) {
  const raw = Object.fromEntries(['rdkit', 'crux', 'rdkit2026'].map((e) => [e, loadRaw('fingerprints', e)]));
  const mols = molecules(rows);
  const n = raw.crux.n;
  const pyPrefix = rawPrefix('fingerprints', 'py-fcfp');
  const pyMeta = existsSync(pyPrefix + '.json') ? readJson(pyPrefix + '.json') : null;
  const pyBits = pyMeta ? readFileSync(pyPrefix + '.bin') : null;
  const pyItems = pyMeta ? pyMeta.read.map((ok, i) => ok ? hash(pyBits.subarray(256 * i, 256 * i + 256)) : '') : null;
  let fpOne = null;
  for (const [name, type, options] of FP_VARIANTS) {
    const v = Object.fromEntries(Object.entries(raw).map(([e, r]) => [e, r.variants[name]]));
    const fcfp = !!options?.useFeatures;
    const cruxOnly = cruxOnlyFp(options);
    const refs = fcfp ? (pyItems ? {py: pyItems} : {}) : {rdkit2026: v.rdkit2026.items};
    const cmp = compareItems({crux: v.crux.items, base: cruxOnly ? null : v.rdkit.items, refs, n,
      refName: fcfp ? 'py' : 'rdkit2026', fallback: new Set(v.crux.fallbacks), traps: new Set(v.crux.traps)});
    const show = async (i) => {
      fpOne ??= Object.fromEntries(await Promise.all(['rdkit', 'crux', 'rdkit2026'].map(async (e) => [e, await ONE.fingerprints(e)])));
      const fps = {};
      for (const e of ['rdkit', 'crux', 'rdkit2026'])
        fps[e] = (await fpOne[e](type, options, [mols[i]])).fps[0];
      const bitsOf = (fp) => fp ? new Set([...fp].flatMap((byte, k) => [0, 1, 2, 3, 4, 5, 6, 7].filter((b) => byte & (1 << b)).map((b) => 8 * k + b))) : new Set();
      const diff = (a, b) => {
        const [x, y] = [bitsOf(a), bitsOf(b)];
        return `${[...x].filter((b) => !y.has(b)).length} bits only Crux, ${[...y].filter((b) => !x.has(b)).length} only the other`;
      };
      return `${clip(mols[i], 90)}: vs RDKit 2024.09 ${diff(fps.crux, fps.rdkit)}; vs 2026.03.6 ${diff(fps.crux, fps.rdkit2026)}`;
    };
    const engines = ['rdkit', 'crux', 'rdkit2026'].filter((e) => v[e]).map((e) => engineStats(e, v[e], n));
    if (fcfp && pyMeta)
      engines.push(engineStats('py', {ms: pyMeta.ms, counts: tally(pyMeta.read.map((r) => r ? 'computed' : 'none'))}, n));
    pushRow(baseRow('fingerprints', name, n, 'molecule', engines, cmp, {
      causes: await withExamples(cmp.causes, show), refCauses: await withExamples(cmp.refCauses, show),
      notes: [cruxOnly ? `Crux's in both modes (Chem's RDKit.js has no ${fcfp ? 'FCFP' : 'chirality'}): no RDKit-mode ` +
        `engine; reference ${fcfp ? 'Python RDKit (GetMorganFeatureAtomInvGen; RDKit.js ignores useFeatures)' :
        'RDKit.js 2026.03.6'}.` : 'RDKit mode: Chem\'s RDKit worker on RDKit.js 2024.09.'].concat(fcfp ?
        ['RDKit.js 2026.03.6 gives plain ECFP for useFeatures (it reads no such option): ' +
          `${v.rdkit2026.items.filter((x, i) => x === raw.rdkit2026.variants.Morgan.items[i]).length} of ${n} ` +
          'identical to its Morgan r2.'] : [], trapNote({traps: v.crux.traps ?? [], panics: v.crux.traps?.length ?
          raw.crux.panics : null}, 'worker'))}));
  }
}};

for (const fn of ['alerts', 'pharmacophore']) {
  FUNCTIONS[fn] = {engines: ['rdkit', 'crux', 'rdkit2026'], async compare(rows) {
    const raw = Object.fromEntries(['rdkit', 'crux', 'rdkit2026'].map((e) => [e, loadRaw(fn, e)]));
    const mols = molecules(rows);
    const n = raw.crux.n;
    const ruleDiff = (a, b) => {
      const [x, y] = [new Set(a.split(',').filter(Boolean)), new Set(b.split(',').filter(Boolean))];
      return [...new Set([...x, ...y])].filter((r) => x.has(r) !== y.has(r)).join(',');
    };
    const cmp = compareItems({crux: raw.crux.items, base: raw.rdkit.items, refs: {rdkit2026: raw.rdkit2026.items}, n,
      fallback: new Set(raw.crux.fallbacks), traps: new Set(raw.crux.traps),
      classify: (i) => `rules ${ruleDiff(raw.crux.items[i], raw.rdkit.items[i])}`});
    const show = (i) => `${clip(mols[i], 90)}: Crux [${raw.crux.items[i]}] RDKit 2024.09 [${raw.rdkit.items[i]}] ` +
      `2026.03.6 [${raw.rdkit2026.items[i]}]`;
    pushRow(baseRow(fn === 'alerts' ? 'structural alerts' : 'pharmacophore features',
      `${Object.keys(raw.crux.rules).length} sets, ${Object.values(raw.crux.rules).reduce((a, b) => a + b, 0)} SMARTS`,
      n, 'molecule', ['rdkit', 'crux', 'rdkit2026'].map((e) => engineStats(e, raw[e], n)), cmp,
      {causes: await withExamples(cmp.causes, show), refCauses: await withExamples(cmp.refCauses, show),
        extra: {perRule: Object.fromEntries(Object.keys(raw.crux.perRule).map((r) => [r,
          {crux: raw.crux.perRule[r], rdkit: raw.rdkit.perRule[r], rdkit2026: raw.rdkit2026.perRule[r]}]))},
        notes: ['Output per molecule: the rule sets (families) with a match.', ...trapNote(raw.crux, 'worker')]}));
  }};
}

FUNCTIONS.search = {engines: ['rdkit+crux (parity.mjs, one process)'], async run(rows) {
  const sets = 'smiles,molblock,self,filtercatalog,bdb-1k';
  const r = spawnSync(process.execPath, ['--max-old-space-size=5500', join(HERE, 'parity.mjs'), String(rows), sets],
    {env: {...process.env, VERBOSE: '1', CRUX_BENCH, DATASET}, encoding: 'utf8', maxBuffer: 1 << 28});
  mkdirSync(join(RAW, 'search'), {recursive: true});
  writeFileSync(join(RAW, 'search', 'parity.log'), r.stdout + r.stderr);
  console.log(r.stdout.split('\n').slice(-3).join('\n'));
  const diffs = [];
  for (const line of r.stdout.split('\n')) {
    const m = line.match(/^!! \[(\S+)\] (.{36}) rdkit +(\S+) crux +(\S+)/);
    if (m)
      diffs.push({set: m[1], label: m[2].trim(), rdkit: +m[3], crux: +m[4]});
  }
  writeJson(inputs(`search-diffs-${rows}.json`), diffs);
  if (diffs.length)
    child('search', 'rdkit2026', rows);
}, async compare(rows) {
  const log = readFileSync(join(RAW, 'search', 'parity.log'), 'utf8');
  const prepared = +log.match(/RDKit: prepared in (\d+) ms/)[1];
  const [, indexed, notParsed, forRdkit] = log.match(/crux: indexed in (\d+) ms, (\d+) rows not parsed, (\d+) of them/);
  const sets = {};
  const differing = [];
  for (const line of log.split('\n')) {
    const m = line.match(/^(!!|RD| {2}) \[(\S+)\] (.{36}) rdkit +(\S+) crux +(\S+) +(\d+) \/ +([\d.]+) ms/);
    if (!m)
      continue;
    const [, flag, set, label, rdHits, cxHits, rdMs, cxMs] = m;
    const s = sets[set] ??= {queries: 0, identical: 0, onRdkit: 0, rdkitMs: 0, cruxMs: 0, rdkitHits: 0, cruxHits: 0};
    s.queries++;
    s.identical += flag !== '!!' ? 1 : 0;
    s.onRdkit += flag === 'RD' ? 1 : 0;
    s.rdkitMs += +rdMs;
    s.cruxMs += +cxMs;
    s.rdkitHits += rdHits === 'undefined' ? 0 : +rdHits;
    s.cruxHits += cxHits === 'undefined' ? 0 : +cxHits;
    if (flag === '!!')
      differing.push(`[${set}] ${label.trim()} RDKit ${rdHits} Crux ${cxHits}`);
  }
  const q = Object.values(sets).reduce((a, s) => a + s.queries, 0);
  const same = Object.values(sets).reduce((a, s) => a + s.identical, 0);
  const rdMs = Object.values(sets).reduce((a, s) => a + s.rdkitMs, 0);
  const cxMs = Object.values(sets).reduce((a, s) => a + s.cruxMs, 0);
  pushRow({function: 'Contains / Not contains search', variant: Object.keys(sets).join(', '), n: q,
    unit: `query over ${rows} molecules`, engines: [
      {engine: 'rdkit', label: LABEL.rdkit, ms: rdMs, usPerItem: us(rdMs, q), prepareMs: prepared,
        prepareUsPerMolecule: us(prepared, rows)},
      {engine: 'crux', label: LABEL.crux, ms: cxMs, usPerItem: us(cxMs, q), prepareMs: +indexed,
        prepareUsPerMolecule: us(+indexed, rows), counts: {'rows Crux did not index': +notParsed,
          'of them RDKit reads (searched on RDKit)': +forRdkit}}],
    identical: {vs: 'rdkit', count: same, of: q, pct: pct(same, q)}, refs: [],
    causes: (() => {
      const ref = loadRaw('search', 'rdkit2026');
      if (!differing.length)
        return [];
      if (!ref)
        return [{cause: 'hit sets differ', count: differing.length, examples: differing.slice(0, 5)}];
      const byCause = new Map();
      const ringOrValence = (q) => /\[(?:[^\]]*?[;&,!\w])?[RrvX]\d/.test(q) || /V2000/.test(q) &&
        q.split('\n').slice(4).some((l) => /^ {3,}-?\d+\.\d{4}/.test(l) && l.length >= 51 && +l.slice(48, 51) !== 0);
      for (const d of ref.items) {
        const cause = (d.rdkit2026 === null ? 'RDKit.js 2026.03.6 reads no query' : d.rdkit2026 === d.crux ?
          'RDKit 2024.09 vs 2026.03 (Crux\'s hit count = RDKit.js 2026.03.6\'s)' : d.rdkit2026 === d.rdkit ?
            'Crux differs from both RDKit releases' : 'all three differ') + (d.rdkit2026 !== null &&
          d.rdkit2026 !== d.crux && ringOrValence(d.query ?? '') ? ': a ring-count/ring-size/valence query (RDKit ' +
            'searches unsanitized SMILES: FastFindRings ring info, no valences; Crux\'s index: the sanitized ' +
            'molecule\'s SSSR)' : '');
        const c = byCause.get(cause) ?? {cause, count: 0, examples: []};
        c.count++;
        if (c.examples.length < 3)
          c.examples.push(`[${d.set}] ${clip(d.set === 'molblock' ? `${d.label} (as a MOL block)` : d.query ?? d.label,
            120).replace(/\n/g, ' ')}: RDKit 2024.09 ${d.rdkit}, Crux ${d.crux}, ` +
            `RDKit.js 2026.03.6 ${d.rdkit2026} hits`);
        byCause.set(cause, c);
      }
      return [...byCause.values()];
    })(),
    extra: {sets}, notes: ['parity.mjs (Chem\'s RDKit search pipeline: pattern fingerprints + unsanitized canonical ' +
      'SMILES; Crux: the column index and the query getQueryMolSafe builds on the Crux module, RDKit for rows Crux ' +
      'does not read); one process runs both engines (parity.mjs\'s design), one thread. Time: the searches\' sum ' +
      '(per query); "prepare": RDKit\'s fingerprints/targets, Crux\'s index build.']});
}};

/** The search types' comparison (searchtypes: the dataset; searchtypesfiles: Chem's test datasets). */
function compareSearchTypes(fn, source) {
  return async (rows) => {
    const raw = Object.fromEntries(['rdkit', 'crux', 'rdkit2026'].map((e) => [e, loadRaw(fn, e)]));
    const items = raw.crux.items;
    const molsAll = source === 'files' ? chemFileMolecules() : molecules(rows);
    const byType = new Map();
    items.forEach((it, i) => {
      const type = it.type.startsWith('Similar') ? it.type.split(' ').slice(0, 2).join(' ') : it.type;
      byType.set(type, [...byType.get(type) ?? [], i]);
    });
    for (const [type, idx] of byType) {
      // a query Crux does not read is searched by RDKit in Crux mode: RDKit mode's hits
      const key = (e) => idx.map((i) => raw[e].items[i].key ?? raw.rdkit.items[i].key);
      const onRdkit = idx.filter((i) => items[i].onRdkit).length;
      const cmp = compareItems({crux: key('crux'), base: key('rdkit'), refs: {rdkit2026: key('rdkit2026')}, n: idx.length,
        classify: (k) => items[idx[k]].set});
      const show = (k) => {
        const i = idx[k];
        const it = items[i];
        const rowsOf = (e) => new Set(raw[e].items[i].rows ?? raw[e].items[i].outside ?? []);
        const [cx, rd, r6] = ['crux', 'rdkit', 'rdkit2026'].map(rowsOf);
        const not = !!it.outside;
        const only = (a, b) => [...a].filter((x) => !b.has(x)).slice(0, 3).map((r) => `${r} ${clip(molsAll[r], 80)}`);
        return `${it.type} [${it.set}] ${clip(it.label, 70)}: hits Crux ${it.key}, RDKit 2024.09 ${raw.rdkit.items[i].key}, ` +
          `2026.03.6 ${raw.rdkit2026.items[i].key}; ${not ? 'outside the hits' : 'hits'} only in Crux: ` +
          `${JSON.stringify(only(cx, rd))}, only in RDKit 2024.09: ${JSON.stringify(only(rd, cx))}`;
      };
      const ms = (e) => idx.reduce((t, i) => t + (raw[e].items[i].ms ?? 0), 0);
      pushRow(baseRow(fn, `${type}${source === 'files' ? ' (Chem\'s test datasets)' : ''}`, idx.length,
        `query over ${raw.crux.rows} molecules`, ['rdkit', 'crux', 'rdkit2026'].map((e) => engineStats(e,
          {ms: ms(e) || null, counts: e === 'crux' ? {...raw.crux.counts, 'queries on RDKit here': onRdkit} : null},
          idx.length, {prepareMs: raw[e].prepareMs, prepareUsPerMolecule: us(raw[e].prepareMs, raw[e].rows)})), cmp,
        {causes: await withExamples(cmp.causes, show), refCauses: await withExamples(cmp.refCauses, show),
          notes: ['Output per query: the hit set (its size and hash). RDKit mode: Chem\'s RDKit search as its pool ' +
            'runs it on a column whose rows its workers cache (getMolSafe\'s sanitized molecules; ≤10,000 per worker), ' +
            'its pattern screen; Crux mode: crux/jobs/search-types.ts, RDKit for the rows Crux does not read and a ' +
            'query Crux does not read. One thread each; "prepare": the pattern fingerprints and canonical SMILES ' +
            '(RDKit: getFingerprints; Crux: the segment keys). Similar times: the Tanimoto only.']}));
    }
  };
}
FUNCTIONS.searchtypes = {engines: ['rdkit', 'crux', 'rdkit2026'], compare: compareSearchTypes('searchtypes', 'dataset')};
FUNCTIONS.searchtypesfiles = {engines: ['rdkit', 'crux', 'rdkit2026'],
  compare: compareSearchTypes('searchtypesfiles', 'files')};

FUNCTIONS.descriptors = {engines: ['crux', 'chempy', 'container'], async run(rows) {
  child('descriptors', 'crux', rows);
  const names = loadRaw('descriptors', 'crux').names;
  writeJson(inputs('descriptor-names-2d.json'), names);
  python(['descriptors', inputs(`molecules-${rows}.json`), inputs('descriptor-names-2d.json'),
    rawPrefix('descriptors', 'chempy')]);
  // 3D: a slice (chem.py embeds each molecule: ~0.1 s each)
  child('descriptors3d', 'crux', rows);
  const slice = parseInt(process.env.SLICE3D ?? '1000');
  writeJson(inputs(`molecules-3d-${slice}.json`), molecules(rows).slice(0, slice));
  writeJson(inputs('descriptor-names-3d.json'), loadRaw('descriptors3d', 'crux').names);
  python(['descriptors', inputs(`molecules-3d-${slice}.json`), inputs('descriptor-names-3d.json'),
    rawPrefix('descriptors3d', 'chempy')]);
  // the container itself, on a slice
  if (process.env.CONTAINER) {
    try {
      const mols = molecules(rows).slice(0, 2000);
      const body = JSON.stringify({molecules: mols, descriptors: names});
      await fetch(`${process.env.CONTAINER}/chem/descriptors`, {method: 'POST',
        body: JSON.stringify({molecules: mols.slice(0, 64), descriptors: names})});
      const t0 = performance.now();
      const res = await (await fetch(`${process.env.CONTAINER}/chem/descriptors`, {method: 'POST', body})).json();
      saveRaw('descriptors', 'container', {n: mols.length, ms: performance.now() - t0, values: res});
    } catch (e) {
      console.log(`  the container (${process.env.CONTAINER}) did not answer: ${e?.message ?? e}`);
    }
  }
}, async compare(rows) {
  const meta = loadRaw('descriptors', 'crux');
  const n = meta.n;
  const names = meta.names;
  const cx = readF64(rawFile('descriptors', 'crux', 'f64'));
  const pyPrefix = rawPrefix('descriptors', 'chempy');
  const pyMeta = readJson(pyPrefix + '.json');
  const py = readF64(pyPrefix + '.f64');
  const mols = molecules(rows);
  const tol = (a, b) => Math.abs(a - b) <= 1e-12 * Math.max(1, Math.abs(b));
  let bitSame = 0;
  let within = 0;
  let nullsBoth = 0;
  const per = {};
  // the values not bit-identical: their distance in units in the last place, and the largest relative difference
  const cxBits = new BigInt64Array(cx.buffer, cx.byteOffset, cx.length);
  const pyBits = new BigInt64Array(py.buffer, py.byteOffset, py.length);
  const ulpBuckets = {'1': 0, '2': 0, '3-4': 0, '5-16': 0, '17-256': 0, '>256': 0};
  let maxRel = 0;
  let maxRelOne = 0;
  const ulpPer = {};
  for (let k = 0; k < names.length; k++) {
    const p = per[names[k]] = {beyond: 0, cruxNullOnly: 0, pyNullOnly: 0, maxAbs: 0, rows: []};
    for (let i = 0; i < n; i++) {
      const a = cx[k * n + i];
      const b = py[k * n + i];
      if (Number.isNaN(a) && Number.isNaN(b)) {
        nullsBoth++;
        bitSame++;
        within++;
        continue;
      }
      if (Number.isNaN(a) || Number.isNaN(b)) {
        Number.isNaN(a) ? p.cruxNullOnly++ : p.pyNullOnly++;
        if (p.rows.length < 2)
          p.rows.push(i);
        continue;
      }
      if (Object.is(a, b) || a === b)
        bitSame++;
      else {
        const d = cxBits[k * n + i] - pyBits[k * n + i];
        const ulp = Math.sign(a) === Math.sign(b) ? Number(d < 0n ? -d : d) : Infinity;
        ulpBuckets[ulp <= 1 ? '1' : ulp <= 2 ? '2' : ulp <= 4 ? '3-4' : ulp <= 16 ? '5-16' : ulp <= 256 ? '17-256' :
          '>256']++;
        maxRel = Math.max(maxRel, Math.abs(a - b) / Math.max(Math.abs(a), Math.abs(b)));
        maxRelOne = Math.max(maxRelOne, Math.abs(a - b) / Math.max(1, Math.abs(a), Math.abs(b)));
        ulpPer[names[k]] = (ulpPer[names[k]] ?? 0) + 1;
      }
      if (tol(a, b))
        within++;
      else {
        p.beyond++;
        p.maxAbs = Math.max(p.maxAbs, Math.abs(a - b));
        if (p.rows.length < 2)
          p.rows.push(i);
      }
    }
  }
  const total = names.length * n;
  const typeDiff = names.filter((nm) => meta.types[nm] !== pyMeta.types[nm]);
  const differing = Object.entries(per).filter(([, p]) => p.beyond + p.cruxNullOnly + p.pyNullOnly)
    .sort((a, b) => (b[1].beyond + b[1].cruxNullOnly + b[1].pyNullOnly) - (a[1].beyond + a[1].cruxNullOnly + a[1].pyNullOnly));
  const causes = differing.map(([name, p]) => ({cause: `${name}: ${p.beyond} beyond 1e-12 (max |Δ| ${p.maxAbs.toPrecision(3)}), ` +
    `${p.cruxNullOnly} null in Crux only, ${p.pyNullOnly} null in chem.py only`, count: p.beyond + p.cruxNullOnly + p.pyNullOnly,
  examples: p.rows.map((i) => `${clip(mols[i], 80)}: Crux ${cx[names.indexOf(name) * n + i]} chem.py ${py[names.indexOf(name) * n + i]}`)}));
  const rowsWith = new Set();
  for (const [name, p] of differing)
    p.rows.forEach((i) => rowsWith.add(i));
  pushRow({function: 'descriptors (2D)', variant: `${names.length} descriptors`, n, unit: 'molecule (all names)',
    engines: [engineStats('chempy', pyMeta, n), engineStats('crux', meta, n, {fallbacks: meta.fallbacks.length})],
    identical: {vs: 'chempy', count: bitSame, of: total, pct: pct(bitSame, total), unit: 'values bit-identical to ' +
      'chem.py on Python RDKit 2026.3.6 (RDKit mode\'s code; the container itself not measured here)'},
    refs: [{engine: 'chempy (|Δ| ≤ 1e-12·max(1,|v|))', count: within, of: total, pct: pct(within, total)}],
    causes, extra: {valuesNullInBoth: nullsBoth, typesDiffer: typeDiff, descriptorsWithDifferences: differing.length,
      notBitIdentical: {ulpBuckets, maxRelDiff: maxRel, descriptors: Object.keys(ulpPer).length,
        top: Object.entries(ulpPer).sort((x, y) => y[1] - x[1]).slice(0, 12)}},
    notes: [`Not bit-identical: ${total - bitSame} values of ${Object.keys(ulpPer).length} descriptors, all floating-` +
      `point rounding (largest |Δ|/max(1,|v|) ${maxRelOne.toPrecision(2)}, largest relative difference ` +
      `${maxRel.toPrecision(2)} on values near 0); distance in units in the last place: ` +
      `${Object.entries(ulpBuckets).map(([b, c]) => `${b}: ${c}`).join(', ')}; most: ${Object.entries(ulpPer)
        .sort((x, y) => y[1] - x[1]).slice(0, 6).map(([nm, c]) => `${nm} ${c}`).join(', ')}.`,
    'RDKit mode computes these in the chem-chem container (chem.py, RDKit 2026.03.5); here chem.py\'s own ' +
      'get_descriptors runs single-process on Python RDKit 2026.3.6 in chunks of 1000 (the container splits a call over ' +
      'its processes). Crux: calcDescriptors in chunks of 1000 (as the Crux workers get them), one thread. ' +
      `Crux mode asks the container for a molecule whose Crux values are all null: ${meta.fallbacks.length} rows.`,
      ...trapNote(meta, 'worker')]});
  // 3D
  const m3 = loadRaw('descriptors3d', 'crux');
  if (m3) {
    const p3prefix = rawPrefix('descriptors3d', 'chempy');
    const p3meta = readJson(p3prefix + '.json');
    const c3 = readF64(rawFile('descriptors3d', 'crux', 'f64'));
    const p3 = readF64(p3prefix + '.f64');
    const n3 = m3.n;
    const stats = {};
    for (let k = 0; k < m3.names.length; k++) {
      const rel = [];
      let nullC = 0;
      let nullP = 0;
      for (let i = 0; i < n3; i++) {
        const a = c3[k * n3 + i];
        const b = p3[k * n3 + i];
        nullC += Number.isNaN(a) ? 1 : 0;
        nullP += Number.isNaN(b) ? 1 : 0;
        if (!Number.isNaN(a) && !Number.isNaN(b))
          rel.push(Math.abs(a - b) / Math.max(1e-9, Math.abs(b)));
      }
      rel.sort((x, y) => x - y);
      stats[m3.names[k]] = {within1e6: rel.filter((x) => x <= 1e-6).length, of: n3,
        medianRelDiff: rel[Math.floor(rel.length / 2)], p90RelDiff: rel[Math.floor(rel.length * 0.9)],
        nullCrux: nullC, nullChemPy: nullP};
    }
    pushRow({function: 'descriptors (3D)', variant: `${m3.names.length} descriptors, first ${n3} molecules`, n: n3,
      unit: 'molecule (all names, with the embedding)',
      engines: [engineStats('chempy', p3meta, n3), engineStats('crux', m3, n3)], identical: null, refs: [], causes: [],
      extra: {perDescriptor: stats}, notes: ['Not comparable value by value: chem.py embeds unseeded (ETKDG; RDKit\'s ' +
        'unseeded source starts as seed 42 in a fresh process and moves on with every embedding), Crux embeds each ' +
        'molecule as a fresh process\'s first embedding; "within 1e-6": relative difference.']});
  }
  const ct = loadRaw('descriptors', 'container');
  if (ct) {
    let same = 0;
    let tot = 0;
    for (let k = 0; k < names.length; k++) {
      const vals = ct.values[names[k]]?.value ?? [];
      for (let i = 0; i < ct.n; i++) {
        tot++;
        const a = vals[i] === null || vals[i] === undefined ? NaN : Number(vals[i]);
        const b = py[k * n + i];
        if ((Number.isNaN(a) && Number.isNaN(b)) || tol(a, b))
          same++;
      }
    }
    pushRow({function: 'descriptors (2D), the container itself', variant: `first ${ct.n} molecules`, n: ct.n,
      unit: 'molecule (all names)', engines: [engineStats('container', ct, ct.n)],
      identical: null, refs: [{engine: 'chempy (|Δ| ≤ 1e-12·max(1,|v|))', count: same, of: tot, pct: pct(same, tot)}],
      causes: [], notes: [`POST ${process.env.CONTAINER}/chem/descriptors, one request (the container splits it over ` +
        'its processes), timed after a 64-molecule request; includes HTTP and JSON.']});
  }
}};

/** A scaffold tree's MOL blocks without their layout: atom coordinates and the bonds' stereo field (the wedges drawn
 * from them), as crux-core tools/chem-wasm/chem_engine_parity.mjs compares them. */
function sansLayout(block) {
  const l = block.split('\n');
  const na = Number(l[3]?.slice(0, 3));
  const nb = Number(l[3]?.slice(3, 6));
  for (let i = 4; i < 4 + na && i < l.length; i++)
    l[i] = l[i].slice(30);
  for (let i = 4 + na; i < 4 + na + nb && i < l.length; i++)
    l[i] = l[i].slice(0, 9) + l[i].slice(12);
  return l.join('\n');
}
function treeSummary(nodes, layout) {
  const out = [];
  const walk = (node, depth) => {
    out.push(`${depth}|${node.smiles}|${layout ? node.scaffold : sansLayout(node.scaffold)}|${node.child_nodes.length}`);
    node.child_nodes.forEach((c) => walk(c, depth + 1));
  };
  nodes.forEach((n) => walk(n, 0));
  return out;
}
FUNCTIONS.scaffoldtree = {engines: ['crux', 'sg'], async run(rows) {
  child('scaffoldtree', 'crux', rows);
  if (!process.env.SG_PYTHON)
    return;
  const budgetMs = parseFloat(process.env.SCAFFOLD_BUDGET_MIN ?? '45') * 60000;
  for (const size of scaffoldSizes(molecules(rows).length)) {
    const file = inputs(`scaffold-input-${size}.json`);
    writeJson(file, molecules(rows).slice(0, size));
    if (!python(['scaffoldtree', file, rawFile('scaffoldtree', `sg-${size}`)], {py: process.env.SG_PYTHON,
      timeoutMs: budgetMs})) {
      saveRaw('scaffoldtree', `sg-${size}-timeout`, {size, budgetMs});
      break;
    }
  }
}, async compare(rows) {
  const meta = loadRaw('scaffoldtree', 'crux');
  if (!meta)
    throw new Error('no Crux scaffold tree run (raw/scaffoldtree/crux.json)');
  const outcomes = [];
  for (const size of Object.keys(meta.sizes).map(Number)) {
    const cx = meta.sizes[size];
    const sgFile = rawFile('scaffoldtree', `sg-${size}`);
    const sgMeta = existsSync(sgFile + '.json') ? readJson(sgFile + '.json') : null;
    const timedOut = loadRaw('scaffoldtree', `sg-${size}-timeout`);
    const engines = [engineStats('crux', cx, size, {counts: cx.raised ? {raised: 1} : {tree: 1}})];
    const base = {function: 'scaffold tree', variant: `first ${size} molecules`, n: size, unit: 'molecule', engines};
    const notes = ['ringCutoff 10, dischargeAndDeradicalize false (getScaffoldTree\'s defaults); the script is Chem\'s ' +
      'scripts/scaffold_tree.py exec\'d as Datagrok runs it.', ...trapNote(meta, 'worker')];
    if (!sgMeta || timedOut) {
      const outcome = cx.raised ? 'Crux raises (the script not run at this size)' : 'Crux alone (the script not run)';
      outcomes.push({size, outcome});
      pushRow({...base, identical: null, refs: [], causes: cx.raised ? [{cause: 'Crux raises', count: 1,
        examples: [cx.raised]}] : [], notes: [timedOut ? `scaffold_tree.py did not finish within ` +
        `${timedOut.budgetMs / 60000} min at this size (killed)${timedOut.note ? `; ${timedOut.note}` : ''}: Crux timed alone.` : 'scaffold_tree.py not run at ' +
        'this size: Crux timed alone.', ...notes]});
      continue;
    }
    engines.unshift(engineStats('sg', sgMeta, size, {counts: sgMeta.raised ? {raised: 1} : {tree: 1}}));
    if (cx.raised || sgMeta.raised) {
      // the same exception: Crux's message names the script function the script raised in (its traceback's last)
      const where = /\(in (?:.*, )?(\w+):\d+\)$/.exec(sgMeta.raised ?? '')?.[1];
      const sameType = !!(cx.raised && sgMeta.raised && where && cx.raised.includes(where));
      const outcome = cx.raised && sgMeta.raised ? `both raise${sameType ? ' (the same exception)' : ''}` :
        cx.raised ? 'Crux raises, the script does not' : 'the script raises, Crux does not';
      outcomes.push({size, outcome});
      pushRow({...base, identical: {vs: 'sg', count: cx.raised && sgMeta.raised ? 1 : 0, of: 1,
        pct: cx.raised && sgMeta.raised ? '100 %' : '0 %', unit: 'outcome (both raise)'}, refs: [],
      causes: [{cause: outcome, count: 1, examples: [`script: ${sgMeta.raised ?? 'a tree'}`, `Crux: ${cx.raised ??
        'a tree'}`]}], notes});
      continue;
    }
    const cxText = readFileSync(rawFile('scaffoldtree', `crux-${size}`), 'utf8');
    const sgText = readFileSync(sgFile, 'utf8');
    const [a, b] = [JSON.parse(cxText), JSON.parse(sgText)];
    const [sa, sb] = [treeSummary(a, false), treeSummary(b, false)];
    let sameNodes = 0;
    for (let k = 0; k < Math.min(sa.length, sb.length); k++)
      sameNodes += sa[k] === sb[k] ? 1 : 0;
    const diffs = sa.flatMap((x, k) => x !== sb[k] ? [k] : []);
    const firstDiff = diffs[0] ?? -1;
    outcomes.push({size, outcome: cxText === sgText ? 'identical JSON' : firstDiff < 0 && sa.length === sb.length ?
      'identical sans layout' : 'differ'});
    pushRow({...base,
      identical: {vs: 'sg', count: cxText === sgText ? 1 : 0, of: 1, pct: cxText === sgText ? '100 %' : '0 %',
        unit: 'JSON text'},
      refs: [{engine: `sg, nodes in order with MOL blocks sans layout (${sb.length} nodes in the script's tree, ` +
        `${sa.length} in Crux's)`, count: sameNodes, of: Math.max(sa.length, sb.length),
      pct: pct(sameNodes, Math.max(sa.length, sb.length))}],
      causes: firstDiff >= 0 ? [{cause: 'nodes that differ (depth|SMILES|MOL block sans layout|children)',
        count: Math.max(sa.length, sb.length) - sameNodes,
        examples: diffs.slice(0, 3).map((k) => `node ${k} (${clip(sb[k]?.split('|')[1], 80)}): Crux vs script ` +
          firstDifference(sa[k], sb[k] ?? '').replace(/\n/g, '⏎'))}] :
        cxText === sgText ? [] : [{cause: 'MOL block coordinates (Crux\'s layout) only', count: 1, examples: []}],
      notes});
  }
  ROWS_OUT[ROWS_OUT.length - 1].extra = {outcomes};
}};

FUNCTIONS.properties = {engines: ['ocl', 'crux', 'py'], async run(rows) {
  child('properties', 'ocl', rows);
  child('properties', 'crux', rows);
  python(['props', inputs(`molecules-${rows}.json`), rawFile('properties', 'py')]);
}, async compare(rows) {
  const [ocl, cx] = ['ocl', 'crux'].map((e) => loadRaw('properties', e));
  const py = loadRaw('properties', 'py');
  const n = cx.n;
  const mols = molecules(rows);
  for (const p of PROPERTY_NAMES) {
    const c = cx.values[p];
    const o = ocl.values[p];
    const y = py.values[p];
    const isNum = typeof c.find((v) => v !== null) === 'number' && p !== 'Molecular formula';
    const float = ['MW', 'LogP', 'LogS', 'PSA'].includes(p);
    let same = 0;
    let twoDec = 0;
    let sumAbs = 0;
    let maxAbs = 0;
    let both = 0;
    let samePy = 0;
    const diffRows = [];
    for (let i = 0; i < n; i++) {
      const [a, b] = [c[i] ?? null, o[i] ?? null];
      if (a === b || (isNum && a !== null && b !== null && Object.is(a, b)))
        same++;
      else if (diffRows.length < 3)
        diffRows.push(i);
      if (float && a !== null && b !== null) {
        both++;
        const d = Math.abs(a - b);
        sumAbs += d;
        maxAbs = Math.max(maxAbs, d);
        twoDec += d < 0.005 ? 1 : 0;
      }
      const z = y[i] ?? null;
      if (a === z || (float && a !== null && z !== null && Math.abs(a - z) <= 1e-12 * Math.max(1, Math.abs(z))))
        samePy++;
    }
    pushRow({function: 'property table', variant: p, n, unit: 'molecule (all ten properties timed together)',
      engines: [engineStats('ocl', ocl, n), engineStats('crux', cx, n), engineStats('py', py, n)],
      identical: {vs: 'ocl', count: same, of: n, pct: pct(same, n)},
      refs: [{engine: 'py (RDKit\'s definitions' + (float ? '; |Δ| ≤ 1e-12·max(1,|v|)' : '') + ')', count: samePy, of: n,
        pct: pct(samePy, n)}],
      causes: same < n ? [{cause: float ? `OpenChemLib's method (agree to ±0.005: ${twoDec} of ${both}; mean |Δ| ` +
        `${(sumAbs / Math.max(1, both)).toFixed(4)}, max |Δ| ${maxAbs.toFixed(3)})` : 'OpenChemLib\'s definition',
      count: n - same, examples: diffRows.map((i) => `${clip(mols[i], 80)}: Crux ${c[i]} OCL ${o[i]} RDKit ${y[i]}`)}] : [],
      notes: trapNote(cx, 'worker')});
  }
}};

FUNCTIONS.reactions = {engines: ['rdkit', 'crux', 'rdkit2026'], async compare(rows) {
  const ENGINES = ['rdkit', 'crux', 'rdkit2026'];
  const raw = Object.fromEntries(ENGINES.map((e) => [e, loadRaw('reactions', e)]));
  if (!raw.crux)
    throw new Error('no Crux reactions run (raw/reactions/crux.json)');
  const missing = ENGINES.filter((e) => !raw[e]);
  for (const e of missing)
    raw[e] = {transformation: {}, twoComponent: {}, deprotect: {}, enumeration: {}};
  const prepared = readJson(inputs(`reaction-inputs-${rows}.json`));
  const mols = prepared.molecules;
  const canon = await loadRdkit2026();
  const canonical = (s) => {
    if (!s || s === '=')
      return s;
    const m = canon.get_mol(s);
    try {
      return m ? m.get_smiles() : null;
    } finally {
      m?.delete();
    }
  };
  const ENGINE_NAME = {rdkit: 'RDKit.js 2024.09', crux: 'Crux', rdkit2026: 'RDKit.js 2026.03.6'};
  for (const part of ['transformation', 'twoComponent', 'deprotect', 'enumeration']) {
    const keys = Object.keys(raw.crux[part]);
    let n = 0;
    let identical = 0;
    let ref = 0;
    let fb = 0;
    const causes = new Map();
    const refCauses = new Map();
    const ms = {rdkit: 0, crux: 0, rdkit2026: 0};
    const detail = {};
    let trapCount = 0;
    const aborted = {rdkit: 0, crux: 0, rdkit2026: 0}; // items an engine aborted on (RDKit.js abort contexts)
    const abortCalls = {rdkit: 0, crux: 0, rdkit2026: 0}; // the aborts themselves
    const crashedItems = {rdkit: 0, crux: 0, rdkit2026: 0};
    const abortExamples = {rdkit: [], crux: [], rdkit2026: []};
    const excCalls = {rdkit: 0, crux: 0, rdkit2026: 0}; // RDKit.js run_reactants exceptions (C++), and their items
    const excItems = {rdkit: 0, crux: 0, rdkit2026: 0};
    const excExamples = {rdkit: [], crux: [], rdkit2026: []};
    let cruxTrapCalls = 0;
    for (const key of keys) {
      const units = Object.fromEntries(ENGINES.map((e) => [e, raw[e][part][key] ?? null]));
      const [rd, cx, r6] = ENGINES.map((e) => units[e] ?? {items: []});
      for (const e of ENGINES)
        ms[e] += units[e]?.ms ?? 0;
      if (cx.threw || rd.threw) {
        detail[key] = {threw: {crux: cx.threw, rdkit: rd.threw}};
        continue;
      }
      const len = Math.max(cx.items.length, rd.items.length, r6.items.length);
      const pair = prepared.twoComponent[key];
      const clean = (x) => x?.split('.').filter((f) => !['O', 'Cl', 'Br', '[Na+]', '[K+]'].includes(f)).join('.');
      const inputsOf = (i) => part === 'transformation' ? [prepared.cleaned[i]] : part === 'deprotect' ? [mols[i]] :
        part === 'twoComponent' && pair ? [clean(pair.c1[i]), clean(pair.c2[i]), pair.c1[i], pair.c2[i]] : [];
      const itemsOf = (contexts) => {
        const set = new Set((contexts ?? []).filter((c) => c && c.part === part && c.key === key)
          .flatMap((c) => c.inputs ?? []));
        return new Set(set.size ? [...Array(len).keys()].filter((i) => inputsOf(i).some((x) => set.has(x))) : []);
      };
      // the items whose molecule(s) Crux panicked on (the contexts the child recorded)
      const traps = itemsOf(raw.crux.trapContexts);
      trapCount += traps.size;
      const abortedOn = Object.fromEntries(ENGINES.map((e) => [e, itemsOf(units[e]?.abortContexts)]));
      const died = Object.fromEntries(ENGINES.map((e) => [e, new Set((units[e]?.items ?? []).flatMap((x, i) =>
        x === 'CRASHED' || x === 'NOT RUN' ? [i] : []))]));
      cruxTrapCalls += cx.traps ?? 0;
      for (const e of ENGINES) {
        const thrownOn = itemsOf(units[e]?.exceptionContexts);
        excCalls[e] += units[e]?.exceptions ?? 0;
        excItems[e] += thrownOn.size;
        for (const i of thrownOn) {
          if (excExamples[e].length < 3)
            excExamples[e].push(`${key}: ${clip(inputsOf(i)[0] ?? `row ${i}`, 100)}`);
        }
        aborted[e] += abortedOn[e].size;
        abortCalls[e] += units[e]?.aborts ?? 0;
        crashedItems[e] += died[e].size;
        for (const i of abortedOn[e]) {
          if (abortExamples[e].length < 3)
            abortExamples[e].push(`${key}: ${clip(inputsOf(i)[0] ?? `row ${i}`, 100)}`);
        }
      }
      const override = (i) => {
        for (const e of ['crux', 'rdkit', 'rdkit2026']) {
          if (died[e].has(i))
            return `${ENGINE_NAME[e]}'s process died on it (no output)`;
        }
        for (const e of ['rdkit', 'rdkit2026', 'crux']) {
          if (abortedOn[e].has(i))
            return `${ENGINE_NAME[e]} aborted on it (Aborted(): an RDKit.js trap; Chem counts an error, the input kept)`;
        }
        return null;
      };
      const forced = new Set(ENGINES.flatMap((e) => [...died[e]]));
      const cmp = compareItems({crux: cx.items, base: missing.includes('rdkit') ? null : rd.items,
        refs: missing.includes('rdkit2026') ? {} : {rdkit2026: r6.items}, n: len, traps, override, forced,
        classify: (i) => traps.has(i) ? (['=', ''].includes(r6.items[i]) ? 'RDKit.js 2026.03.6 gives no product there ' +
          'either' : 'RDKit.js 2026.03.6 gives a product') :
          part === 'deprotect' && canonical(cx.items[i]) === canonical(rd.items[i]) ?
            'the MOL block text only (the same canonical SMILES)' :
            cx.items[i] && rd.items[i] && stripStereo(cx.items[i]) === stripStereo(rd.items[i]) ? 'stereo marks only' :
              null});
      n += len;
      identical += cmp.identical ?? 0;
      ref += cmp.identicalRefs.rdkit2026 ?? 0;
      fb += cx.fallbacks ?? 0;
      for (const [map, from] of [[causes, cmp.causes], [refCauses, cmp.refCauses]]) {
        for (const [cause, {count, rows: rr}] of from) {
          const c = map.get(cause) ?? {count: 0, rows: []};
          c.count += count;
          for (const i of rr) {
            if (c.rows.length < 3)
              c.rows.push([key, i]);
          }
          map.set(cause, c);
        }
      }
      detail[key] = {n: len, identical: cmp.identical, identicalRdkit2026: cmp.identicalRefs.rdkit2026,
        reacted: Object.fromEntries(ENGINES.map((e) => [e, (units[e]?.items ?? []).filter((x) => x !== '=').length])),
        ms: Object.fromEntries(ENGINES.map((e) => [e, units[e]?.ms])), fallbacks: cx.fallbacks,
        ...Object.fromEntries(ENGINES.filter((e) => units[e]?.aborts).map((e) => [`aborts ${e}`, units[e].aborts]))};
    }
    const show = async ([key, i]) => {
      const [rd, cx, r6] = ENGINES.map((e) => raw[e][part][key]?.items?.[i]);
      const input = part === 'transformation' || part === 'deprotect' ? mols[i] :
        part === 'twoComponent' ? (prepared.twoComponent[key] ? `${prepared.twoComponent[key].r1[i]} + ${prepared.twoComponent[key].r2[i]}` : '') : `row ${i}`;
      const same = part === 'deprotect' ? ` (as canonical SMILES: ${canonical(cx) === canonical(rd) ? 'the same' : 'different'})` : '';
      return `${key}: ${clip(input, 70)} → Crux ${clip(cx, 120)} | RDKit 2024.09 ${clip(rd, 120)} | 2026.03.6 ${clip(r6, 120)}${same}`;
    };
    const unit = {transformation: 'molecule × reaction', twoComponent: 'reactant pair', deprotect: 'molecule × fragment',
      enumeration: 'output row'}[part];
    const timesN = part === 'enumeration' ? Math.max(1, n) : n;
    const engineNotes = ENGINES.flatMap((e) => [
      missing.includes(e) ? `${ENGINE_NAME[e]}: no run (its child failed).` : null,
      abortCalls[e] ? `${ENGINE_NAME[e]} aborted ${abortCalls[e]} times (Aborted(), each a trap Chem's loop catches ` +
        `as an error: the input kept) on ${aborted[e]} items${abortExamples[e].length ? `, e.g. ${abortExamples[e]
          .join('; ')}` : ''}.` : null,
      crashedItems[e] ? `${ENGINE_NAME[e]}'s process died on ${crashedItems[e]} items (outputs 'CRASHED'/'NOT RUN').` :
        null,
      excCalls[e] ? `${ENGINE_NAME[e]}'s run_reactants threw (a C++ exception, which Chem counts as an error: the input ` +
        `kept) ${excCalls[e]} times, on ${excItems[e]} items${excExamples[e].length ? `, e.g. ${excExamples[e].join('; ')}` :
        ''}.` : null]).filter(Boolean);
    const restartNote = part === 'transformation' ? ENGINES.flatMap((e) => {
      const rs = raw[e].restarts ?? [];
      const heap = rs.filter((r) => !r.afterPanics);
      const pan = rs.filter((r) => r.afterPanics);
      const died = raw[e].crashed ?? [];
      return [rs.length ? `${ENGINE_NAME[e]}: ${rs.length} fresh processes between two units (${heap.length ?
        `${heap.length} with its wasm heap past ${HEAP_RESTART_BYTES / 2 ** 30} GB: ${heap.slice(0, 4).map((r) =>
          `${r.heapGb} GB before ${r.before}`).join(', ')}${heap.length > 4 ? ', …' : ''}` : ''}${heap.length &&
        pan.length ? '; ' : ''}${pan.length ? `${pan.length} after a unit with ${e === 'crux' ? 'Crux panics' :
          'run_reactants exceptions'}` : ''}).` : null,
      died.length ? `${ENGINE_NAME[e]}'s process died ${died.length} times (${Object.entries(tally(died.map((c) =>
        c.error || `exit ${c.exit}`))).map(([m, k]) => `${k}× ${clip(m, 90)}`).join('; ')}), in ${[...new Set(died.map(
        (c) => c.key.split('#')[0]))].length} units: each such unit was run again in chunks of 1000 in fresh processes ` +
        '(a chunk it died in again, item by item); items lost: see the counts above.' : null].filter(Boolean);
    }) : [];
    pushRow({function: 'reactions', variant: {transformation: `transformation, ${keys.length} one-reactant default reactions`,
      twoComponent: `two-component pairwise, ${keys.length - 1} default reactions + the demo matrix`,
      deprotect: `Deprotect (${keys.join(', ')})`, enumeration: `Reaction Enumeration (${keys.join('; ')})`}[part], n,
    unit, engines: ENGINES.map((e) => ({engine: e, label: LABEL[e], ms: missing.includes(e) ? null : ms[e],
      usPerItem: missing.includes(e) ? null : us(ms[e], timesN), fallbacks: e === 'crux' ? fb : null,
      counts: {...(abortCalls[e] ? {aborts: abortCalls[e], 'items aborted on': aborted[e]} : {}),
        ...(crashedItems[e] ? {'items the process died on': crashedItems[e]} : {})}})),
    identical: missing.includes('rdkit') ? null : {vs: 'rdkit', count: identical, of: n, pct: pct(identical, n)},
    refs: missing.includes('rdkit2026') ? [] : [{engine: 'rdkit2026', count: ref, of: n, pct: pct(ref, n)}],
    causes: await withExamples(causes, show), refCauses: await withExamples(refCauses, show), extra: {detail,
      canonicalizeMs: part === 'transformation' ? prepared.canonicalMs : undefined},
    notes: [...engineNotes, ...restartNote, ...(raw.crux.trapContexts ?? []).some((c) => c.part === part) ? [`Crux panicked on ` +
      `${(raw.crux.trapContexts ?? []).filter((c) => c.part === part).length} reaction runs (${trapCount} items): ` +
      `${[...new Set((raw.crux.trapContexts ?? []).filter((c) => c.part === part).map((c) => `${c.key} on ` +
        `${c.inputs[c.inputs.length - 1]}`))].slice(0, 5).join('; ')}. Chem's EngineReaction gives such a molecule ` +
      'RDKit\'s product and marks the Crux module failed (crux-module.ts noteCruxError): RDKit for the rest of the ' +
      'session. Here Crux ran on (the stub keeps the module, the panicked reaction object is made again), so every ' +
      'other item is Crux\'s.'] : [],
    ...part === 'transformation' ? ['Chem\'s runTransformationReaction per reaction over the dataset; RDKit\'s ' +
      `canonical SMILES of the molecules (convertMolNotation, RDKit in both modes) computed once beforehand: ` +
      `${prepared.canonicalMs.toFixed(0)} ms, not in the times.`] : part === 'deprotect' ? ['Output: a MOL block for a ' +
      'deprotected molecule, the input where none.'] : []]});
  }
}};

FUNCTIONS.mmpa = {engines: ['rdkit', 'crux', 'rdkit2026'], async compare(rows) {
  const raw = Object.fromEntries(['rdkit', 'crux', 'rdkit2026'].map((e) => [e, loadRaw('mmpa', e)]));
  const mols = molecules(rows);
  const n = raw.crux.n;
  const classify = (i) => {
    const [a, b] = [JSON.parse(raw.crux.items[i]), JSON.parse(raw.rdkit.items[i])];
    if (a[0] !== b[0])
      return `the molecule's SMILES differs${stripStereo(a[0]) === stripStereo(b[0]) ? ' (stereo marks only)' : ''}`;
    // the cuts as a multiset: each pair's pieces (stereo marks stripped from each piece before sorting)
    const sorted = (f, strip = false) => JSON.stringify((f ?? []).map((p) => (strip ? p.map(stripStereo) : p).join('.'))
      .sort());
    if (sorted(a[1]) === sorted(b[1]))
      return 'the same cuts in another order';
    if (sorted(a[1], true) === sorted(b[1], true))
      return 'stereo marks in the fragments';
    return (a[1]?.length ?? 0) !== (b[1]?.length ?? 0) ? 'another number of cuts' : 'fragments differ';
  };
  const cmp = compareItems({crux: raw.crux.items, base: raw.rdkit.items, refs: {rdkit2026: raw.rdkit2026.items}, n,
    fallback: new Set(raw.crux.fallbacks), traps: new Set(raw.crux.traps), classify});
  const show = (i) => `${clip(mols[i], 80)}: Crux ${clip(raw.crux.items[i], 150)} | RDKit 2024.09 ` +
    `${clip(raw.rdkit.items[i], 150)}`;
  pushRow(baseRow('MMPA fragmentation', 'get_mmpa_frags(1, 1, 20) via mmpGetFragments', n, 'molecule',
    ['rdkit', 'crux', 'rdkit2026'].map((e) => engineStats(e, raw[e], n)), cmp,
    {causes: await withExamples(cmp.causes, show), refCauses: await withExamples(cmp.refCauses, show),
      notes: ['Output per molecule: its canonical SMILES and its single-cut fragment pairs, in order.',
        ...trapNote(raw.crux, 'worker')]}));
}};

FUNCTIONS.rgroups = {engines: ['rdkit', 'crux', 'rdkit2026'], async compare(rows) {
  const raw = Object.fromEntries(['rdkit', 'crux', 'rdkit2026'].map((e) => [e, loadRaw('rgroups', e)]));
  const series = readJson(inputs(`rgroup-series-${rows}.json`));
  const n = raw.crux.n;
  // a series' output by column name: {name: [smiles, atoms, bonds] per molecule} (atoms and bonds: R columns only)
  const parsed = (e, k) => raw[e].items[k].startsWith('THROWS') ? null : JSON.parse(raw[e].items[k]);
  const byName = (r) => Object.fromEntries(r.colNames.map((c, j) => [c, r.smiles[j].map((s, i) =>
    JSON.stringify([s, j ? r.atoms[j - 1][i] ?? null : null, j ? r.bonds[j - 1][i] ?? null : null]))]));
  // the decomposition with its columns matched by name (RDKit orders the labels as strings, R1 R10 R11 R2 ...; Crux
  // by number, R1 R2 ... R10: Chem's columns, the trellis's default Y axis and the highlight colours follow the order)
  const sameByName = (a, b) => {
    if (!a || !b || a.colNames.length !== b.colNames.length)
      return false;
    const [x, y] = [byName(a), byName(b)];
    return a.colNames.every((c) => y[c] && JSON.stringify(x[c]) === JSON.stringify(y[c]));
  };
  // cells: each series' columns matched by name, molecule by molecule (a column one engine lacks: every cell differs)
  let cells = 0;
  const sameCells = {rdkit: 0, rdkit2026: 0};
  const cellCauses = new Map();
  for (let k = 0; k < n; k++) {
    const a = parsed('crux', k);
    for (const ref of ['rdkit', 'rdkit2026']) {
      const b = parsed(ref, k);
      if (!a || !b)
        continue;
      const [x, y] = [byName(a), byName(b)];
      for (const col of new Set([...a.colNames, ...b.colNames])) {
        for (let i = 0; i < series[k].rows.length; i++) {
          if (ref === 'rdkit')
            cells++;
          const [cs, rs] = [x[col] ? JSON.parse(x[col][i])[0] : undefined, y[col] ? JSON.parse(y[col][i])[0] : undefined];
          if (cs === rs) {
            sameCells[ref]++;
            continue;
          }
          if (ref !== 'rdkit2026')
            continue;
          const cause = !x[col] || !y[col] ? 'a column one engine lacks' : cs && rs && stripStereo(cs) === stripStereo(rs) ?
            'stereo marks only' : 'other';
          const e = cellCauses.get(cause) ?? {count: 0, examples: []};
          e.count++;
          if (e.examples.length < 3)
            e.examples.push(`series ${k} (${series[k].core}) ${col} row ${i}: Crux ${cs} 2026.03.6 ${rs}`);
          cellCauses.set(cause, e);
        }
      }
    }
  }
  const classify = (k) => {
    const [a, b, c] = ['crux', 'rdkit', 'rdkit2026'].map((e) => parsed(e, k));
    if (sameByName(a, c) || sameByName(a, b))
      return 'the same decomposition (every cell and highlight), its columns in another order: RDKit orders the ' +
        'labels as strings (Core, R1, R10, R2, ...), Crux by number (R1, ..., R9, R10)';
    return a && c && a.colNames.length !== c.colNames.length ? 'other R groups (another set of labels)' :
      'other R groups or highlights';
  };
  const decompSame = {rdkit: 0, rdkit2026: 0};
  for (let k = 0; k < n; k++) {
    for (const ref of ['rdkit', 'rdkit2026'])
      decompSame[ref] += sameByName(parsed('crux', k), parsed(ref, k)) ? 1 : 0;
  }
  const cmp = compareItems({crux: raw.crux.items, base: raw.rdkit.items, refs: {rdkit2026: raw.rdkit2026.items}, n,
    fallback: new Set(raw.crux.fallbacks), traps: new Set(raw.crux.traps), classify});
  const show = (k) => `series ${k} (core ${series[k].core}, ${series[k].rows.length} molecules): Crux ` +
    `${clip(raw.crux.items[k], 140)} | RDKit 2024.09 ${clip(raw.rdkit.items[k], 140)}`;
  const mN = raw.crux.molecules;
  const row = baseRow('R-group decomposition', `${n} series, ${mN} molecules (Greedy, the dialog's options)`, n,
    'series (time per molecule in µs/item)', ['rdkit', 'crux', 'rdkit2026'].map((e) => engineStats(e, raw[e], mN)), cmp,
    {causes: await withExamples(cmp.causes, show), refCauses: await withExamples(cmp.refCauses, show),
      extra: {cells: {total: cells, identicalRdkit: sameCells.rdkit, identicalRdkit2026: sameCells.rdkit2026,
        causes: Object.fromEntries(cellCauses)}, sameDecomposition: decompSame},
      notes: [`Series identical as text (column order included): ${cmp.identical} to RDKit mode, ` +
        `${cmp.identicalRefs.rdkit2026} to RDKit.js 2026.03.6; the same decomposition with the columns matched by name ` +
        `(every R-group SMILES and highlight): ${decompSame.rdkit} of ${n} to RDKit mode, ${decompSame.rdkit2026} to ` +
        'RDKit.js 2026.03.6. The rest differ in column order only, every series with 10 or more R groups: RDKit ' +
        '(MinimalLib\'s std::map, both releases and Python 2026.3.6) orders the labels as strings (Core, Mol, R1, R10, ' +
        'R2, ...), Crux (crux-wasm rgroup.rs rdkit_label_order) by number with Mol last (Core, R1, ..., R10, Mol). In ' +
        'Chem the R-group columns are added in that order, the trellis\'s Y axis defaults to the second R column ' +
        '(R10 in RDKit mode, R2 in Crux mode) and the highlight colours go by column position.',
      `Cells (R-group SMILES per molecule, columns matched by name): ${sameCells.rdkit} of ${cells} identical to RDKit ` +
        `mode (${pct(sameCells.rdkit, cells)}), ${sameCells.rdkit2026} to RDKit.js 2026.03.6 ` +
        `(${pct(sameCells.rdkit2026, cells)}). Output per series: column names, R-group SMILES, the highlight atoms ` +
        'and bonds.', ...trapNote(raw.crux, 'worker')]});
  row.refs.push({engine: 'rdkit2026, columns matched by name', count: decompSame.rdkit2026, of: n,
    pct: pct(decompSame.rdkit2026, n)}, {engine: 'rdkit (RDKit mode), columns matched by name', count: decompSame.rdkit,
    of: n, pct: pct(decompSame.rdkit, n)});
  pushRow(row);
}};

FUNCTIONS.fmcs = {engines: ['rdkit', 'crux', 'rdkit2026'], async compare(rows) {
  const raw = Object.fromEntries(['rdkit', 'crux', 'rdkit2026'].map((e) => [e, loadRaw('fmcs', e)]));
  const clusters = readJson(inputs(`fmcs-clusters-${rows}.json`));
  const n = raw.crux.n;
  const r6 = await loadRdkit2026();
  const sameQuery = (a, b) => {
    const [qa, qb] = [a, b].map((s) => {
      try {
        return r6.get_qmol(s);
      } catch {
        return null;
      }
    });
    try {
      return !!qa && !!qb && qa.get_num_atoms() === qb.get_num_atoms() && qa.get_num_bonds() === qb.get_num_bonds() &&
        qa.get_substruct_match(qb) !== '{}' && qb.get_substruct_match(qa) !== '{}';
    } finally {
      qa?.delete();
      qb?.delete();
    }
  };
  const smarts = (e) => raw[e].items.map((s) => JSON.parse(s)[0]);
  const smiles = (e) => raw[e].items.map((s) => JSON.parse(s)[1]);
  const [cs, rs] = [smarts('crux'), smarts('rdkit')];
  const cmp = compareItems({crux: cs, base: rs, refs: {rdkit2026: smarts('rdkit2026')}, n,
    fallback: new Set(raw.crux.fallbacks), traps: new Set(raw.crux.traps), classify: (k) => !rs[k] && cs[k] ?
      'RDKit timed out (45 s: no MCS), Crux finished within the limit' : sameQuery(cs[k], rs[k]) ?
        'the same query, its atoms written in another order' : 'another MCS'});
  const cmpSmiles = compareItems({crux: smiles('crux'), base: smiles('rdkit'), refs: {rdkit2026: smiles('rdkit2026')}, n});
  const show = (k) => `cluster ${k} (${clusters[k].rows.length} molecules, scaffold ${clip(clusters[k].scaffold, 60)}): ` +
    `Crux ${clip(cs[k], 150)} | RDKit 2024.09 ${clip(rs[k], 150)}`;
  const mN = raw.crux.molecules;
  pushRow(baseRow('FMCS (Cluster MCS)', `${n} clusters, ${mN} molecules, exact atoms and bonds, 45 s per cluster`, n,
    'cluster (µs/item per cluster)', ['rdkit', 'crux', 'rdkit2026'].map((e) => engineStats(e, raw[e], n,
      {maxClusterMs: Math.max(...raw[e].times)})), cmp,
    {causes: await withExamples(cmp.causes, show), refCauses: await withExamples(cmp.refCauses, show),
      extra: {smiles: {identicalRdkit: cmpSmiles.identical, identicalRdkit2026: cmpSmiles.identicalRefs.rdkit2026}},
      notes: [`Output: FMCS's SMARTS per cluster (compared), and the SMILES Cluster MCS shows (clusterMcsSmiles): ` +
        `${cmpSmiles.identical} of ${n} identical to RDKit mode, ${cmpSmiles.identicalRefs.rdkit2026} to RDKit.js ` +
        '2026.03.6.', ...trapNote(raw.crux, 'worker')]}));
}};

const quantile = (v, p) => {
  if (!v.length)
    return null;
  const s = [...v].sort((a, b) => a - b);
  return s[Math.min(s.length - 1, Math.floor(s.length * p))];
};
FUNCTIONS.mces = {engines: ['crux', 'py', 'rdkit (FMCS)', 'crux (FMCS)'], async run(rows) {
  child('mces', 'crux', rows);
  writeJson(inputs('mces-options.json'), loadRaw('mces', 'crux').options);
  python(['mces', inputs(`mces-pairs-${rows}.json`), inputs('mces-options.json'), rawFile('mces', 'py')]);
  for (const engine of ['rdkit', 'crux']) {
    try {
      child('mces-fmcs', engine, rows);
    } catch (e) {
      console.log(`  ${e.message}: the other engines go on`);
    }
  }
}, async compare(rows) {
  const cx = loadRaw('mces', 'crux');
  const py = loadRaw('mces', 'py');
  const pairs = readJson(inputs(`mces-pairs-${rows}.json`));
  const n = cx.n;
  const strip = (o) => JSON.stringify(o.error ? {error: o.error} : o.results);
  let same = 0;
  const causes = new Map();
  for (let i = 0; i < n; i++) {
    const [a, b] = [cx.items[i], py.items[i]];
    if (strip(a) === strip(b)) {
      same++;
      continue;
    }
    const ra = a.results?.[0];
    const rb = b.results?.[0];
    const cause = !ra || !rb ? 'one gives no MCES (below the threshold)' : ra.timedOut || rb.timedOut ?
      'a run timed out (5 s)' : ra.smarts !== rb.smarts ? 'another MCES' : 'the same SMARTS, other matches or numbers';
    const c = causes.get(cause) ?? {count: 0, rows: []};
    c.count++;
    if (c.rows.length < 3)
      c.rows.push(i);
    causes.set(cause, c);
  }
  const show = (i) => `${clip(pairs[i].smiles1, 60)} / ${clip(pairs[i].smiles2, 60)}: Crux ${clip(strip(cx.items[i]), 150)} ` +
    `| Python ${clip(strip(py.items[i]), 150)}`;
  const tCx = cx.items.map((o) => o.s).filter((x) => x !== undefined);
  const tPy = py.items.map((o) => o.s).filter((x) => x !== undefined);
  const [fr, fc] = ['rdkit', 'crux'].map((e) => loadRaw('mces-fmcs', e) ?? {ms: null, times: [], counts: {
    'no run (its child failed)': n}, items: []});
  const fmcsEngine = (e, raw, label) => ({engine: e, label, ms: raw.ms, usPerItem: raw.times.length ?
    us(raw.ms, raw.times.length) : null, medianMs: quantile(raw.times, 0.5) ?? 0, p95Ms: quantile(raw.times, 0.95) ?? 0,
  maxMs: raw.times.length ? Math.max(...raw.times) : 0, counts: raw.counts});
  const fmcsNotes = [['RDKit.js 2024.09', fr], ['Crux', fc]].flatMap(([name, raw]) => [
    raw.pairsRun !== undefined && raw.pairsRun < n ? `${name}'s FMCS ran on the first ${raw.pairsRun} pairs only (a time ` +
      'budget; its time and per-pair figures are over those).' : null,
    raw.threw?.length ? `${name}'s FMCS threw on ${raw.threw.length} pairs (time and per-pair figures: the other pairs): ` +
      `${raw.threw.slice(0, 3).map((t) => `pair ${t.pair} (${clip(pairs[t.pair].smiles1, 50)} / ` +
        `${clip(pairs[t.pair].smiles2, 50)}): ${t.error}, heap ${t.heapGb} GB`).join('; ')}.` : null,
    raw.crashed?.length ? `${name}'s process died on ${raw.crashed.length} pairs (${raw.crashed.slice(0, 3).map((c) =>
      `${c.key}: exit ${c.exit}${c.signal ? ` ${c.signal}` : ''}, heap ${c.heapGb} GB${c.error ? `, ${clip(c.error, 120)}` :
        ''}`).join('; ')}).` : null,
    raw.restarts?.length ? `${name}: ${raw.restarts.length} fresh processes (its wasm heap past ` +
      `${HEAP_RESTART_BYTES / 2 ** 30} GB between two pairs: ${raw.restarts.slice(0, 4).map((r) => `${r.heapGb} GB ` +
        `before ${r.before}`).join(', ')}${raw.restarts.length > 4 ? ', …' : ''}).` : null,
    raw.aborts?.count ? `${name}: ${raw.aborts.count} RDKit.js aborts.` : null]).filter(Boolean);
  const outcome = (o) => o.error ? `error ${o.error}` : !o.results?.length ? 'no MCES (below the threshold)' :
    o.results[0].timedOut ? 'MCES, timed out (5 s: the best so far)' : 'MCES';
  py.counts = tally(py.items.map(outcome));
  cx.counts = {...tally(cx.items.map(outcome)), ...Object.fromEntries(Object.entries(cx.counts).map(([k, v]) =>
    [`Chem job: ${k}`, v]))};
  pushRow({function: 'Rascal MCES (activity cliff MCS)', variant: `${n} similar pairs (Morgan Tanimoto ≥ 0.6), ` +
    `options ${JSON.stringify(cx.options)}`, n, unit: 'pair', engines: [
    engineStats('py', py, n, {medianMs: 1000 * quantile(tPy, 0.5), p95Ms: 1000 * quantile(tPy, 0.95)}),
    engineStats('crux', cx, n, {medianMs: 1000 * quantile(tCx, 0.5), p95Ms: 1000 * quantile(tCx, 0.95),
      chemJobMs: cx.jobMs, chemJobUsPerPair: us(cx.jobMs, n), counts: cx.counts}),
    fmcsEngine('rdkit', fr, `${LABEL.rdkit}: FMCS of the pair (RDKit mode's cliff MCS)`),
    fmcsEngine('crux', fc, `${LABEL.crux}: FMCS of the pair (Crux mode's fallback)`)],
  identical: null, refs: [{engine: 'py (Rascal FindMCES, every field)', count: same, of: n, pct: pct(same, n)}],
  causes: await withExamples(causes, show), extra: {fmcsIdenticalCruxRdkit: fc.items.filter((s, i) => s === fr.items[i] &&
    s !== 'THREW' && s !== 'CRASHED' && s !== 'NOT RUN').length},
  notes: ['RDKit mode gives an activity cliff FMCS (one connected MCS); Crux mode Rascal\'s MCES (the most bonds ' +
    'shared, possibly in fragments): not the same algorithm. Rascal is compared with Python RDKit field by field ' +
    '(SMARTS, atom and bond matches, similarities bit for bit); FMCS is timed on the same pairs (60 s guard per pair, ' +
    'RDKit mode has none). Crux\'s "Chem job" time: crux.worker.ts pairUncommonAtomsAndBonds (read, get_mces, the ' +
    'atoms and bonds outside it).', ...fmcsNotes, ...trapNote(cx, 'worker')]});
}};

FUNCTIONS.molzip = {engines: ['crux', 'rdkit2026', 'py'], async run(rows) {
  child('molzip', 'crux', rows);
  child('molzip', 'rdkit2026', rows);
  python(['molzip', inputs(`molzip-inputs-${rows}.json`), rawFile('molzip', 'py')]);
}, async compare(rows) {
  const [cx, r6, py] = ['crux', 'rdkit2026', 'py'].map((e) => loadRaw('molzip', e));
  const list = readJson(inputs(`molzip-inputs-${rows}.json`));
  const n = cx.n;
  const trapSet = new Set(cx.traps);
  const cmp = compareItems({crux: cx.items, base: null, refs: {rdkit2026: r6.items, py: py.smiles}, n,
    traps: trapSet, classify: (i) => trapSet.has(i) ?
      (/=[CN][\/\\]\[\*:\d+\]\)[\/\\]\[\*:\d+\]|\[\*:\d+\][\/\\][CN]=[CN][\/\\]\[\*:\d+\]/.test(list[i].smiles) ?
        'a stereo double bond between two attachment points' : null) :
      /\[N@/.test(cx.items[i] ?? '') && !/\[N@/.test(r6.items[i] ?? '') ? 'a stereo N Crux keeps (as Python RDKit)' :
        cx.items[i] && r6.items[i] && stripStereo(cx.items[i]) === stripStereo(r6.items[i]) ? 'stereo marks only' : null});
  const parse = loadRaw('parse', 'crux');
  const roundTrip = parse ? list.filter((x, i) => cx.items[i] && cx.items[i] === parse.items[x.row]).length : null;
  const show = (i) => `${clip(list[i].smiles, 90)} → Crux ${clip(cx.items[i], 90)} | RDKit.js 2026.03.6 ` +
    `${clip(r6.items[i], 90)} | Python ${clip(py.smiles[i], 90)}`;
  pushRow(baseRow('molzip', `${n} inputs (each molecule cut at 1 or 2 acyclic bonds, [*:n] labels)`, n, 'input',
    [engineStats('crux', cx, n), engineStats('rdkit2026', r6, n), engineStats('py', py, n,
      {counts: tally(py.smiles.map((s) => s ? 'joined' : 'refused or unread'))})], cmp,
    {causes: [], refCauses: await withExamples(cmp.refCauses, show),
      notes: ['RDKit mode = Crux (Chem\'s RDKit.js 2024.09 has no molzip). Crux: Chem\'s molzip job (one molecule of ' +
        'fragments); RDKit.js 2026.03.6: molzip(first piece, the others); Python: Chem.molzip(mol); each result read ' +
        'again and written as canonical SMILES. ' + (roundTrip !== null ? `Round trip: ${roundTrip} of ${n} give back ` +
        'the molecule\'s own canonical SMILES (Crux).' : ''), ...trapNote(cx, 'molzip')]}));
}};

// -- report

function mdTable() {
  const fmtMs = (e) => e.ms === null || e.ms === undefined ? '-' :
    `${e.ms.toFixed(0)} ms (${e.usPerItem === null ? '-' : e.usPerItem >= 100 ? e.usPerItem.toFixed(0) : e.usPerItem.toFixed(1)} µs)`;
  const counts = (e) => e.counts ? Object.entries(e.counts).map(([k, v]) => `${k} ${v}`).join(', ') : '-';
  const lines = ['| Function | n | Engine: total (per item) | Read / failed / threw | Identical to RDKit mode | ' +
    'Identical to reference | Differences by cause |', '|---|---|---|---|---|---|---|'];
  for (const r of ROWS_OUT) {
    const more = (e) => [e.label?.includes(': ') ? e.label.split(': ').slice(1).join(': ') : null,
      e.fallbacks ? `${e.fallbacks} fallbacks` : null, e.traps ? `${e.traps} traps` : null,
      e.medianMs !== undefined ? `median ${e.medianMs.toFixed(1)} ms, p95 ${e.p95Ms.toFixed(1)} ms` +
        (e.maxMs !== undefined ? `, max ${e.maxMs.toFixed(0)} ms` : '') : null,
      e.maxClusterMs !== undefined ? `slowest ${e.maxClusterMs.toFixed(0)} ms` : null,
      e.prepareMs !== undefined ? `prepare ${e.prepareMs.toFixed(0)} ms (${e.prepareUsPerMolecule.toFixed(1)} µs/molecule)` :
        null, e.chemJobMs !== undefined ? `Chem's job ${e.chemJobMs.toFixed(0)} ms (${e.chemJobUsPerPair.toFixed(0)} µs)` :
        null].filter(Boolean).join(', ');
    const engines = r.engines.map((e) => `**${e.engine}** ${fmtMs(e)}${more(e) ? ` (${more(e)})` : ''}`).join('<br>');
    const read = r.engines.map((e) => `${e.engine}: ${counts(e)}`).join('<br>');
    const ident = r.identical ? `${r.identical.count} / ${r.identical.of} (${r.identical.pct})${r.identical.unit ?
      ` ${r.identical.unit}` : ''}` : 'n/a';
    const refs = r.refs.map((x) => `${x.engine}: ${x.count} / ${x.of} (${x.pct})`).join('<br>') || '-';
    const causes = [...(r.causes ?? []), ...(r.refCauses ?? [])].map((c) => `${c.cause}: ${c.count}`).join('<br>') || '-';
    lines.push(`| ${r.function}${r.variant ? ` — ${r.variant}` : ''} | ${r.n} ${r.unit ?? ''} | ${engines} | ${read} | ` +
      `${ident} | ${refs} | ${causes} |`);
  }
  const details = ['', '## Examples and notes', ''];
  for (const r of ROWS_OUT) {
    details.push(`### ${r.function}${r.variant ? ` — ${r.variant}` : ''}`, '');
    for (const note of r.notes ?? [])
      details.push(`- ${note}`);
    for (const c of [...(r.causes ?? []), ...(r.refCauses ?? [])]) {
      details.push(`- **${c.cause}** (${c.count})`);
      for (const ex of c.examples ?? [])
        details.push(`  - \`${String(ex).replace(/`/g, '\'')}\``);
    }
    if (r.extra)
      details.push('', '```', JSON.stringify(r.extra, null, 1).slice(0, 6000), '```');
    details.push('');
  }
  return lines.join('\n') + '\n' + details.join('\n');
}

// ---- entry ---------------------------------------------------------------------------------------------------------

const isMain = process.argv[1] && resolvePath(process.argv[1]) === fileURLToPath(import.meta.url);
if (!isMain) {
  // imported (a probe of one call): nothing runs
} else if (process.argv[2] === '--child') {
  const [fn, engine, rows] = process.argv.slice(3);
  // a Rust panic in the Crux file (console_error_panic_hook writes it, then the call traps): counted, not printed
  const consoleError = console.error;
  console.error = (...a) => {
    const text = a.map(String).join(' ');
    if (text.includes('panicked at')) {
      PANICS.push(text.split('\n').slice(0, 2).join(' ').replace(/\s+/g, ' ').slice(0, 200));
      PANIC_CTX.push(globalThis.__panicContext?.() ?? null);
    } else if (text.startsWith('Aborted(')) { // RDKit.js's abort(): counted with the run's context, the first few shown
      ABORTS.push({ctx: globalThis.__panicContext?.() ?? null, heap: heapBytes(), what: text.slice(0, 80)});
      if (ABORTS.length <= 3)
        consoleError(`${text} (heap ${gb(heapBytes())} GB, ${JSON.stringify(globalThis.__panicContext?.() ?? null)})`);
    } else
      consoleError(...a);
  };
  // wasm-bindgen frees a Crux object when its JS wrapper is collected (FinalizationRegistry): an object a panic left
  // borrowed (a reaction whose run_reactants panicked, made again here) traps there, outside any call ("memory access
  // out of bounds" in ReactionFinalization), which would end the process: counted instead (the object is not reused)
  const finalizerTraps = [];
  process.on('uncaughtException', (e) => {
    if (e instanceof WebAssembly.RuntimeError && /Finalization/.test(e.stack ?? '')) {
      finalizerTraps.push(String(e.stack).split('\n').find((l) => /Finalization/.test(l))?.trim().slice(0, 80));
      return;
    }
    consoleError(e?.stack ?? e);
    process.exit(1);
  });
  const res = await RUN[fn](engine, parseInt(rows));
  res.panics = {count: PANICS.length, messages: tally(PANICS)};
  if (finalizerTraps.length)
    res.finalizerTraps = {count: finalizerTraps.length, where: tally(finalizerTraps)};
  res.aborts ??= {count: ABORTS.length, heapGb: ABORTS.length ? gb(Math.max(...ABORTS.map((x) => x.heap))) : null};
  res.heapGb = gb(heapBytes());
  if (PANICS.length)
    console.log(`  ${fn}/${engine}: ${PANICS.length} Crux panics: ${JSON.stringify(tally(PANICS))}`);
  if (ABORTS.length)
    console.log(`  ${fn}/${engine}: ${ABORTS.length} RDKit.js aborts in this process`);
  saveRaw(fn, engine, res);
  process.exit(0);
} else {
  const args = process.argv.slice(2);
  const rows = parseInt(args.find((a) => /^\d+$/.test(a)) ?? '100000');
  const only = args.find((a) => a.startsWith('only='))?.slice(5).split(',') ?? Object.keys(FUNCTIONS);
  const compareOnly = args.includes('--compare-only');
  mkdirSync(RAW, {recursive: true});
  if (only.some((fn) => !fn.startsWith('searchtypes'))) // the search types make their inputs themselves
    await prepareInputs(rows);
  const resultsFile = join(OUT, 'results.json');
  const previous = existsSync(resultsFile) ? readJson(resultsFile) : {rows: []};
  const failed = new Set();
  for (const fn of only) {
    const f = FUNCTIONS[fn];
    if (!f)
      throw new Error(`unknown function ${fn}`);
    console.log(`== ${fn}`);
    const t0 = performance.now();
    // OUT/skip-run.txt: functions whose run is skipped (their outputs made by another invocation), compared only
    const skipRun = !process.env.FORCE_RUN && existsSync(join(OUT, 'skip-run.txt')) && readFileSync(join(OUT,
      'skip-run.txt'), 'utf8').split(/\s+/).includes(fn);
    if (!compareOnly && !skipRun) {
      if (f.run)
        await f.run(rows);
      else {
        for (const engine of f.engines) {
          try {
            child(fn, engine, rows);
          } catch (e) { // one engine's failure ends only its own run
            console.log(`  ${e.message}: the other engines go on`);
          }
        }
      }
    }
    const before = ROWS_OUT.length;
    try {
      await f.compare(rows);
    } catch (e) { // a comparison that fails keeps the function's previous rows (and says so)
      ROWS_OUT.splice(before);
      failed.add(fn);
      console.log(`== ${fn}: the comparison failed (${String(e?.stack ?? e).split('\n').slice(0, 3).join(' ')}): ` +
        'its previous rows kept');
      continue;
    }
    for (const r of ROWS_OUT.slice(before))
      r.key = fn;
    console.log(`== ${fn} done in ${((performance.now() - t0) / 1000).toFixed(0)} s`);
  }
  // results accumulate: a function run again replaces its rows
  const done = new Set(only.filter((fn) => !failed.has(fn)));
  const kept = previous.rows.filter((r) => !done.has(r.key));
  const order = Object.keys(FUNCTIONS);
  ROWS_OUT.splice(0, ROWS_OUT.length, ...[...kept, ...ROWS_OUT].sort((a, b) => order.indexOf(a.key) - order.indexOf(b.key)));
  const header = {dataset: DATASET, rows, node: process.version, date: new Date().toISOString(),
    crux: readFileSync(join(CHEM_DIR, 'vendor', 'chem-wasm', 'VERSION.chem-wasm'), 'utf8').match(/crux_core_rev: (\w+)/)[1]};
  writeFileSync(resultsFile, JSON.stringify({...header, rows: ROWS_OUT}, null, 1));
  writeFileSync(join(OUT, 'results.md'), `# Chem: Crux vs RDKit mode on ${rows} molecules\n\n` +
    `Dataset ${DATASET}, ${header.date}, Node ${header.node}, Crux ${header.crux}. Times: one thread, warm ` +
    `(module load and compile excluded; the first ${WARMUP} items run once before timing, excluded).\n\n` + mdTable() +
    // OUT/investigation.md: the causes followed up by hand (probes), appended as they are
    (existsSync(join(OUT, 'investigation.md')) ? '\n' + readFileSync(join(OUT, 'investigation.md'), 'utf8') : ''));
  console.log(`results: ${resultsFile}, ${join(OUT, 'results.md')}`);
}
