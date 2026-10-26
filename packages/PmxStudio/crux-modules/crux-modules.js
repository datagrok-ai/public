// Crux module loader ([decision:wasm-modules]): crux-wasm is built once per domain (chem, ml, qm,
// pkpd, ...; crates/crux-wasm/src/modules.rs) and a page downloads only the modules it calls.
// `call(name, table, params, opts)` looks the function up in manifest.json (generated from the
// crux-fn registry, checked by `cargo test -p crux-wasm --test modules`), instantiates its module
// once (in a dedicated Web Worker where there are workers, else in this thread; a module of the
// page's thread, manifest `thread: main` (chem), in this thread on a page) and calls it
// through the leaf ABI ([decision:wasm-leaf-abi]): columns in (typed arrays, plain arrays,
// Datagrok raw columns, `{columns, meta}`: ./crux-leaf.js), the call's result out as columns
// (`{fn, rows, main: {rows, columns, col}, named, meta}`), no Arrow IPC on either side.
// `callFn(name, ipc, params)` (Arrow IPC bytes in, the call's outputs envelope out: a CXF1 frame,
// ./crux-outputs.js decodes it) needs a validation build ([decision:wasm-leaf-abi] P1c). Modules never share memory: resources (`fnResource`)
// are replayed into every module that loads. An isolated function (manifest `isolated`: crux-fn
// `Footprint::Isolated`, e.g. sql.query; [decision:sql-placement]) never runs on a page's main
// thread: where its module would be in-thread (`worker: false`, no Worker), it gets a worker of
// its own, or the call is refused (an in-thread `load`, or no Worker on a page).
//
//   import { CruxModules } from './crux-modules/crux-modules.js';
//   import { mainIpc } from './crux-modules/crux-outputs.js';
//   const cm = await CruxModules.fromBase(new URL('crux-modules/', location.href));
//   const r = await cm.call('stats.ols', { y, x }, { formula: 'y ~ x' });   // y, x: Float64Array
//   r.main.col.estimate                                                    // a Float64Array
//
// Node / custom hosts pass `load(module)` returning the initialised wasm-bindgen exports (in-thread):
//   new CruxModules(manifest, { load: async (m) => require(`${dir}/${m}/crux_wasm.js`) })

// (crux-leaf.js ships beside this file; in the source tree it is crates/crux-leaf/js/)
const { Leaf, table: leafTable, decodeResult, transferables } = await import('./crux-leaf.js').catch(() => import('../../crux-leaf/js/crux-leaf.js'))

const enc = (v) => (v === undefined || v === null ? '' : typeof v === 'string' ? v : JSON.stringify(v))

/** The module a function runs in (`undefined`: no module has it). Case-insensitive, as the registry. */
export function moduleFor(manifest, name) {
  const f = manifest.functions
  if (f[name] !== undefined) return f[name]
  if (!manifest._lower) manifest._lower = Object.fromEntries(Object.entries(f).map(([k, m]) => [k.toLowerCase(), m]))
  return manifest._lower[String(name).trim().toLowerCase()]
}

/** The URL of a module's file (`kind`: 'js' or 'wasm') under `base`, with its content hash
 *  (manifest `hash`, written by modules.py) as `?v=`: an immutable URL a server may let the browser
 *  cache for good (web/dev/static-delivery.ts). The glue fetches its wasm by its own hash. */
export function moduleUrl(manifest, module, base, kind = 'js') {
  const m = manifest.modules[module]
  const u = new URL(m[kind], base)
  const v = m.hash?.[kind]
  if (v) u.searchParams.set('v', v)
  return u.href
}

/** In-thread handle over a module's wasm-bindgen exports. */
function directHandle(exports) {
  const leaf = typeof exports.leafAbi === 'function' ? new Leaf(exports) : null
  return {
    exports,
    leaf,
    call: async (op, args) => {
      switch (op) {
        case 'leafCall':
          return leaf.call(args.name, args.table, args.params, args.opts)
        case 'leafCallAsync':
          return leaf.callAsync(args.name, args.table, args.params, args.opts)
        case 'callFn':
          if (typeof exports.callFn !== 'function') throw Object.assign(new Error(`callFn is in validation builds only ([decision:wasm-leaf-abi] P1c: modules.py build --validation; node packages have it): call the module through the leaf ABI (call)`), { code: 'no_callfn' })
          return args.opts ? exports.callFn(args.name, args.ipc, args.params, args.opts) : exports.callFn(args.name, args.ipc, args.params)
        case 'fnDescribe':
          return exports.fnDescribe(args.name)
        case 'fnList':
          return exports.fnList()
        case 'fnResource':
          return exports.fnResource(args.name, args.id, args.bytes)
        case 'fnDataset':
          return exports.fnDataset(args.name, args.manifest, args.base)
        case 'cruxModules':
          return exports.cruxModules()
        default:
          throw new Error(`unknown op ${op}`)
      }
    },
    terminate: () => {},
  }
}

/** Worker handle: crux-module-worker.js instantiates the module in its own thread. */
function workerHandle(workerUrl, jsUrl) {
  const w = new Worker(workerUrl, { type: 'module' })
  let next = 0
  const pending = new Map()
  w.onmessage = (e) => {
    const { id, ok, result, error, code, progress } = e.data
    const p = pending.get(id)
    if (!p) return
    if (progress) return p.onProgress?.(...progress)
    pending.delete(id)
    ok ? p.resolve(result) : p.reject(Object.assign(new Error(error), { code }))
  }
  w.onerror = (e) => {
    const err = new Error(`crux module worker (${jsUrl}): ${e.message || 'failed'}`)
    for (const p of pending.values()) p.reject(err)
    pending.clear()
  }
  const call = (op, args = {}, transfer = [], onProgress) =>
    new Promise((resolve, reject) => {
      const id = next++
      pending.set(id, { resolve, reject, onProgress })
      w.postMessage({ id, op, args }, transfer)
    })
  return { call, terminate: () => w.terminate(), ready: call('init', { js: jsUrl }) }
}

export class CruxModules {
  /**
   * `manifest`: manifest.json. Options: `base` (URL of the directory holding manifest.json, the
   * worker script and <module>/crux_wasm.js), `load(module)` (in-thread loader, overrides the rest),
   * `worker` (false: never use a Worker).
   */
  constructor(manifest, options = {}) {
    this.manifest = manifest
    this.options = options
    this.handles = new Map()
    this.resources = []
  }

  /** Fetches `<base>/manifest.json`. */
  static async fromBase(base, options = {}) {
    const b = new URL(String(base).endsWith('/') ? String(base) : `${base}/`, globalThis.location?.href)
    // (revalidated every time: it names the modules' current hashes)
    const manifest = await (await fetch(new URL('manifest.json', b), { cache: 'no-cache' })).json()
    return new CruxModules(manifest, { ...options, base: b })
  }

  /** Every function name the modules serve, with its module. */
  functions() {
    return { ...this.manifest.functions }
  }

  moduleFor(name) {
    return moduleFor(this.manifest, name)
  }

  /** The modules loaded (or loading) so far. */
  loaded() {
    return [...this.handles.keys()]
  }

  /** The handle of `module`, instantiated on first use (one per module, kept). `isolated`: for a
   *  function that never runs on a page's main thread (a worker of its own where the module's
   *  handle would be in-thread there). */
  handle(module, isolated = false) {
    if (!this.manifest.modules[module]) return Promise.reject(new Error(`unknown crux module ${module}`))
    const own = isolated && this.inThread()
    const key = own ? `${module}#isolated` : module
    let h = this.handles.get(key)
    if (!h) {
      h = this._open(module, own)
      this.handles.set(key, h)
      h.catch(() => this.handles.delete(key))
    }
    return h
  }

  /** Whether this loader's modules run on a page's main thread. */
  inThread() {
    const { load, worker } = this.options
    return typeof document !== 'undefined' && (!!load || worker === false || typeof Worker === 'undefined')
  }

  /** Whether `name` is isolated (never on a page's main thread; manifest `isolated`). */
  isolated(name) {
    const iso = this.manifest.isolated
    if (!Array.isArray(iso) || !iso.length) return false
    if (!this.manifest._isolated) this.manifest._isolated = new Set(iso.map((n) => n.toLowerCase()))
    return this.manifest._isolated.has(String(name).trim().toLowerCase())
  }

  async _open(module, ownWorker = false) {
    const { load, base, worker } = this.options
    let h
    if (ownWorker && (load || !base || typeof Worker === 'undefined')) {
      throw new Error(`crux module ${module}: its isolated functions never run on the page's main thread, and this loader has no worker for them (pass \`base\` where Worker exists)`)
    }
    if (load && !ownWorker) {
      h = directHandle(await load(module))
      ;(this.loadedSync ??= new Map()).set(module, h)
    } else {
      if (!base) throw new Error('CruxModules: pass `base` (or use CruxModules.fromBase) or `load`')
      const jsUrl = moduleUrl(this.manifest, module, base)
      // a module of the page's thread (manifest `thread: main`: chem alone, [decision:wasm-leaf-abi]
      // step 3) loads in-thread on a page: its calls are synchronous (`callSync`), its exports direct
      const mainThread = this.manifest.modules[module]?.thread === 'main' && typeof document !== 'undefined'
      if (ownWorker || (!mainThread && worker !== false && typeof Worker !== 'undefined')) {
        h = workerHandle(new URL('crux-module-worker.js', base).href, jsUrl)
        await h.ready
      } else {
        const exports = await import(/* @vite-ignore */ jsUrl)
        if (typeof exports.default === 'function') await exports.default()
        h = directHandle(exports)
        ;(this.loadedSync ??= new Map()).set(module, h)
      }
    }
    for (const r of this.resources) await h.call(r.op || 'fnResource', r)
    return h
  }

  /** Loads `module` without calling anything (warm-up / prefetch). */
  async load(module) {
    await this.handle(module)
  }

  /** The in-thread wasm-bindgen exports of `module` (its own API: `CollectionBuilder`, `get_mol`,
   *  ...); needs `load` or `worker: false`. */
  async exports(module) {
    const h = await this.handle(module)
    if (!h.exports) throw new Error(`crux module ${module} runs in a worker: construct CruxModules with { worker: false } for its exports`)
    return h.exports
  }

  /** Runs `name` through its module's leaf ABI: `table` columns (crux-leaf.js's input forms),
   *  `params` an object or JSON text, `opts` the call's (`bind`, `strict`, `error_column`,
   *  `variant`, `exact`; `onProgress(done, total, unit, phase, message)`; `transfer: true` moves
   *  the input's buffers to a module Worker instead of copying them). Resolves to the result
   *  (crux-leaf.js decodeResult: `{fn, rows, main: {rows, columns, col}, named, meta}`); a call
   *  that awaits the device (a WebGPU variant) runs as `callAsync`. Errors carry `code`. */
  async call(name, table, params, opts) {
    return this._leaf('leafCall', name, table, params, opts)
  }

  /** `call` that may await the device: `variant: auto` may pick a WebGPU variant. */
  async callAsync(name, table, params, opts) {
    return this._leaf('leafCallAsync', name, table, params, opts)
  }

  /** `call` returning synchronously, for a module loaded in this thread (`load` / `worker:
   *  false`; the main-thread chemistry): throws when it is not loaded yet or runs in a Worker. */
  callSync(name, table, params, opts) {
    const module = this.moduleFor(name)
    const h = this.loadedSync?.get(module)
    if (!h?.leaf) throw new Error(`crux module ${module} is not loaded in this thread (await load(${JSON.stringify(module)}) with an in-thread loader first)`)
    return h.leaf.call(name, table, params, opts)
  }

  async _leaf(op, name, table, params, opts = {}) {
    const module = this.moduleFor(name)
    if (module === undefined) throw Object.assign(new Error(`unknown function ${JSON.stringify(name)} (no crux module serves it)`), { code: 'unknown_fn' })
    const h = await this.handle(module, this.isolated(name))
    if (h.leaf) {
      try {
        return await h.call(op, { name, table, params, opts })
      } catch (e) {
        // a call that awaits the device: as callAsync (the result is the same)
        if (op === 'leafCall' && e?.code === 'needs_async') return h.call('leafCallAsync', { name, table, params, opts })
        throw e
      }
    }
    if (h.exports) throw new Error(`crux module ${module} has no leaf ABI (rebuild it: tools/parity/wasm-modules/modules.py build)`)
    const { onProgress, onLog, transfer, ...rest } = opts || {}
    const t = leafTable(table)
    const send = (o) => h.call(o, { name, table: t, params: params ?? null, opts: rest, progress: !!onProgress }, transfer ? transferables(t) : [], onProgress)
    try {
      return decodeResult(await send(op))
    } catch (e) {
      // (as in-thread: a call that awaits the device runs async; a transferred input is gone)
      if (op === 'leafCall' && e?.code === 'needs_async' && !transfer) return decodeResult(await send('leafCallAsync'))
      throw e
    }
  }

  /** Runs `name` on Arrow IPC bytes; `params` / `opts` objects or JSON text. Resolves to the call's outputs envelope (crux-outputs.js).
   *  Validation builds only (the module's `callFn`, [decision:wasm-leaf-abi] P1c): what ships is called by `call`. */
  async callFn(name, ipc, params, opts) {
    const module = this.moduleFor(name)
    if (module === undefined) throw new Error(`unknown function ${JSON.stringify(name)} (no crux module serves it)`)
    const h = await this.handle(module, this.isolated(name))
    const o = enc(opts)
    return h.call('callFn', { name, ipc, params: enc(params), opts: o || undefined })
  }

  /** One function's spec (JSON text), from its module. */
  async fnDescribe(name) {
    const module = this.moduleFor(name)
    if (module === undefined) throw new Error(`unknown function ${JSON.stringify(name)}`)
    return (await this.handle(module)).call('fnDescribe', { name })
  }

  /** The registry of `module` (JSON text; it also lists the functions it carries as dependencies). */
  async fnList(module = this.manifest.core) {
    return (await this.handle(module)).call('fnList', {})
  }

  /** A resource (weights, ...) for later calls, in every module (loaded now or later). */
  async fnResource(name, id, bytes) {
    const r = { name, id, bytes }
    this.resources.push(r)
    await Promise.all([...this.handles.values()].map(async (h) => (await h).call('fnResource', r)))
  }

  /** A dataset read in place by byte ranges ([decision:dataset-resource]), in every module:
   *  `manifest` its crux-dataset.json text, `base` a base URL (crux-server `/data/files/<path>`) or
   *  `{path: File | Blob | bytes | url}`. Resolves to its id (`sha256:…`, from the first module).
   *  A module in this thread (no workers) can only read it by URL off the main thread. */
  async fnDataset(name, manifest, base) {
    const r = { op: 'fnDataset', name, manifest, base }
    this.resources.push(r)
    const ids = await Promise.all([...this.handles.values()].map(async (h) => (await h).call('fnDataset', r)))
    return ids[0]
  }

  /** Stops every module's worker; the next call reloads. */
  async terminate() {
    for (const h of this.handles.values()) (await h.catch(() => null))?.terminate()
    this.handles.clear()
  }
}
