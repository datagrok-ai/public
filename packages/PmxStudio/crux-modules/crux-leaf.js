// The leaf ABI, JS side ([decision:wasm-leaf-abi], P0). Plain JS, no dependencies.
//
// A leaf module exports leafAbi / leafInfo / leafCall / leafCallAsync / leafRead / leafDescribe /
// leafFree / leafCancel / leafResource (crates/crux-leaf-wasm). This file turns what a caller has
// into the ABI's canonical columns and the results back into typed arrays:
//
//   import { Leaf } from './crux-leaf.js'
//   const leaf = new Leaf(await import('./cardiac/crux_leaf_wasm.js').then(async (m) => (await m.default(), m)))
//   const r = leaf.call('cardiac.block_hill', { conc: new Float64Array([10]), ic50: [10], h: [1] })
//   r.main.col.block                                      // Float64Array [0.5]
//   const k = leaf.call('cardiac.block_hill', t, {}, { keep: true })      // stays in the leaf
//   leaf.call('cardiac.block_hill', { conc: { handle: k.handle, column: 'factor' }, ic50: [1], h: [1] })
//   leaf.read(k.handle); leaf.free(k.handle)
//
// Column inputs (each becomes {name, format, length, nullCount, offset, buffers: [Uint8Array|null],
// children, dictionary}: the Arrow C Data Interface layout):
//   - a typed array (Float64Array 'g', Float32Array 'f', Int32Array 'i', BigInt64Array 'l', ...);
//   - a plain array of numbers / strings / booleans / bigints / number arrays (null = null);
//   - {format, values, offsets?, validity?, children?, dictionary?, nullValue?}: explicit layout;
//     `nullValue` marks nulls by a sentinel in `values` (Datagrok's raw data: DG.FLOAT_NULL,
//     DG.INT_NULL), `dictionary` (an array of strings or a column) makes `values` category indexes
//     (Datagrok string columns: getRawData() indexes + categories);
//   - {handle, column, of?, offset?, length?}: a column of a kept result of the same leaf (no copy,
//     no transfer; `of: 'input'`: the column sent with that call; offset / length: a zero-copy
//     slice, a chunk of a stashed table, P1a);
//   - already canonical ({format, length, buffers}).
// Results: {fn, handle?, rows, main: {rows, columns, col: {name: values}}, named, meta}: each column
// {name, format, length, nullCount, validity, offsets, values, children, dictionary, meta} with
// typed arrays (`values` of a utf8 column: its UTF-8 bytes; `toArray(column)` gives JS values).

const TYPED = {
  g: Float64Array, f: Float32Array, l: BigInt64Array, L: BigUint64Array, i: Int32Array, I: Uint32Array,
  s: Int16Array, S: Uint16Array, c: Int8Array, C: Uint8Array, tdD: Int32Array, tdm: BigInt64Array, e: Uint16Array,
}
const FORMAT_OF = new Map([
  [Float64Array, 'g'], [Float32Array, 'f'], [BigInt64Array, 'l'], [BigUint64Array, 'L'], [Int32Array, 'i'], [Uint32Array, 'I'],
  [Int16Array, 's'], [Uint16Array, 'S'], [Int8Array, 'c'], [Uint8Array, 'C'], [Uint8ClampedArray, 'C'],
])
const enc = new TextEncoder()
const dec = new TextDecoder()

const bytes = (ta) => (ta instanceof Uint8Array ? ta : new Uint8Array(ta.buffer, ta.byteOffset, ta.byteLength))
const typedOf = (format) => TYPED[format] || (format.startsWith('ts') || format.startsWith('tD') ? BigInt64Array : null)
const view = (u, T) => {
  if (!u) return null
  if (u.byteOffset % T.BYTES_PER_ELEMENT) u = u.slice()     // a copy only for a misaligned buffer
  return new T(u.buffer, u.byteOffset, Math.floor(u.byteLength / T.BYTES_PER_ELEMENT))
}

/** An LSB-first bitmap of `n` bits (set where `bit(i)`), always materialized: a Boolean column's
 *  value bitmap (an all-true column has every bit set, not a null buffer). */
function bitsOf(n, bit) {
  const b = new Uint8Array((n + 7) >> 3)
  for (let i = 0; i < n; i++) if (bit(i)) b[i >> 3] |= 1 << (i & 7)
  return b
}

/** A validity bitmap of `n` rows: null (with 0 nulls) when every row is valid. Not for values. */
function bitmapOf(n, isValid) {
  const b = new Uint8Array((n + 7) >> 3)
  let nulls = 0
  for (let i = 0; i < n; i++) {
    if (isValid(i)) b[i >> 3] |= 1 << (i & 7)
    else nulls++
  }
  return nulls ? { validity: b, nulls } : { validity: null, nulls: 0 }
}

function utf8Column(name, strs) {
  const parts = new Array(strs.length)
  let total = 0
  for (let i = 0; i < strs.length; i++) {
    const s = strs[i]
    if (s === null || s === undefined) continue
    parts[i] = enc.encode(String(s))
    total += parts[i].length
  }
  const offsets = new Int32Array(strs.length + 1)
  const data = new Uint8Array(total)
  let at = 0
  for (let i = 0; i < strs.length; i++) {
    if (parts[i]) { data.set(parts[i], at); at += parts[i].length }
    offsets[i + 1] = at
  }
  const { validity, nulls } = bitmapOf(strs.length, (i) => strs[i] !== null && strs[i] !== undefined)
  return { name, format: 'u', length: strs.length, nullCount: nulls, offset: 0, buffers: [validity, bytes(offsets), data] }
}

function fromArray(name, a) {
  const first = a.find((x) => x !== null && x !== undefined)
  const n = a.length
  if (typeof first === 'string' || first === undefined) return utf8Column(name, a)
  const { validity, nulls } = bitmapOf(n, (i) => a[i] !== null && a[i] !== undefined)
  if (typeof first === 'number') {
    const v = new Float64Array(n)
    for (let i = 0; i < n; i++) v[i] = a[i] ?? 0
    return { name, format: 'g', length: n, nullCount: nulls, offset: 0, buffers: [validity, bytes(v)] }
  }
  if (typeof first === 'bigint') {
    const v = new BigInt64Array(n)
    for (let i = 0; i < n; i++) v[i] = a[i] ?? 0n
    return { name, format: 'l', length: n, nullCount: nulls, offset: 0, buffers: [validity, bytes(v)] }
  }
  if (typeof first === 'boolean') {
    return { name, format: 'b', length: n, nullCount: nulls, offset: 0, buffers: [validity, bitsOf(n, (i) => a[i] === true)] }
  }
  if (Array.isArray(first) || ArrayBuffer.isView(first)) {
    // a list of numbers per row
    const offsets = new Int32Array(n + 1)
    let total = 0
    for (let i = 0; i < n; i++) { total += a[i] ? a[i].length : 0; offsets[i + 1] = total }
    const flat = new Float64Array(total)
    for (let i = 0; i < n; i++) if (a[i]) flat.set(a[i], offsets[i])
    const child = { name: 'item', format: 'g', length: total, nullCount: 0, offset: 0, buffers: [null, bytes(flat)] }
    return { name, format: '+l', length: n, nullCount: nulls, offset: 0, buffers: [validity, bytes(offsets)], children: [child] }
  }
  throw new Error(`column ${name}: unsupported values (${typeof first})`)
}

function fromLayout(name, c) {
  const f = c.format ?? (ArrayBuffer.isView(c.values) ? FORMAT_OF.get(c.values.constructor) : undefined)
  if (!f) throw new Error(`column ${name}: no format (and values not a typed array)`)
  const T = typedOf(f)
  let values = c.values
  if (values && !ArrayBuffer.isView(values) && Array.isArray(values)) values = T ? T.from(values) : values
  const n = c.length ?? (f === 'u' || f === 'z' || f === '+l' ? c.offsets.length - 1 : values?.length ?? 0)
  let validity = c.validity ? bytes(c.validity) : null
  let nullCount = c.nullCount
  if (c.nullValue !== undefined && values) {
    const nv = c.nullValue
    const isNull = nv !== nv ? (x) => x !== x : (x) => x === nv
    const v = bitmapOf(n, (i) => !isNull(values[i]) && (!validity || (validity[i >> 3] >> (i & 7)) & 1))
    validity = v.validity
    nullCount = v.nulls
  }
  const out = { name, format: f, length: n, offset: c.offset || 0, buffers: [validity] }
  if (nullCount !== undefined) out.nullCount = nullCount
  if (f === 'b') out.buffers.push(Array.isArray(values) ? bitsOf(n, (i) => values[i] === true) : values ? bytes(values) : new Uint8Array((n + 7) >> 3))
  else if (f === 'u' || f === 'z' || f === 'U' || f === 'Z') out.buffers.push(bytes(c.offsets), bytes(values))
  else if (f === '+l' || f === '+L') out.buffers.push(bytes(c.offsets))
  else if (f === '+s' || f.startsWith('+w:')) {}
  else out.buffers.push(bytes(values))
  if (c.children) out.children = c.children.map((ch, i) => column(ch.name ?? `f${i}`, ch))
  if (c.dictionary) out.dictionary = column('', c.dictionary)
  if (c.meta) out.meta = c.meta
  return out
}

/** A binary column (`z`) of byte arrays (null: a null), as a canonical column. */
export function bytesColumn(name, arrays) {
  const n = arrays.length
  const offsets = new Int32Array(n + 1)
  for (let i = 0; i < n; i++) offsets[i + 1] = offsets[i] + (arrays[i] ? arrays[i].length : 0)
  const data = new Uint8Array(offsets[n])
  for (let i = 0; i < n; i++) if (arrays[i]) data.set(arrays[i], offsets[i])
  const { validity, nulls } = bitmapOf(n, (i) => arrays[i] !== null && arrays[i] !== undefined)
  return { name, format: 'z', length: n, nullCount: nulls, offset: 0, buffers: [validity, bytes(offsets), data] }
}

/** One input as a canonical column (or a {name, handle, column} reference). */
export function column(name, v) {
  if (v && typeof v === 'object' && !ArrayBuffer.isView(v) && !Array.isArray(v)) {
    if (v.handle !== undefined) {
      const r = { name, handle: v.handle, column: v.column ?? name }
      if (v.of) r.of = v.of
      if (v.offset !== undefined || v.length !== undefined) { r.offset = v.offset ?? 0; r.length = v.length }
      return r
    }
    if (v.buffers && v.format) return { ...v, name }
    if (v.format || ArrayBuffer.isView(v.values)) return fromLayout(name, v)
    // {values: plain array, meta}: a column with its field metadata (semType, units)
    if (Array.isArray(v.values)) return v.meta ? { ...fromArray(name, v.values), meta: v.meta } : fromArray(name, v.values)
  }
  if (ArrayBuffer.isView(v)) {
    const f = FORMAT_OF.get(v.constructor)
    if (!f) throw new Error(`column ${name}: unsupported typed array ${v.constructor.name}`)
    return { name, format: f, length: v.length, nullCount: 0, offset: 0, buffers: [null, bytes(v)] }
  }
  if (Array.isArray(v)) return fromArray(name, v)
  throw new Error(`column ${name}: a typed array, an array, a layout object or {handle, column} expected`)
}

/** A table input: `{name: column input}`, `[columns]`, or `{rows?, columns: [...], meta?}`
 *  (`meta`: the table's metadata, an Arrow schema's; a column's `meta` its field's: what binding
 *  by semantic type and the facets read). */
export function table(t) {
  if (!t) return { columns: [] }
  if (Array.isArray(t)) return { columns: t.map((c) => column(c.name, c)) }
  if (Array.isArray(t.columns)) {
    const o = { rows: t.rows, columns: t.columns.map((c) => column(c.name, c)) }
    if (t.meta) o.meta = t.meta
    return o
  }
  return { columns: Object.entries(t).map(([k, v]) => column(k, v)) }
}

/** A table's columns without their buffers (names, formats, children, metadata): `leafPlan`'s
 *  input, what a check needs. */
export function describe(t) {
  const strip = (c) => {
    const o = { name: c.name, format: c.format }
    if (c.children) o.children = c.children.map(strip)
    if (c.dictionary) o.dictionary = strip(c.dictionary)
    if (c.meta) o.meta = c.meta
    return o
  }
  const tt = t && Array.isArray(t.columns) && t.columns.every((c) => c.format && (c.buffers || c.handle === undefined)) ? t : table(t)
  const o = { columns: tt.columns.map(strip) }
  if (tt.meta) o.meta = tt.meta
  return o
}

/** The ArrayBuffers of a canonical table or result (unique; for postMessage transfer lists). */
export function transferables(t, memory) {
  const seen = new Set()
  const walk = (c) => {
    for (const b of c.buffers || []) if (b && b.buffer !== memory) seen.add(b.buffer)
    for (const ch of c.children || []) walk(ch)
    if (c.dictionary) walk(c.dictionary)
  }
  for (const c of t?.columns || []) walk(c)
  return [...seen].filter((b) => !(typeof SharedArrayBuffer !== 'undefined' && b instanceof SharedArrayBuffer))
}

/** A canonical result column -> typed arrays. */
export function decodeColumn(c) {
  const f = c.format
  const b = c.buffers || []
  const out = { name: c.name, format: f, length: c.length, nullCount: c.nullCount, offset: c.offset || 0, validity: b[0] || null }
  if (c.meta) out.meta = c.meta
  if (f === 'u' || f === 'z') { out.offsets = view(b[1], Int32Array); out.values = b[2] || new Uint8Array(0) }
  else if (f === 'U' || f === 'Z') { out.offsets = view(b[1], BigInt64Array); out.values = b[2] || new Uint8Array(0) }
  else if (f === '+l') out.offsets = view(b[1], Int32Array)
  else if (f === '+L') out.offsets = view(b[1], BigInt64Array)
  else if (f === 'b') out.values = b[1] || null
  else if (typedOf(f)) out.values = view(b[1], typedOf(f))
  if (c.children) out.children = c.children.map(decodeColumn)
  if (c.dictionary) out.dictionary = decodeColumn(c.dictionary)
  return out
}

const valid = (c, i) => !c.validity || (c.validity[(c.offset + i) >> 3] >> ((c.offset + i) & 7)) & 1

/** A decoded column as JS values (numbers, strings, booleans, arrays, null). */
export function toArray(c) {
  const n = c.length, o = c.offset || 0, out = new Array(n)
  for (let i = 0; i < n; i++) {
    if (!valid(c, i)) { out[i] = null; continue }
    const k = o + i
    switch (c.format) {
      case 'u': out[i] = dec.decode(c.values.subarray(c.offsets[k], c.offsets[k + 1])); break
      case 'z': out[i] = c.values.slice(c.offsets[k], c.offsets[k + 1]); break
      case 'b': out[i] = !!((c.values[k >> 3] >> (k & 7)) & 1); break
      case '+l': {
        const ch = toArray(c.children[0])
        out[i] = ch.slice(c.offsets[k], c.offsets[k + 1])
        break
      }
      default: out[i] = c.values ? (typeof c.values[k] === 'bigint' ? Number(c.values[k]) : c.values[k]) : null
    }
  }
  if (c.dictionary) {
    const d = toArray(c.dictionary)
    return out.map((x) => (x === null ? null : d[x]))
  }
  return out
}

export function decodeTable(t) {
  if (!t) return t
  const columns = (t.columns || []).map((c) => (c.buffers ? decodeColumn(c) : c))
  const col = {}
  for (const c of columns) col[c.name] = c.values ?? c
  const out = { rows: t.rows, columns, col }
  // (the table's own metadata: an Arrow schema's)
  if (t.meta) out.meta = t.meta
  return out
}

/** A raw leaf result with typed arrays. */
export function decodeResult(r) {
  const named = {}
  for (const [k, v] of Object.entries(r.named || {})) named[k] = v.kind === 'table' ? { ...v, table: decodeTable(v.table) } : v
  return { ...r, main: decodeTable(r.main), named }
}

/** The JSON options and the JS-side ones of a call. */
// `keep`: false (read out), true (kept, evictable past the cap) or 'pin' (kept, never evicted
// until freed or unpinned: a run's live data)
function splitOpts(opts = {}) {
  const { keep, callId, onProgress, onLog, transfer, ...rest } = opts
  return { keep: keep === 'pin' ? 'pin' : !!keep, callId, onProgress, onLog, json: Object.keys(rest).length ? JSON.stringify(rest) : undefined }
}

// Leaf instance numbers of this page (a handle's high bits): each leaf Worker, and each respawn
// of one, gets the next, so a handle of an instance that is gone is refused by its successor
// instead of naming the successor's objects. A random start keeps tabs apart (29 bits).
let nextInstance = 1 + Math.floor(Math.random() * 2 ** 28)
const newInstance = () => (nextInstance = (nextInstance % (2 ** 29 - 2)) + 1)

/** An in-thread leaf over its wasm-bindgen exports (synchronous `call`). */
export class Leaf {
  constructor(exports) {
    this.x = exports
    if (exports.leafAbi() !== 1) throw new Error(`leaf ABI ${exports.leafAbi()}: this loader speaks 1`)
  }
  info() { return this.x.leafInfo() }
  /** Synchronous; throws an Error with `code` (`needs_async`: use callAsync). */
  call(name, t, params, opts) {
    const o = splitOpts(opts)
    return decodeResult(this.x.leafCall(name, table(t), params ?? null, o.json, o.keep, o.callId, o.onProgress, o.onLog))
  }
  async callAsync(name, t, params, opts) {
    const o = splitOpts(opts)
    return decodeResult(await this.x.leafCallAsync(name, table(t), params ?? null, o.json, o.keep, o.callId, o.onProgress, o.onLog))
  }
  /** Copies a kept result out (fresh, transferable buffers). `which`: `'main'`, `'named'` or
   *  one named output's name (absent: all of them). */
  read(handle, which) { return decodeResult(this.x.leafRead(handle, false, which)) }
  /** Views into wasm memory: read them before calling into the leaf again. */
  view(handle, which) { return decodeResult(this.x.leafRead(handle, true, which)) }
  /** Checks a call without running it: `{fn, want: {columns, meta}, shape}`, or throws its
   *  error (`resources`: names of resources the call will be given besides this leaf's own). */
  plan(name, t, params, opts, resources = [], relaxVariant = false) {
    const o = splitOpts(opts)
    return this.x.leafPlan(name, describe(t), params ?? null, o.json, resources, relaxVariant)
  }
  /** A dataset for range reads (crux-wasm's fnDataset; resolves its id). */
  dataset(name, manifest, base) { return this.x.fnDataset(name, manifest, base) }
  describe(handle) { return this.x.leafDescribe(handle) }
  /** Stashes a table in the leaf (no call): `{handle, rows, main}` described; later calls take its
   *  columns by `{handle, column, offset?, length?}`. */
  put(t, keep = true) { return this.x.leafPut(table(t), keep) }
  /** Pins or unpins a kept result (a pinned one is never evicted). */
  pin(handle, on = true) { return this.x.leafPin(handle, on) }
  /** `{instance, cap}` of the handle table (set from `cfg` where given). */
  config(cfg = {}) { return this.x.leafConfig(cfg) }
  /** Columns of kept results, each the concatenation of its parts ({handle, column, of?, offset?,
   *  length?}), compacted, one copy out: `{columns: [{name, parts}], free?: [handle]}`. */
  gather(spec) { return decodeTable(this.x.leafGather(spec)) }
  free(handle) { return this.x.leafFree(handle) }
  cancel(callId) { return this.x.leafCancel(callId) }
  resource(name, id, bytes) { return this.x.leafResource(name, id, bytes) }
  /** Resource bytes stored by identity: a call sees them only bound by `resources: {name: {id}}`. */
  store(id, bytes, lease) { return this.x.leafStore(id, bytes, lease) }
  /** Drops stored resource `id`; true when it was stored. */
  unstore(id) { return this.x.leafUnstore(id) }
  /** Which of `ids` are not stored (evicted under memory pressure, unstored, never stored); with
   *  `lease`, the others are pinned under it until a call naming it (`opts.lease`) ends. */
  stored(ids, lease) { return this.x.leafStored(ids, lease) }
  /** Ends lease `lease` (its resources evictable again). */
  release(lease) { return this.x.leafRelease(lease) }
  stats() { return this.x.leafStats() }
}

/** Instantiates a leaf's wasm-bindgen web glue (`js`: its URL; `module`: a compiled
 *  WebAssembly.Module to reuse, e.g. on a respawn). */
export async function loadLeaf(js, module) {
  const m = await import(/* @vite-ignore */ js)
  await m.default(module ? { module_or_path: module } : undefined)
  return new Leaf(m)
}

/** A leaf in a Web Worker (crux-leaf-worker.js): the same calls, as Promises. Inputs are
 *  copied to the Worker (a structured clone: the caller keeps its arrays); `transfer: true`
 *  (`call` / `callAsync` opts, `put`'s second argument) moves them instead: every ArrayBuffer
 *  behind an input view is detached here, whole, so every other view over the same buffer (a
 *  subarray, another column in one allocation, the Arrow table it came from) reads as empty
 *  after. Results always move by transfer (they are fresh). `respawn()`: terminate (the cancel of
 *  a synchronous call) and start again from the compiled module, replaying resources. */
export class LeafWorker {
  /** `onFail(error)`: the Worker failed (its script, an import, a trap that ends it): every
   *  request it had is rejected (`code: 'worker_failed'`), its handles are gone; `respawn()`
   *  starts a new one. */
  constructor(workerUrl, js, module, { onFail } = {}) {
    this.workerUrl = workerUrl
    this.js = js
    this.module = module
    this.resources = []
    this.onFail = onFail
    this._start()
  }
  _start() {
    this.w = new Worker(this.workerUrl, { type: 'module' })
    this.next = 1
    this.pending = new Map()
    this.failed = null
    // this instance's number (its handles' high bits); a respawn takes a new one
    this.instance = newInstance()
    const w = this.w
    const fail = (why) => {
      if (this.w !== w || this.failed) return
      this.failed = Object.assign(new Error(`leaf Worker (${this.js}) failed: ${why}`), { code: 'worker_failed' })
      for (const p of this.pending.values()) p.reject(this.failed)
      this.pending.clear()
      this.onFail?.(this.failed)
    }
    w.onerror = (e) => {
      e.preventDefault?.()
      fail(e.message || 'error')
    }
    w.onmessageerror = () => fail('a message could not be read')
    w.onmessage = (e) => {
      if (e.data?.dead) return fail(e.data.dead)
      const { id, ok, result, error, code, progress, log } = e.data
      const p = this.pending.get(id)
      if (!p) return
      if (progress) return p.onProgress?.(...progress)
      if (log) return p.onLog?.(...log)
      this.pending.delete(id)
      if (ok) p.resolve(result)
      else p.reject(Object.assign(new Error(error), { code }))
    }
    this.ready = this._post('init', { js: this.js, module: this.module, instance: this.instance }).then(async (info) => {
      for (const r of this.resources) await this._post(r.op || 'resource', r)
      return info
    })
    // (a failed start rejects `ready`: callers see it, nothing waits for ever)
    this.ready.catch(() => {})
  }
  _post(op, args, transfer = [], onProgress, onLog) {
    return new Promise((resolve, reject) => {
      if (this.failed) return reject(this.failed)
      const id = this.next++
      this.pending.set(id, { resolve, reject, onProgress, onLog })
      this.w.postMessage({ id, op, args }, transfer)
    })
  }
  /** The port another thread (the orchestrator) calls this leaf through. */
  async port() {
    const ch = new MessageChannel()
    await this.ready
    await this._post('port', { port: ch.port1 }, [ch.port1])
    return ch.port2
  }
  async call(name, t, params, opts = {}) {
    await this.ready
    const tt = table(t)
    const { onProgress, onLog, transfer = false, ...rest } = opts
    return decodeResult(await this._post('call', { name, table: tt, params, opts: rest, progress: !!onProgress, log: !!onLog }, transfer ? transferables(tt) : [], onProgress, onLog))
  }
  async callAsync(name, t, params, opts = {}) {
    await this.ready
    const tt = table(t)
    const { onProgress, onLog, transfer = false, ...rest } = opts
    return decodeResult(await this._post('callAsync', { name, table: tt, params, opts: rest, progress: !!onProgress, log: !!onLog }, transfer ? transferables(tt) : [], onProgress, onLog))
  }
  async read(handle, which) { await this.ready; return decodeResult(await this._post('read', { handle, which })) }
  async plan(name, t, params, opts = {}, resources = [], relaxVariant = false) {
    await this.ready
    const { onProgress, onLog, keep, callId, ...rest } = opts
    return this._post('plan', { name, table: describe(t), params, opts: rest, resources, relaxVariant })
  }
  /** A dataset for range reads in this Worker (crux-wasm's fnDataset: sync reads are allowed
   *  here); replayed on a respawn. Resolves its id. */
  async dataset(name, manifest, base) {
    this.resources.push({ op: 'dataset', name, manifest, base })
    await this.ready
    return this._post('dataset', { name, manifest, base })
  }
  /** `keep`: true, or 'pin' (a run's input: never evicted until freed). */
  async put(t, transfer = false, keep = true) { await this.ready; const tt = table(t); return this._post('put', { table: tt, keep }, transfer ? transferables(tt) : []) }
  async pin(handle, on = true) { await this.ready; return this._post('pin', { handle, on }) }
  /** `{instance, cap}` of the leaf's handle table, its cap set from `cfg.cap` where given (not
   *  replayed on a respawn: a new instance starts at the default). */
  async config(cfg = {}) { await this.ready; return this._post('config', cfg) }
  async gather(spec) { await this.ready; return decodeTable(await this._post('gather', { spec })) }
  async free(handle) { await this.ready; return this._post('free', { handle }) }
  async cancel(callId) { return this._post('cancel', { callId }) }
  /** Resource bytes by identity (bound per call by `resources: {name: {id}}`); replayed on a respawn. */
  async store(id, bytes, lease) {
    // (replayed on a respawn, without its lease: one entry per identity)
    this.resources = this.resources.filter((r) => !(r.op === 'store' && r.id === id))
    this.resources.push({ op: 'store', id, bytes })
    await this.ready
    return this._post('store', { id, bytes, lease })
  }
  /** Ends lease `lease` (its stored resources evictable again). */
  async release(lease) { await this.ready; return this._post('release', { lease }) }
  /** Drops stored resource `id` here (and from the respawn replay). */
  async unstore(id) {
    this.resources = this.resources.filter((r) => !(r.op === 'store' && r.id === id))
    await this.ready
    return this._post('unstore', { id })
  }
  /** Which of `ids` are not stored in the leaf (the others marked used). */
  async stored(ids, lease) { await this.ready; return this._post('stored', { ids, lease }) }
  async stats() { await this.ready; return this._post('stats', {}) }
  async resource(name, id, bytes) {
    this.resources.push({ name, id, bytes })
    await this.ready
    return this._post('resource', { name, id, bytes })
  }
  /** Ends a running synchronous call: terminate the Worker, fail its calls (`cancelled`), start
   *  a new one from the compiled module (resources replayed). Resolves when the new one is ready. */
  async respawn(why = 'cancelled') {
    this.w.terminate()
    for (const p of this.pending.values()) p.reject(Object.assign(new Error(why), { code: why === 'cancelled' ? 'cancelled' : 'worker_failed' }))
    this._start()
    return this.ready
  }
  terminate() { this.w.terminate() }
}
