// One leaf in its own Web Worker ([decision:wasm-leaf-abi]; crux-leaf.js `LeafWorker` starts it):
// `init` instantiates the leaf's web glue (from a compiled WebAssembly.Module when given: a respawn
// compiles nothing), then call / callAsync / put / gather / read / describe / free / cancel /
// resource run here. Inputs arrive transferred; results (fresh buffers: the leaf's copy out) are
// transferred back, raw (the receiver decodes). Progress is posted at most every 50 ms.
//
// Ports (P1a): `port` adds a MessagePort that speaks the same messages (the orchestrator's);
// `link {peer, port}` a port to another leaf Worker (made by the orchestrator: leaf to leaf, the
// data never passes through it); `send {to, spec}` gathers columns here and puts them on leaf
// `to` over its link, answering with the handle there. A port carries requests ({id, op, args})
// and answers ({id, ok, ...}) both ways.

import { loadLeaf, transferables } from './crux-leaf.js'

let leaf = null
const peers = new Map() // leaf id -> requester over its link port

function requester(port) {
  let next = 1
  const pending = new Map()
  const call = (op, args, transfer = []) => new Promise((resolve, reject) => {
    const id = next++
    pending.set(id, { resolve, reject })
    port.postMessage({ id, op, args }, transfer)
  })
  call.answer = (m) => {
    const p = pending.get(m.id)
    if (!p || m.progress) return
    pending.delete(m.id)
    m.ok ? p.resolve(m.result) : p.reject(Object.assign(new Error(m.error), { code: m.code }))
  }
  return call
}

function serve(port) {
  port.onmessage = (e) => (e.data.op ? handle((m, t) => port.postMessage(m, t), e.data) : port.requester?.answer(e.data))
}

function resultTransfer(result) {
  if (!result || typeof result !== 'object') return []
  const transfer = result.main ? transferables(result.main) : result.columns ? transferables(result) : []
  for (const v of Object.values(result.named || {})) {
    if (v.table) transfer.push(...transferables(v.table))
    if (v.bytes) transfer.push(v.bytes.buffer)
  }
  return transfer
}

async function handle(post, { id, op, args }) {
  try {
    let result = null
    if (op === 'init') {
      leaf = await loadLeaf(args.js, args.module)
      // the host's number for this instance: its handles carry it
      if (args.instance) leaf.config({ instance: args.instance })
      result = { ...leaf.info(), instance: leaf.config().instance }
    } else if (op === 'port') {
      serve(args.port)
      result = true
    } else if (op === 'link') {
      args.port.requester = requester(args.port)
      peers.set(args.peer, args.port.requester)
      serve(args.port)
      result = true
    } else if (!leaf) {
      throw new Error('leaf not initialised')
    } else if (op === 'call' || op === 'callAsync') {
      const opts = { ...args.opts }
      if (args.progress) {
        let last = 0
        opts.onProgress = (...p) => {
          const now = performance.now()
          if (now - last >= 50) { last = now; post({ id, progress: p }) }
        }
      }
      // the call's log lines, as they come (a caller that asked: the orchestrator's task log)
      if (args.log) opts.onLog = (level, msg) => post({ id, log: [level, msg] })
      const { keep, callId, onProgress, onLog, ...rest } = opts
      const json = Object.keys(rest).length ? JSON.stringify(rest) : undefined
      const kept = args.keep ?? keep ?? false
      result = op === 'call'
        ? leaf.x.leafCall(args.name, args.table, args.params ?? null, json, kept, args.callId ?? callId, onProgress, onLog)
        : await leaf.x.leafCallAsync(args.name, args.table, args.params ?? null, json, kept, args.callId ?? callId, onProgress, onLog)
    } else if (op === 'put') {
      result = leaf.x.leafPut(args.table, args.keep ?? true)
    } else if (op === 'pin') {
      result = leaf.pin(args.handle, args.on ?? true)
    } else if (op === 'gather') {
      result = leaf.x.leafGather(args.spec)
    } else if (op === 'send') {
      // columns of this leaf onto another, over their link: gathered here (one copy out), put there
      const to = peers.get(args.to)
      if (!to) throw new Error(`no link to leaf ${args.to}`)
      const t = leaf.x.leafGather(args.spec)
      // (pinned there: a run's live columns)
      result = await to('put', { table: t, keep: 'pin' }, transferables(t))
    } else if (op === 'read') {
      result = leaf.x.leafRead(args.handle, false, args.which ?? undefined)
    } else if (op === 'plan') {
      const { keep, callId, ...rest } = args.opts || {}
      const json = Object.keys(rest).length ? JSON.stringify(rest) : undefined
      result = leaf.x.leafPlan(args.name, args.table, args.params ?? null, json, args.resources || [], !!args.relaxVariant)
    } else if (op === 'source') {
      // one object read in place by range (crux-wasm's fnSource: a URL the server serves)
      if (typeof leaf.x.fnSource !== 'function') throw new Error('this leaf takes no sources (no fnSource)')
      result = leaf.x.fnSource(args.name, args.id, args.src) ?? true
    } else if (op === 'dataset') {
      // reads block (sync XHR / FileReaderSync): allowed here, in the Worker
      if (typeof leaf.x.fnDataset !== 'function') throw new Error('this leaf takes no datasets (no fnDataset)')
      result = leaf.x.fnDataset(args.name, args.manifest, args.base)
    } else if (op === 'describe') {
      result = leaf.x.leafDescribe(args.handle)
    } else if (op === 'free') {
      if (args.handles) for (const h of args.handles) leaf.free(h)
      else result = leaf.free(args.handle)
    } else if (op === 'cancel') {
      result = leaf.cancel(args.callId)
    } else if (op === 'config') {
      // the handle table's cap (`{cap}`; the instance is init's)
      const { instance, ...rest } = args || {}
      result = leaf.config(rest)
    } else if (op === 'store') {
      // resource bytes by identity: only a call that binds them by {id} sees them (counted against
      // the handle cap; evicted when no call holds them and room is needed: `stored` tells)
      // `lease`: pinned under it in the same step, until the call naming it ends (or `release`)
      result = leaf.x.leafStore(args.id, args.bytes, args.lease ?? undefined) ?? true
    } else if (op === 'unstore') {
      result = leaf.x.leafUnstore(args.id)
    } else if (op === 'stored') {
      // the identities of `ids` not stored here (the host stores them again); the others leased
      result = typeof leaf.x.leafStored === 'function' ? Array.from(leaf.x.leafStored(args.ids || [], args.lease ?? undefined)) : args.ids || []
    } else if (op === 'release') {
      result = typeof leaf.x.leafRelease === 'function' ? leaf.x.leafRelease(args.lease) : false
    } else if (op === 'resource') {
      result = leaf.resource(args.name, args.id, args.bytes)
    } else if (op === 'stats') {
      result = leaf.stats()
    } else {
      throw new Error(`unknown op ${op}`)
    }
    post({ id, ok: true, result }, resultTransfer(result))
  } catch (e) {
    // a trap (a Rust panic, a wasm fault): this instance may be left inconsistent; the call fails
    // and the host is told the Worker is done (it restarts it from the compiled module)
    const trap = typeof WebAssembly !== 'undefined' && e instanceof WebAssembly.RuntimeError
    post({ id, ok: false, error: String(e?.message ?? e), code: trap ? 'worker_failed' : e?.code })
    if (trap) self.postMessage({ dead: `the leaf trapped: ${e.message}` })
  }
}

self.onmessage = (e) => handle((m, t) => self.postMessage(m, t), e.data)
