// One crux module in its own Web Worker ([decision:wasm-modules]; crux-modules.js starts it):
// `init` imports the module's wasm-bindgen glue (target web) and instantiates it, then
// leafCall / leafCallAsync (the leaf ABI, [decision:wasm-leaf-abi]: columns in, the result's
// columns out, their fresh buffers transferred back) / callFn / fnDescribe / fnList /
// fnResource run in this thread. The callFn result's buffer is transferred back (awaited first:
// a GPU call's callFn returns a Promise).
// (a static import: this Worker runs from the shipped layout only, crux-leaf.js beside it; its
// message handler must be set when the module is evaluated)
import { transferables } from './crux-leaf.js'

let mod = null

function resultTransfer(r) {
  const t = r?.main ? transferables(r.main) : []
  for (const v of Object.values(r?.named || {})) {
    if (v.table) t.push(...transferables(v.table))
    if (v.bytes) t.push(v.bytes.buffer)
  }
  return t
}

self.onmessage = async (e) => {
  const { id, op, args } = e.data
  try {
    let result = null
    const transfer = []
    if (op === 'init') {
      mod = await import(args.js)
      if (typeof mod.default === 'function') await mod.default()
      result = typeof mod.cruxModules === 'function' ? mod.cruxModules() : null
    } else if (!mod) {
      throw new Error('crux module not initialised')
    } else if (op === 'callFn') {
      if (typeof mod.callFn !== 'function') throw Object.assign(new Error(`callFn is in validation builds only ([decision:wasm-leaf-abi] P1c: modules.py build --validation; node packages have it): call the module through the leaf ABI (call)`), { code: 'no_callfn' })
      // callFn returns the outputs envelope's bytes, or a Promise of them for a call that awaits the device
      // (crux_fn::needs_async: device = gpu)
      result = await (args.opts ? mod.callFn(args.name, args.ipc, args.params, args.opts) : mod.callFn(args.name, args.ipc, args.params))
      transfer.push(result.buffer)
    } else if (op === 'leafCall' || op === 'leafCallAsync') {
      const json = args.opts && Object.keys(args.opts).length ? JSON.stringify(args.opts) : undefined
      let last = 0
      const onProgress = args.progress ? (...p) => { const now = performance.now(); if (now - last >= 50) { last = now; self.postMessage({ id, progress: p }) } } : undefined
      result = op === 'leafCall'
        ? mod.leafCall(args.name, args.table, args.params ?? null, json, false, undefined, onProgress, undefined)
        : await mod.leafCallAsync(args.name, args.table, args.params ?? null, json, false, undefined, onProgress, undefined)
      transfer.push(...resultTransfer(result))
    } else if (op === 'fnDescribe') {
      result = mod.fnDescribe(args.name)
    } else if (op === 'fnList') {
      result = mod.fnList()
    } else if (op === 'fnResource') {
      mod.fnResource(args.name, args.id, args.bytes)
    } else if (op === 'fnDataset') {
      // reads block (sync XHR / FileReaderSync): allowed here, in the worker
      result = mod.fnDataset(args.name, args.manifest, args.base)
    } else if (op === 'cruxModules') {
      result = mod.cruxModules()
    } else {
      throw new Error(`unknown op ${op}`)
    }
    self.postMessage({ id, ok: true, result }, transfer)
  } catch (err) {
    self.postMessage({ id, ok: false, error: String((err && err.message) || err), code: err?.code })
  }
}
