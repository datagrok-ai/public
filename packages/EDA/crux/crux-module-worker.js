// One crux module in its own Web Worker ([decision:wasm-modules]; crux-modules.mjs starts it):
// `init` imports the module's wasm-bindgen glue (target web) and instantiates it, then
// callFn / callFnRow / fnDescribe / fnList / fnResource run in this thread. The callFn result's buffer is
// transferred back (awaited first: a GPU call's callFn returns a Promise).

let mod = null

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
      // callFn returns the IPC bytes, or a Promise of them for a call that awaits the device
      // (crux_fn::needs_async: device = gpu)
      result = await (args.opts ? mod.callFn(args.name, args.ipc, args.params, args.opts) : mod.callFn(args.name, args.ipc, args.params))
      transfer.push(result.buffer)
    } else if (op === 'callFnRow') {
      result = mod.callFnRow(args.name, args.inputs, args.params)
    } else if (op === 'fnDescribe') {
      result = mod.fnDescribe(args.name)
    } else if (op === 'fnList') {
      result = mod.fnList()
    } else if (op === 'fnResource') {
      mod.fnResource(args.name, args.id, args.bytes)
    } else if (op === 'cruxModules') {
      result = mod.cruxModules()
    } else {
      throw new Error(`unknown op ${op}`)
    }
    self.postMessage({ id, ok: true, result }, transfer)
  } catch (err) {
    self.postMessage({ id, ok: false, error: String((err && err.message) || err) })
  }
}
