// A crux-fn call's outputs envelope ([decision:crux-fn-outputs]): what crux-wasm's callFn /
// callFnAsync / runShape (and crux-server's Arrow mode, `crux fn run --out O.cxf`) return. One CXF1
// frame (crux-fn-core src/frame.rs, src/outputs.rs): `CXF1`, the JSON header's length (u32 LE), the
// header, then each attachment as its length (u64 LE) and bytes. The header:
//
//   {"v": 1, "fn": "<name>", "main": "<main output name>",
//    "outputs": {"<name>": {"kind": "table", "att": n}            // an Arrow IPC stream attachment
//                       | {"kind": "value", "value": …}            // inline JSON
//                       | {"kind": "binary", "media": "…", "att": n}}}
//
//   import { decodeOutputs, mainIpc } from './crux-outputs.js'
//   tableFromIPC(mainIpc(callFn('mol.qed', ipc, '{}')))      // the main table (single-table functions)
//   const o = decodeOutputs(callFn('t.fit', ipc, '{}'))
//   o.outputs.coef.ipc; o.outputs.objf.value; o.outputs.model.bytes; o.outputs.model.media
//
// No dependencies (no apache-arrow): a table stays its IPC bytes, a view into the envelope (no copy).

const MAGIC = [0x43, 0x58, 0x46, 0x31] // "CXF1"
const VERSION = 1

/** Whether `b` starts like an envelope (a CXF1 frame). */
export function isOutputs(b) {
  return b instanceof Uint8Array && b.length >= 8 && MAGIC.every((m, i) => b[i] === m)
}

/** The frame's header (parsed JSON) and its attachments (views into `b`). */
function frame(b) {
  if (!isOutputs(b)) throw new Error('crux outputs: not an outputs envelope (no CXF1 header)')
  const dv = new DataView(b.buffer, b.byteOffset, b.byteLength)
  const n = dv.getUint32(4, true)
  if (8 + n > b.length) throw new Error('crux outputs: truncated envelope')
  const header = JSON.parse(new TextDecoder().decode(b.subarray(8, 8 + n)))
  const atts = []
  let at = 8 + n
  while (at < b.length) {
    if (at + 8 > b.length) throw new Error('crux outputs: truncated envelope')
    const len = dv.getUint32(at, true) + dv.getUint32(at + 4, true) * 2 ** 32
    at += 8
    if (at + len > b.length) throw new Error('crux outputs: truncated envelope')
    atts.push(b.subarray(at, at + len))
    at += len
  }
  return { header, atts }
}

/**
 * The outputs of an envelope: `{fn, main, outputs, names, mainIpc()}`. `outputs[name]` is
 * `{kind: 'table', ipc}` (Arrow IPC stream bytes), `{kind: 'value', value}` or `{kind: 'binary',
 * media, bytes}`; `main` names the main table (undefined: the function has none); `names`: every
 * output's name, the main table first.
 */
export function decodeOutputs(b) {
  const { header: h, atts } = frame(b)
  if (h.v !== VERSION) throw new Error(`crux outputs: envelope version ${h.v} (this side reads ${VERSION})`)
  const att = (o, name) => {
    if (o.blob !== undefined) throw new Error(`crux outputs: output ${name} is a blob reference (resolved by Crux Compute)`)
    const a = atts[o.att]
    if (a === undefined) throw new Error(`crux outputs: output ${name}: attachment ${o.att} missing`)
    return a
  }
  const outputs = {}
  for (const [name, o] of Object.entries(h.outputs ?? {})) {
    if (o.kind === 'table') outputs[name] = { kind: 'table', ipc: att(o, name) }
    else if (o.kind === 'value') outputs[name] = { kind: 'value', value: o.value }
    else if (o.kind === 'binary') outputs[name] = { kind: 'binary', media: o.media, bytes: att(o, name) }
    else throw new Error(`crux outputs: output ${name}: unknown kind ${o.kind}`)
  }
  const main = h.main ?? undefined
  if (main !== undefined && outputs[main]?.kind !== 'table') throw new Error(`crux outputs: main output ${main} missing`)
  const names = Object.keys(outputs)
  if (main !== undefined) names.sort((a, c) => (a === main ? -1 : c === main ? 1 : 0))
  return {
    fn: h.fn,
    main,
    outputs,
    names,
    /** The main table's IPC stream bytes (an error when the function has none). */
    mainIpc() {
      if (main === undefined) throw new Error(`${h.fn}: has no main output (only named outputs: ${names.join(', ')})`)
      return outputs[main].ipc
    },
  }
}

/** The main table's Arrow IPC stream bytes of an envelope (single-table functions). */
export function mainIpc(b) {
  return decodeOutputs(b).mainIpc()
}

/**
 * An envelope of `outputs` (`{name: {kind: 'table', ipc} | {kind: 'value', value} | {kind: 'binary',
 * media, bytes}}`) of function `fn` whose main table is `main` (omit: none). The JS side's writer:
 * tests, and hosts that build results themselves.
 */
export function encodeOutputs(fn, outputs, main) {
  const atts = []
  const entries = {}
  const ordered = Object.keys(outputs)
  if (main !== undefined) {
    if (outputs[main]?.kind !== 'table') throw new Error(`crux outputs: main output ${main} is not a table`)
    ordered.sort((a, c) => (a === main ? -1 : c === main ? 1 : 0))
  }
  for (const name of ordered) {
    const o = outputs[name]
    if (o.kind === 'table') entries[name] = { kind: 'table', att: atts.push(o.ipc) - 1 }
    else if (o.kind === 'value') entries[name] = { kind: 'value', value: o.value }
    else if (o.kind === 'binary') entries[name] = { kind: 'binary', media: o.media, att: atts.push(o.bytes) - 1 }
    else throw new Error(`crux outputs: output ${name}: unknown kind ${o.kind}`)
  }
  const head = { v: VERSION, fn, ...(main !== undefined ? { main } : {}), outputs: entries }
  const json = new TextEncoder().encode(JSON.stringify(head))
  const out = new Uint8Array(8 + json.length + atts.reduce((s, a) => s + 8 + a.length, 0))
  const dv = new DataView(out.buffer)
  out.set(MAGIC, 0)
  dv.setUint32(4, json.length, true)
  out.set(json, 8)
  let at = 8 + json.length
  for (const a of atts) {
    dv.setUint32(at, a.length % 2 ** 32, true)
    dv.setUint32(at + 4, Math.floor(a.length / 2 ** 32), true)
    out.set(a, at + 8)
    at += 8 + a.length
  }
  return out
}
