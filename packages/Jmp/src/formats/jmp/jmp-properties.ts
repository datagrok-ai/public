// The column properties a .jmp file keeps (notes, units, formulas, value labels and the rest). Not in JMPReader.jl.
// The container is jmp-convert's (github.com/davidmathers/jmp-convert, MIT); what the entries hold was read
// from public files (JMP_PLUGIN.md §4.6). Provisional until confirmed on more files.
//
// After a column's six descriptor bytes, 5 bytes for a numeric column (dt1 0x01, 0x0A) or 4 for any other,
// then u16 count, u64 total length, and count × (u16 id, u64 length, payload).
//   0x07 = JSL properties, each a u32 length and JSL text such as `Notes("...")`, ended by 0xFFFFFFFF
//   0x0E = value labels, u32 n, u32 length, n keys (f64 for a numeric column, NUL-ended text otherwise),
//         u32 n, u32 length, n NUL-ended labels

import {ByteReader} from './jmp-reader';
import {JslNode, JslSyntaxError, parseJsl} from './jsl';

export interface JmpColumnProperty {
  /** The property's name as JMP writes it, e.g. "Notes" or "Value Ordering". */
  name: string;
  /** The JSL text, e.g. `Notes("Was the offer accepted?")`. */
  text: string;
  /** The parsed property; undefined when it couldn't be read. */
  node?: Extract<JslNode, {type: 'call'}>;
}

export interface JmpValueLabels {
  keys: (number | string)[];
  labels: string[];
}

export interface JmpColumnProperties {
  properties: JmpColumnProperty[];
  valueLabels?: JmpValueLabels;
}

const PROP_JSL = 0x07;
const PROP_VALUE_LABELS = 0x0e;
const END_OF_LIST = 0xffffffff;
const MAX_PROPERTIES = 1000;
const NUMERIC_STORAGE = new Set([0x01, 0x0a, 0x03, 0xff, 0xfe, 0xfc]);

const utf8 = new TextDecoder('utf-8');

/** Reads the properties of the column whose descriptor bytes end at `afterDescriptor`. A column whose
 * header doesn't have the shape above has none. */
export function readColumnProperties(bytes: Uint8Array, afterDescriptor: number, dt1: number): JmpColumnProperties {
  const gaps = dt1 === 0x01 || dt1 === 0x0a ? [5, 4] : [4, 5];
  for (const gap of gaps) {
    const entries = propertyEntries(bytes, afterDescriptor + gap);
    if (entries)
      return decode(entries, dt1);
  }
  return {properties: []};
}

function propertyEntries(bytes: Uint8Array, at: number): [number, Uint8Array][] | null {
  try {
    const r = new ByteReader(bytes);
    r.seek(at);
    const n = r.u16();
    const total = r.u64();
    const start = r.pos;
    if (n > MAX_PROPERTIES || total > r.length - start)
      return null;
    const out: [number, Uint8Array][] = [];
    for (let k = 0; k < n; k++) {
      const id = r.u16();
      const len = r.u64();
      out.push([id, r.bytesAt(len)]);
    }
    return r.pos - start === total ? out : null;
  } catch {
    return null;
  }
}

function decode(entries: [number, Uint8Array][], dt1: number): JmpColumnProperties {
  const out: JmpColumnProperties = {properties: []};
  const numeric = NUMERIC_STORAGE.has(dt1);
  for (const [id, payload] of entries) {
    if (id === PROP_JSL)
      out.properties.push(...jslProperties(payload));
    else if (id === PROP_VALUE_LABELS)
      out.valueLabels = valueLabels(payload, numeric);
  }
  return out;
}

function jslProperties(payload: Uint8Array): JmpColumnProperty[] {
  const r = new ByteReader(payload);
  const out: JmpColumnProperty[] = [];
  while (r.pos + 4 <= r.length) {
    const len = r.u32();
    if (len === END_OF_LIST || len > r.length - r.pos)
      break;
    const text = utf8.decode(r.bytesAt(len));
    const name = /^\s*(?:"([^"]*)"n?|([A-Za-z_][\w ]*?))\s*\(/.exec(text);
    if (!name)
      continue;
    let node: JmpColumnProperty['node'];
    try {
      const parsed = parseJsl(text);
      if (parsed.type === 'call')
        node = parsed;
    } catch (e) {
      if (!(e instanceof JslSyntaxError))
        throw e;
    }
    out.push({name: name[1] ?? name[2], text, node});
  }
  return out;
}

function valueLabels(payload: Uint8Array, numeric: boolean): JmpValueLabels | undefined {
  try {
    const r = new ByteReader(payload);
    const n = r.u32();
    const keyBytes = r.bytesAt(r.u32());
    const n2 = r.u32();
    const labels = nulSeparated(r.bytesAt(r.u32()), n);
    if (n2 !== n || labels.length !== n || r.pos !== r.length)
      return undefined;
    if (numeric) {
      if (keyBytes.length !== 8 * n)
        return undefined;
      const v = new DataView(keyBytes.buffer, keyBytes.byteOffset, keyBytes.byteLength);
      return {keys: Array.from({length: n}, (_, i) => v.getFloat64(8 * i, true)), labels};
    }
    const keys = nulSeparated(keyBytes, n);
    return keys.length === n ? {keys, labels} : undefined;
  } catch {
    return undefined;
  }
}

function nulSeparated(bytes: Uint8Array, n: number): string[] {
  const out: string[] = [];
  let start = 0;
  for (let i = 0; i < bytes.length && out.length < n; i++) {
    if (bytes[i] === 0) {
      out.push(utf8.decode(bytes.subarray(start, i)));
      start = i + 1;
    }
  }
  return out;
}
