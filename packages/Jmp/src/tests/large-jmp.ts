// Makes a large .jmp file for performance tests (JMP_PLUGIN.md §8), from JMPReader.jl's compressed.jmp with its
// seven columns (number, date-time, two text columns, duration, date, long text) refilled with `nrows`
// rows of varied values, so gzip can't shrink them away. Every block is rebuilt in the layout the parser
// reads and the metadata's row count and column offsets are patched; bytes the parser doesn't read are
// copied from the template. The file is valid for this reader; JMP itself was never asked to open one.
// It has no platform imports, so the Node tool and the browser test share it.

import {ByteReader} from '../formats/jmp/jmp-reader';
import {readHeader, readMetadata} from '../formats/jmp/jmp-parser';
import {GZIP_SECTION_START} from '../formats/jmp/jmp-types';

const METADATA_ANCHOR = [0x07, 0x00, 0x08, 0x00, 0x00, 0x00];
const JMP_2000 = 3029529600; // 2000-01-01 in seconds since 1904
const WORDS = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'eta', 'theta', 'iota', 'kappa', 'lambda', 'mu'];

async function gzip(data: Uint8Array): Promise<Uint8Array> {
  const stream = new Blob([data as BlobPart]).stream().pipeThrough(new CompressionStream('gzip'));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}

/** A seeded generator (mulberry32), so the same size gives the same file. 32-bit integer arithmetic
 * throughout: a multiply in doubles would lose bits and repeat after a short cycle. */
function random(seed: number): () => number {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function concat(parts: Uint8Array[]): Uint8Array {
  const out = new Uint8Array(parts.reduce((n, p) => n + p.length, 0));
  let at = 0;
  for (const p of parts) {
    out.set(p, at);
    at += p.length;
  }
  return out;
}

function indexOf(bytes: Uint8Array, seq: number[] | Uint8Array, from = 0): number {
  outer: for (let i = from; i <= bytes.length - seq.length; i++) {
    for (let j = 0; j < seq.length; j++) {
      if (bytes[i + j] !== seq[j])
        continue outer;
    }
    return i;
  }
  return -1;
}

function u64(n: number): Uint8Array {
  const b = new Uint8Array(8);
  new DataView(b.buffer).setBigUint64(0, BigInt(n), true);
  return b;
}

/** Eight-byte numbers, 5% missing (NaN). */
function numbers(nrows: number, rnd: () => number, value: () => number): Uint8Array {
  const out = new Uint8Array(nrows * 8);
  const v = new DataView(out.buffer);
  for (let k = 0; k < nrows; k++)
    v.setFloat64(k * 8, rnd() < 0.05 ? NaN : value(), true);
  return out;
}

/** The compressed variable-width text layout: i64 (length − 8), width byte size, u32 max width, widths, text. */
function variableStrings(values: string[]): Uint8Array {
  const enc = new TextEncoder();
  const bytes = values.map((s) => enc.encode(s));
  const max = bytes.reduce((m, b) => Math.max(m, b.length), 0);
  const wb = max < 128 ? 1 : 2;
  const total = 8 + 1 + 4 + values.length * wb + bytes.reduce((n, b) => n + b.length, 0);
  const out = new Uint8Array(total);
  const v = new DataView(out.buffer);
  v.setBigInt64(0, BigInt(total - 8), true);
  out[8] = wb;
  v.setUint32(9, max, true);
  let p = 13;
  for (const b of bytes) {
    if (wb === 1)
      out[p] = b.length;
    else
      v.setInt16(p, b.length, true);
    p += wb;
  }
  for (const b of bytes) {
    out.set(b, p);
    p += b.length;
  }
  return out;
}

/** Fixed-width text: `width` bytes per row, zero-padded. */
function fixedStrings(nrows: number, width: number, rnd: () => number): Uint8Array {
  const out = new Uint8Array(nrows * width);
  for (let k = 0; k < nrows; k++)
    out[k * width] = 0x41 + Math.floor(rnd() * 26);
  return out;
}

/** compressed.jmp's columns, in order, and how to fill each. */
function columnData(name: string, nrows: number, rnd: () => number): Uint8Array {
  switch (name) {
  case 'numeric': return numbers(nrows, rnd, () => (rnd() - 0.5) * 1e4);
  case 'y-m-d h:m:s': return numbers(nrows, rnd, () => JMP_2000 + Math.floor(rnd() * 8e8));
  case 'yyyy-mm-dd': return numbers(nrows, rnd, () => JMP_2000 + Math.floor(rnd() * 9000) * 86400);
  case 'min:s': return numbers(nrows, rnd, () => Math.floor(rnd() * 36000));
  case 'character1': return fixedStrings(nrows, 2, rnd);
  case 'character11': return variableStrings(Array.from({length: nrows}, () =>
    `${WORDS[Math.floor(rnd() * WORDS.length)]}-${Math.floor(rnd() * 500)}`));
  case 'character130': return variableStrings(Array.from({length: nrows}, () => {
    const n = 3 + Math.floor(rnd() * 15);
    return Array.from({length: n}, () => WORDS[Math.floor(rnd() * WORDS.length)] + Math.floor(rnd() * 1e4)).join(' ');
  }));
  default: throw new Error(`large-jmp: the template has an unexpected column "${name}"`);
  }
}

/** A .jmp file of `nrows` rows built from `template` (the bytes of compressed.jmp). */
export async function makeLargeJmp(template: Uint8Array, nrows: number, seed = 42): Promise<Uint8Array> {
  const r = new ByteReader(template);
  const info = readMetadata(r, readHeader(r));
  const offsets = info.columns.map((c) => c.offset);
  const rnd = random(seed);

  const blocks: Uint8Array[] = [];
  for (let i = 0; i < info.ncols; i++) {
    const start = offsets[i];
    const end = i + 1 < info.ncols ? offsets[i + 1] : template.length;
    const marker = indexOf(template.subarray(start, end), GZIP_SECTION_START);
    if (marker < 0)
      throw new Error(`large-jmp: column "${info.columns[i].name}" isn't compressed in the template`);
    const header = template.subarray(start, start + marker + GZIP_SECTION_START.length);
    const raw = columnData(info.columns[i].name, nrows, rnd);
    const gz = await gzip(raw);
    blocks.push(concat([header, u64(gz.length), u64(raw.length), gz]));
  }

  // metadata: the row count sits 38 bytes before the end of the anchor (readMetadata); the column offsets
  // follow the column count
  const meta = template.slice(0, offsets[0]);
  const anchor = indexOf(meta, METADATA_ANCHOR);
  new DataView(meta.buffer).setBigInt64(anchor + METADATA_ANCHOR.length - 38, BigInt(nrows), true);
  const table = indexOf(meta, concat([u64(offsets[0]), u64(offsets[1])]));
  if (anchor < 0 || table < 0)
    throw new Error('large-jmp: the template\'s metadata isn\'t laid out as expected');
  let at = offsets[0];
  for (let i = 0; i < info.ncols; i++) {
    new DataView(meta.buffer).setBigUint64(table + i * 8, BigInt(at), true);
    at += blocks[i].length;
  }
  return concat([meta, ...blocks]);
}

/** Rows that make a file of about `bytes` (measured: about 53 bytes a row). */
export function rowsForSize(bytes: number): number {
  return Math.round(bytes / 53);
}

