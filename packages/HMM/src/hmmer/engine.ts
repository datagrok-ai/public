// Synced from Rusty-HMMER web/engine.ts by web/sync-datagrok.mjs — do not edit here.
/** Typed wrapper of the lean HMMER WebAssembly engine (crates/hmmer-web).
 *
 * The module has no imports, so it instantiates the same way in browsers,
 * Web Workers and Node. Sequences go in as one packed buffer; results come
 * back as packed little-endian records decoded here (layout: `pack_query` in
 * crates/hmmer-web/src/lib.rs). Keep `ABI_VERSION` in sync with the crate. */

export const ABI_VERSION = 2;
const RESULT_MAGIC = 0x31524d48;

/** Option bits of the packed options record. */
export const Options = {
  alignments: 1,
  allHits: 2,
  cutGa: 4,
  cutTc: 8,
  cutNc: 16,
  max: 32,
  noBias: 64,
  noNull2: 128,
  t: 256,
  domT: 512,
  incT: 1024,
  incdomT: 2048,
  seed: 4096,
} as const;

/** Per-query status codes. */
export const QueryStatus = {
  ok: 0,
  invalidResidue: 1,
  empty: 2,
  failed: 3,
} as const;

/** `p7_IS_REPORTED`, `p7_IS_INCLUDED` hit flags. */
export const HitFlags = {included: 1, reported: 2} as const;

export type Alphabet = 'amino' | 'dna' | 'rna';

export interface ModelInfo {
  index: number;
  length: number;
  name: string;
  accession: string;
  description: string;
  /** `M + 1` characters (position 0 unused), or empty. */
  consensus: string;
  reference: string;
  gathering: [number, number] | null;
  trusted: [number, number] | null;
  noise: [number, number] | null;
  maxLength: number;
}

export interface DatabaseInfo {
  alphabet: Alphabet;
  models: ModelInfo[];
}

export interface Alignment {
  model: string;
  match: string;
  target: string;
  /** Posterior probability line (`0-9`, `*`, `.` at deletions), if any. */
  posterior: string | null;
  /** The model's RF line (`.` at insertions), if the model has one. */
  reference: string | null;
}

export interface DomainHit {
  reported: boolean;
  included: boolean;
  envFrom: number;
  envTo: number;
  aliFrom: number;
  aliTo: number;
  hmmFrom: number;
  hmmTo: number;
  bitscore: number;
  bias: number;
  accuracy: number;
  lnP: number;
  /** `exp(lnP)` as the engine computes it (E-values are this times Z or domZ). */
  p: number;
  cEvalue: number;
  iEvalue: number;
  /** Values as HMMER's text output prints them (`%.1f`, `%.2g`), read back. */
  printed: {bitscore: number; bias: number; cEvalue: number; iEvalue: number};
  alignment: Alignment | null;
}

export interface Hit {
  /** Model index (scan) or sequence index (search). */
  target: number;
  flags: number;
  score: number;
  preScore: number;
  sumScore: number;
  lnP: number;
  /** `exp(lnP)` as the engine computes it (the E-value is this times Z). */
  p: number;
  evalue: number;
  expectedDomains: number;
  regions: number;
  clustered: number;
  overlaps: number;
  envelopes: number;
  bestDomain: number;
  domains: DomainHit[];
}

export interface QueryResult {
  status: number;
  z: number;
  domZ: number;
  hits: Hit[];
}

/** HMMER's search options (hmmsearch/hmmscan command-line equivalents). */
export interface SearchOptions {
  /** Return alignment display lines (needed for numbering and display). */
  alignments?: boolean;
  /** Return every registered hit and domain, not only reported ones. */
  allHits?: boolean;
  /** `--cut_ga`, `--cut_tc`, `--cut_nc`. */
  cutoff?: 'ga' | 'tc' | 'nc';
  /** `--max`: disable the acceleration filters. */
  max?: boolean;
  /** `--nobias`, `--nonull2`. */
  noBias?: boolean;
  noNull2?: boolean;
  /** `-E`, `--domE`, `--incE`, `--incdomE` (defaults 10, 10, 0.01, 0.01). */
  evalue?: number;
  domEvalue?: number;
  incEvalue?: number;
  incdomEvalue?: number;
  /** `-T`, `--domT`, `--incT`, `--incdomT`: bit score thresholds instead of E-values. */
  score?: number;
  domScore?: number;
  incScore?: number;
  incdomScore?: number;
  /** `-Z`, `--domZ`: search space sizes. */
  z?: number;
  domZ?: number;
  /** `--F1`, `--F2`, `--F3`. */
  f1?: number;
  f2?: number;
  f3?: number;
  /** `--seed` (default 42). */
  seed?: number;
}

/** A sequence to compare: residues only (named by its index) or name + residues. */
export type SequenceInput = string | {name: string; residues: string};

interface Exports {
  memory: WebAssembly.Memory;
  hw_version(): number;
  hw_alloc(length: number): number;
  hw_free(pointer: number): void;
  hw_len(pointer: number): number;
  hw_error(): number;
  hw_db_load(pointer: number): number;
  hw_db_free(handle: number): void;
  hw_db_info(handle: number): number;
  hw_scan(handle: number, sequences: number, options: number): number;
  hw_search(handle: number, model: number, sequences: number, options: number): number;
}

const latin1 = new TextDecoder('latin1');
const utf8 = new TextDecoder();
const encoder = new TextEncoder();

class Reader {
  private at = 0;
  private readonly bytes: Uint8Array;
  private readonly view: DataView;
  constructor(bytes: Uint8Array) {
    this.bytes = bytes;
    this.view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  }
  u32(): number {
    const value = this.view.getUint32(this.at, true);
    this.at += 4;
    return value;
  }
  f32(): number {
    const value = this.view.getFloat32(this.at, true);
    this.at += 4;
    return value;
  }
  f64(): number {
    const value = this.view.getFloat64(this.at, true);
    this.at += 8;
    return value;
  }
  ascii(length: number): string {
    const text = latin1.decode(this.bytes.subarray(this.at, this.at + length));
    this.at += length;
    return text;
  }
  text(): string {
    const length = this.u32();
    const text = utf8.decode(this.bytes.subarray(this.at, this.at + length));
    this.at += length;
    return text;
  }
}

function packSequences(sequences: SequenceInput[]): Uint8Array {
  const names: Uint8Array[] = new Array(sequences.length);
  const residues: Uint8Array[] = new Array(sequences.length);
  let size = 4;
  for (let i = 0; i < sequences.length; i++) {
    const s = sequences[i];
    names[i] = encoder.encode(typeof s === 'string' ? String(i) : s.name);
    residues[i] = encoder.encode(typeof s === 'string' ? s : s.residues);
    size += 8 + names[i].length + residues[i].length;
  }
  const bytes = new Uint8Array(size);
  const view = new DataView(bytes.buffer);
  view.setUint32(0, sequences.length, true);
  let at = 4;
  for (let i = 0; i < sequences.length; i++) {
    for (const part of [names[i], residues[i]]) {
      view.setUint32(at, part.length, true);
      bytes.set(part, at + 4);
      at += 4 + part.length;
    }
  }
  return bytes;
}

/** The 128-byte options record (`Options` in crates/hmmer-web). */
function packOptions(o: SearchOptions): Uint8Array {
  let bits = 0;
  if (o.alignments) bits |= Options.alignments;
  if (o.allHits) bits |= Options.allHits;
  if (o.cutoff === 'ga') bits |= Options.cutGa;
  if (o.cutoff === 'tc') bits |= Options.cutTc;
  if (o.cutoff === 'nc') bits |= Options.cutNc;
  if (o.max) bits |= Options.max;
  if (o.noBias) bits |= Options.noBias;
  if (o.noNull2) bits |= Options.noNull2;
  if (o.score !== undefined) bits |= Options.t;
  if (o.domScore !== undefined) bits |= Options.domT;
  if (o.incScore !== undefined) bits |= Options.incT;
  if (o.incdomScore !== undefined) bits |= Options.incdomT;
  if (o.seed !== undefined) bits |= Options.seed;
  const bytes = new Uint8Array(128);
  const view = new DataView(bytes.buffer);
  view.setUint32(0, bits, true);
  view.setUint32(4, o.seed ?? 0, true);
  const reals = [o.evalue, o.domEvalue, o.incEvalue, o.incdomEvalue, o.score, o.domScore, o.incScore,
    o.incdomScore, o.z, o.domZ, o.f1, o.f2, o.f3];
  reals.forEach((value, i) => view.setFloat64(8 + 8 * i, value ?? 0, true));
  return bytes;
}

function readDomain(r: Reader): DomainHit {
  const flags = r.u32();
  const [envFrom, envTo, aliFrom, aliTo, hmmFrom, hmmTo] = [r.u32(), r.u32(), r.u32(), r.u32(), r.u32(), r.u32()];
  const bitscore = r.f32();
  const bias = r.f64();
  const accuracy = r.f32();
  const lnP = r.f64();
  const p = r.f64();
  const cEvalue = r.f64();
  const iEvalue = r.f64();
  const printed = {bitscore: r.f64(), bias: r.f64(), cEvalue: r.f64(), iEvalue: r.f64()};
  const columns = r.u32();
  let alignment: Alignment | null = null;
  if (columns > 0) {
    const model = r.ascii(columns);
    const match = r.ascii(columns);
    const target = r.ascii(columns);
    const posterior = flags & 4 ? r.ascii(columns) : null;
    const reference = flags & 8 ? r.ascii(columns) : null;
    alignment = {model, match, target, posterior, reference};
  }
  return {
    reported: (flags & 1) !== 0, included: (flags & 2) !== 0,
    envFrom, envTo, aliFrom, aliTo, hmmFrom, hmmTo,
    bitscore, bias, accuracy, lnP, p, cEvalue, iEvalue, printed, alignment,
  };
}

function readQuery(r: Reader): QueryResult {
  const status = r.u32();
  const z = r.f64();
  const domZ = r.f64();
  const count = r.u32();
  const hits: Hit[] = new Array(count);
  for (let h = 0; h < count; h++) {
    const target = r.u32();
    const flags = r.u32();
    const score = r.f32();
    const preScore = r.f32();
    const sumScore = r.f32();
    const lnP = r.f64();
    const p = r.f64();
    const evalue = r.f64();
    const expectedDomains = r.f32();
    const [regions, clustered, overlaps, envelopes, bestDomain] = [r.u32(), r.u32(), r.u32(), r.u32(), r.u32()];
    const domains: DomainHit[] = new Array(r.u32());
    for (let d = 0; d < domains.length; d++) domains[d] = readDomain(r);
    hits[h] = {target, flags, score, preScore, sumScore, lnP, p, evalue, expectedDomains,
      regions, clustered, overlaps, envelopes, bestDomain, domains};
  }
  return {status, z, domZ, hits};
}

/** An instantiated engine. One instance is single-threaded; use one per worker. */
export class HmmerEngine {
  private readonly x: Exports;
  private constructor(exports: Exports) {
    this.x = exports;
  }

  static async instantiate(source: BufferSource | WebAssembly.Module): Promise<HmmerEngine> {
    const instance = source instanceof WebAssembly.Module ?
      await WebAssembly.instantiate(source, {}) :
      (await WebAssembly.instantiate(source, {})).instance;
    const engine = new HmmerEngine(instance.exports as unknown as Exports);
    const version = engine.x.hw_version();
    if (version !== ABI_VERSION)
      throw new Error(`HMMER engine ABI ${version}, expected ${ABI_VERSION}`);
    return engine;
  }

  private put(bytes: Uint8Array): number {
    const pointer = this.x.hw_alloc(bytes.length);
    if (pointer === 0) throw new Error('HMMER engine: out of memory');
    new Uint8Array(this.x.memory.buffer, pointer, bytes.length).set(bytes);
    return pointer;
  }

  /** Copy a returned buffer out of linear memory and free it. */
  private take(pointer: number): Uint8Array {
    if (pointer === 0) throw new Error(`HMMER engine: ${this.lastError()}`);
    const length = this.x.hw_len(pointer);
    const bytes = new Uint8Array(this.x.memory.buffer, pointer, length).slice();
    this.x.hw_free(pointer);
    return bytes;
  }

  private lastError(): string {
    const pointer = this.x.hw_error();
    const length = this.x.hw_len(pointer);
    const text = utf8.decode(new Uint8Array(this.x.memory.buffer, pointer, length));
    this.x.hw_free(pointer);
    return text;
  }

  /** Load a model database: `.h3m` (hmmpress) or ASCII HMMER text bytes. */
  loadDatabase(bytes: Uint8Array): HmmDatabase {
    const input = this.put(bytes);
    const handle = this.x.hw_db_load(input);
    this.x.hw_free(input);
    if (handle === 0) throw new Error(`HMMER engine: ${this.lastError()}`);
    const r = new Reader(this.take(this.x.hw_db_info(handle)));
    const alphabet = (['amino', 'dna', 'rna'] as const)[r.u32()];
    const models: ModelInfo[] = new Array(r.u32());
    for (let i = 0; i < models.length; i++) {
      const length = r.u32();
      const [name, accession, description, consensus, reference] = [r.text(), r.text(), r.text(), r.text(), r.text()];
      const pair = (): [number, number] | null => {
        const a = r.f32();
        const b = r.f32();
        return Number.isNaN(a) ? null : [a, b];
      };
      const [gathering, trusted, noise] = [pair(), pair(), pair()];
      models[i] = {index: i, length, name, accession, description, consensus, reference,
        gathering, trusted, noise, maxLength: r.u32()};
    }
    return new HmmDatabase(this, this.x, handle, {alphabet, models});
  }

  /** @internal */
  runScan(handle: number, sequences: SequenceInput[], options: SearchOptions): QueryResult[] {
    const input = this.put(packSequences(sequences));
    const record = this.put(packOptions(options));
    const result = this.x.hw_scan(handle, input, record);
    this.x.hw_free(input);
    this.x.hw_free(record);
    const r = new Reader(this.take(result));
    if (r.u32() !== RESULT_MAGIC) throw new Error('HMMER engine: bad result buffer');
    const results: QueryResult[] = new Array(r.u32());
    for (let q = 0; q < results.length; q++) results[q] = readQuery(r);
    return results;
  }

  /** @internal */
  runSearch(handle: number, model: number, sequences: SequenceInput[], options: SearchOptions): QueryResult {
    const input = this.put(packSequences(sequences));
    const record = this.put(packOptions(options));
    const result = this.x.hw_search(handle, model, input, record);
    this.x.hw_free(input);
    this.x.hw_free(record);
    const r = new Reader(this.take(result));
    if (r.u32() !== RESULT_MAGIC || r.u32() !== 1) throw new Error('HMMER engine: bad result buffer');
    return readQuery(r);
  }
}

/** A loaded model database. */
export class HmmDatabase {
  readonly info: DatabaseInfo;
  private readonly engine: HmmerEngine;
  private readonly x: Exports;
  private handle: number;

  /** @internal */
  constructor(engine: HmmerEngine, exports: Exports, handle: number, info: DatabaseInfo) {
    this.engine = engine;
    this.x = exports;
    this.handle = handle;
    this.info = info;
  }

  /** `hmmscan`: each sequence against every model (Z = number of models). */
  scan(sequences: SequenceInput[], options: SearchOptions = {}): QueryResult[] {
    return this.engine.runScan(this.live(), sequences, options);
  }

  /** `hmmsearch`: one model against every sequence (Z = number of sequences,
   * or `options.z` when one search is split across workers). */
  search(model: number, sequences: SequenceInput[], options: SearchOptions = {}): QueryResult {
    return this.engine.runSearch(this.live(), model, sequences, options);
  }

  free(): void {
    if (this.handle !== 0) this.x.hw_db_free(this.handle);
    this.handle = 0;
  }

  private live(): number {
    if (this.handle === 0) throw new Error('HMMER database was freed');
    return this.handle;
  }
}
