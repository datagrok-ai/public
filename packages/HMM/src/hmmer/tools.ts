// Synced from crux-core crates/hmmer-web/web/tools.ts by crates/hmmer-web/web/sync-datagrok.mjs — do not edit here.
/** Typed wrapper of the lazily loaded tools module (crates/hmmer-web-tools,
 * web/dist/hmmer_tools.wasm): phmmer query models, hmmalign (traces per chunk
 * of sequences, then the alignment from all traces) and hmmlogo values. Like
 * engine.ts it runs in any JS runtime without imports; one instance is
 * single-threaded. */

export const TOOLS_ABI_VERSION = 1;
const RESULT_MAGIC = 0x31544d48; // "HMT1"

interface Exports {
  memory: WebAssembly.Memory;
  ht_version(): number;
  ht_alloc(length: number): number;
  ht_free(pointer: number): void;
  ht_len(pointer: number): number;
  ht_error(): number;
  ht_seqmodel(request: number): number;
  ht_logo(models: number, mode: number): number;
  ht_align_traces(request: number): number;
  ht_align_rows(request: number): number;
  ht_extract(models: number, index: number): number;
  ht_build(request: number): number;
}

/** A profile built by hmmbuild from aligned rows. */
export interface BuiltProfile {
  /** The HMMER3/f text C hmmbuild writes (DATE: `Thu Jan  1 00:00:00 1970`). */
  text: string;
  /** hmmbuild's summary output (header and table). */
  summary: string;
}

/** hmmalign's alignment as `--outformat afa` prints it (unwrapped rows). */
export interface ProfileAlignment {
  /** One aligned row per input sequence: upper-case residues in consensus
   * columns, lower-case insertions, `-` deletions, `.` insert padding. */
  rows: string[];
  /** Per column: true for consensus (match) columns, the j-th holding model node j. */
  consensus: boolean[];
}

/** A request builder for the tools ABI (little-endian words and length-prefixed bytes). */
class Request {
  private readonly parts: Uint8Array[] = [];
  word(n: number): this {
    const bytes = new Uint8Array(4);
    new DataView(bytes.buffer).setUint32(0, n, true);
    this.parts.push(bytes);
    return this;
  }
  real(x: number): this {
    const bytes = new Uint8Array(8);
    new DataView(bytes.buffer).setFloat64(0, x, true);
    this.parts.push(bytes);
    return this;
  }
  raw(bytes: Uint8Array): this {
    this.parts.push(bytes);
    return this;
  }
  bytes(bytes: Uint8Array): this {
    return this.word(bytes.length).raw(bytes);
  }
  text(s: string): this {
    return this.bytes(encoder.encode(s));
  }
  sequences(sequences: {name: string; residues: string}[]): this {
    this.word(sequences.length);
    for (const s of sequences) this.text(s.name).text(s.residues);
    return this;
  }
  finish(): Uint8Array {
    const out = new Uint8Array(this.parts.reduce((n, p) => n + p.length, 0));
    let at = 0;
    for (const p of this.parts) {
      out.set(p, at);
      at += p.length;
    }
    return out;
  }
}

/** phmmer's query model options (`--popen`, `--pextend`, `--mx`, `--seed`). */
export interface QueryModelOptions {
  popen?: number;
  pextend?: number;
  /** Built-in score matrix (BLOSUM45/50/62/80/90, PAM30/70/120/240); default BLOSUM62. */
  matrix?: string;
  /** Builder RNG seed; default 42 as phmmer. */
  seed?: number;
}

/** `hmmlogo --height_relent_all` (default), `--height_relent_abovebg`, `--height_score`. */
export type LogoMode = 'relent_all' | 'relent_abovebg' | 'score';
const LOGO_MODES: LogoMode[] = ['relent_all', 'relent_abovebg', 'score'];

/** hmmlogo values of one model; arrays are indexed by match position 0..M-1. */
export interface Logo {
  /** Maximum expected height (bits), or NaN in score mode. */
  maxHeight: number;
  /** `heights[position][residue]` in alphabet order. */
  heights: Float32Array[];
  relativeEntropy: Float32Array;
  insertProbability: Float32Array;
  insertLength: Float32Array;
  occupancy: Float32Array;
}

const utf8 = new TextDecoder();
const latin1 = new TextDecoder('latin1');
const encoder = new TextEncoder();

export class HmmerTools {
  private readonly x: Exports;
  private constructor(exports: Exports) {
    this.x = exports;
  }

  static async instantiate(source: BufferSource | WebAssembly.Module): Promise<HmmerTools> {
    const instance = source instanceof WebAssembly.Module ?
      await WebAssembly.instantiate(source, {}) :
      (await WebAssembly.instantiate(source, {})).instance;
    const tools = new HmmerTools(instance.exports as unknown as Exports);
    const version = tools.x.ht_version();
    if (version !== TOOLS_ABI_VERSION)
      throw new Error(`HMMER tools ABI ${version}, expected ${TOOLS_ABI_VERSION}`);
    return tools;
  }

  private put(bytes: Uint8Array): number {
    const pointer = this.x.ht_alloc(bytes.length);
    if (pointer === 0) throw new Error('HMMER tools: out of memory');
    new Uint8Array(this.x.memory.buffer, pointer, bytes.length).set(bytes);
    return pointer;
  }

  /** Grow linear memory to hold `bytes` more (allocated and freed: the allocator keeps it). */
  private reserve(bytes: number): void {
    const pointer = this.x.ht_alloc(bytes);
    if (pointer !== 0) this.x.ht_free(pointer);
  }

  /** Copy a returned buffer out of linear memory and free it. */
  private take(pointer: number): Uint8Array {
    if (pointer === 0) {
      const error = this.x.ht_error();
      const text = utf8.decode(new Uint8Array(this.x.memory.buffer, error, this.x.ht_len(error)));
      this.x.ht_free(error);
      throw new Error(`HMMER tools: ${text}`);
    }
    const bytes = new Uint8Array(this.x.memory.buffer, pointer, this.x.ht_len(pointer)).slice();
    this.x.ht_free(pointer);
    return bytes;
  }

  private call(input: Uint8Array, f: (pointer: number) => number): Uint8Array {
    const pointer = this.put(input);
    try {
      return this.take(f(pointer));
    } finally {
      this.x.ht_free(pointer);
    }
  }

  /** phmmer's calibrated query models of `queries`, as one `.h3m` database
   * (load it into the search engine and search: identical to phmmer). */
  queryModels(queries: {name: string; residues: string}[], options: QueryModelOptions = {}): Uint8Array {
    const request = new Request().word(options.seed ?? 0).real(options.popen ?? -1).real(options.pextend ?? -1)
      .text(options.matrix ?? '').sequences(queries).finish();
    return this.call(request, (pointer) => this.x.ht_seqmodel(pointer));
  }

  /** hmmalign traces of `sequences` against model `index` of `models` (opaque
   * packed bytes; chunks of one list may be traced separately, in any worker). */
  alignTraces(models: Uint8Array, index: number, sequences: {name: string; residues: string}[]): Uint8Array {
    const request = new Request().bytes(models).word(index).sequences(sequences).finish();
    return this.call(request, (pointer) => this.x.ht_align_traces(pointer));
  }

  /** The alignment of `sequences` from their traces (chunks in order), as
   * `hmmalign [--trim] --outformat afa` prints it. */
  alignRows(models: Uint8Array, index: number, sequences: {name: string; residues: string}[],
    traceChunks: Uint8Array[], trim = false): ProfileAlignment {
    const request = new Request().bytes(models).word(index).word(trim ? 1 : 0).sequences(sequences);
    let count = 0;
    for (const chunk of traceChunks) {
      const view = new DataView(chunk.buffer, chunk.byteOffset, chunk.byteLength);
      if (view.getUint32(0, true) !== RESULT_MAGIC) throw new Error('HMMER tools: bad trace chunk');
      count += view.getUint32(4, true);
    }
    request.word(count);
    for (const chunk of traceChunks) request.raw(chunk.subarray(8));
    const bytes = this.call(request.finish(), (pointer) => this.x.ht_align_rows(pointer));
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    if (view.getUint32(0, true) !== RESULT_MAGIC) throw new Error('HMMER tools: bad alignment result');
    const n = view.getUint32(4, true);
    const alen = view.getUint32(8, true);
    const consensus = Array.from(bytes.subarray(12, 12 + alen), (b) => b === 1);
    const rows: string[] = [];
    for (let i = 0; i < n; i++) {
      const at = 12 + alen * (i + 1);
      rows.push(latin1.decode(bytes.subarray(at, at + alen)));
    }
    return {rows, consensus};
  }

  /** `hmmbuild [--amino|--dna|--rna] --informat afa -n name` of aligned rows
   * (equal lengths; `-`/`.` gaps). Search with the returned text to reproduce
   * `hmmbuild` followed by `hmmsearch`. */
  build(rows: {name: string; aligned: string}[], name: string, alphabet: 'amino' | 'dna' | 'rna' | null = 'amino'):
    BuiltProfile {
    const afa = rows.map((r) => `>${r.name}\n${r.aligned}\n`).join('');
    // A large alignment grows linear memory once (~21 bytes per cell needed), not in
    // many small steps: 100,000 rows x 120 columns built in 3.3 s on a fresh module, 0.6 s so.
    if (afa.length > 4_000_000) this.reserve(afa.length * 24);
    const code = alphabet === null ? 0 : ['amino', 'dna', 'rna'].indexOf(alphabet) + 1;
    const bytes = this.call(new Request().word(code).text(name).text(afa).finish(),
      (pointer) => this.x.ht_build(pointer));
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    if (view.getUint32(0, true) !== RESULT_MAGIC) throw new Error('HMMER tools: bad build result');
    const textLength = view.getUint32(4, true);
    const text = utf8.decode(bytes.subarray(8, 8 + textLength));
    const summaryLength = view.getUint32(8 + textLength, true);
    const summary = utf8.decode(bytes.subarray(12 + textLength, 12 + textLength + summaryLength));
    return {text, summary};
  }

  /** Model `index` of `models` (`.h3m` or HMMER text) as its own `.h3m` bytes. */
  extract(models: Uint8Array, index: number): Uint8Array {
    return this.call(models, (pointer) => this.x.ht_extract(pointer, index));
  }

  /** hmmlogo values of every model in `models` (`.h3m` or HMMER text). */
  logos(models: Uint8Array, mode: LogoMode = 'relent_all'): Logo[] {
    const bytes = this.call(models, (pointer) => this.x.ht_logo(pointer, LOGO_MODES.indexOf(mode)));
    const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    let at = 0;
    const u32 = () => (at += 4, view.getUint32(at - 4, true));
    const f32 = () => (at += 4, view.getFloat32(at - 4, true));
    if (u32() !== RESULT_MAGIC) throw new Error('HMMER tools: bad logo result');
    const logos: Logo[] = [];
    for (let n = u32(); n > 0; n--) {
      const m = u32();
      const k = u32();
      const logo: Logo = {maxHeight: f32(), heights: [], relativeEntropy: new Float32Array(m),
        insertProbability: new Float32Array(m), insertLength: new Float32Array(m), occupancy: new Float32Array(m)};
      for (let i = 0; i < m; i++) {
        const row = new Float32Array(k);
        for (let r = 0; r < k; r++) row[r] = f32();
        logo.heights.push(row);
        logo.relativeEntropy[i] = f32();
        logo.insertProbability[i] = f32();
        logo.insertLength[i] = f32();
        logo.occupancy[i] = f32();
      }
      logos.push(logo);
    }
    return logos;
  }
}
