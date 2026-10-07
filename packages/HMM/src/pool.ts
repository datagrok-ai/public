/** A pool of HMMER workers sharing one compiled WebAssembly module.
 *
 * Comparisons of different query sequences are independent, so a batch is
 * split into chunks that idle workers pull in turn; results are identical for
 * any number of workers. The pool starts on first use and terminates its
 * workers after a minute without work, releasing their memory. */
import * as DG from 'datagrok-api/dg';
import type {DatabaseInfo, QueryResult, SearchOptions, SequenceInput} from './hmmer/engine.ts';
import {mergeSearch, partOptions, type MergedSearch, type SearchPart} from './hmmer/search.ts';
import type {AnarciOptions, SequenceResult} from './hmmer/anarci/anarci.ts';
import type {WorkerRequest, WorkerResponse} from './workers/protocol';

const IDLE_MS = 60_000;
/** Below this many sequences per worker, extra workers cost more than they save. */
const MIN_CHUNK = 16;

interface PoolWorker {
  worker: Worker;
  pending: Map<number, {resolve: (value: unknown) => void; reject: (error: Error) => void}>;
  ready: Promise<unknown>;
  databases: Set<string>;
}

export class HmmerPool {
  private static instance: HmmerPool | null = null;
  private readonly base: string;
  private module: Promise<WebAssembly.Module> | null = null;
  private workers: PoolWorker[] = [];
  private nextId = 1;
  private active = 0;
  private idleTimer: ReturnType<typeof setTimeout> | null = null;
  /** Model databases by key, kept to load them into workers started later. */
  private readonly databases = new Map<string, ArrayBuffer>();

  private constructor(base: string) {
    this.base = base;
  }

  /** The shared pool; `base` is the package's dist/ URL. */
  static get(pkg: DG.Package): HmmerPool {
    return HmmerPool.instance ??= new HmmerPool(`${pkg.webRoot}dist/`);
  }

  /** Workers used for a batch: one per spare hardware thread, at most 8. */
  static get maxWorkers(): number {
    const threads = typeof navigator !== 'undefined' ? navigator.hardwareConcurrency || 2 : 2;
    return Math.max(1, Math.min(8, threads - 1));
  }

  private compile(): Promise<WebAssembly.Module> {
    return this.module ??= (async () => {
      const response = await fetch(this.base + 'hmmer_web.wasm');
      if (!response.ok) throw new Error(`Failed to load the HMMER engine (HTTP ${response.status})`);
      return WebAssembly.compile(await response.arrayBuffer());
    })().catch((e) => {
      this.module = null;
      throw e;
    });
  }

  private call<T>(w: PoolWorker, request: WorkerRequest, transfer: Transferable[] = []): Promise<T> {
    const id = this.nextId++;
    return new Promise<T>((resolve, reject) => {
      w.pending.set(id, {resolve: resolve as (value: unknown) => void, reject});
      w.worker.postMessage({id, request}, transfer);
    });
  }

  private async spawn(): Promise<PoolWorker> {
    const module = await this.compile();
    const worker = new Worker(new URL('./workers/hmmer.worker', import.meta.url));
    const w: PoolWorker = {worker, pending: new Map(), ready: Promise.resolve(), databases: new Set()};
    worker.onmessage = (event: MessageEvent<WorkerResponse>) => {
      const message = event.data;
      const entry = w.pending.get(message.id);
      if (!entry) return;
      w.pending.delete(message.id);
      if (message.ok) entry.resolve(message.result);
      else entry.reject(new Error(message.error));
    };
    worker.onerror = (event) => {
      for (const {reject} of w.pending.values()) reject(new Error(event.message || 'HMMER worker failed'));
      w.pending.clear();
    };
    w.ready = this.call(w, {op: 'init', module, base: this.base});
    await w.ready;
    return w;
  }

  /** Make sure `count` workers are running, with every registered database loaded. */
  private async ensure(count: number): Promise<PoolWorker[]> {
    const missing = count - this.workers.length;
    if (missing > 0)
      this.workers.push(...await Promise.all(Array.from({length: missing}, () => this.spawn())));
    const used = this.workers.slice(0, count);
    await Promise.all(used.map(async (w) => {
      for (const [key, bytes] of this.databases) {
        if (w.databases.has(key)) continue;
        await this.call(w, {op: 'load', key, bytes: bytes.slice(0)});
        w.databases.add(key);
      }
    }));
    return used;
  }

  private begin(): void {
    this.active++;
    if (this.idleTimer) clearTimeout(this.idleTimer);
    this.idleTimer = null;
  }

  private end(): void {
    if (--this.active > 0) return;
    this.idleTimer = setTimeout(() => this.terminate(), IDLE_MS);
  }

  /** Stop all workers (they restart on the next request). */
  terminate(): void {
    for (const w of this.workers) w.worker.terminate();
    this.workers = [];
  }

  /** Run `make(chunk, offset)` for consecutive chunks of `items` on up to `workers` workers. */
  async map<T, R>(items: T[], make: (chunk: T[], offset: number) => WorkerRequest, workers = HmmerPool.maxWorkers):
    Promise<R[]> {
    this.begin();
    try {
      const count = Math.max(1, Math.min(workers, Math.ceil(items.length / MIN_CHUNK)));
      const used = await this.ensure(count);
      // Several chunks per worker keep all of them busy until the end.
      const size = Math.max(1, Math.ceil(items.length / (count * 4)));
      const chunks: T[][] = [];
      const offsets: number[] = [];
      for (let i = 0; i < items.length; i += size) {
        chunks.push(items.slice(i, i + size));
        offsets.push(i);
      }
      const results: R[][] = new Array(chunks.length);
      let next = 0;
      await Promise.all(used.map(async (w) => {
        while (next < chunks.length) {
          const i = next++;
          results[i] = await this.call<R[]>(w, make(chunks[i], offsets[i]));
        }
      }));
      return results.flat();
    } finally {
      this.end();
    }
  }

  /** ANARCI numbering of `(name, sequence)` pairs. */
  anarci(sequences: [string, string][], options: AnarciOptions, workers?: number): Promise<SequenceResult[]> {
    return this.map(sequences, (chunk) => ({op: 'anarci', sequences: chunk, options}), workers);
  }

  /** Register a model database (HMMER text or .h3m bytes); returns its description. */
  async loadDatabase(key: string, bytes: Uint8Array): Promise<DatabaseInfo> {
    this.begin();
    try {
      const copy = bytes.slice().buffer;
      this.databases.set(key, copy);
      const [first] = await this.ensure(1);
      for (const w of this.workers) w.databases.delete(key);
      const info = await this.call<DatabaseInfo>(first, {op: 'load', key, bytes: copy.slice(0)});
      first.databases.add(key);
      return info;
    } catch (e) {
      this.databases.delete(key);
      throw e;
    } finally {
      this.end();
    }
  }

  hasDatabase(key: string): boolean {
    return this.databases.has(key);
  }

  /** `hmmscan` of every sequence against a registered database. */
  scan(key: string, sequences: SequenceInput[], options: SearchOptions = {}, workers?: number): Promise<QueryResult[]> {
    return this.map(sequences, (chunk) => ({op: 'scan', key, sequences: chunk, options}), workers);
  }

  /** `hmmsearch` of one database model against `sequences` (unique names),
   * split across workers and merged; identical to one search (hmmer/search.ts). */
  async search(key: string, model: number, sequences: {name: string; residues: string}[],
    options: SearchOptions = {}, workers?: number): Promise<MergedSearch> {
    const total = sequences.length;
    const parts = await this.map<{name: string; residues: string}, SearchPart>(sequences, (chunk, offset) =>
      ({op: 'search', key, model, sequences: chunk, options: partOptions(options, total), offset}), workers);
    return mergeSearch(parts, options, total, sequences.map((s) => s.name));
  }
}
