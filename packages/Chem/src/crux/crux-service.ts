import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {WorkerMessageBusClient} from '../worker-message-bus-client';
import {getRdKitService} from '../utils/chem-common-rdkit';
import {Fingerprint, hasNewLines} from '../utils/chem-common';
import {compileCrux, loadCrux} from '../engine/crux-module';
import {CruxCollection, CruxWorkerApi, DescriptorValues} from './crux-api';
import type {IFpResult} from '../rdkit-service/rdkit-service-worker-similarity';
import type {RuleId} from '../panels/structural-alerts';

const CRUX_INDEX_TAG = '.crux-index';
/** Rows per segment in the first round over the workers: small, so the first hits come early; doubles each round. */
const FIRST_SEGMENT_SIZE = 1000;
const MIN_SEGMENT_SIZE = 500;
const MAX_SEGMENT_SIZE = 50000;
/** Molecules per job (descriptors, fingerprints): a column's molecules are split across the workers in such chunks. */
const CHUNK_SIZE = 1000;
/**
 * A job's time limit: this much, plus its kind's time per molecule (`JOB_MS_PER_MOLECULE`); past it the workers
 * restart and the caller's RDKit path (the container, the script, the RDKit pool) runs.
 */
const JOB_BASE_MS = 30000;
/** How a job's timeout starts its message. */
const JOB_TIMEOUT = 'a Crux job took longer than';
const JOB_MS_PER_MOLECULE = {descriptors: 500, fingerprints: 20, alerts: 50, scaffoldTree: 200};
/**
 * crux keeps ~380 bytes per indexed molecule (~1.1 GB at the limit, spread over the workers); least recently used
 * column indexes are dropped past this many rows.
 */
const MAX_INDEXED_ROWS = 3000000;
/**
 * At most this many crux workers. Each holds a crux WebAssembly instance, and a page has room for only so many
 * WebAssembly memories, ~120 in Chrome (shared by its workers): sized by `hardwareConcurrency` alone, crux's pool and
 * RDKit's (`RdKitService`) filled it on a 32-thread machine, and past ~64 threads they would not both start. Past
 * this many, the searches gain little.
 */
const MAX_WORKERS = 16;

class CruxWorkerClient extends WorkerMessageBusClient {
  constructor() {
    super(new Worker(new URL('./crux.worker', import.meta.url)));
  }

  moduleInit(module: WebAssembly.Module): Promise<unknown> {
    const init = this.call('module::init', [module], false);
    this._ready = init.catch((e) => console.error('Chem | crux worker init failed:', e));
    return init;
  }

  build = (key: string, smiles: string[]): Promise<Uint32Array> => this.call('build', [key, smiles]);

  search = (key: string, smarts: string): Promise<Uint32Array> => this.call('search', [key, smarts]);

  drop = (keys: string[]): Promise<void> => this.call('drop', [keys]);

  calcDescriptors = (molecules: (string | null)[], names: string[]): Promise<{[name: string]: DescriptorValues}> =>
    this.call('calcDescriptors', [molecules, names]);

  scaffoldTree = (smiles: string[], ringCutoff: number, discharge: boolean): Promise<string> =>
    this.call('scaffoldTree', [smiles, ringCutoff, discharge]);

  getFingerprints = (type: Fingerprint, molecules: string[], smiles: boolean):
    Promise<IFpResult & {unread: number[]}> => this.call('getFingerprints', [type, molecules, smiles]);

  getStructuralAlerts = (alerts: {[rule in RuleId]?: string[]}, molecules: string[]):
    Promise<{alerts: {[rule in RuleId]?: boolean[]}, unread: number[]}> =>
    this.call('getStructuralAlerts', [alerts, molecules]);

  failForTests = (message: string): Promise<unknown> => this.call('fail::test', [message]);

  busyForTests = (ms: number): Promise<unknown> => this.call('busy::test', [ms]);
}

/** Contiguous rows of a column indexed by one crux collection that lives in one worker. */
export interface CruxSegment {
  key: string;
  start: number;
  end: number;
  worker: number;
  /** Of the column values the segment was built from. */
  hash: number;
  /** Molecules to build the collection from, released once the worker has built it. */
  smiles: string[] | null;
  /** Segment-local rows crux could not parse, and their molecules ('' once RDKit could not parse one either). */
  failed: Uint32Array | null;
  unparsed: string[] | null;
  built: Promise<void> | null;
}

export interface CruxColumnIndex {
  id: number;
  length: number;
  version: number;
  segments: CruxSegment[];
  lastUsed: number;
  /** Running searches and builds; segments replaced meanwhile are freed once it drops to zero. */
  users: number;
  released: CruxSegment[];
  /** Forgotten by a {@link CruxService.restart}: the searches still running on it fail. */
  lost?: boolean;
}

/**
 * Crux substructure indexes of molecule columns, split into segments over a pool of web workers.
 * A segment is built lazily by the first search (or warm-up) that reaches it, so the first hits come
 * while the rest of the column is still being indexed.
 */
export class CruxService {
  private module: WebAssembly.Module | null = null;
  private workers: CruxWorkerClient[] = [];
  /** Counts {@link restart}s: a search that failed on older workers does not restart the new ones. */
  generation = 0;
  private validator: CruxCollection | null = null;
  private indexes = new Map<number, CruxColumnIndex>();
  private lastId = 0;
  private indexLock: Promise<unknown> = Promise.resolve();

  static async create(): Promise<CruxService> {
    const service = new CruxService();
    // Chem's one Crux file, compiled once (engine/crux-module.ts) and shared, so the workers only instantiate it
    const crux = await loadCrux('search') as unknown as CruxWorkerApi;
    const module = await compileCrux();
    service.validator = new crux.CollectionBuilder(true).finish();
    service.module = module;
    try {
      await service.startWorkers();
    } catch (e) {
      // The workers that did start are let go: the next search creates the service again, and each failed attempt
      // used to leave a pool running (hundreds of workers once the page had no room for another WebAssembly memory)
      for (const worker of service.workers)
        worker.terminate();
      service.validator?.free();
      throw e;
    }
    grok.events.onTableRemoved.subscribe((e) => {
      for (const col of e.args.dataFrame.columns) {
        const index = service.indexes.get(col.temp[CRUX_INDEX_TAG]);
        if (index)
          service.dropIndex(index);
      }
    });
    return service;
  }

  private startWorkers(): Promise<unknown> {
    const workerCount = Math.max(1, Math.min(navigator.hardwareConcurrency - 2, MAX_WORKERS));
    this.workers = Array.from({length: workerCount}, () => new CruxWorkerClient());
    return Promise.all(this.workers.map((w) => w.moduleInit(this.module!)));
  }

  /**
   * Replaces the workers of `generation`, which frees all their memory, and forgets every index: the searches still
   * running on them fail, the next ones index their columns again.
   */
  restart(generation: number): void {
    if (generation !== this.generation)
      return;
    this.generation++;
    for (const worker of this.workers)
      worker.terminate();
    for (const index of this.indexes.values())
      index.lost = true;
    this.indexes.clear();
    // a failed worker init is logged by the worker client
    this.startWorkers().catch(() => {});
  }

  /** The workers indexing and searching: `hardwareConcurrency - 2`, at most {@link MAX_WORKERS}. */
  get workerCount(): number {
    return this.workers.length;
  }

  isValidSmarts(smarts: string): boolean {
    try {
      this.validator!.substructureSearch(smarts, 0);
      return true;
    } catch {
      return false;
    }
  }

  /** The column's index, created or refreshed (only the segments whose molecules changed) as needed. */
  getIndex(col: DG.Column): Promise<CruxColumnIndex> {
    const index = this.indexLock.then(() => this.prepareIndex(col));
    this.indexLock = index.catch(() => {});
    return index;
  }

  private async prepareIndex(col: DG.Column): Promise<CruxColumnIndex> {
    const version = col.version;
    let index = this.indexes.get(col.temp[CRUX_INDEX_TAG]);
    if (index?.length !== col.length || index.version !== version) {
      const values = readColumn(col);
      if (index?.length === col.length)
        await this.replaceChangedSegments(index, values);
      else {
        if (index)
          this.dropIndex(index);
        index = await this.createIndex(values);
        col.temp[CRUX_INDEX_TAG] = index.id;
        this.dropLeastRecentlyUsed(index);
      }
      index.version = version;
    }
    index.lastUsed = performance.now();
    return index;
  }

  /** Searches every segment, calling `onSegment` with its segment-local hits as soon as the segment is done. */
  async search(index: CruxColumnIndex, smarts: string, onSegment: (segment: CruxSegment, hits: Uint32Array) => void,
    isCancelled: () => boolean): Promise<void> {
    await this.use(index, (segment, worker) => this.ensureBuilt(segment)
      .then(() => isCancelled() ? null : worker.search(segment.key, smarts))
      .then((hits) => {
        if (hits && !isCancelled())
          onSegment(segment, hits);
      }), isCancelled);
  }

  /** Builds the segments not built yet. */
  async build(index: CruxColumnIndex, isCancelled: () => boolean): Promise<void> {
    await this.use(index, (segment) => this.ensureBuilt(segment), isCancelled);
  }

  /**
   * Runs `action` over the index segments, one at a time per worker, in row order, until the index is dropped; fails
   * if a restart forgets it.
   */
  private async use(index: CruxColumnIndex, action: (segment: CruxSegment, worker: CruxWorkerClient) => Promise<void>,
    isCancelled: () => boolean): Promise<void> {
    const segments = index.segments;
    index.users++;
    try {
      await Promise.all(this.workers.map(async (worker, w) => {
        for (const segment of segments) {
          if (segment.worker !== w)
            continue;
          if (index.lost)
            throw new Error('crux restarted');
          if (isCancelled() || !this.indexes.has(index.id))
            return;
          await action(segment, worker);
        }
      }));
    } finally {
      if (--index.users === 0 && !index.lost)
        this.release(index);
    }
  }

  /**
   * Descriptor values as the chem-chem container's `POST /chem/descriptors` gives them (`{name: {type, value}}`, group
   * names expanded, NaN and infinities null), the molecules split across the workers in chunks.
   */
  async calcDescriptors(molecules: (string | null)[], names: string[]): Promise<{[name: string]: DescriptorValues}> {
    const parts = await this.inChunks(molecules, JOB_MS_PER_MOLECULE.descriptors,
      (worker, chunk) => worker.calcDescriptors(chunk, names));
    const result: {[name: string]: DescriptorValues} = {};
    for (const name of Object.keys(parts[0])) {
      // Python's type of the first value ('double' for a column without one): the first chunk with a value says it
      const typed = parts.find((part) => part[name].value.some((v) => v !== null)) ?? parts[0];
      result[name] = {type: typed[name].type, value: parts.flatMap((part) => part[name].value)};
    }
    return result;
  }

  /**
   * Fingerprints as the RDKit pool's `getFingerprints` gives them (the same worker code, run on Crux): one per molecule
   * (null where none), the canonical SMILES when asked, and the rows whose molecule Crux cannot read (`unread`).
   */
  async getFingerprints(type: Fingerprint, molecules: string[], smiles: boolean):
    Promise<IFpResult & {unread: number[]}> {
    const parts = await this.inChunks(molecules, JOB_MS_PER_MOLECULE.fingerprints,
      (worker, chunk) => worker.getFingerprints(type, chunk, smiles));
    return {fps: parts.flatMap((part) => part.fps), smiles: smiles ? parts.flatMap((part) => part.smiles!) : null,
      unread: parts.flatMap((part, k) => part.unread.map((row) => k * CHUNK_SIZE + row))};
  }

  /**
   * Which molecules match each group's SMARTS, as the RDKit pool's `getStructuralAlerts` gives it (the same worker
   * code, run on Crux): `[group, a flag per molecule]` per group, and the rows whose molecule Crux cannot read.
   */
  async getStructuralAlerts(alerts: {[rule in RuleId]?: string[]}, molecules: string[]):
    Promise<{alerts: [RuleId, boolean[]][], unread: number[]}> {
    const parts = await this.inChunks(molecules, JOB_MS_PER_MOLECULE.alerts,
      (worker, chunk) => worker.getStructuralAlerts(alerts, chunk));
    const result: {[rule in RuleId]?: boolean[]} = {};
    for (const part of parts) {
      for (const rule of Object.keys(part.alerts) as RuleId[])
        result[rule] = (result[rule] ?? []).concat(part.alerts[rule]!);
    }
    return {alerts: Object.entries(result) as [RuleId, boolean[]][],
      unread: parts.flatMap((part, k) => part.unread.map((row) => k * CHUNK_SIZE + row))};
  }

  /**
   * `job` over the molecules in chunks of {@link CHUNK_SIZE}, the results in chunk order. A worker runs its calls one
   * after another, so each worker is handed its next chunk when the last one is done, and a chunk's time limit runs
   * from its own start; a failed chunk ends the call (the chunks not started are not run).
   */
  private async inChunks<M, R>(molecules: M[], msPerMolecule: number,
    job: (worker: CruxWorkerClient, chunk: M[]) => Promise<R>, chunkSize = CHUNK_SIZE, baseMs = JOB_BASE_MS):
    Promise<R[]> {
    const chunks: M[][] = [];
    for (let start = 0; start < molecules.length; start += chunkSize)
      chunks.push(molecules.slice(start, start + chunkSize));
    if (chunks.length === 0)
      chunks.push([]);
    const results: R[] = new Array(chunks.length);
    let next = 0;
    let failed = false;
    await Promise.all(this.workers.slice(0, chunks.length).map(async (worker) => {
      while (!failed && next < chunks.length) {
        const i = next++;
        try {
          results[i] = await this.job(job(worker, chunks[i]), baseMs + chunks[i].length * msPerMolecule);
        } catch (e) {
          failed = true;
          throw e;
        }
      }
    }));
    return results;
  }

  /** scripts/scaffold_tree.py's JSON of the molecules, from one worker. */
  scaffoldTree(smiles: string[], ringCutoff: number, discharge: boolean): Promise<string> {
    return this.job(this.workers[0].scaffoldTree(smiles, ringCutoff, discharge),
      JOB_BASE_MS + smiles.length * JOB_MS_PER_MOLECULE.scaffoldTree);
  }

  /**
   * A job's result within `timeoutMs`. A wasm trap (the worker reports `RuntimeError: ...`; out of memory), a timeout
   * or a worker gone restarts the workers: they may hold a broken instance or a call that never ends; the next jobs and
   * searches get new ones. Another failure (a molecule Crux refuses, a "not supported by Crux" call) keeps them and
   * their search indexes.
   */
  private async job<R>(work: Promise<R>, timeoutMs: number): Promise<R> {
    const generation = this.generation;
    let timer: ReturnType<typeof setTimeout> | undefined;
    const timeout = new Promise<never>((_, reject) =>
      timer = setTimeout(() => reject(new Error(`${JOB_TIMEOUT} ${timeoutMs} ms`)), timeoutMs));
    try {
      return await Promise.race([work, timeout]);
    } catch (e: any) {
      if (/RuntimeError|unreachable|out of bounds|memory|Worker terminated/i.test(String(e?.message ?? e)) ||
        String(e?.message ?? e).startsWith(JOB_TIMEOUT))
        this.restart(generation);
      throw e;
    } finally {
      clearTimeout(timer);
    }
  }

  /** Tests: lets the workers go (the service is not used after). */
  terminateForTests(): void {
    for (const worker of this.workers)
      worker.terminate();
    for (const index of this.indexes.values())
      index.lost = true;
    this.indexes.clear();
  }

  /** Tests: a job that fails in its worker with `message` (by default as a trapped instance would: a restart). */
  failJobForTests(message = 'simulated: RuntimeError: unreachable'): Promise<unknown> {
    return this.job(this.workers[0].failForTests(message), JOB_BASE_MS);
  }

  /**
   * Tests: `chunks` chunks of one molecule each, every one a job that keeps its worker busy `ms` (as a long Crux call
   * does), with a limit of `limitMs` per chunk; the number of chunks done.
   */
  async busyChunksForTests(chunks: number, ms: number, limitMs: number): Promise<number> {
    const done = await this.inChunks(Array.from({length: chunks}, (_, i) => i), 0, (worker) => worker.busyForTests(ms),
      1, limitMs);
    return done.length;
  }

  private ensureBuilt(segment: CruxSegment): Promise<void> {
    segment.built ??= this.workers[segment.worker].build(segment.key, segment.smiles!).then((failed) => {
      segment.failed = failed;
      segment.unparsed = Array.from(failed, (row) => segment.smiles![row]);
      segment.smiles = null;
    }, (e) => {
      segment.built = null;
      throw e;
    });
    return segment.built;
  }

  private async createIndex(values: string[]): Promise<CruxColumnIndex> {
    const smiles = await toCruxSmiles(values);
    const segments: CruxSegment[] = [];
    const workerCount = this.workers.length;
    for (let start = 0, size = FIRST_SEGMENT_SIZE; start < values.length; size = Math.min(MAX_SEGMENT_SIZE, size * 2)) {
      // the last round splits what is left evenly over the workers
      const roundSize = Math.max(MIN_SEGMENT_SIZE, Math.min(size, Math.ceil((values.length - start) / workerCount)));
      for (let worker = 0; worker < workerCount && start < values.length; worker++) {
        const end = Math.min(values.length, start + roundSize);
        segments.push(this.createSegment(values, smiles.slice(start, end), start, worker));
        start = end;
      }
    }
    const index: CruxColumnIndex = {id: ++this.lastId, length: values.length, version: -1, segments, lastUsed: 0,
      users: 0, released: []};
    this.indexes.set(index.id, index);
    return index;
  }

  private createSegment(values: string[], smiles: string[], start: number, worker: number): CruxSegment {
    const end = start + smiles.length;
    return {key: `${++this.lastId}`, start, end, worker, hash: hashStrings(values, start, end), smiles,
      failed: null, unparsed: null, built: null};
  }

  private async replaceChangedSegments(index: CruxColumnIndex, values: string[]): Promise<void> {
    const changed = index.segments.filter((s) => hashStrings(values, s.start, s.end) !== s.hash);
    const smiles = await toCruxSmiles(changed.flatMap((s) => values.slice(s.start, s.end)));
    let offset = 0;
    // running searches keep iterating the segments they started with, the replaced ones are freed after them
    index.segments = index.segments.map((segment) => {
      if (!changed.includes(segment))
        return segment;
      index.released.push(segment);
      const length = segment.end - segment.start;
      offset += length;
      return this.createSegment(values, smiles.slice(offset - length, offset), segment.start, segment.worker);
    });
    if (index.users === 0)
      this.release(index);
  }

  private dropLeastRecentlyUsed(current: CruxColumnIndex): void {
    let rows = 0;
    for (const index of this.indexes.values())
      rows += index.length;
    const byAge = Array.from(this.indexes.values())
      .filter((i) => i !== current && i.users === 0)
      .sort((a, b) => a.lastUsed - b.lastUsed);
    for (const index of byAge) {
      if (rows <= MAX_INDEXED_ROWS)
        break;
      rows -= index.length;
      this.dropIndex(index);
    }
  }

  private dropIndex(index: CruxColumnIndex): void {
    this.indexes.delete(index.id);
    index.released.push(...index.segments);
    if (index.users === 0)
      this.release(index);
  }

  private release(index: CruxColumnIndex): void {
    const keys: string[][] = this.workers.map(() => []);
    for (const segment of index.released.splice(0)) {
      if (segment.built)
        keys[segment.worker].push(segment.key);
    }
    for (let w = 0; w < this.workers.length; w++) {
      if (keys[w].length > 0)
        this.workers[w].drop(keys[w]).catch((e) => console.warn('Chem | crux segments were not freed:', e));
    }
  }
}

/** Column values, '' for empty cells and with CXSMILES extensions cut off. */
function readColumn(col: DG.Column): string[] {
  const values = Array.from(col.toList() as (string | null)[], (v) => v ?? '');
  for (let i = 0; i < values.length; i++) {
    const extension = hasNewLines(values[i]) ? -1 : values[i].indexOf(' |');
    if (extension > 0)
      values[i] = values[i].substring(0, extension);
  }
  return values;
}

/** The molecules as SMILES crux can parse: molblocks are converted by RDKit. */
async function toCruxSmiles(values: string[]): Promise<string[]> {
  const molblockRows: number[] = [];
  for (let i = 0; i < values.length; i++) {
    if (hasNewLines(values[i]))
      molblockRows.push(i);
  }
  if (molblockRows.length === 0)
    return values;
  const converted = await (await getRdKitService())
    .convertMolNotation(molblockRows.map((i) => values[i]), DG.chem.Notation.Smiles);
  const smiles = values.slice();
  for (let k = 0; k < molblockRows.length; k++)
    smiles[molblockRows[k]] = converted[k];
  return smiles;
}

function hashStrings(values: string[], start: number, end: number): number {
  let hash = 0x811c9dc5;
  for (let row = start; row < end; row++) {
    const value = values[row];
    for (let i = 0; i < value.length; i++)
      hash = Math.imul(hash ^ value.charCodeAt(i), 0x01000193);
    hash = Math.imul(hash ^ 0x0a, 0x01000193);
  }
  return hash >>> 0;
}
