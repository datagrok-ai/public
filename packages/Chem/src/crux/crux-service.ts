import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {WorkerMessageBusClient} from '../worker-message-bus-client';
import {getRdKitService} from '../utils/chem-common-rdkit';
import {hasNewLines} from '../utils/chem-common';
//@ts-ignore
import initCrux, {Collection, CollectionBuilder} from './crux_wasm.js';

const CRUX_INDEX_TAG = '.crux-index';
/** Rows per segment in the first round over the workers: small, so the first hits come early; doubles each round. */
const FIRST_SEGMENT_SIZE = 1000;
const MIN_SEGMENT_SIZE = 500;
const MAX_SEGMENT_SIZE = 50000;
/** crux keeps ~1 KB per indexed molecule; least recently used column indexes are dropped past this many rows. */
const MAX_INDEXED_ROWS = 3000000;

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
}

/**
 * Crux substructure indexes of molecule columns, split into segments over a pool of web workers.
 * A segment is built lazily by the first search (or warm-up) that reaches it, so the first hits come
 * while the rest of the column is still being indexed.
 */
export class CruxService {
  private workers: CruxWorkerClient[] = [];
  private validator: Collection | null = null;
  private indexes = new Map<number, CruxColumnIndex>();
  private lastId = 0;
  private indexLock: Promise<unknown> = Promise.resolve();

  static async create(webRoot: string): Promise<CruxService> {
    const service = new CruxService();
    const response = await fetch(`${webRoot}/dist/crux_wasm_bg.wasm`);
    if (!response.ok)
      throw new Error(`crux wasm is not available: ${response.status} ${response.statusText}`);
    // compiled once and shared, so the workers only instantiate it
    const module = await WebAssembly.compile(await response.arrayBuffer());
    await initCrux({module_or_path: module});
    service.validator = new CollectionBuilder(true).finish();
    const workerCount = Math.max(1, navigator.hardwareConcurrency - 2);
    service.workers = Array.from({length: workerCount}, () => new CruxWorkerClient());
    await Promise.all(service.workers.map((w) => w.moduleInit(module)));
    grok.events.onTableRemoved.subscribe((e) => {
      for (const col of e.args.dataFrame.columns) {
        const index = service.indexes.get(col.temp[CRUX_INDEX_TAG]);
        if (index)
          service.dropIndex(index);
      }
    });
    return service;
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

  /** Runs `action` over the index segments, one at a time per worker, in row order, until the index is dropped. */
  private async use(index: CruxColumnIndex, action: (segment: CruxSegment, worker: CruxWorkerClient) => Promise<void>,
    isCancelled: () => boolean): Promise<void> {
    const segments = index.segments;
    index.users++;
    try {
      await Promise.all(this.workers.map(async (worker, w) => {
        for (const segment of segments) {
          if (segment.worker !== w)
            continue;
          if (isCancelled() || !this.indexes.has(index.id))
            return;
          await action(segment, worker);
        }
      }));
    } finally {
      if (--index.users === 0)
        this.release(index);
    }
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
