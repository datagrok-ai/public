import * as DG from 'datagrok-api/dg';
import * as grok from 'datagrok-api/grok';
import {WorkerMessageBusClient} from '../worker-message-bus-client';
import {Fingerprint, hasNewLines} from '../utils/chem-common';
import {compileCrux, loadCrux, noteCruxError} from '../engine/crux-module';
import {DescriptorTiming, INT_NONE, PackedDescriptorValues, isNoneValue} from './crux-api';
import {CLUSTER_MCS_TIMEOUT_S} from './crux-api';
import {getQueryMolSafe} from '../utils/mol-creation_rdkit';
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';
import type {IFpResult} from '../rdkit-service/rdkit-service-worker-similarity';
import type {RuleId} from '../panels/structural-alerts';
import type {FingerprintOptions} from '../utils/fingerprint-options';
import type {IMmpFragmentsResult, IRGroupAnalysisResult} from '../rdkit-service/rdkit-service-worker-substructure';
import type {RGroupCluster} from './crux-analysis-jobs';

const CRUX_INDEX_TAG = '.crux-index';
/** Rows per segment in the first round over the workers: small, so the first hits come early; doubles each round. */
const FIRST_SEGMENT_SIZE = 1000;
const MIN_SEGMENT_SIZE = 500;
const MAX_SEGMENT_SIZE = 50000;
/** The most molecules per job (descriptors, fingerprints, alerts): a column's molecules are split across the
 * workers in chunks of this size, or smaller ones when there are fewer than this many per worker, so that every worker
 * works. */
const CHUNK_SIZE = 1000;
/**
 * A job's time limit: this much, plus its kind's time per molecule (`JOB_MS_PER_MOLECULE`); past it the workers
 * restart and the caller's RDKit path (the container, the script, the RDKit pool) runs.
 */
const JOB_BASE_MS = 30000;
/** How a job's timeout starts its message. */
const JOB_TIMEOUT = 'a Crux job took longer than';
const JOB_MS_PER_MOLECULE = {descriptors: 500, fingerprints: 20, alerts: 50, scaffoldTree: 200, properties: 20};
// Rascal's options for an activity cliff's MCS (CruxService.pairUncommonAtomsAndBonds): RDKit's defaults but these
// three, measured on the cliffs of Chem's tests and demo at their cutoffs (2026-10-21)
/** get_mces' 0.5 gave one cliff of the cutoff-20 test no MCES, 0.4 none. */
const PAIR_MCS_SIMILARITY_THRESHOLD = 0.3;
/** RDKit's 1000 gave 22 of 24 triterpenoid and 185 of 185 peptide cliffs no MCES ("too many matching bond pairs"),
 * 5000 still 156 peptide ones, 10000 none (with 16 MB of wasm memory). */
const PAIR_MCS_MAX_BOND_MATCH_PAIRS = 10000;
/** Seconds; RDKit's 60 is over a Crux job's limit, which restarts the pool. Past it Rascal gives its best MCES so far
 * (36 of the 185 peptide cliffs, 6.9 s at most with the setup; every other cliff took under 0.1 s). */
const PAIR_MCS_TIMEOUT_S = 5;
/**
 * The analyses' time limits (CruxService.mmpGetFragments and the rest): per molecule (`mmpa`, `rGroups`, `link`,
 * `molzip`) or per cluster, a cluster being one job (`clusterMCS`: FMCS's Timeout on Crux, crux-analysis-jobs.ts;
 * `clusterRGroups`: RdKitService.clusterRGroups' minute; both plus JOB_BASE_MS); `mcs` for one FMCS: RDKit's default
 * FMCS Timeout, an hour.
 */
const ANALYSIS_MS = {mmpa: 50, rGroups: 100, link: 20, molzip: 20, clusterMCS: CLUSTER_MCS_TIMEOUT_S * 1000,
  clusterRGroups: 60000, mcs: 3600000};
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
/**
 * The fields of a CXSMILES extension the Crux module reads as RDKit does (checked against RDKit 2026.03 and Chem's
 * RDKit.js on Chem's datasets, 2026-10-20): coordinates, stereo groups (&n, on, a, r), wedges (w, wU, wD), radicals
 * (^n) and atom labels (checked apart: {@link CX_QUERY_LABEL}). Its reader skips the query features (positional
 * variation m:, ring bond count rb:, unsaturation u:, substitution s:, c:/t: bond lists, link nodes LN:, ...), so a
 * query holding any other field runs on RDKit.
 */
const CX_FIELDS_CRUX_READS = [/^\([^)]*\)/, /^\$(?!_AV:)[^$]*\$/, /^(?:&|o)\d+:\d+(?:,\d+)*/, /^a:\d+(?:,\d+)*/,
  /^r(?=,|$)/, /^w[UD]?:\d+\.\d+(?:,\d+\.\d+)*/, /^\^\d:\d+(?:,\d+)*/];
/** CXSMILES atom labels of query atoms (RDKit's Q_e, QH_p, AH_p, X_p, XH_p, M_p, MH_p, star_e), which the Crux module
 * reads as plain dummy atoms (matching nothing). */
const CX_QUERY_LABEL = /_[pe]$/;

/** Whether the Crux module reads the CXSMILES extension of `molString`, if it has one, as RDKit does. */
function cruxReadsCxExtension(molString: string): boolean {
  const extension = hasNewLines(molString) ? null : /\s\|([^|]*)\|\s*$/.exec(molString);
  let rest = extension?.[1] ?? '';
  while (rest !== '') {
    const field = CX_FIELDS_CRUX_READS.map((re) => re.exec(rest)?.[0]).find((f) => f !== undefined);
    if (field === undefined ||
      (field.startsWith('$') && field.slice(1, -1).split(';').some((label) => CX_QUERY_LABEL.test(label))))
      return false;
    rest = rest.slice(field.length);
    if (rest.startsWith(','))
      rest = rest.slice(1);
    else if (rest !== '')
      return false;
  }
  return true;
}

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

  search = (key: string, molString: string, molBlockFailover: string): Promise<Uint32Array> =>
    this.call('search', [key, molString, molBlockFailover]);

  drop = (keys: string[]): Promise<void> => this.call('drop', [keys]);

  calcDescriptors = (molecules: (string | null)[], names: string[]):
    Promise<{[name: string]: PackedDescriptorValues}> => this.call('calcDescriptors', [molecules, names]);

  calcDescriptorsTimed = (molecules: (string | null)[], names: string[]):
    Promise<{values: {[name: string]: PackedDescriptorValues}, wasmMs: number, parseMs: number}> =>
    this.call('calcDescriptors::timed', [molecules, names]);

  scaffoldTree = (smiles: string[], ringCutoff: number, discharge: boolean): Promise<string> =>
    this.call('scaffoldTree', [smiles, ringCutoff, discharge]);

  getFingerprints = (type: Fingerprint, molecules: string[], smiles: boolean, options?: FingerprintOptions | null):
    Promise<IFpResult & {unread: number[]}> => this.call('getFingerprints', [type, molecules, smiles, options]);

  getStructuralAlerts = (alerts: {[rule in RuleId]?: string[]}, molecules: string[]):
    Promise<{alerts: {[rule in RuleId]?: boolean[]}, unread: number[]}> =>
    this.call('getStructuralAlerts', [alerts, molecules]);

  chemProperties = (molecules: (string | null)[], names: string[]):
    Promise<{[name: string]: (number | string | null)[]}> => this.call('chemProperties', [molecules, names]);

  pairUncommonAtomsAndBonds = (molecule1: string, molecule2: string, details: string): Promise<ISubstruct[] | null> =>
    this.call('pairUncommonAtomsAndBonds', [molecule1, molecule2, details]);

  /** An analysis job of crux-analysis-jobs.ts by name. */
  analysis = <R>(op: string, args: unknown[]): Promise<R> => this.call(op, args);

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
  /** The main thread's Crux module: whether Crux reads a query ({@link readsQuery}). */
  private crux: RDModule | null = null;
  private indexes = new Map<number, CruxColumnIndex>();
  private lastId = 0;
  private indexLock: Promise<unknown> = Promise.resolve();

  static async create(): Promise<CruxService> {
    const service = new CruxService();
    // Chem's one Crux file, compiled once (engine/crux-module.ts) and shared, so the workers only instantiate it
    service.crux = await loadCrux('search') as unknown as RDModule;
    service.module = await compileCrux();
    try {
      await service.startWorkers();
    } catch (e) {
      // The workers that did start are let go: the next search creates the service again, and each failed attempt
      // used to leave a pool running (hundreds of workers once the page had no room for another WebAssembly memory)
      for (const worker of service.workers)
        worker.terminate();
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

  /**
   * Whether Crux reads the query the way Chem's RDKit search builds it (getQueryMolSafe gives a molecule on the Crux
   * module): the workers then search with that molecule; any query they read they search as RDKit matches it. A
   * CXSMILES whose extension holds query features the Crux module does not read yet ({@link CX_FIELDS_CRUX_READS}) is
   * not read.
   */
  readsQuery(molString: string, molBlockFailover: string): boolean {
    if (!cruxReadsCxExtension(molString) || !cruxReadsCxExtension(molBlockFailover))
      return false;
    try {
      const query = getQueryMolSafe(molString, molBlockFailover, this.crux!);
      query?.delete();
      return !!query;
    } catch (e) {
      noteCruxError(e);
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

  /**
   * Searches every segment for the query (see {@link readsQuery}), calling `onSegment` with its segment-local hits as
   * soon as the segment is done.
   */
  async search(index: CruxColumnIndex, molString: string, molBlockFailover: string,
    onSegment: (segment: CruxSegment, hits: Uint32Array) => void, isCancelled: () => boolean): Promise<void> {
    await this.use(index, (segment, worker) => this.ensureBuilt(segment)
      .then(() => isCancelled() ? null : worker.search(segment.key, molString, molBlockFailover))
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
  async calcDescriptors(molecules: (string | null)[], names: string[], timing?: DescriptorTiming):
    Promise<{[name: string]: PackedDescriptorValues}> {
    const started = performance.now();
    const parts = await this.inChunks(molecules, JOB_MS_PER_MOLECULE.descriptors, !timing ?
      (worker, chunk) => worker.calcDescriptors(chunk, names) :
      async (worker, chunk) => {
        const sent = performance.now();
        const {values, wasmMs, parseMs} = await worker.calcDescriptorsTimed(chunk, names);
        timing.chunks.push({worker: this.workers.indexOf(worker), size: chunk.length, ms: performance.now() - sent,
          wasmMs, parseMs});
        return values;
      });
    const merging = performance.now();
    if (timing) {
      timing.workers = this.workers.length;
      timing.phases['chunks'] = merging - started;
    }
    const result: {[name: string]: PackedDescriptorValues} = {};
    for (const name of Object.keys(parts[0])) {
      // Python's type of the first value ('double' for a column without one): the first chunk with a value says it
      const typed = parts.find((part) => part[name].value.some((_, i) => !isNoneValue(part[name], i))) ?? parts[0];
      const type = typed[name].type;
      const value = type === 'int' ? new Int32Array(molecules.length) : new Float64Array(molecules.length);
      let offset = 0;
      for (const part of parts) {
        const chunk = part[name];
        if (chunk.type === type)
          value.set(chunk.value, offset);
        else { // a chunk of another type (only one without a value): its missing values as this type's
          for (let i = 0; i < chunk.value.length; i++)
            value[offset + i] = isNoneValue(chunk, i) ? (type === 'int' ? INT_NONE : NaN) : chunk.value[i];
        }
        offset += chunk.value.length;
      }
      result[name] = {type, value};
    }
    if (timing)
      timing.phases['merge'] = performance.now() - merging;
    return result;
  }

  /**
   * Fingerprints as the RDKit pool's `getFingerprints` gives them (the same worker code, run on Crux): one per molecule
   * (null where none), the canonical SMILES when asked, and the rows whose molecule Crux cannot read (`unread`).
   * `options`: radius, size, chirality, features (utils/fingerprint-options.ts).
   */
  async getFingerprints(type: Fingerprint, molecules: string[], smiles: boolean, options?: FingerprintOptions | null):
    Promise<IFpResult & {unread: number[]}> {
    const starts: number[] = [];
    const parts = await this.inChunks(molecules, JOB_MS_PER_MOLECULE.fingerprints,
      (worker, chunk) => worker.getFingerprints(type, chunk, smiles, options), CHUNK_SIZE, JOB_BASE_MS, starts);
    return {fps: parts.flatMap((part) => part.fps), smiles: smiles ? parts.flatMap((part) => part.smiles!) : null,
      unread: parts.flatMap((part, k) => part.unread.map((row) => starts[k] + row))};
  }

  /**
   * Which molecules match each group's SMARTS, as the RDKit pool's `getStructuralAlerts` gives it (the same worker
   * code, run on Crux): `[group, a flag per molecule]` per group, and the rows whose molecule Crux cannot read.
   */
  async getStructuralAlerts(alerts: {[rule in RuleId]?: string[]}, molecules: string[]):
    Promise<{alerts: [RuleId, boolean[]][], unread: number[]}> {
    const starts: number[] = [];
    const parts = await this.inChunks(molecules, JOB_MS_PER_MOLECULE.alerts,
      (worker, chunk) => worker.getStructuralAlerts(alerts, chunk), CHUNK_SIZE, JOB_BASE_MS, starts);
    const result: {[rule in RuleId]?: boolean[]} = {};
    for (const part of parts) {
      for (const rule of Object.keys(part.alerts) as RuleId[])
        result[rule] = (result[rule] ?? []).concat(part.alerts[rule]!);
    }
    return {alerts: Object.entries(result) as [RuleId, boolean[]][],
      unread: parts.flatMap((part, k) => part.unread.map((row) => starts[k] + row))};
  }

  /**
   * Chem's property table (widgets/properties.ts: MW, HBA, ..., `Molecular formula`) by crux's `chemProperties`
   * (RDKit's definitions), the molecules in chunks across the workers: a value or null per molecule for each name.
   */
  async chemProperties(molecules: (string | null)[], names: string[]):
    Promise<{[name: string]: (number | string | null)[]}> {
    const parts = await this.inChunks(molecules, JOB_MS_PER_MOLECULE.properties,
      (worker, chunk) => worker.chemProperties(chunk, names));
    return Object.fromEntries(names.map((name) => [name, parts.flatMap((part) => part[name])]));
  }

  /**
   * `job` over the molecules in chunks of at most `chunkSize` ({@link CHUNK_SIZE}), the results in chunk order, each
   * chunk's first row in `starts`. Fewer than `chunkSize` molecules per worker are split evenly over the workers. A
   * worker runs its calls one after another, so each worker is handed its next chunk when the last one is done, and a
   * chunk's time limit runs from its own start; a failed chunk ends the call (the chunks not started are not run).
   */
  private async inChunks<M, R>(molecules: M[], msPerMolecule: number,
    job: (worker: CruxWorkerClient, chunk: M[]) => Promise<R>, chunkSize = CHUNK_SIZE, baseMs = JOB_BASE_MS,
    starts: number[] = []): Promise<R[]> {
    const size = Math.max(1, Math.min(chunkSize, Math.ceil(molecules.length / Math.max(1, this.workers.length))));
    const chunks: M[][] = [];
    for (let start = 0; start < molecules.length; start += size) {
      starts.push(start);
      chunks.push(molecules.slice(start, start + size));
    }
    if (chunks.length === 0) {
      starts.push(0);
      chunks.push([]);
    }
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
   * Each molecule's atoms and bonds outside the MCS of the two (what an activity cliff's tooltip highlights) as RDKit's
   * Rascal MCES finds it (crux's `get_mces`, on one worker): the largest set of bonds the two share, possibly in
   * several fragments, where FMCS finds one connected substructure. Null for a pair whose similarity estimate is under
   * the threshold.
   */
  pairUncommonAtomsAndBonds(molecule1: string, molecule2: string): Promise<ISubstruct[] | null> {
    const options = {similarityThreshold: PAIR_MCS_SIMILARITY_THRESHOLD,
      maxBondMatchPairs: PAIR_MCS_MAX_BOND_MATCH_PAIRS, timeout: PAIR_MCS_TIMEOUT_S};
    return this.job(this.workers[0].pairUncommonAtomsAndBonds(molecule1, molecule2, JSON.stringify(options)),
      JOB_BASE_MS);
  }

  // Analyses (engine/analysis-jobs.ts sends them, crux-analysis-jobs.ts runs them): the RDKit worker's own code on
  // Crux, so the results have the RDKit pool's shapes; molzip, which Chem's RDKit.js lacks.

  /**
   * MMPA fragments as the RDKit pool's `mmpGetFragments` gives them, the molecules in chunks across the workers, and
   * the rows whose molecule Crux did not read (`unread`).
   */
  async mmpGetFragments(molecules: string[]): Promise<IMmpFragmentsResult & {unread: number[]}> {
    const starts: number[] = [];
    const parts = await this.inChunks(molecules, ANALYSIS_MS.mmpa, (worker, chunk) =>
      worker.analysis<IMmpFragmentsResult & {unread: number[]}>('mmpGetFragments', [chunk]), CHUNK_SIZE, JOB_BASE_MS,
    starts);
    // concat, not flatMap: a molecule RDKit's code does not read leaves its row empty, and the rows stay in place
    return {frags: ([] as [string, string][][]).concat(...parts.map((part) => part.frags)),
      smiles: ([] as string[]).concat(...parts.map((part) => part.smiles)),
      unread: parts.flatMap((part, k) => part.unread.map((row) => starts[k] + row))};
  }

  /** The RDKit pool's `getRGroups` (one decomposition, on the least busy worker); null where Crux cannot read a
   * molecule. */
  rGroupAnalysis(molecules: string[], core: string, coreIsQMol: boolean, options?: string):
    Promise<IRGroupAnalysisResult | null> {
    return this.job(this.leastBusy().analysis('rGroupAnalysis', [molecules, core, coreIsQMol, options]),
      JOB_BASE_MS + molecules.length * ANALYSIS_MS.rGroups);
  }

  /** The RDKit pool's `getMCS` (FMCS's SMARTS, on the least busy worker); null where Crux cannot read a molecule. */
  mostCommonStructure(molecules: string[], exactAtoms: boolean, exactBonds: boolean): Promise<string | null> {
    return this.job(this.leastBusy().analysis('mostCommonStructure', [molecules, exactAtoms, exactBonds]),
      JOB_BASE_MS + ANALYSIS_MS.mcs);
  }

  /**
   * Each cluster's FMCS SMARTS, as the RDKit pool's `clusterMCS(clusters, ..., rawSmarts)` gives it ('' for one whose
   * FMCS ran past the RDKit pool's 45 s), each cluster a job ({@link eachCluster}): null for a cluster Crux cannot do
   * (the RDKit pool does it); `timedOut`, the clusters past their time limit (no MCS).
   */
  clusterMCS(clusters: string[][], exactAtoms: boolean, exactBonds: boolean):
    Promise<{results: (string | null)[], timedOut: number[]}> {
    return this.eachCluster<string>(clusters, JOB_BASE_MS + ANALYSIS_MS.clusterMCS, 'clusterMCS',
      [exactAtoms, exactBonds]);
  }

  /** Each cluster's R-group decomposition, as the RDKit pool's `clusterRGroups` gives it, each cluster a job
   * ({@link eachCluster}): null for a cluster Crux cannot do; `timedOut`, the clusters past their time limit (none). */
  clusterRGroups(clusters: RGroupCluster[]): Promise<{results: (IRGroupAnalysisResult | null)[], timedOut: number[]}> {
    return this.eachCluster<IRGroupAnalysisResult>(clusters, JOB_BASE_MS + ANALYSIS_MS.clusterRGroups,
      'clusterRGroups', []);
  }

  /**
   * The analysis job `op` per cluster, as the RDKit pool's processClusterQueue runs it: each worker takes the next
   * cluster when its last is done; a cluster past `limitMs` is listed in `timedOut` and its worker replaced
   * ({@link restartWorker}), and the other clusters go on. A cluster whose call fails (a trap, a lost worker) is null.
   */
  private async eachCluster<R>(clusters: unknown[], limitMs: number, op: string, args: unknown[]):
    Promise<{results: (R | null)[], timedOut: number[]}> {
    const results = new Array<R | null>(clusters.length).fill(null);
    const timedOut: number[] = [];
    let next = 0;
    await Promise.all(this.workers.map(async (_, w) => {
      while (next < clusters.length) {
        const i = next++;
        let timer: ReturnType<typeof setTimeout> | undefined;
        const late = new Promise<'late'>((resolve) => timer = setTimeout(() => resolve('late'), limitMs));
        try {
          const done = await Promise.race([this.workers[w].analysis<(R | null)[]>(op, [[clusters[i]], ...args]),
            late]);
          if (done === 'late') {
            timedOut.push(i);
            this.restartWorker(w);
          } else
            results[i] = done[0];
        } catch (e: any) {
          if (/RuntimeError|unreachable|out of bounds|memory/i.test(String(e?.message ?? e)))
            this.restartWorker(w); // its instance may be broken
        } finally {
          clearTimeout(timer);
        }
      }
    }));
    return {results, timedOut};
  }

  /** The worker with the fewest calls in flight (on a tie a later one: worker 0 takes the one-worker jobs). */
  private leastBusy(): CruxWorkerClient {
    return this.workers.reduceRight((best, worker) => worker.pendingCalls < best.pendingCalls ? worker : best);
  }

  /**
   * Replaces worker `w`, stuck in a call past its limit, as the RDKit pool replaces one (RdKitService.restartWorker):
   * the column indexes with a segment on it are dropped (a search running on one fails and falls back to RDKit; the
   * next search indexes the column again); the other workers and the other indexes stay.
   */
  private restartWorker(w: number): void {
    this.workers[w].terminate();
    this.workers[w] = new CruxWorkerClient();
    this.workers[w].moduleInit(this.module!).catch(() => {}); // a failed init is logged by the worker client
    for (const index of [...this.indexes.values()]) {
      if (index.segments.some((segment) => segment.worker === w)) {
        index.lost = true;
        this.dropIndex(index);
      }
    }
  }

  /** The RDKit pool's `linkRGroupFragments`, joined by molzip, the rows in chunks across the workers. */
  async linkRGroupFragments(cores: string[], fragmentColumns: string[][], attachIdx: (number | number[])[]):
    Promise<string[]> {
    const rows = cores.map((core, i) => [core, ...fragmentColumns.map((column) => column[i])]);
    return (await this.inChunks(rows, ANALYSIS_MS.link, (worker, chunk) =>
      worker.analysis<string[]>('linkRGroupFragments', [chunk, attachIdx]))).flat();
  }

  /** RDKit's molzip of each molecule's pieces (`details`: MolzipParams as JSON), the molecules in chunks across the
   * workers: the joined molecules' SMILES ('' where one is not zipped) and the first refusal's message. */
  async molzip(molecules: string[], details: string): Promise<{smiles: string[], error: string | null}> {
    const parts = await this.inChunks(molecules, ANALYSIS_MS.molzip, (worker, chunk) =>
      worker.analysis<{smiles: string[], error: string | null}>('molzip', [chunk, details]));
    return {smiles: parts.flatMap((part) => part.smiles), error: parts.find((part) => part.error)?.error ?? null};
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
    const smiles = values; // a molblock is read by the worker that indexes it (crux.worker.ts build)
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
    const smiles = changed.flatMap((s) => values.slice(s.start, s.end));
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
