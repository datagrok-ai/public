/**
 * The Crux worker: Chem's one Crux file (vendor/chem-wasm), instantiated from the Module the main thread compiled
 * (engine/crux-module.ts). It keeps the substructure search's column segments (`build`, `search`, `drop`) and runs the
 * jobs of the job tables by name (`JOBS`: crux/jobs/*.ts, crux-analysis-jobs.ts), each given the worker's
 * `WorkerJobContext` (crux-api.ts). A new job is an entry in a table; CruxService sends it (`jobInChunks` for a column,
 * `idleWorkerJob` for one long call, or a method of its own).
 *
 * The module recovers from a trap inside the call (crux-recover.ts: a fresh instance, the objects rebuilt, the call
 * tried once more; a second trap refuses that call alone, a plain error the job reports as a refusal), so the worker
 * goes on. The search collections keep no journal (their molecules are many): a build the engine was renewed in is made
 * again, the molecule that trapped twice a row Crux could not read (the main thread's RDKit asks about it, as about any
 * such row); a collection made before a renewal is lost, and a search of it says so ({@link INDEX_LOST}): the service
 * drops that column's index and indexes it again at the next search.
 */
import initChemWasm, * as chemWasm from '../../vendor/chem-wasm/chem_wasm.js';
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {CruxCollection, CruxWorkerApi, WorkerJob, WorkerJobContext} from './crux-api';
import {RdKitServiceWorker} from '../rdkit-service/rdkit-service-worker';
import {MAX_MOL_CACHE_SIZE} from '../rdkit-service/rdkit-service-worker-base';
import {getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {INDEX_LOST, isWasmTrap} from './crux-api';
import {recoverable} from './crux-recover';
import {ANALYSIS_JOBS} from './crux-analysis-jobs';
import {BATCH_JOBS} from './jobs/batch';
import {BIOPOLYMER_JOBS} from './jobs/biopolymer';
import {CALCULATION_JOBS} from './jobs/calculations';
import {MOLECULE_JOBS} from './jobs/molecule-jobs';
import {SEARCH_TYPE_JOBS, dropSearchTypeKeys} from './jobs/search-types';
import {SCRIPT_JOBS} from './jobs/python-free';
import {SYNTHON_JOBS} from './jobs/synthon-search';
import {INCHI_JOBS} from './jobs/inchi';
import {REACTION_JOBS} from './jobs/reactions';

/** The module, recovering (see the header). */
const recovering = recoverable(chemWasm, {
  settings: (g: typeof chemWasm) => g.use_legacy_stereo_perception(false), // as module::init sets it
  log: 'Chem | Crux worker |',
  unjournaled: ['Collection', 'CollectionBuilder'],
});
const glue = recovering.crux;
const crux = glue as unknown as CruxWorkerApi;
const module = glue as unknown as RDModule;
const ctx: Worker = self as any;
/** Each segment's collection, with the module's renewal count when it was made (a later renewal lost it). */
const collections = new Map<string, {collection: CruxCollection, renewals: number}>();
/** An input Crux's SMILES reader refuses: what a build adds for a molecule that trapped twice, so that the row keeps
 * its place among the collection's inputs and is reported as one Crux could not read. */
const UNREADABLE = ')';
/** Tests: the row of the next builds before whose read the module's test trap is armed (each attempt), or -1. */
let buildTrapRow = -1;
/** The RDKit worker's code over the Crux module, which answers RDKit.js's calls (conversions, the analyses). Made on
 * first use, without its molecule cache (its budget spent from the start): each job parses its molecules and frees
 * them; and without its pauses every 50 molecules ({@link NO_PAUSES}). */
let rdkitShaped: RdKitServiceWorker | null = null;
/**
 * The RDKit worker's code pauses (a `setTimeout`) every `_terminationCheckDelay` molecules to read the RDKit pool's
 * terminate flag, which nothing sets here (CruxService replaces a worker instead), and from the sixth pause on the
 * browser holds each one 4 ms (nested timers): a conversion of 1,000 molecules slept ~80 ms. One pause per call.
 */
const NO_PAUSES = Number.MAX_SAFE_INTEGER;
/** The pause interval of {@link rdkitShaped} (tests compare it with the RDKit pool's 50: `pauses::test`). */
let pauseEvery = NO_PAUSES;
function rdkitShapedWorker(): RdKitServiceWorker {
  if (!rdkitShaped) {
    rdkitShaped = new RdKitServiceWorker(module, '');
    rdkitShaped._cacheCounter = MAX_MOL_CACHE_SIZE;
    rdkitShaped._terminationCheckDelay = pauseEvery;
  }
  return rdkitShaped;
}

/** Every table's jobs by name; a name in two tables is an error. */
function jobTable(...tables: {[op: string]: WorkerJob}[]): {[op: string]: WorkerJob} {
  const jobs: {[op: string]: WorkerJob} = {};
  for (const table of tables) {
    for (const [op, job] of Object.entries(table)) {
      if (Object.hasOwn(jobs, op))
        throw new Error(`Chem | two Crux worker jobs are named ${op}`);
      jobs[op] = job;
    }
  }
  return jobs;
}

const JOBS = jobTable(CALCULATION_JOBS, BATCH_JOBS, MOLECULE_JOBS, ANALYSIS_JOBS, SEARCH_TYPE_JOBS, SCRIPT_JOBS,
  SYNTHON_JOBS, BIOPOLYMER_JOBS, INCHI_JOBS, REACTION_JOBS);

/**
 * Builds the collection of one index segment from its SMILES and MOL blocks (a MOL block read as Chem's getMolSafe
 * reads it, its H atoms kept); returns the (segment-local) rows Crux could not read: the search asks RDKit on the main
 * thread about them (a MOL block with query features, say).
 */
function build(key: string, molecules: string[]): Uint32Array {
  const fresh = () => {
    const b = new crux.CollectionBuilder(true);
    // what RDKit's default read rejects, the way RDKit reads it without Kekulize (Chem's getMolSafe retry)
    b.setLenient(true);
    return b;
  };
  let builder = fresh();
  let failed: number[] = [];
  /** Rows the engine was renewed in once, and twice (these are added as {@link UNREADABLE}). */
  const retried = new Set<number>();
  const refused = new Set<number>();
  // one molecule per call: the builder only counts failures, the NOT CONTAINS complement needs their rows
  for (let i = 0; i < molecules.length; i++) {
    const renewals = recovering.renewals();
    if (i === buildTrapRow && !refused.has(i))
      recovering.armTestTrap(1); // tests: a trap in this row's read, on every attempt (a molecule that always traps)
    try {
      if (builder.addMany([refused.has(i) ? UNREADABLE : molecules[i]]) === 0)
        failed.push(i);
    } catch (e) {
      if (recovering.renewals() === renewals)
        throw e;
      // the engine was renewed in this call (a trap): the builder is lost. It is made again from row 0, row i tried
      // once more there; renewed in again, row i is added as one Crux cannot read.
      (retried.has(i) ? refused : retried).add(i);
      builder = fresh();
      failed = [];
      i = -1;
    }
  }
  collections.set(key, {collection: builder.finish(), renewals: recovering.renewals()});
  return new Uint32Array(failed);
}

/**
 * The segment's molecules that match the query Chem's RDKit search builds (getQueryMolSafe), read by Crux: they are
 * the molecules RDKit's get_substruct_match matches (stereo ignored, as there). The main thread checked that Crux
 * reads the query (CruxService.readsQuery). A query of no atoms (a MOL block of none; a query that does not read, with
 * no failover) matches no molecule, as in RDKit.
 */
function search(key: string, molString: string, molBlockFailover: string): Uint32Array {
  const query = getQueryMolSafe(molString, molBlockFailover, module);
  if (!query)
    throw new Error(`Crux reads no query from ${molString}`);
  try {
    if (query.get_num_atoms() === 0)
      return new Uint32Array(0);
    const held = collections.get(key)!;
    if (held.renewals !== recovering.renewals())
      throw new Error(`${INDEX_LOST} ${key} (made again at the next search)`);
    return held.collection.substructureSearchMol(query, 0);
  } finally {
    query.delete();
  }
}

function drop(keys: string[]): void {
  for (const key of keys) {
    collections.get(key)?.collection.free(); // a lost one's free does nothing (the glue: old objects are inert)
    collections.delete(key);
  }
  dropSearchTypeKeys(keys);
}

ctx.addEventListener('message', async (e: MessageEvent) => {
  const {op, args} = e.data;
  const port = e.ports[0];
  try {
    let result;
    const transfer: Transferable[] = [];
    if (op === 'module::init') {
      await initChemWasm({module_or_path: args[0]});
      chemWasm.use_legacy_stereo_perception(false); // as the RDKit workers (rdkit.worker.ts)
    } else if (op === 'build')
      result = build(args[0], args[1]);
    else if (op === 'search')
      result = search(args[0], args[1], args[2]);
    else if (op === 'drop')
      drop(args[0]);
    else if (Object.hasOwn(JOBS, op)) {
      const context: WorkerJobContext = {crux: glue, module, rdkitWorker: rdkitShapedWorker, transfer,
        renewals: () => recovering.renewals()};
      result = await JOBS[op](context, ...args);
    } else if (op === 'fail::test') // tests: a job that fails with this message
      throw new Error(args[0]);
    else if (op === 'trap::test') // tests: crux-core's real test trap (calls, instances), crux-recover.ts
      recovering.armTestTrap(args[0], args[1]);
    else if (op === 'buildTrap::test') // tests: the row of the next builds whose read always traps (-1: none)
      buildTrapRow = args[0];
    else if (op === 'pauses::test') { // tests: the RDKit-shaped worker's pause interval (null: none, as by default)
      pauseEvery = args[0] ?? NO_PAUSES;
      if (rdkitShaped)
        rdkitShaped._terminationCheckDelay = pauseEvery;
    } else if (op === 'busy::test') { // tests: a job that keeps the worker busy, as a long Crux call does
      const end = performance.now() + args[0];
      let spins = 0;
      while (performance.now() < end)
        spins++;
      result = spins;
    } else
      throw new Error(`Unknown operation: ${op}`);
    const dead = recovering.failed(); // no fresh instance could be made: the instance answers nothing more
    if (dead !== null)
      throw new WebAssembly.RuntimeError(`a Crux call trapped, ${op} fails (${(dead as Error)?.message ?? dead})`);
    port.postMessage({op: op, retval: result}, transfer);
  } catch (err: any) {
    // a trap reported as one: the service restarts the workers for it (and not for a refusal)
    const trapped = isWasmTrap(err) || recovering.failed() !== null;
    port.postMessage({error: trapped ? `RuntimeError: ${err?.message ?? err}` :
      err instanceof Error ? err.message : String(err)});
  }
});
