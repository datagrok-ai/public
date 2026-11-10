/**
 * The Crux worker: Chem's one Crux file (vendor/chem-wasm), instantiated from the Module the main thread compiled
 * (engine/crux-module.ts). It keeps the substructure search's column segments (`build`, `search`, `drop`) and runs the
 * jobs of the job tables by name (`JOBS`: crux/jobs/*.ts, crux-analysis-jobs.ts), each given the worker's
 * `WorkerJobContext` (crux-api.ts). A new job is an entry in a table; CruxService sends it (`jobInChunks` for a column,
 * `idleWorkerJob` for one long call, or a method of its own).
 */
import initChemWasm, * as chemWasm from '../../vendor/chem-wasm/chem_wasm.js';
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {CruxCollection, CruxWorkerApi, WorkerJob, WorkerJobContext} from './crux-api';
import {RdKitServiceWorker} from '../rdkit-service/rdkit-service-worker';
import {MAX_MOL_CACHE_SIZE} from '../rdkit-service/rdkit-service-worker-base';
import {getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {isWasmTrap} from './crux-api';
import {trapNoted, watchTraps} from './trap-watch';
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

/** The module, every call watched for a trap its caller swallows (trap-watch.ts): later jobs fail as a trap. */
const glue = watchTraps(chemWasm);
const crux = glue as unknown as CruxWorkerApi;
const module = glue as unknown as RDModule;
const ctx: Worker = self as any;
const collections = new Map<string, CruxCollection>();
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
  const builder = new crux.CollectionBuilder(true);
  // what RDKit's default read rejects, the way RDKit reads it without Kekulize (Chem's getMolSafe retry)
  builder.setLenient(true);
  const failed: number[] = [];
  // one molecule per call: the builder only counts failures, the NOT CONTAINS complement needs their rows
  for (let i = 0; i < molecules.length; i++) {
    if (builder.addMany([molecules[i]]) === 0)
      failed.push(i);
  }
  collections.set(key, builder.finish());
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
    return query.get_num_atoms() === 0 ? new Uint32Array(0) : collections.get(key)!.substructureSearchMol(query, 0);
  } finally {
    query.delete();
  }
}

function drop(keys: string[]): void {
  for (const key of keys) {
    collections.get(key)?.free();
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
      const context: WorkerJobContext = {crux: glue, module, rdkitWorker: rdkitShapedWorker, transfer};
      result = await JOBS[op](context, ...args);
    } else if (op === 'fail::test') // tests: a job that fails with this message
      throw new Error(args[0]);
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
    const swallowed = trapNoted(); // a trap a molecule's `catch` let pass: the instance answers nothing more
    if (swallowed !== null)
      throw new WebAssembly.RuntimeError(`a Crux call trapped, ${op} fails (${swallowed})`);
    port.postMessage({op: op, retval: result}, transfer);
  } catch (err: any) {
    // a trap reported as one: the service restarts the workers for it (and not for a refusal)
    const trapped = isWasmTrap(err) || trapNoted() !== null;
    port.postMessage({error: trapped ? `RuntimeError: ${err?.message ?? err}` :
      err instanceof Error ? err.message : String(err)});
  }
});
