import type {RDModule, RDMol} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {ChemWasm} from '../engine/crux-module';
import type {AssembleOptions} from '@datagrok-libraries/chem-meta/src/macromolecule-assembler';
import type {RdKitServiceWorker} from '../rdkit-service/rdkit-service-worker';

/**
 * What Chem calls in the one Crux file (vendor/chem-wasm) beyond RDKit.js's API: the search Collection (the API of the
 * search-only build it replaces), descriptors and the scaffold tree; and the shape of the Crux worker's jobs. Rendering
 * uses the RDKit.js-shaped calls (`RDModule`, chem-meta).
 */

/** What every job of the Crux worker (crux.worker.ts `JOBS`: crux/jobs/, crux-analysis-jobs.ts) gets first. */
export interface WorkerJobContext {
  /** The worker's Crux module: chem-wasm's own exports beside RDKit.js's calls. */
  crux: ChemWasm;
  /** The same module as RDKit.js's (`get_mol`, ...), for the helpers that take either engine's module (getMolSafe). */
  module: RDModule;
  /** The RDKit worker's code (rdkit-service/) over the Crux module, made on first use, without a molecule cache. */
  rdkitWorker: () => RdKitServiceWorker;
  /** Buffers of the result that go to the main thread without a copy. */
  transfer: Transferable[];
}

/** A job of the Crux worker, run by name (CruxService: `jobInChunks`, `idleWorkerJob`, a typed method). */
export type WorkerJob = (ctx: WorkerJobContext, ...args: any[]) => unknown;

/** Whether `e` leaves a Crux instance broken: a wasm trap (`WebAssembly.RuntimeError`), or a stack overflow out of a
 * wasm call (`RangeError: Maximum call stack size exceeded` in Chrome and Safari, `InternalError: too much recursion`
 * in Firefox; e.g. a recursive writer on a very large molecule), whose unwind does not restore the wasm stack pointer:
 * after a few the instance reads and writes out of bounds. */
export function isWasmTrap(e: unknown): e is Error {
  return e instanceof WebAssembly.RuntimeError || (e instanceof Error &&
    (e instanceof RangeError || e.name === 'InternalError') && /call stack|too much recursion/i.test(e.message));
}

/** Rethrows a wasm trap ({@link isWasmTrap}), which a job must not turn into "no result": the call fails, so the
 * service replaces the worker (its instance may be broken) and RDKit does the work. Other errors are the job's to
 * handle. */
export function throwIfTrap(e: unknown): void {
  if (isWasmTrap(e))
    throw e;
}

export interface CruxCollection {
  /**
   * The input positions of the molecules RDKit's `get_substruct_match(query)` matches (0 = no limit): `query` is a
   * molecule of the same Crux module (Chem's getQueryMolSafe builds it); stereo is ignored, as there.
   */
  substructureSearchMol(query: RDMol, limit: number): Uint32Array;
  free(): void;
}

export interface CruxCollectionBuilder {
  /** Read what RDKit's default read rejects the way RDKit reads it without Kekulize (Chem's getMolSafe retry). */
  setLenient(lenient: boolean): void;
  /** The number of molecules added: SMILES, or MOL blocks read as Chem's getMolSafe reads them (H atoms kept). */
  addMany(molecules: string[]): number;
  finish(): CruxCollection;
}

/** One FilterCatalog alert crux's `chemAlerts` found in a molecule: its name, and per match the molecule's atoms. */
export interface ToxicityAlert {
  name: string;
  atoms: number[][];
}

/** A descriptor's values as the chem-chem container returns them: Python's type of its first value, nulls kept. */
export interface DescriptorValues {
  type: 'int' | 'double';
  value: (number | null)[];
}

/**
 * A descriptor's values as the Crux workers send them and Chem merges them: a typed array, a missing value as
 * {@link INT_NONE} (ints) or NaN (doubles; the container never gives NaN: it gives NaN and infinities as missing).
 */
export interface PackedDescriptorValues {
  type: 'int' | 'double';
  value: Int32Array | Float64Array;
}

/** A missing int descriptor value in {@link PackedDescriptorValues} (Datagrok's int null). */
export const INT_NONE = -2147483648;

/** The seconds a cluster's FMCS may take on Crux (CruxService.clusterMCS, crux-analysis-jobs.ts): '' past it, as
 * RdKitService.clusterMCS gives a cluster whose worker it restarts after 45 s. */
export const CLUSTER_MCS_TIMEOUT_S = 45;

/** Whether `value[i]` of a packed descriptor is missing. */
export function isNoneValue(values: PackedDescriptorValues, i: number): boolean {
  const v = values.value[i];
  return values.type === 'int' ? v === INT_NONE : Number.isNaN(v);
}

/**
 * Diagnosis of a Crux descriptors call (setDescriptorTimingForTests): its phases on the main thread in ms, the pool
 * size, and each chunk's round trip (`ms`), the worker's Crux call (`wasmMs`) and its JSON.parse (`parseMs`).
 */
export interface DescriptorTiming {
  rows: number;
  names: number;
  workers: number;
  phases: {[phase: string]: number};
  chunks: {worker: number, size: number, ms: number, wasmMs: number, parseMs: number}[];
}

// The biopolymer calls (Bio-on-Crux PLAN §4.3, crux-core's chem-wasm `chem_biopolymer.rs`; types from the generated
// vendor/chem-wasm/chem_wasm.d.ts): a monomer library compiled once and assembled from per batch (the workers,
// jobs/biopolymer.ts) or per row (the main thread, engine/biopolymer.ts); the stereo pass and the V2000 writer.

/** A Pistoia monomer library read by Crux (`new ChemMonomerLibrary(json)`: throws for text that is no JSON array;
 * each template is compiled at its first use, every one by `errors()`). */
export type CruxMonomerLibrary = InstanceType<ChemWasm['ChemMonomerLibrary']>;

/** A batch's result (`ChemMonomerLibrary.assemble`): per row a V3000 block (null where the row has an error, '' for a
 * row of no monomers), an AsmLayout (passed through: a value Chem does not list is kept) and an error. */
export type CruxAssembledBatch = ReturnType<CruxMonomerLibrary['assemble']>;

/** One row's result (`ChemMonomerLibrary.assembleOne`), with `atomPosition` / `bondPosition` under the option `map`. */
export type CruxAssembledRow = ReturnType<CruxMonomerLibrary['assembleOne']>;

/** The biopolymer exports of Chem's Crux file. */
export type CruxBiopolymerApi = Pick<ChemWasm, 'ChemMonomerLibrary' | 'chemMarkStereoAbsolute' | 'chemToV2000'>;

/** The biopolymer calls of `crux` (a worker's module or the main thread's); throws "not supported by Crux" for a Crux
 * file older than them. */
export function cruxBiopolymer(crux: unknown): CruxBiopolymerApi {
  const api = crux as Partial<CruxBiopolymerApi>;
  if (typeof api.ChemMonomerLibrary !== 'function' || typeof api.chemMarkStereoAbsolute !== 'function' ||
    typeof api.chemToV2000 !== 'function')
    throw new Error('not supported by Crux: biopolymer assembly (a Crux file older than its biopolymer calls)');
  return api as CruxBiopolymerApi;
}

/**
 * AssembleOptions (chem-meta) as the binding reads them, the same for a batch on the workers and a row on the main
 * thread: its keys only (another key is "not supported by Crux"), `map` only where `map` is set (a row's call;
 * `assemble` refuses it). The depictor by default, as Bio's beautify was: Crux's own default sends big rows to
 * CoordGen, about 8x slower (an antibody VH chain 1.6 s instead of 0.2 s), so only `nonlinearEngine: 'auto'` asks for
 * that switch (`autoCoordgen`).
 */
export function assembleOptionsJson(options: AssembleOptions, map: boolean): string {
  const engine = options.nonlinearEngine ?? 'depictor';
  return JSON.stringify({layout: options.layout, autoCoordgen: engine === 'auto',
    ...(engine === 'auto' ? {} : {nonlinearEngine: engine}), ...(map && options.map ? {map: true} : {})});
}

/** How a Crux worker says it does not hold the monomer library a batch names (jobs/biopolymer.ts): the service then
 * sends the library and asks once more. Not a trap's message (CruxService's `isTrap`). */
export const BIOPOLYMER_LIBRARY_MISSING = 'biopolymer: the worker holds no monomer library';

/** The monomer libraries a Crux worker keeps compiled (jobs/biopolymer.ts), the least recently used freed first. */
export const BIOPOLYMER_LIBRARIES = 4;

export interface CruxWorkerApi {
  CollectionBuilder: new (buildIndex: boolean) => CruxCollectionBuilder;
  /** The container's `GET /chem/descriptors/tree` body (JSON text): groups by RDKit module, the 3D group too. */
  descriptorTree(): string;
  /**
   * The container's `POST /chem/descriptors` body (JSON text, `{name: DescriptorValues}`): molecules as chem.py reads
   * them (null or unreadable: nulls), descriptor or group names, in request order.
   */
  calcDescriptors(molecules: (string | null)[], names: string[]): string;
  /** scripts/scaffold_tree.py's JSON text (ScaffoldGraph's Scaffold Tree); throws where the script fails. */
  scaffoldTree(smiles: string[], ringCutoff: number, dischargeAndDeradicalize: boolean): string;
  /**
   * RDKit's Rascal MCES of two molecules of the module (`rascalMCES`, Python's `FindMCES`): a JSON array, `[]` below
   * `details.similarityThreshold` (0.5 by default), else one result `{smartsString, atomMatches, bondMatches,
   * similarity, timedOut, numFragments, ...}`. `details`: RDKit's `RascalOptions` fields by name.
   */
  get_mces(mol1: RDMol, mol2: RDMol, details?: string): string;
}
