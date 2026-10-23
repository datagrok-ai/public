// The Crux worker: Chem's one Crux file (vendor/chem-wasm), instantiated from the Module the main thread compiled
// (engine/crux-module.ts). Substructure search over column segments; in Crux mode ("Chemical engine") descriptors, the
// scaffold tree, an activity cliff's MCS, and RDKit's own worker code run on Crux (fingerprints, structural alerts;
// MMPA, R-groups, FMCS and fragment linking: crux-analysis-jobs.ts).
import initChemWasm, * as chemWasm from '../../vendor/chem-wasm/chem_wasm.js';
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {ISubstruct} from '@datagrok-libraries/chem-meta/src/types';
import type {CruxCollection, CruxWorkerApi, DescriptorValues, PackedDescriptorValues} from './crux-api';
import {INT_NONE} from './crux-api';
import {RdKitServiceWorker} from '../rdkit-service/rdkit-service-worker';
import {notConverted} from '../rdkit-service/rdkit-service-worker-substructure';
import {MAX_MOL_CACHE_SIZE} from '../rdkit-service/rdkit-service-worker-base';
import {getMolSafe, getQueryMolSafe} from '../utils/mol-creation_rdkit';
import {getUncommonAtomsAndBondsOfMatches} from '../utils/chem-common';
import {assertFingerprintOptions} from '../utils/fingerprint-options';
import {ANALYSIS_JOBS} from './crux-analysis-jobs';

const crux = chemWasm as unknown as CruxWorkerApi;
const ctx: Worker = self as any;
const collections = new Map<string, CruxCollection>();
/** The RDKit worker's code over the Crux module, which answers RDKit.js's calls: the same fingerprints (Morgan r2/2048
 * after remove_hs, pattern, ...) as the RDKit pool's (the alerts: {@link structuralAlerts}, the same matching). Made on
 * first use, without its molecule cache (its budget spent from the start): each job parses its molecules and frees
 * them. */
let rdkitShaped: RdKitServiceWorker | null = null;
function rdkitShapedWorker(): RdKitServiceWorker {
  if (!rdkitShaped) {
    rdkitShaped = new RdKitServiceWorker(chemWasm as unknown as RDModule, '');
    rdkitShaped._cacheCounter = MAX_MOL_CACHE_SIZE;
  }
  return rdkitShaped;
}

/** The rows of `molecules` with a molecule Crux cannot read: the RDKit pool does them (engine/worker-jobs.ts). */
function unread(molecules: string[], done: (i: number) => boolean): number[] {
  const rows: number[] = [];
  for (let i = 0; i < molecules.length; i++) {
    if (!molecules[i] || done(i))
      continue;
    const mol = getMolSafe(molecules[i], {}, chemWasm as unknown as RDModule).mol;
    if (mol)
      mol.delete();
    else
      rows.push(i);
  }
  return rows;
}

/**
 * Which molecules match each rule's SMARTS, as the RDKit worker's `getStructuralAlerts` matches them (a rule's flag is
 * set by its first matching SMARTS), each molecule parsed once; and the rows with a molecule Crux cannot read.
 */
function structuralAlerts(alerts: {[rule: string]: string[]}, molecules: string[]):
  {alerts: {[rule: string]: boolean[]}, unread: number[]} {
  const module = chemWasm as unknown as RDModule;
  const queries = Object.entries(alerts).map(([rule, smarts]) => [rule, smarts.map((s) => {
    try {
      return module.get_qmol(s);
    } catch {
      return null;
    }
  })] as const);
  const found: {[rule: string]: boolean[]} = {};
  for (const [rule] of queries)
    found[rule] = new Array<boolean>(molecules.length).fill(false);
  const unread: number[] = [];
  try {
    for (let i = 0; i < molecules.length; i++) {
      const mol = getMolSafe(molecules[i], {}, module).mol;
      if (!mol) {
        if (molecules[i])
          unread.push(i);
        continue;
      }
      try {
        for (const [rule, smarts] of queries) {
          for (const query of smarts) {
            try {
              if (query && mol.get_substruct_match(query) !== '{}') {
                found[rule][i] = true;
                break;
              }
            } catch {}
          }
        }
      } finally {
        mol.delete();
      }
    }
  } finally {
    for (const [, smarts] of queries)
      smarts.forEach((query) => query?.delete());
  }
  return {alerts: found, unread};
}

/** The container-shaped descriptor values as typed arrays (the main thread gets their buffers, not copies). */
function packDescriptors(values: {[name: string]: DescriptorValues}, transfer: Transferable[]):
  {[name: string]: PackedDescriptorValues} {
  const packed: {[name: string]: PackedDescriptorValues} = {};
  for (const [name, {type, value}] of Object.entries(values)) {
    const array = type === 'int' ? new Int32Array(value.length) : new Float64Array(value.length);
    const none = type === 'int' ? INT_NONE : NaN;
    for (let i = 0; i < value.length; i++)
      array[i] = value[i] ?? none;
    packed[name] = {type, value: array};
    transfer.push(array.buffer);
  }
  return packed;
}

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
 * no failover) matches no molecule, as in RDKit (Crux's search matched every one before crux-core d71ad16e7).
 */
function search(key: string, molString: string, molBlockFailover: string): Uint32Array {
  const query = getQueryMolSafe(molString, molBlockFailover, chemWasm as unknown as RDModule);
  if (!query)
    throw new Error(`Crux reads no query from ${molString}`);
  try {
    return query.get_num_atoms() === 0 ? new Uint32Array(0) : collections.get(key)!.substructureSearchMol(query, 0);
  } finally {
    query.delete();
  }
}

/**
 * Each molecule's atoms and bonds outside the RDKit Rascal MCES of the two (`get_mces`; `details`: its options), each
 * molecule read as Chem's getMolSafe reads it; null when the pair is below the similarity threshold. Throws for a
 * molecule Crux cannot read.
 */
function pairUncommonAtomsAndBonds(molecule1: string, molecule2: string, details: string): ISubstruct[] | null {
  const read = [molecule1, molecule2].map((m) => getMolSafe(m, {}, chemWasm as unknown as RDModule));
  try {
    const unreadable = read.findIndex((m) => !m.mol || m.isQMol);
    if (unreadable >= 0)
      throw new Error(`Crux reads no molecule from ${[molecule1, molecule2][unreadable]}`);
    const [mces] = JSON.parse(crux.get_mces(read[0].mol!, read[1].mol!, details));
    return !mces ? null : read.map(({mol}, k) => getUncommonAtomsAndBondsOfMatches(mol!,
      mces.atomMatches.map((match: number[]) => match[k]), mces.bondMatches.map((match: number[]) => match[k])));
  } finally {
    read.forEach((m) => m.mol?.delete());
  }
}

function drop(keys: string[]): void {
  for (const key of keys) {
    collections.get(key)?.free();
    collections.delete(key);
  }
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
    else if (op === 'calcDescriptors')
      result = packDescriptors(JSON.parse(crux.calcDescriptors(args[0], args[1])), transfer);
    else if (op === 'calcDescriptors::timed') { // diagnosis: the Crux call and the parse (and packing) timed apart
      const started = performance.now();
      const text = crux.calcDescriptors(args[0], args[1]);
      const computed = performance.now();
      result = {values: packDescriptors(JSON.parse(text), transfer), wasmMs: computed - started,
        parseMs: performance.now() - computed};
    } else if (op === 'scaffoldTree')
      result = crux.scaffoldTree(args[0], args[1], args[2]);
    else if (op === 'getFingerprints') {
      assertFingerprintOptions(chemWasm as unknown as RDModule, args[0], args[3]); // a refusal: the RDKit pool's call
      const fps = await rdkitShapedWorker().getFingerprints(args[0], args[1], args[2], args[3]);
      result = {...fps, unread: unread(args[1], (i) => fps.fps[i] !== null)};
    } else if (op === 'getStructuralAlerts')
      result = structuralAlerts(args[0], args[1]);
    else if (op === 'convertMolNotation') { // the column conversions (engine/worker-jobs.ts)
      const results = await rdkitShapedWorker().convertMolNotation(args[0], args[1], args[2]);
      const none = notConverted(args[1]);
      result = {results, none, unconverted: (args[0] as string[]).flatMap((m, i) => m && results[i] === none ? [i] : [])};
    } else if (op === 'chemProperties') // Chem's property table (widgets/properties.ts)
      result = JSON.parse(chemWasm.chemProperties(args[0], args[1]));
    else if (op === 'chemQed') // Crux mode's drug-likeness: QED and its properties (widgets/drug-likeness.ts)
      result = JSON.parse(chemWasm.chemQed(args[0]));
    else if (op === 'chemAlerts') // Crux mode's toxicity: FilterCatalog's Brenk and NIH alerts (widgets/toxicity.ts)
      result = JSON.parse(chemWasm.chemAlerts(args[0], args[1], args[2]));
    else if (op === 'pairUncommonAtomsAndBonds')
      result = pairUncommonAtomsAndBonds(args[0], args[1], args[2]);
    else if (Object.hasOwn(ANALYSIS_JOBS, op)) // MMPA, R-groups, FMCS, fragment linking, molzip (crux-analysis-jobs.ts)
      result = await ANALYSIS_JOBS[op](rdkitShapedWorker(), chemWasm as unknown as RDModule, ...args);
    else if (op === 'fail::test') // tests: a job that fails with this message
      throw new Error(args[0]);
    else if (op === 'busy::test') { // tests: a job that keeps the worker busy, as a long Crux call does
      const end = performance.now() + args[0];
      let spins = 0;
      while (performance.now() < end)
        spins++;
      result = spins;
    } else
      throw new Error(`Unknown operation: ${op}`);
    port.postMessage({op: op, retval: result}, transfer);
  } catch (err: any) {
    // a trap reported as one: the service restarts the workers for it (and not for a refusal)
    port.postMessage({error: err instanceof WebAssembly.RuntimeError ? `RuntimeError: ${err.message}` :
      err instanceof Error ? err.message : String(err)});
  }
});
