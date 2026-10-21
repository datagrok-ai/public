// The Crux worker's analysis jobs (crux.worker.ts runs them; engine/analysis-jobs.ts sends them): MMPA fragmentation,
// R-group decomposition, FMCS and fragment linking run the RDKit worker's own code
// (rdkit-service-worker-substructure.ts) on the Crux module, which answers RDKit.js's calls; molzip is RDKit's molzip,
// which Chem's RDKit.js lacks.
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {RdKitServiceWorker} from '../rdkit-service/rdkit-service-worker';
import {IMmpFragmentsResult, IRGroupAnalysisResult, hasMolzip, zippedSmiles}
  from '../rdkit-service/rdkit-service-worker-substructure';
import {getMolSafe} from '../utils/mol-creation_rdkit';
import {CLUSTER_MCS_TIMEOUT_S} from './crux-api';

/** A cluster's R-group decomposition, as RdKitService.clusterRGroups takes it. */
export interface RGroupCluster {molecules: string[], core: string, coreIsQMol: boolean, options?: string}

/** The rows of `molecules` Crux reads no molecule from (an empty one is no such row). */
function unreadRows(molecules: string[], module: RDModule): number[] {
  const rows: number[] = [];
  for (let i = 0; i < molecules.length; i++) {
    if (!molecules[i])
      continue;
    const mol = getMolSafe(molecules[i], {}, module).mol;
    if (mol)
      mol.delete();
    else
      rows.push(i);
  }
  return rows;
}

/** Whether get_mmpa_frags reads the molecule and returns, finding no fragmentation (its lists null, as RDKit.js gives
 * them where fragmentMol returns false), where RdKitServiceWorker.mmpGetFragments gave no SMILES. */
function noFragmentation(worker: RdKitServiceWorker, molecule: string): boolean {
  const mol = worker.getMolWithSmilesCheck(molecule);
  try {
    const frags = mol?.get_mmpa_frags(1, 1, 20);
    return !!frags && !frags.cores && !frags.sidechains;
  } catch {
    return false;
  } finally {
    mol?.delete();
  }
}

/** Each cluster's result, one after another; null where Crux cannot read one of its molecules or the call fails (the
 * RDKit pool does that cluster). */
async function perCluster<C, R>(clusters: C[], molecules: (cluster: C) => string[], module: RDModule,
  run: (cluster: C) => R | Promise<R>): Promise<(R | null)[]> {
  const results: (R | null)[] = [];
  for (const cluster of clusters) {
    try {
      results.push(unreadRows(molecules(cluster), module).length > 0 ? null : await run(cluster));
    } catch {
      results.push(null);
    }
  }
  return results;
}

/**
 * The pieces of each molecule joined by RDKit's molzip (`details`: MolzipParams as JSON), each molecule read as Chem's
 * getMolSafe reads it, the result's canonical SMILES ({@link zippedSmiles}); '' for a molecule not read or one molzip
 * refuses (a label used twice, a valence over the limit), and the first refusal's message.
 */
function molzip(module: RDModule, molecules: string[], details: string): {smiles: string[], error: string | null} {
  if (!hasMolzip(module))
    throw new Error('not supported by Crux: molzip'); // a Crux file older than crux-core's molzip
  let error: string | null = null;
  const smiles = molecules.map((molecule) => {
    const read = molecule ? getMolSafe(molecule, {}, module) : null;
    try {
      return !read?.mol || read.isQMol ? '' : zippedSmiles(module.molzip(read.mol, details), module);
    } catch (e: any) {
      error ??= `${molecule}: ${e?.message ?? e}`;
      return '';
    } finally {
      read?.mol?.delete();
    }
  });
  return {smiles, error};
}

type Job = (worker: RdKitServiceWorker, module: RDModule, ...args: any[]) => unknown;

/** The jobs by name (the Crux worker's `op`), each with the RDKit worker over the Crux module and the module. */
export const ANALYSIS_JOBS: {[op: string]: Job} = {
  /**
   * MMPA fragments as RdKitServiceWorker.mmpGetFragments gives them, and the rows the RDKit pool does (`unread`): a
   * molecule Crux did not read, or whose fragmentation threw (a "not supported by Crux" refusal, say). RDKit's code
   * gives both no SMILES and no fragments; a molecule RDKit itself does not fragment (get_mmpa_frags gives no lists:
   * fragmentMol found no cut) gets the same from Crux and stays.
   */
  mmpGetFragments: (worker, module, molecules: string[]): IMmpFragmentsResult & {unread: number[]} => {
    const result = worker.mmpGetFragments(molecules);
    const unread: number[] = [];
    for (let i = 0; i < molecules.length; i++) {
      if (molecules[i] && !result.smiles[i] && !noFragmentation(worker, molecules[i]))
        unread.push(i);
    }
    return {...result, unread};
  },
  /** One R-group decomposition (RdKitServiceWorker.rGroupAnalysis); null where Crux cannot read a molecule. */
  rGroupAnalysis: async (worker, module, molecules: string[], core: string, coreIsQMol: boolean, options?: string):
    Promise<IRGroupAnalysisResult | null> => unreadRows(molecules, module).length > 0 ? null :
    worker.rGroupAnalysis(molecules, core, coreIsQMol, options),
  /** The FMCS SMARTS of the molecules (RdKitServiceWorker.mostCommonStructure); null where Crux cannot read one. */
  mostCommonStructure: (worker, module, molecules: string[], exactAtoms: boolean, exactBonds: boolean): string | null =>
    unreadRows(molecules, module).length > 0 ? null : worker.mostCommonStructure(molecules, exactAtoms, exactBonds),
  /** Each cluster's FMCS SMARTS within {@link CLUSTER_MCS_TIMEOUT_S} ('' past it); null for one the RDKit pool does. */
  clusterMCS: (worker, module, clusters: string[][], exactAtoms: boolean, exactBonds: boolean) =>
    perCluster(clusters, (molecules) => molecules, module,
      (molecules) => worker.mostCommonStructure(molecules, exactAtoms, exactBonds, CLUSTER_MCS_TIMEOUT_S)),
  /** Each cluster's R-group decomposition; null for one the RDKit pool does. */
  clusterRGroups: (worker, module, clusters: RGroupCluster[]) => perCluster(clusters, (c) => c.molecules, module,
    (c) => worker.rGroupAnalysis(c.molecules, c.core, c.coreIsQMol, c.options)),
  /** RdKitServiceWorker.linkRGroupFragments of rows `[core, a fragment per column...]`: by molzip on Crux. */
  linkRGroupFragments: (worker, module, rows: string[][], attachIdx: (number | number[])[]): string[] => {
    if (!hasMolzip(module))
      throw new Error('not supported by Crux: molzip');
    return worker.linkRGroupFragments(rows.map((row) => row[0]),
      attachIdx.map((_, p) => rows.map((row) => row[p + 1])), attachIdx);
  },
  molzip: (_worker, module, molecules: string[], details: string) => molzip(module, molecules, details),
};
