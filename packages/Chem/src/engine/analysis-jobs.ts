/**
 * MMPA fragmentation, R-group decomposition and FMCS by the "Chemical engine" setting (README.md beside this file). In
 * Crux mode the Crux pool runs the RDKit worker's own code on the Crux module (crux/crux-analysis-jobs.ts), so the
 * results have the RDKit pool's shapes; the RDKit pool does the molecules or clusters Crux cannot, and the whole call
 * when Crux fails. In RDKit mode, the RDKit pool, as before.
 */
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {IMmpFragmentsResult, IRGroupAnalysisResult} from '../rdkit-service/rdkit-service-worker-substructure';
import type {RGroupCluster} from '../crux/crux-analysis-jobs';
import {clusterMcsSmiles} from '../rdkit-service/rdkit-service';
import {getRdKitModule, getRdKitService} from '../utils/chem-common-rdkit';
import {getCruxService} from '../crux/crux-searches';
import {isCruxRequested} from './chem-engine';
import {cruxModuleSync} from './crux-module';
import {cruxOrRdKitAsync, warnCruxFallback} from './fallback';
import {forCrux} from './worker-jobs';

/** `RdKitService.mmpGetFragments(molecules)` (MMPA's single cuts) by the engine in effect. */
export async function mmpFragmentsByEngine(molecules: string[]): Promise<IMmpFragmentsResult> {
  const rdkit = async () => (await getRdKitService()).mmpGetFragments(molecules);
  if (!isCruxRequested())
    return rdkit();
  return cruxOrRdKitAsync('MMPA fragments', async () => {
    const {frags, smiles, unread} = await (await getCruxService()).mmpGetFragments(forCrux(molecules));
    if (unread.length) { // the molecules Crux did not read: the RDKit pool's fragments for them
      const patch = await (await getRdKitService()).mmpGetFragments(unread.map((i) => molecules[i]));
      unread.forEach((row, k) => {
        frags[row] = patch.frags[k];
        smiles[row] = patch.smiles[k];
      });
      const done = unread.filter((_, k) => patch.smiles[k]);
      if (done.length) {
        warnCruxFallback('MMPA fragments', `RDKit fragmented ${done.length} molecules Crux could not, e.g. ` +
          molecules[done[0]]);
      }
    }
    return {frags, smiles};
  }, rdkit);
}

/** `RdKitService.getRGroups` (one R-group decomposition) by the engine in effect; the caller holds the Chem critical
 * section, as there. */
export async function rGroupsByEngine(molecules: string[], core: string, coreIsQMol: boolean, options?: string):
  Promise<IRGroupAnalysisResult> {
  const rdkit = async () => (await getRdKitService()).getRGroups(molecules, core, coreIsQMol, options);
  if (!isCruxRequested())
    return rdkit();
  return cruxOrRdKitAsync('R-group decomposition', async () => {
    const result = await (await getCruxService()).rGroupAnalysis(molecules, core, coreIsQMol, options);
    if (!result)
      throw new Error('Crux reads not every molecule');
    return result;
  }, rdkit);
}

/** `RdKitService.getMCS` (FMCS's SMARTS of the molecules) by the engine in effect. */
export async function mcsByEngine(molecules: string[], exactAtoms: boolean, exactBonds: boolean): Promise<string> {
  const rdkit = async () => (await getRdKitService()).getMCS(molecules, exactAtoms, exactBonds);
  if (!isCruxRequested())
    return rdkit();
  return cruxOrRdKitAsync('MCS', async () => {
    const mcs = await (await getCruxService()).mostCommonStructure(molecules, exactAtoms, exactBonds);
    if (mcs === null)
      throw new Error('Crux reads not every molecule');
    return mcs;
  }, rdkit);
}

/**
 * `RdKitService.clusterMCS` (each cluster's MCS: SMILES, or FMCS's SMARTS with `rawSmarts`) by the engine in effect; in
 * Crux mode the SMILES are written by the Crux module, by RDKit's for an MCS (a one-molecule cluster's molecule) Crux
 * does not read.
 */
export async function clusterMcsByEngine(clusters: string[][], exactAtoms: boolean, exactBonds: boolean,
  rawSmarts = false): Promise<string[]> {
  const rdkit = async () => (await getRdKitService()).clusterMCS(clusters, exactAtoms, exactBonds, rawSmarts);
  if (!isCruxRequested())
    return rdkit();
  return cruxOrRdKitAsync('cluster MCS', async () => {
    // as RdKitService.clusterMCS: a one-molecule cluster is its own MCS, an empty one has none
    const mcs = clusters.map((molecules) => molecules.length === 1 ? molecules[0] : '');
    const multi = clusters.flatMap((molecules, i) => molecules.length > 1 ? [i] : []);
    const found = await (await getCruxService()).clusterMCS(multi.map((i) => forCrux(clusters[i])), exactAtoms,
      exactBonds);
    multi.forEach((i, k) => mcs[i] = found.results[k] ?? '');
    // the clusters Crux could not do: the RDKit pool's MCS of them (one past its time limit has none, as on RDKit)
    const rest = multi.filter((_, k) => found.results[k] === null && !found.timedOut.includes(k));
    if (rest.length) {
      const patch = await (await getRdKitService()).clusterMCS(rest.map((i) => clusters[i]), exactAtoms, exactBonds,
        true);
      rest.forEach((i, k) => mcs[i] = patch[k]);
      warnCruxFallback('cluster MCS', `RDKit found the MCS of ${rest.length} clusters, e.g. ${clusters[rest[0]][0]}`);
    }
    if (rawSmarts)
      return mcs;
    const smiles = clusterMcsSmiles(forCrux(mcs), clusters, cruxModuleSync() as RDModule);
    const unwritten = mcs.flatMap((it, i) => it && !smiles[i] ? [i] : []);
    if (unwritten.length) { // an MCS the Crux module does not read: RDKit's SMILES of it
      const patch = clusterMcsSmiles(unwritten.map((i) => mcs[i]), unwritten.map((i) => clusters[i]), getRdKitModule());
      unwritten.forEach((i, k) => smiles[i] = patch[k]);
      if (patch.some((s) => s))
        warnCruxFallback('cluster MCS', `RDKit wrote the SMILES of ${unwritten.length} MCS, e.g. ${mcs[unwritten[0]]}`);
    }
    return smiles;
  }, rdkit);
}

/** `RdKitService.clusterRGroups` (each cluster's R-group decomposition, null for none) by the engine in effect. */
export async function clusterRGroupsByEngine(clusters: RGroupCluster[]): Promise<(IRGroupAnalysisResult | null)[]> {
  const rdkit = async (some: RGroupCluster[]) => (await getRdKitService()).clusterRGroups(some);
  if (!isCruxRequested())
    return rdkit(clusters);
  return cruxOrRdKitAsync('cluster R-groups', async () => {
    const {results: result, timedOut} = await (await getCruxService()).clusterRGroups(clusters);
    // the clusters Crux could not do: the RDKit pool's decomposition of them (one past its time limit has none, as on
    // RDKit)
    const rest = clusters.flatMap((_, i) => result[i] === null && !timedOut.includes(i) ? [i] : []);
    if (rest.length) {
      const patch = await rdkit(rest.map((i) => clusters[i]));
      rest.forEach((i, k) => result[i] = patch[k]);
      const done = rest.filter((_, k) => patch[k] !== null);
      if (done.length) {
        warnCruxFallback('cluster R-groups', `RDKit decomposed ${done.length} clusters Crux could not, e.g. ` +
          clusters[done[0]].molecules[0]);
      }
    }
    return result;
  }, () => rdkit(clusters));
}
