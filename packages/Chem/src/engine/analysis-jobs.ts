/**
 * MMPA fragmentation, R-group decomposition, FMCS and fragment linking by the "Chemical engine" setting (README.md
 * beside this file), and molzip. In Crux mode the Crux pool runs the RDKit worker's own code on the Crux module
 * (crux/crux-analysis-jobs.ts), so the results have the RDKit pool's shapes; the RDKit pool does the molecules, rows or
 * clusters Crux cannot, and the whole call when Crux fails. In RDKit mode, the RDKit pool, as before. molzip runs on
 * Crux in both modes: Chem's RDKit.js (2024.09) has none.
 */
import type {RDModule} from '@datagrok-libraries/chem-meta/src/rdkit-api';
import type {IMmpFragmentsResult, IRGroupAnalysisResult} from '../rdkit-service/rdkit-service-worker-substructure';
import type {RGroupCluster} from '../crux/crux-analysis-jobs';
import {clusterMcsSmiles} from '../utils/mol-creation_rdkit';
import {getRdKitModule, getRdKitService} from '../utils/chem-common-rdkit';
import {getCruxService} from '../crux/crux-searches';
import {byEngineAsync} from './chem-engine';
import {cruxModuleSync} from './crux-module';
import {rdkitForRows} from './fallback';
import {forCrux} from './worker-jobs';

/** `RdKitService.mmpGetFragments(molecules)` (MMPA's single cuts) by the engine in effect. */
export async function mmpFragmentsByEngine(molecules: string[]): Promise<IMmpFragmentsResult> {
  const rdkit = async () => (await getRdKitService()).mmpGetFragments(molecules);
  return byEngineAsync('MMPA fragments', async () => {
    const {frags, smiles, unread} = await (await getCruxService()).mmpGetFragments(forCrux(molecules));
    await rdkitForRows('MMPA fragments', unread,
      async (rows) => (await getRdKitService()).mmpGetFragments(rows.map((i) => molecules[i])), (patch, k, row) => {
        frags[row] = patch.frags[k];
        smiles[row] = patch.smiles[k];
        return !!patch.smiles[k];
      }, (row) => molecules[row]);
    return {frags, smiles};
  }, rdkit);
}

/** `RdKitService.getRGroups` (one R-group decomposition) by the engine in effect; the caller holds the Chem critical
 * section, as there. */
export async function rGroupsByEngine(molecules: string[], core: string, coreIsQMol: boolean, options?: string):
  Promise<IRGroupAnalysisResult> {
  const rdkit = async () => (await getRdKitService()).getRGroups(molecules, core, coreIsQMol, options);
  return byEngineAsync('R-group decomposition', async () => {
    const result = await (await getCruxService()).rGroupAnalysis(molecules, core, coreIsQMol, options);
    if (!result)
      throw new Error('Crux reads not every molecule');
    return result;
  }, rdkit);
}

/** `RdKitService.getMCS` (FMCS's SMARTS of the molecules) by the engine in effect. */
export async function mcsByEngine(molecules: string[], exactAtoms: boolean, exactBonds: boolean): Promise<string> {
  const rdkit = async () => (await getRdKitService()).getMCS(molecules, exactAtoms, exactBonds);
  return byEngineAsync('MCS', async () => {
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
  return byEngineAsync('cluster MCS', async () => {
    // as RdKitService.clusterMCS: a one-molecule cluster is its own MCS, an empty one has none
    const mcs = clusters.map((molecules) => molecules.length === 1 ? molecules[0] : '');
    const multi = clusters.flatMap((molecules, i) => molecules.length > 1 ? [i] : []);
    const found = await (await getCruxService()).clusterMCS(multi.map((i) => forCrux(clusters[i])), exactAtoms,
      exactBonds);
    multi.forEach((i, k) => mcs[i] = found.results[k] ?? '');
    // the clusters Crux could not do: the RDKit pool's MCS of them (one past its time limit has none, as on RDKit);
    // each is a fallback, an MCS RDKit does not find too
    const rest = multi.filter((_, k) => found.results[k] === null && !found.timedOut.includes(k));
    await rdkitForRows('cluster MCS', rest,
      async (rows) => (await getRdKitService()).clusterMCS(rows.map((i) => clusters[i]), exactAtoms, exactBonds, true),
      (patch, k, i) => {
        mcs[i] = patch[k];
        return true;
      }, (i) => clusters[i][0]);
    if (rawSmarts)
      return mcs;
    const smiles = clusterMcsSmiles(forCrux(mcs), clusters, cruxModuleSync() as RDModule);
    // an MCS the Crux module does not read: RDKit's SMILES of it
    await rdkitForRows('cluster MCS', mcs.flatMap((it, i) => it && !smiles[i] ? [i] : []),
      async (rows) => clusterMcsSmiles(rows.map((i) => mcs[i]), rows.map((i) => clusters[i]), getRdKitModule()),
      (patch, k, i) => !!(smiles[i] = patch[k]), (i) => mcs[i]);
    return smiles;
  }, rdkit);
}

/** `RdKitService.clusterRGroups` (each cluster's R-group decomposition, null for none) by the engine in effect. */
export async function clusterRGroupsByEngine(clusters: RGroupCluster[]): Promise<(IRGroupAnalysisResult | null)[]> {
  const rdkit = async (some: RGroupCluster[]) => (await getRdKitService()).clusterRGroups(some);
  return byEngineAsync('cluster R-groups', async () => {
    const {results: result, timedOut} = await (await getCruxService()).clusterRGroups(clusters);
    // the clusters Crux could not do: the RDKit pool's decomposition of them (one past its time limit has none, as on
    // RDKit)
    const rest = clusters.flatMap((_, i) => result[i] === null && !timedOut.includes(i) ? [i] : []);
    await rdkitForRows('cluster R-groups', rest, (rows) => rdkit(rows.map((i) => clusters[i])),
      (patch, k, i) => (result[i] = patch[k]) !== null, (i) => clusters[i].molecules[0]);
    return result;
  }, () => rdkit(clusters));
}

/**
 * `RdKitService.linkRGroupFragments` (each row's core joined with its fragments, canonical SMILES, '' for none) by the
 * engine in effect: in Crux mode by RDKit's molzip, the RDKit pool's string assembly for a row molzip does not join.
 */
export async function linkRGroupFragmentsByEngine(cores: string[], fragmentColumns: string[][],
  attachIdx: (number | number[])[]): Promise<string[]> {
  const rdkit = async (rows: number[] | null) => (await getRdKitService()).linkRGroupFragments(
    rows ? rows.map((i) => cores[i]) : cores,
    rows ? fragmentColumns.map((c) => rows.map((i) => c[i])) : fragmentColumns, attachIdx);
  return byEngineAsync('fragment linking', async () => {
    const smiles = await (await getCruxService()).linkRGroupFragments(forCrux(cores), fragmentColumns, attachIdx);
    await rdkitForRows('fragment linking', smiles.flatMap((s, i) => s ? [] : [i]), rdkit,
      (patch, k, i) => !!(smiles[i] = patch[k]), (i) => cores[i]);
    return smiles;
  }, () => rdkit(null));
}

/**
 * `RdKitService.mmpLinkFragments` (each core joined with its fragment at `[*:1]`, canonical SMILES, '' for none: the
 * MMP viewer's generations) by the engine in effect: in Crux mode the RDKit worker's
 * code on the Crux workers, the RDKit pool for a row Crux does not join.
 */
export async function mmpLinkFragmentsByEngine(cores: string[], fragments: string[]): Promise<string[]> {
  const rdkit = async (rows: number[] | null) => (await getRdKitService()).mmpLinkFragments(
    rows ? rows.map((i) => cores[i]) : cores, rows ? rows.map((i) => fragments[i]) : fragments);
  return byEngineAsync('mmpLinkFragments', async () => {
    const smiles = await (await getCruxService()).mmpLinkFragments(forCrux(cores), fragments);
    await rdkitForRows('mmpLinkFragments', smiles.flatMap((s, i) => s || !cores[i] ? [] : [i]), rdkit,
      (patch, k, i) => !!(smiles[i] = patch[k]), (i) => cores[i]);
    return smiles;
  }, () => rdkit(null));
}

/**
 * RDKit's MolzipParams by name (`RDKit::MolzipParams`, Code/GraphMol/ChemTransforms/MolFragmenter.h), each optional
 * with RDKit's default: `Label` (how pieces name their attachment points: `AtomMapNumber`, `[*:1]`, the
 * default; `Isotope`, `[1*]`; `FragmentOnBonds`, the labels of RDKit's FragmentOnBonds; `AtomType`, the atoms of
 * `AtomSymbols`; `AtomProperty`, the unsigned-int atom property `AtomProperty`); `EnforceValenceRules` (true: a join
 * over an atom's valence is refused).
 */
export interface MolzipParams {
  Label?: 'AtomMapNumber' | 'Isotope' | 'FragmentOnBonds' | 'AtomType' | 'AtomProperty';
  AtomSymbols?: string[];
  AtomProperty?: string;
  EnforceValenceRules?: boolean;
}

/**
 * RDKit's molzip (`molzip(mol, params)`, crux-core's port, run on the Crux pool in both engines): each molecule holds
 * its pieces as fragments (SMILES `C[*:1].N[*:1]`, or a MOL block of several fragments), joined where two of them
 * carry the same attachment point (a label without a partner stays a dummy atom). The joined molecules' canonical
 * SMILES, '' for a molecule that does not read or that molzip refuses (a label used twice, a valence over the limit;
 * the first refusal is logged). Throws when the Crux module does not load.
 */
export async function molzip(molecules: string[], params: MolzipParams = {}): Promise<string[]> {
  const {smiles, error} = await (await getCruxService()).molzip(molecules, JSON.stringify(params));
  if (error)
    console.warn(`Chem | molzip did not join every molecule: ${error}`);
  return smiles;
}
