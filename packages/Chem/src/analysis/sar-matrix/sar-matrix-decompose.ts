import {getRdKitModule, getRdKitService} from '../../utils/chem-common-rdkit';
import {IRGroupAnalysisResult} from '../../rdkit-service/rdkit-service-worker-substructure';
import {logSarTime} from './sar-matrix-types';

const MIN_ANCHOR_HEAVY_ATOMS = 4;
export const MAX_SAR_CLUSTER_SIZE = 300;
/** Permissive floor: anchoring on the ring scaffold makes the core a small fraction of a large
 *  molecule by design, so this only rejects a degenerate anchor (median coverage across cluster). */
const MIN_MCS_COVERAGE = 0.2;
const MIN_CLEAN_FRACTION = 0.5;

const R_GROUP_OPTIONS = JSON.stringify({
  matchingStrategy: 'Greedy',
  includeTargetMolInResults: true,
  onlyMatchAtRGroups: false,
});

export interface PositionRecord {
  molIdx: number;
  coreSmiles: string;
  values: {[position: string]: string};
}

export interface ClusterDecomposition {
  records: PositionRecord[];
  positions: string[];
}

function positionNumber(position: string): number {
  return Number.parseInt(position.replace(/^\D+/, ''), 10);
}

/** Clean = single connected fragment with exactly one attachment point that is this position's own
 *  dummy (empty value = H, also clean). Rejects bad-anchor artifacts: multiple/foreign dummies or
 *  disconnected merges. */
function isCleanSubstituent(fragment: string, expectedNumber: number): boolean {
  if (!fragment)
    return true;
  if (fragment.includes('.'))
    return false;
  const dummies = [...fragment.matchAll(/\[\*:(\d+)\]/g)];
  return dummies.length === 1 && Number.parseInt(dummies[0][1], 10) === expectedNumber;
}

function median(values: number[]): number {
  if (values.length === 0)
    return 0;
  const sorted = [...values].sort((a, b) => a - b);
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
}

// Each miss is a main-thread RDKit parse and keys recur heavily across clusters. Cleared wholesale
// past the cap so neither map grows unbounded across a session.
const heavyAtomCache = new Map<string, number>();
const canonicalCoreCache = new Map<string, string>();
const DECOMPOSE_CACHE_MAX = 10000;

function heavyAtomCount(rdkit: ReturnType<typeof getRdKitModule>, smiles: string): number {
  const cached = heavyAtomCache.get(smiles);
  if (cached !== undefined)
    return cached;
  let mol = null;
  let count = 0;
  try {
    mol = rdkit.get_mol(smiles);
    count = mol && mol.is_valid() ? mol.get_num_atoms() : 0;
  } catch {
    count = 0;
  } finally {
    mol?.delete();
  }
  if (heavyAtomCache.size >= DECOMPOSE_CACHE_MAX)
    heavyAtomCache.clear();
  heavyAtomCache.set(smiles, count);
  return count;
}

/** Canonical SMILES of a matched core (the matrix row key), or the raw string when unparseable. */
function canonicalCore(rdkit: ReturnType<typeof getRdKitModule>, coreRaw: string): string {
  const cached = canonicalCoreCache.get(coreRaw);
  if (cached !== undefined)
    return cached;
  let coreSmiles = coreRaw;
  let coreMol = null;
  try {
    coreMol = rdkit.get_mol(coreRaw);
    if (coreMol?.is_valid())
      coreSmiles = coreMol.get_smiles();
  } catch {
    // keep the raw string
  } finally {
    coreMol?.delete();
  }
  if (canonicalCoreCache.size >= DECOMPOSE_CACHE_MAX)
    canonicalCoreCache.clear();
  canonicalCoreCache.set(coreRaw, coreSmiles);
  return coreSmiles;
}

interface ClusterPrep {
  /** Index into the caller's cluster array, where this cluster's result lands. */
  idx: number;
  molIdx: number[];
  smiles: string[];
  anchor: string | null;
  /** True when the anchor is a query SMARTS (MCS with generic atoms), parsed by the worker via
   *  `get_qmol` rather than `get_mol`. */
  anchorIsQuery: boolean;
}

/** Reject an anchor too small to carry a shared multi-position core, or covering too little of the
 *  molecules — the leftover would become oversized R-groups that decompose badly. */
function anchorIsUsable(prep: ClusterPrep, rdkit: ReturnType<typeof getRdKitModule>): boolean {
  if (!prep.anchor)
    return false;
  let anchorHeavyAtoms = 0;
  let anchorMol = null;
  try {
    anchorMol = rdkit.get_qmol(prep.anchor);
    anchorHeavyAtoms = anchorMol?.get_num_atoms() ?? 0;
  } catch {
    return false;
  } finally {
    anchorMol?.delete();
  }
  if (anchorHeavyAtoms < MIN_ANCHOR_HEAVY_ATOMS)
    return false;
  const coverages = prep.smiles
    .map((s) => heavyAtomCount(rdkit, s))
    .filter((n) => n > 0)
    .map((n) => anchorHeavyAtoms / n);
  return median(coverages) >= MIN_MCS_COVERAGE;
}

/** Build records from one cluster's RGD output, dropping molecules with any malformed substituent
 *  or a multi-component core. Returns null (caller falls back) when the RGD failed, assigned no R
 *  positions, or too few molecules decomposed cleanly. */
function buildDecomposition(prep: ClusterPrep, res: IRGroupAnalysisResult | null,
  rdkit: ReturnType<typeof getRdKitModule>): ClusterDecomposition | null {
  // res.smiles[0] is the Core column; the rest are R1, R2, ...
  if (!res || res.colNames.length < 2)
    return null;
  const positions = res.colNames.slice(1);
  const positionNums = positions.map(positionNumber);
  const records: PositionRecord[] = [];
  let matched = 0;
  for (let i = 0; i < prep.molIdx.length; i++) {
    const coreRaw = res.smiles[0][i];
    if (!coreRaw)
      continue; // molecule didn't match the anchor
    matched++;

    const values: {[position: string]: string} = {};
    let clean = true;
    for (let p = 0; p < positions.length; p++) {
      const value = res.smiles[p + 1][i] ?? '';
      if (!isCleanSubstituent(value, positionNums[p])) {
        clean = false;
        break;
      }
      values[positions[p]] = value;
    }
    if (!clean || coreRaw.includes('.'))
      continue;

    records.push({molIdx: prep.molIdx[i], coreSmiles: canonicalCore(rdkit, coreRaw), values});
  }

  if (records.length < 2 || records.length < matched * MIN_CLEAN_FRACTION)
    return null;
  return {records, positions};
}

/**
 * Decomposes every cluster against one shared anchor scaffold each, so R1, R2, ... come out aligned
 * by construction. Anchor priority per cluster: a single shared Bemis-Murcko ring scaffold; else the
 * generic MCS of several related scaffolds; else the MCS of the whole molecules. Returns one entry
 * per input cluster, or null where no useful clean anchor exists (caller falls back to the
 * single-position construction). Batched through a worker-per-cluster queue so independent clusters
 * run in parallel and stuck MCS/RGD workers are killed after a timeout rather than hanging the run.
 */
export async function decomposeClusters(clusterMembers: number[][], molecules: string[],
  scaffolds: string[] = []): Promise<(ClusterDecomposition | null)[]> {
  const results = new Array<ClusterDecomposition | null>(clusterMembers.length).fill(null);
  const rdkit = getRdKitModule();
  const service = await getRdKitService();

  // Phase 1 — pick each cluster's anchor source, queueing the ones that need an MCS.
  const preps: ClusterPrep[] = [];
  const scaffoldMcs: {input: string[], prep: ClusterPrep}[] = [];
  const wholeMolMcs: {input: string[], prep: ClusterPrep}[] = [];
  for (let i = 0; i < clusterMembers.length; i++) {
    const molIdx = clusterMembers[i];
    if (molIdx.length < 2 || molIdx.length > MAX_SAR_CLUSTER_SIZE)
      continue;
    const prep: ClusterPrep = {
      idx: i, molIdx, smiles: molIdx.map((k) => molecules[k]), anchor: null, anchorIsQuery: false,
    };
    preps.push(prep);
    const clusterScaffolds = scaffolds.length ?
      [...new Set(molIdx.map((k) => scaffolds[k]).filter((s) => s))] : [];
    if (clusterScaffolds.length === 1)
      prep.anchor = clusterScaffolds[0]; // concrete substructure, no MCS needed
    else if (clusterScaffolds.length >= 2) {
      prep.anchorIsQuery = true;
      scaffoldMcs.push({input: clusterScaffolds, prep});
    } else {
      prep.anchorIsQuery = true;
      wholeMolMcs.push({input: prep.smiles, prep});
    }
  }

  // Phase 2 — MCS passes. `rawSmarts` (last arg): under Any-atom comparison the MCS is a query
  // SMARTS whose generic atoms a SMILES round-trip would destroy.
  if (scaffoldMcs.length > 0) {
    const t = performance.now();
    const mcs = await service.clusterMCS(scaffoldMcs.map((q) => q.input), false, true, true);
    logSarTime(`MCS anchors over scaffolds (${scaffoldMcs.length} clusters)`, t);
    scaffoldMcs.forEach((q, j) => {
      if (mcs[j])
        q.prep.anchor = mcs[j];
      else // scaffold MCS failed — retry over the whole molecules
        wholeMolMcs.push({input: q.prep.smiles, prep: q.prep});
    });
  }
  if (wholeMolMcs.length > 0) {
    const t = performance.now();
    const mcs = await service.clusterMCS(wholeMolMcs.map((q) => q.input), false, true, true);
    logSarTime(`MCS anchors over whole molecules (${wholeMolMcs.length} clusters)`, t);
    wholeMolMcs.forEach((q, j) => q.prep.anchor = mcs[j] || null);
  }

  // Phase 3 — validate anchors on the main thread, then batch-decompose the survivors.
  const decomposable = preps.filter((prep) => anchorIsUsable(prep, rdkit));
  if (decomposable.length === 0)
    return results;

  const t = performance.now();
  const rgd = await service.clusterRGroups(decomposable.map((prep) => ({
    molecules: prep.smiles, core: prep.anchor!, coreIsQMol: prep.anchorIsQuery, options: R_GROUP_OPTIONS,
  })));
  logSarTime(`R-group decomposition (${decomposable.length} clusters)`, t);

  decomposable.forEach((prep, j) => results[prep.idx] = buildDecomposition(prep, rgd[j], rdkit));
  return results;
}
