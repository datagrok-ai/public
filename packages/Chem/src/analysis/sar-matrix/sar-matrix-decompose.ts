import {getRdKitModule, getRdKitService} from '../../utils/chem-common-rdkit';
import {IRGroupAnalysisResult} from '../../rdkit-service/rdkit-service-worker-substructure';
import {logSarTime} from './sar-matrix-types';

/** An anchor scaffold with fewer heavy atoms than this can't carry a useful shared multi-position
 *  core; the cluster falls back to the single-position construction instead. */
const MIN_ANCHOR_HEAVY_ATOMS = 4;
/** R-Group Decomposition over a cluster larger than this is slow relative to what it adds; beyond
 *  this size the cluster falls back to the single-position construction. */
export const MAX_SAR_CLUSTER_SIZE = 300;
/** The anchor must cover at least this fraction of a molecule's heavy atoms (median across the
 *  cluster). Anchoring on the ring scaffold makes the core a legitimately small fraction of a large
 *  molecule (big R-groups are expected), so this is a permissive floor that only rejects a
 *  degenerate anchor (e.g. a lone ring shared across otherwise-unrelated molecules). */
const MIN_MCS_COVERAGE = 0.2;
/** A cluster is rejected (falls back to single-position) unless at least this fraction of the
 *  molecules that matched the anchor decompose into clean, single-attachment substituents. */
const MIN_CLEAN_FRACTION = 0.5;

/** Pre-serialized: the worker call takes the RGD options as a JSON string. */
const R_GROUP_OPTIONS = JSON.stringify({
  matchingStrategy: 'Greedy',
  includeTargetMolInResults: true,
  onlyMatchAtRGroups: false,
});

/** Per-molecule outcome of one cluster's anchor-scaffold R-Group Decomposition. */
export interface PositionRecord {
  molIdx: number;
  /** Canonical SMILES of the concrete matched core (row key), carrying `[*:N]` dummies. */
  coreSmiles: string;
  /** Substituent SMILES (each carrying its own `[*:N]` dummy) keyed by position name ('R1', ...). */
  values: {[position: string]: string};
}

export interface ClusterDecomposition {
  records: PositionRecord[];
  /** Position names R-Group Decomposition assigned, in its own order. */
  positions: string[];
}

/** Position number embedded in an R-group column name ('R7' -> 7), or NaN. */
function positionNumber(position: string): number {
  return Number.parseInt(position.replace(/^\D+/, ''), 10);
}

/**
 * A substituent is clean only if it is a single connected fragment carrying exactly one attachment
 * point, and that point is this position's own dummy. This rejects the artifacts a bad anchor
 * produces: two attachment points (`Cc(c:[*:7])on:[*:7]`), disconnected merges
 * (`FC(F)(F)CN[*:2].O=[*:2]`), or a foreign position's dummy (`...[*:4])\\[*:7]`).
 * An empty value (H at this position) is clean.
 */
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

/** Every miss here is a main-thread RDKit parse, and both keys recur heavily: the coverage median
 *  parses each cluster member, and clusters share members across the run, while one cluster's records
 *  collapse to a handful of distinct cores. Cleared wholesale past the cap so neither map can grow
 *  unbounded across a session. */
const heavyAtomCache = new Map<string, number>();
const canonicalCoreCache = new Map<string, string>();
const DECOMPOSE_CACHE_MAX = 10000;

/** Heavy-atom count of a SMILES, or 0 if it can't be parsed. */
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

/** Canonical SMILES of a matched core (the matrix row key), or the raw string when RDKit can't parse
 *  it. Memoized: the record loop canonicalizes the same few cores once per matched molecule. */
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

/** Per-cluster state threaded through the batched phases of {@link decomposeClusters}. */
interface ClusterPrep {
  /** Index into the caller's cluster array, where this cluster's result lands. */
  idx: number;
  molIdx: number[];
  smiles: string[];
  /** The scaffold the cluster is decomposed against; null until resolved (or unresolvable). */
  anchor: string | null;
  /** True when the anchor is an MCS (a query SMARTS with generic atoms) rather than a concrete
   *  structure — the RGD worker then parses it with `get_qmol`. */
  anchorIsQuery: boolean;
}

/** #2/#3 — reject an anchor too small to carry a shared multi-position core, or one covering too
 *  little of the molecules: the leftover would become oversized R-groups that decompose badly. */
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

/**
 * #1 — records from one cluster's RGD output: validate every substituent and drop the whole molecule
 * if any is malformed, and reject a multi-component core (a bad match). Returns `null` when the RGD
 * itself failed or its worker was killed as stuck, when it assigned no R positions, or when too few
 * of the matched molecules decomposed cleanly — the anchor is wrong for this cluster, so it falls
 * back rather than building a matrix from a handful of survivors.
 */
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
 * Decomposes every cluster against one shared anchor scaffold each (the same
 * `exactAtomSearch=false, exactBondSearch=true` MCS the interactive R-Groups Analysis "MCS" button
 * uses), so every molecule's R1, R2, ... columns come out aligned by construction.
 *
 * The anchor per cluster: the shared ring **scaffold** (Bemis-Murcko, computed upstream) when there
 * is exactly one — anchoring there rather than on the MCS of whole molecules keeps peripheral
 * variation out of the core. Several related scaffolds (a ring hop) anchor on their generic MCS, so
 * both match with aligned R-labels; no scaffolds (Murcko unavailable) falls back to the MCS of the
 * whole molecules.
 *
 * Batched on purpose: the MCS anchors and the R-group decompositions each go through a
 * worker-per-cluster queue ({RdKitService.clusterMCS} / `clusterRGroups`), so independent
 * clusters run truly in parallel — calling the per-cluster service entries under a `Promise.all`
 * would only serialize them all on the chem critical section, on a single worker. The queue also
 * survives RDKit's known MCS/RGD hangs: a stuck worker is killed after a timeout and its cluster
 * comes back empty.
 *
 * Returns one entry per input cluster; `null` where the cluster has no useful, clean shared
 * multi-position anchor — too small or too large, an anchor covering too little of the molecules,
 * too many malformed substituents, or a stuck/failed worker. Callers then fall back to the
 * single-position construction, so a heterogeneous cluster degrades to a simpler view instead of
 * rendering garbage fragments.
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
      prep.anchor = clusterScaffolds[0]; // one shared scaffold — a concrete substructure, no MCS needed
    else if (clusterScaffolds.length >= 2) {
      prep.anchorIsQuery = true;
      scaffoldMcs.push({input: clusterScaffolds, prep});
    } else {
      prep.anchorIsQuery = true;
      wholeMolMcs.push({input: prep.smiles, prep});
    }
  }

  // Phase 2 — the MCS passes, each one batched call. `rawSmarts`: under Any-atom comparison the MCS
  // is a query SMARTS whose generic atoms a SMILES round-trip would destroy.
  if (scaffoldMcs.length > 0) {
    const t = performance.now();
    const mcs = await service.clusterMCS(scaffoldMcs.map((q) => q.input), false, true, true);
    logSarTime(`MCS anchors over scaffolds (${scaffoldMcs.length} clusters)`, t);
    scaffoldMcs.forEach((q, j) => {
      if (mcs[j])
        q.prep.anchor = mcs[j];
      else // the scaffold MCS failed (or its worker hung) — retry over the whole molecules
        wholeMolMcs.push({input: q.prep.smiles, prep: q.prep});
    });
  }
  if (wholeMolMcs.length > 0) {
    const t = performance.now();
    const mcs = await service.clusterMCS(wholeMolMcs.map((q) => q.input), false, true, true);
    logSarTime(`MCS anchors over whole molecules (${wholeMolMcs.length} clusters)`, t);
    wholeMolMcs.forEach((q, j) => q.prep.anchor = mcs[j] || null);
  }

  // Phase 3 — validate the anchors on the main thread (cheap, memoized RDKit calls)...
  const decomposable = preps.filter((prep) => anchorIsUsable(prep, rdkit));
  if (decomposable.length === 0)
    return results;

  // ...and run every surviving cluster's R-Group Decomposition in one batched, parallel pass.
  const t = performance.now();
  const rgd = await service.clusterRGroups(decomposable.map((prep) => ({
    molecules: prep.smiles, core: prep.anchor!, coreIsQMol: prep.anchorIsQuery, options: R_GROUP_OPTIONS,
  })));
  logSarTime(`R-group decomposition (${decomposable.length} clusters)`, t);

  // Phase 4 — per-cluster post-processing, all synchronous.
  decomposable.forEach((prep, j) => results[prep.idx] = buildDecomposition(prep, rgd[j], rdkit));
  return results;
}
