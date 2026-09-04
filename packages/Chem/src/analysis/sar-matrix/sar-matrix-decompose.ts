import {getRdKitModule, getRdKitService} from '../../utils/chem-common-rdkit';
import {IRGroupAnalysisResult} from '../../rdkit-service/rdkit-service-worker-substructure';
import {logSarTime} from './sar-matrix-types';

const MIN_ANCHOR_HEAVY_ATOMS = 4;
export const MAX_SAR_CLUSTER_SIZE = 300;
/** Counted, not timed: a timeout retires the worker running it, and the clusters queued behind that
 *  worker are then dropped, so a busier machine yields fewer matrices over the same data. The MCS is
 *  a subgraph of the cluster's smallest member, so that member's atom count bounds the search. Past
 *  the ceiling a cluster gets no anchor and falls back to a single position. */
const MCS_SAMPLE_SIZE = 8;
const MAX_MCS_ANCHOR_ATOMS = 30;
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
const anchorAtomCache = new Map<string, number>();
const canonicalCoreCache = new Map<string, string>();
const DECOMPOSE_CACHE_MAX = 10000;

/** Heavy atoms of a structure, memoised: keys recur heavily across clusters and each miss is a
 *  main-thread RDKit parse. `asQuery` reads the string as a SMARTS anchor, which carries no validity
 *  flag of its own; anything that will not parse counts zero. */
function cachedAtomCount(smiles: string, rdkit: ReturnType<typeof getRdKitModule>,
  asQuery = false): number {
  const cache = asQuery ? anchorAtomCache : heavyAtomCache;
  const cached = cache.get(smiles);
  if (cached !== undefined)
    return cached;
  let mol = null;
  let count = 0;
  try {
    mol = asQuery ? rdkit.get_qmol(smiles) : rdkit.get_mol(smiles);
    count = asQuery ? mol?.get_num_atoms() ?? 0 : (mol && mol.is_valid() ? mol.get_num_atoms() : 0);
  } catch {
    count = 0;
  } finally {
    mol?.delete();
  }
  if (cache.size >= DECOMPOSE_CACHE_MAX)
    cache.clear();
  cache.set(smiles, count);
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

/** Distinct structures, smallest first, capped at {@link MCS_SAMPLE_SIZE}: the MCS cannot exceed the
 *  smallest member, so those determine it while larger ones only widen the search. The resulting
 *  anchor still faces `anchorIsUsable` against the whole cluster. Sorting by size then string keeps
 *  the sample identical across runs. */
function mcsSample(smiles: string[], rdkit: ReturnType<typeof getRdKitModule>): string[] {
  // Unparseable input counts zero atoms, which would otherwise sort it to the front and make it the
  // member the affordability gate measures — letting a molecule nothing can read stand in for a small
  // one and admit an unbounded search.
  const ranked = [...new Set(smiles)]
    .map((s) => ({s, atoms: cachedAtomCount(s, rdkit)}))
    .filter((r) => r.atoms > 0);
  ranked.sort((a, b) => a.atoms - b.atoms || (a.s < b.s ? -1 : a.s > b.s ? 1 : 0));
  return ranked.slice(0, MCS_SAMPLE_SIZE).map((r) => r.s);
}

/**
 * A cluster's site key as an anchor substructure, or null when it yields nothing usable.
 *
 * The key's marks are cut points rather than atoms, so they come off and the ring system underneath
 * is the anchor. Leaving them in would make the anchor demand a substituent at each site and drop
 * every molecule carrying a hydrogen there. What remains is re-parsed: removing a mark can leave a
 * fragment RDKit will not read, and an anchor it cannot parse is worse than none.
 */
function siteKeyAnchor(siteKey: string): string | null {
  if (!siteKey)
    return null;
  const bare = siteKey.replace(/\(\[\d*\*(?::\d+)?\]\)/g, '').replace(/\[\d*\*(?::\d+)?\]/g, '');
  if (!bare)
    return null;
  let mol = null;
  try {
    mol = getRdKitModule().get_mol(bare);
    return mol?.is_valid() ? mol.get_smiles() : null;
  } catch {
    return null;
  } finally {
    mol?.delete();
  }
}

/** Whether a sample from {@link mcsSample} is worth attempting: at least two structures to compare,
 *  and a smallest member small enough to bound the search. */
function mcsIsAffordable(sample: string[], rdkit: ReturnType<typeof getRdKitModule>): boolean {
  return sample.length >= 2 && cachedAtomCount(sample[0], rdkit) <= MAX_MCS_ANCHOR_ATOMS;
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
 *  molecules — the leftover would become oversized R-groups that decompose badly. Takes the anchor
 *  rather than the prep, so a caller can weigh one without first writing it onto the prep. */
function anchorIsUsable(anchor: string | null, smiles: string[],
  rdkit: ReturnType<typeof getRdKitModule>): boolean {
  if (!anchor)
    return false;
  const anchorHeavyAtoms = cachedAtomCount(anchor, rdkit, true);
  if (anchorHeavyAtoms < MIN_ANCHOR_HEAVY_ATOMS)
    return false;
  const coverages = smiles
    .map((s) => cachedAtomCount(s, rdkit))
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
 * by construction. Anchor priority per cluster: the site key the grouping already derived, else the
 * MCS of the whole molecules. Returns one entry per input cluster, or null where no useful clean
 * anchor exists (caller falls back to the single-position construction). Batched through a
 * worker-per-cluster queue so independent clusters run in parallel.
 *
 * `siteKeys` keeps the MCS a last resort: a cluster grouped by site reduces to that remainder by
 * construction, so its scaffold is already known exactly. The key still has to clear the same
 * usability bar an MCS anchor faces, so a key too small to anchor falls through to the MCS.
 */
export async function decomposeClusters(clusterMembers: number[][], molecules: string[],
  useMcsAnchors: boolean = true,
  siteKeys: string[] = []): Promise<(ClusterDecomposition | null)[]> {
  const results = new Array<ClusterDecomposition | null>(clusterMembers.length).fill(null);
  const rdkit = getRdKitModule();
  const service = await getRdKitService();

  // Phase 1 — pick each cluster's anchor source, queueing the ones that need an MCS.
  const preps: ClusterPrep[] = [];
  const wholeMolMcs: {input: string[], prep: ClusterPrep}[] = [];
  for (let i = 0; i < clusterMembers.length; i++) {
    const molIdx = clusterMembers[i];
    if (molIdx.length < 2 || molIdx.length > MAX_SAR_CLUSTER_SIZE)
      continue;
    const prep: ClusterPrep = {
      idx: i, molIdx, smiles: molIdx.map((k) => molecules[k]), anchor: null, anchorIsQuery: false,
    };
    preps.push(prep);
    const siteAnchor = siteKeyAnchor(siteKeys[i] ?? '');
    // Exact shared scaffold, already derived by the grouping. Measured against the same bar the MCS
    // anchors face and only taken when it clears it: a key can be a small fraction of the molecules
    // it keys, and a cluster the MCS could still anchor should reach the MCS.
    if (siteAnchor && anchorIsUsable(siteAnchor, prep.smiles, rdkit)) {
      prep.anchor = siteAnchor;
      continue;
    }
    if (!useMcsAnchors)
      continue; // no anchor without an MCS; the assembler builds this cluster from a single position
    // Left without an anchor when the cluster is past the ceiling; the assembler then builds it from
    // a single position rather than a shared core.
    const sample = mcsSample(prep.smiles, rdkit);
    if (mcsIsAffordable(sample, rdkit)) {
      prep.anchorIsQuery = true;
      wholeMolMcs.push({input: sample, prep});
    }
  }

  // Phase 2 — MCS pass. Atoms compare by element (2nd arg), not generically: interchangeable atoms
  // are what make this search explode, and by element a benzene anchors a benzene, not any 6-ring.
  // `rawSmarts` (last arg) keeps RDKit's SMARTS, which a SMILES round-trip would not survive intact.
  if (wholeMolMcs.length > 0) {
    const t = performance.now();
    const mcs = await service.clusterMCS(wholeMolMcs.map((q) => q.input), true, true, true);
    logSarTime(`MCS anchors over whole molecules (${wholeMolMcs.length} clusters)`, t);
    wholeMolMcs.forEach((q, j) => q.prep.anchor = mcs[j] || null);
  }

  // Phase 3 — validate anchors on the main thread, then batch-decompose the survivors.
  const decomposable = preps.filter((prep) => anchorIsUsable(prep.anchor, prep.smiles, rdkit));
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
