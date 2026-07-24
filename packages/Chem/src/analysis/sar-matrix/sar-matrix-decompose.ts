import * as DG from 'datagrok-api/dg';

import {getMCS} from '../../utils/most-common-subs';
import {getRdKitModule} from '../../utils/chem-common-rdkit';
import {rGroupsMinilib} from '../r-group-analysis';

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

const R_GROUP_OPTIONS = {
  matchingStrategy: 'Greedy',
  includeTargetMolInResults: true,
  onlyMatchAtRGroups: false,
};

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

/**
 * The scaffold the cluster is decomposed against. Anchoring on the shared **ring scaffold**
 * (Bemis-Murcko, computed upstream) rather than the MCS of whole molecules keeps peripheral
 * variation out of the core: a varying methyl or an alkyl homolog (Me/Et/iPr) becomes an R-group
 * instead of splitting one scaffold into several near-identical core rows or truncating the R-group.
 *
 * - One shared scaffold in the cluster → use it directly (concrete substructure).
 * - Several related scaffolds (a ring hop) → their generic MCS, so both match with aligned R-labels.
 * - No scaffolds (Murcko unavailable) → fall back to the MCS of whole molecules.
 */
async function computeAnchor(molIdx: number[], molCol: DG.Column<string>, scaffolds: string[]):
  Promise<{anchor: string | null, anchorIsQuery: boolean}> {
  const clusterScaffolds = scaffolds.length ?
    [...new Set(molIdx.map((i) => scaffolds[i]).filter((s) => s))] : [];
  if (clusterScaffolds.length === 1)
    return {anchor: clusterScaffolds[0], anchorIsQuery: false};
  if (clusterScaffolds.length >= 2) {
    const scafCol = DG.DataFrame.fromColumns([DG.Column.fromStrings('scaffold', clusterScaffolds)]).col('scaffold')!;
    scafCol.semType = DG.SEMTYPE.MOLECULE;
    const mcs = await getMCS(scafCol, false, true);
    if (mcs)
      return {anchor: mcs, anchorIsQuery: true};
  }
  return {anchor: (await getMCS(molCol, false, true)) || null, anchorIsQuery: true};
}

/** Heavy-atom count of a SMILES, or 0 if it can't be parsed. */
function heavyAtomCount(rdkit: ReturnType<typeof getRdKitModule>, smiles: string): number {
  let mol = null;
  try {
    mol = rdkit.get_mol(smiles);
    return mol && mol.is_valid() ? mol.get_num_atoms() : 0;
  } catch {
    return 0;
  } finally {
    mol?.delete();
  }
}

/**
 * Computes one generic-atom-query scaffold shared by every molecule in the cluster (the same
 * `getMCS(col, exactAtomSearch=false, exactBondSearch=true)` the interactive R-Groups Analysis
 * "MCS" button uses), then decomposes the whole cluster against it in a single R-Group
 * Decomposition call — every molecule's R1, R2, ... columns come out aligned by construction.
 *
 * Returns `null` when the cluster has no useful, clean shared multi-position anchor — the anchor
 * covers too little of the molecules (they share no real core), or too many molecules decompose
 * into malformed substituents. Callers then fall back to the single-position construction, so a
 * heterogeneous cluster degrades to a simpler view instead of rendering garbage fragments.
 */
export async function decomposeCluster(molIdx: number[], molecules: string[], scaffolds: string[] = []):
  Promise<ClusterDecomposition | null> {
  if (molIdx.length < 2 || molIdx.length > MAX_SAR_CLUSTER_SIZE)
    return null;

  const smiles = molIdx.map((i) => molecules[i]);
  // rGroupsMinilib reads molecules.dataFrame.rowCount, so the column must belong to a DataFrame.
  const molDf = DG.DataFrame.fromColumns([DG.Column.fromStrings('mol', smiles)]);
  const molCol = molDf.col('mol')!;
  molCol.semType = DG.SEMTYPE.MOLECULE;

  const {anchor, anchorIsQuery} = await computeAnchor(molIdx, molCol, scaffolds);
  if (!anchor)
    return null;

  const rdkit = getRdKitModule();
  let anchorHeavyAtoms = 0;
  let anchorMol = null;
  try {
    anchorMol = rdkit.get_qmol(anchor);
    anchorHeavyAtoms = anchorMol?.get_num_atoms() ?? 0;
  } catch {
    return null;
  } finally {
    anchorMol?.delete();
  }
  if (anchorHeavyAtoms < MIN_ANCHOR_HEAVY_ATOMS)
    return null;

  // #3 — reject a scaffold that covers too little of the molecules: the leftover would become
  // oversized R-groups that decompose badly.
  const coverages = smiles
    .map((s) => heavyAtomCount(rdkit, s))
    .filter((n) => n > 0)
    .map((n) => anchorHeavyAtoms / n);
  if (median(coverages) < MIN_MCS_COVERAGE)
    return null;

  const {rGroups} = await rGroupsMinilib(molCol, anchor, anchorIsQuery, 0, R_GROUP_OPTIONS);
  // rGroups[0] is the Core column; the rest are R1, R2, ...
  if (rGroups.length < 2)
    return null;

  const coreCol = rGroups[0];
  const positions = rGroups.slice(1).map((c) => c.name);
  const positionNums = positions.map(positionNumber);
  const records: PositionRecord[] = [];
  let matched = 0;
  for (let i = 0; i < molIdx.length; i++) {
    const coreRaw = coreCol.get(i);
    if (!coreRaw)
      continue; // molecule didn't match the anchor
    matched++;

    // #1 — validate every substituent; drop the whole molecule if any is malformed, and reject a
    // multi-component core (a bad match).
    const values: {[position: string]: string} = {};
    let clean = true;
    for (let p = 0; p < positions.length; p++) {
      const value = rGroups[p + 1].get(i) ?? '';
      if (!isCleanSubstituent(value, positionNums[p])) {
        clean = false;
        break;
      }
      values[positions[p]] = value;
    }
    if (!clean || coreRaw.includes('.'))
      continue;

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
    records.push({molIdx: molIdx[i], coreSmiles, values});
  }

  // #1 — if too few of the matched molecules decomposed cleanly, the anchor is wrong for this
  // cluster: fall back rather than build a matrix from a handful of survivors.
  if (records.length < 2 || records.length < matched * MIN_CLEAN_FRACTION)
    return null;
  return {records, positions};
}
