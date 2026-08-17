// import * as grok from 'datagrok-api/grok'; // only used by the disabled Murcko call (see below)
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {getRdKitService} from '../../utils/chem-common-rdkit';
import {getMmpFrags} from '../molecular-matched-pairs/mmp-analysis/mmpa-fragments';
import {SCALING_METHODS} from '../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {scaleActivity} from '../molecular-matched-pairs/mmp-viewer/mmpa-utils';
import {assembleMultiPositionMatrix} from './sar-matrix-assemble';
import {decomposeClusters} from './sar-matrix-decompose';
import {computeMatrixConfidence} from './sar-matrix-confidence';
import {buildMatchedSeries, buildCoarserLevels, clusterRelatedCores, groupSeriesBySite}
  from './sar-matrix-clustering';
import {rankMatrices, SarRankScheme} from './sar-matrix-ranking';
import {logSarTime, SarMatrix, SarMatrixCell} from './sar-matrix-types';

/**
 * Link each virtual analog's row core to its column substituent in one batched worker call.
 *
 * Single-cut matrices only (empty `refValues`): `mmpLinkFragments` substitutes `[*:1]` alone, so a
 * multi-attachment core would return unresolved `[*:2]…` dummies. Decomposed matrices are linked during
 * assembly instead. `positions` holds only the VARYING positions, so it cannot be the discriminator.
 */
async function linkVirtualStructures(matrices: SarMatrix[]): Promise<void> {
  const cores: string[] = [];
  const fragments: string[] = [];
  const targets: SarMatrixCell[] = [];

  for (const matrix of matrices) {
    if (Object.keys(matrix.refValues).length > 0)
      continue;
    for (let ri = 0; ri < matrix.rows.length; ri++) {
      for (let ci = 0; ci < matrix.columns.length; ci++) {
        const cell = matrix.cells[ri][ci];
        if (cell.kind === 'virtual' && cell.smiles === null) {
          cores.push(matrix.rows[ri].coreSmiles);
          fragments.push(matrix.columns[ci].substSmiles);
          targets.push(cell);
        }
      }
    }
  }
  if (targets.length === 0)
    return;

  const linked = await (await getRdKitService()).mmpLinkFragments(cores, fragments);
  for (let i = 0; i < targets.length; i++)
    targets[i].smiles = linked[i] ? linked[i] : null;
}

/** Row cores, column substituents and observed compounds of an assembled matrix, as comparable sets. */
function matrixContent(matrix: SarMatrix): {rows: Set<string>, cols: Set<string>, mols: Set<number>} {
  const mols = new Set<number>();
  for (const row of matrix.cells) {
    for (const cell of row) {
      if (cell.kind === 'real' && cell.molIdx !== null)
        mols.add(cell.molIdx);
    }
  }
  return {
    rows: new Set(matrix.rows.map((r) => r.keySmiles)),
    cols: new Set(matrix.columns.map((c) => `${c.position}${c.substSmiles}`)),
    mols,
  };
}

/**
 * Drop matrices whose rows, columns and compounds are all subsets of a larger one. All three must
 * match: a different core/substituent split over the same compounds is a distinct view, not a crop.
 * Runs on assembled matrices, since two differently-keyed groups can assemble into the same table.
 */
function dropCroppedMatrices(matrices: SarMatrix[]): SarMatrix[] {
  const content = matrices.map(matrixContent);
  const subset = <T>(a: Set<T>, b: Set<T>): boolean => a.size <= b.size && [...a].every((x) => b.has(x));
  return matrices.filter((_matrix, i) => !matrices.some((_other, j) => {
    if (i === j)
      return false;
    const a = content[i];
    const b = content[j];
    if (matrices[i].level !== matrices[j].level) {
      // A coarser matrix containing finer ones is the hierarchy, not redundancy; drop it only when it
      // restates a finer matrix exactly, keeping the finer one.
      const identical = a.rows.size === b.rows.size && a.cols.size === b.cols.size &&
        a.mols.size === b.mols.size;
      if (!identical || matrices[i].level < matrices[j].level)
        return false;
    }
    const smaller = a.rows.size < b.rows.size || a.cols.size < b.cols.size || a.mols.size < b.mols.size;
    // On identical content keep the earlier index, so the two cannot eliminate each other.
    if (!smaller && !(a.rows.size === b.rows.size && a.cols.size === b.cols.size && j < i))
      return false;
    return subset(a.rows, b.rows) && subset(a.cols, b.cols) && subset(a.mols, b.mols);
  }));
}

/** "Series A", ..., "Series Z", then "Series 27", "Series 28", ... */
function seriesLabel(i: number): string {
  return i < 26 ? `Series ${String.fromCharCode(65 + i)}` : `Series ${i + 1}`;
}

/** How series become the rows of one matrix. */
export enum SarGrouping {
  /** Cut the cores a second time and group the ones that leave the same remainder. */
  Site = 'Site',
  /** Cluster the cores by fingerprint similarity. */
  Similarity = 'Similarity',
}

/** Ceiling on {@link SarMatrixParams.fragmentationLevels}: each level is another fragmentation pass,
 *  and deep tiers mostly restate the ones below (a key reduced to a bare ring cannot strip further). */
export const MAX_SERIES_LEVELS = 5;

export interface SarMatrixParams {
  scaling: SCALING_METHODS;
  fragmentCutoff: number;
  predictVirtual: boolean;
  grouping: SarGrouping;
  /** Tiers of matrices, counted as the UI numbers them (3 => navigator shows L1, L2, L3). Each tier
   *  past 1 relates whole matrices by cutting their keys again, giving fewer and broader ones. */
  fragmentationLevels: number;
  /** Whether higher scaled activity is more potent. Passed in rather than re-derived from `scaling`
   *  so ranking and the cards' best/mean text agree with the direction the user set. */
  higherIsBetter: boolean;
  /** Core-similarity clustering threshold (lower groups more distant cores). Similarity grouping only. */
  threshold: number;
  rankScheme: SarRankScheme;
}

/**
 * Bemis-Murcko scaffolds — DISABLED. The server-side `Chem:MurckoScaffolds` script runs in an
 * on-demand container that costs ~30 s per call and effectively never resolves on large sets, stalling
 * the whole pipeline. It only sharpened the decomposition anchor; assembly falls back to whole-molecule
 * MCS without it. Returns [] so that fallback always runs. To re-enable, restore the body below and the
 * `grok` import (and a timeout so it can never block again).
 */
async function computeMurckoScaffolds(_molList: string[]): Promise<string[]> {
  return [];
  /*
  const MURCKO_TIMEOUT_MS = 10000;
  const compute = (async (): Promise<string[]> => {
    // MurckoScaffolds parses with MolFromSmiles, so convert molblock input to SMILES first.
    const smilesList = await (await getRdKitService()).convertMolNotation(_molList, DG.chem.Notation.Smiles);
    const tmp = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', smilesList)]);
    tmp.col('smiles')!.semType = DG.SEMTYPE.MOLECULE;
    await grok.functions.call('Chem:MurckoScaffolds', {data: tmp, smiles: 'smiles'});
    const scaffolds = tmp.col('scaffolds');
    return scaffolds ? scaffolds.toList().map((s) => s ?? '') : [];
  })().catch(() => [] as string[]);
  const timeout = new Promise<string[]>((resolve) => setTimeout(() => resolve([]), MURCKO_TIMEOUT_MS));
  return Promise.race([compute, timeout]);
  */
}

/**
 * Run the full SAR Matrix pipeline: fragment molecules, build matched series, group related cores,
 * assemble a potency matrix per group (predicting virtual analogs), and return them ranked.
 */
export async function runSarMatrix(molecules: DG.Column, activity: DG.Column<number>,
  params: SarMatrixParams): Promise<SarMatrix[]> {
  const tTotal = performance.now();
  const molList = molecules.toList();
  const scaledCol = params.scaling === SCALING_METHODS.NONE ? activity : scaleActivity(activity, params.scaling);
  // Map missing activities to NaN so assemblers skip them: scaleActivity passes the null sentinel
  // through unchanged, and read as a number it would poison the Free-Wilson fit. Bound the scan by
  // row count, not buffer length — column storage has spare capacity past the last row.
  const scaled = scaledCol.getRawData();
  const activities = new Float32Array(activity.length);
  for (let i = 0; i < activities.length; i++)
    activities[i] = activity.isNone(i) ? NaN : scaled[i];

  const tFrag = performance.now();
  const logAnd = <T>(stage: string) => (r: T): T => {
    logSarTime(stage, tFrag);
    return r;
  };
  const [[frags], scaffolds] = await Promise.all([
    getMmpFrags(molList).then(logAnd(`MMP fragmentation (${molList.length} molecules)`)),
    computeMurckoScaffolds(molList).then(logAnd('Murcko scaffolds')),
  ]);
  let t = performance.now();
  const series = buildMatchedSeries(frags, params.fragmentCutoff);
  logSarTime(`matched series (${series.length} series)`, t);
  t = performance.now();
  const base = params.grouping === SarGrouping.Site ?
    await groupSeriesBySite(series) :
    await clusterRelatedCores(series, params.threshold);
  logSarTime(`grouping by ${params.grouping.toLowerCase()} (${base.length} groups)`, t);
  // Clamped, not trusted: a programmatic caller can pass any number and each level costs a pass.
  const levels = Math.min(MAX_SERIES_LEVELS, Math.max(1, params.fragmentationLevels));
  t = performance.now();
  const clusters = await buildCoarserLevels(base, levels - 1);
  logSarTime(`coarser levels (${clusters.length} clusters)`, t);

  // Decompose all clusters in one batched pass before assembly so they run parallel across workers;
  // per-cluster calls under the Promise.all below would serialize on one worker. A null decomposition
  // falls back to single-position construction inside the assembler.
  const clusterMembers = clusters.map((c) => [...new Set(c.series.flatMap((s) => s.members.map((m) => m.molIdx)))]);
  t = performance.now();
  const decomps = await decomposeClusters(clusterMembers, molList, scaffolds);
  logSarTime(`decomposition total (${decomps.filter(Boolean).length}/${clusters.length} clusters decomposed)`, t);

  // A matrix needs >=2 rows to compare cores; folding can collapse several series onto one row, so
  // this is checked after assembly.
  t = performance.now();
  const assembled = (await Promise.all(clusters
    .map((cluster, i) => assembleMultiPositionMatrix(cluster, molList, activities, params.predictVirtual, decomps[i]))))
    .filter((matrix) => matrix.realCount > 0 && matrix.columns.length > 0 && matrix.rows.length >= 2);
  logSarTime(`assembly (${assembled.length} matrices)`, t);
  const matrices = dropCroppedMatrices(assembled);
  if (matrices.length < assembled.length)
    _package.logger.debug(`SAR Matrix: dropped ${assembled.length - matrices.length} cropped matrices`);

  // Re-link each survivor to its nearest surviving ancestor so dropping a matrix does not orphan the
  // ones below it. The chain is read from clusters, which hold links through groups that never assembled.
  const clusterParent = new Map(clusters.map((cluster) => [cluster.id, cluster.parentId]));
  const survived = new Set(matrices.map((matrix) => matrix.id));
  for (const matrix of matrices) {
    let parent = clusterParent.get(matrix.id);
    while (parent !== undefined && !survived.has(parent))
      parent = clusterParent.get(parent);
    matrix.parentId = parent;
  }

  // Labels assigned before ranking so re-ranking never renames a matrix.
  matrices.forEach((matrix, i) => matrix.label = seriesLabel(i));
  t = performance.now();
  matrices.forEach((matrix) => matrix.confidence = computeMatrixConfidence(matrix));
  logSarTime('confidence (leave-one-out)', t);

  if (params.predictVirtual) {
    t = performance.now();
    await linkVirtualStructures(matrices);
    logSarTime('virtual structure linking', t);
  }

  const ranked = rankMatrices(matrices, params.rankScheme, params.higherIsBetter);
  logSarTime(`TOTAL (${ranked.length} matrices)`, tTotal);
  return ranked;
}
