import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {getRdKitService} from '../../utils/chem-common-rdkit';
import {getMmpFrags} from '../molecular-matched-pairs/mmp-analysis/mmpa-fragments';
import {SCALING_METHODS} from '../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {scaleActivity} from '../molecular-matched-pairs/mmp-viewer/mmpa-utils';
import {assembleMultiPositionMatrix} from './sar-matrix-assemble';
import {computeMatrixConfidence} from './sar-matrix-confidence';
import {buildMatchedSeries, buildCoarserLevels, clusterRelatedCores, groupSeriesBySite}
  from './sar-matrix-clustering';
import {rankMatrices, SarRankScheme} from './sar-matrix-ranking';
import {SarMatrix, SarMatrixCell} from './sar-matrix-types';

/**
 * Build the structure of every virtual analog by linking its row core to its column
 * substituent, reusing the same `mmpLinkFragments` call the MMP Generations tab uses.
 * All matrices are linked in one batched worker call.
 *
 * Single-cut matrices only: `mmpLinkFragments` substitutes `[*:1]` alone, so a core carrying more than
 * one attachment point would come back with unresolved `[*:2]…` dummies — a structurally wrong molecule
 * that then flows into the export and the make-list. Decomposed matrices are linked during assembly
 * (`linkVirtualCellStructures`, multi-attachment); a cell still null there failed that linking and must
 * stay null, which every consumer already handles.
 *
 * The discriminator is `refValues`, not `positions`: `positions` holds only the VARYING positions, so a
 * core decomposed at R1/R2 that varies at R1 alone still carries a `[*:2]`. Only the single-cut
 * assembler leaves `refValues` empty.
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
 * Drop matrices that are a crop of another — every row, every column and every compound already
 * present in a larger one. Such a matrix shows a corner of a table the user can already see, so it
 * costs a navigator entry and buys nothing.
 *
 * Only all three together justify dropping. A matrix over the same compounds that splits them
 * differently between core and substituent introduces rows or columns the larger one lacks, and that
 * is a different view of the chemistry rather than a duplicate of it.
 *
 * The check runs on assembled matrices, not on the groups they came from: decomposition means two
 * groups keyed differently can still assemble into the same table.
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
      // A coarser matrix legitimately contains the finer ones — that is the hierarchy, not redundancy.
      // It earns its place only by being broader, so one that restates a finer matrix exactly is
      // dropped in favour of it; the finer one is the real matrix and the coarser only renames it.
      const identical = a.rows.size === b.rows.size && a.cols.size === b.cols.size &&
        a.mols.size === b.mols.size;
      if (!identical || matrices[i].level < matrices[j].level)
        return false;
    }
    const smaller = a.rows.size < b.rows.size || a.cols.size < b.cols.size || a.mols.size < b.mols.size;
    // Identical content: the earlier index is kept, so the two cannot eliminate each other.
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

/** Ceiling on {@link SarMatrixParams.fragmentationLevels}. Every level is another fragmentation pass
 *  over every key, and the levels stop earning their cost quickly: a key reduced to a bare ring has
 *  nothing left to strip, so the deeper tiers increasingly just restate the ones below them. */
export const MAX_SERIES_LEVELS = 5;

export interface SarMatrixParams {
  scaling: SCALING_METHODS;
  fragmentCutoff: number;
  predictVirtual: boolean;
  grouping: SarGrouping;
  /** How many tiers of matrices to build. 1 gives the matrices themselves — one per site, their rows
   *  being the series. Each tier past that relates whole matrices to each other by cutting their keys
   *  again, giving fewer and broader ones the higher it goes. Counted the way the UI numbers them, so
   *  3 here means the navigator shows L1, L2 and L3. */
  fragmentationLevels: number;
  /** Whether a higher scaled activity is the more potent one. Passed in rather than re-derived from
   *  `scaling`, because the caller also honours an explicit direction the user set: deriving it twice
   *  lets the first ranking sort by one rule while every card's best/mean text reports the other. */
  higherIsBetter: boolean;
  /** Core-similarity clustering threshold (lower groups more distant cores). Similarity grouping only. */
  threshold: number;
  rankScheme: SarRankScheme;
}

/**
 * Bemis-Murcko scaffold per molecule, used to form scaffold-homogeneous clusters. Runs the
 * server-side `Chem:MurckoScaffolds` RDKit script on a throwaway dataframe so the user's table is
 * not modified. Returns [] when the compute environment is unavailable, so clustering falls back
 * to the single-cut core.
 */
async function computeMurckoScaffolds(molList: string[]): Promise<string[]> {
  try {
    // The Python script parses with MolFromSmiles, so molblock input (e.g. an SDF-derived table)
    // yields empty scaffolds. Convert to SMILES in the worker first.
    const smilesList = await (await getRdKitService()).convertMolNotation(molList, DG.chem.Notation.Smiles);
    const tmp = DG.DataFrame.fromColumns([DG.Column.fromStrings('smiles', smilesList)]);
    tmp.col('smiles')!.semType = DG.SEMTYPE.MOLECULE;
    await grok.functions.call('Chem:MurckoScaffolds', {data: tmp, smiles: 'smiles'});
    const scaffolds = tmp.col('scaffolds');
    return scaffolds ? scaffolds.toList().map((s) => s ?? '') : [];
  } catch (e) {
    _package.logger.warning(`SAR Matrix: Murcko scaffolds unavailable, clustering by single-cut core: ${e}`);
    return [];
  }
}

/**
 * Run the full SAR Matrix pipeline: fragment the molecules, build matched series,
 * group related cores (by site or by similarity, per `params.grouping`), assemble a potency matrix
 * per group (predicting virtual analogs), and return the matrices ranked by the chosen scheme.
 *
 * @param molecules Molecule column.
 * @param activity Numerical activity column.
 * @param params Analysis parameters.
 */
export async function runSarMatrix(molecules: DG.Column, activity: DG.Column<number>,
  params: SarMatrixParams): Promise<SarMatrix[]> {
  const molList = molecules.toList();
  const scaledCol = params.scaling === SCALING_METHODS.NONE ? activity : scaleActivity(activity, params.scaling);
  // Missing activities must not become "real" cells: scaleActivity passes the null sentinel through
  // unchanged, and a sentinel read as a number would render as an extreme potency and poison the
  // Free-Wilson fit. Map every missing value to NaN so the assemblers can skip it.
  // The scan is bounded by the row count, never by the raw buffer's length: column storage is
  // allocated with spare capacity, so the buffer runs past the last row and probing a row that does
  // not exist fails the platform's own bounds check.
  const scaled = scaledCol.getRawData();
  const activities = new Float32Array(activity.length);
  for (let i = 0; i < activities.length; i++)
    activities[i] = activity.isNone(i) ? NaN : scaled[i];

  const [[frags], scaffolds] = await Promise.all([getMmpFrags(molList), computeMurckoScaffolds(molList)]);
  const series = buildMatchedSeries(frags, params.fragmentCutoff);
  const base = params.grouping === SarGrouping.Site ?
    await groupSeriesBySite(series) :
    await clusterRelatedCores(series, params.threshold);
  // Level 1 built the series and level 2 the matrices, so anything past 2 adds a coarser matrix above
  // the ones already there, holding their compounds over rows that agree a further cut deeper.
  // Clamped rather than trusted: the inputs that set this are bounded, but a programmatic caller can
  // pass any number, and each extra level costs another fragmentation pass over every key.
  const levels = Math.min(MAX_SERIES_LEVELS, Math.max(1, params.fragmentationLevels));
  const clusters = await buildCoarserLevels(base, levels - 1);

  // A group needs two rows to be a matrix at all: comparing cores is the whole point, and a single
  // row carries no core comparison and no row effect for the additive model to use. Folding can
  // collapse a group of several series onto one row, so this is checked after assembly, not before.
  const assembled = (await Promise.all(clusters
    .map((cluster) => assembleMultiPositionMatrix(cluster, molList, activities, params.predictVirtual, scaffolds))))
    .filter((matrix) => matrix.realCount > 0 && matrix.columns.length > 0 && matrix.rows.length >= 2);
  const matrices = dropCroppedMatrices(assembled);
  if (matrices.length < assembled.length)
    _package.logger.debug(`SAR Matrix: dropped ${assembled.length - matrices.length} cropped matrices`);

  // Dropping a matrix must not orphan the ones below it, so each survivor is re-linked to its nearest
  // surviving ancestor. The chain is read from the clusters, which still hold links through groups
  // that never assembled at all.
  const clusterParent = new Map(clusters.map((cluster) => [cluster.id, cluster.parentId]));
  const survived = new Set(matrices.map((matrix) => matrix.id));
  for (const matrix of matrices) {
    let parent = clusterParent.get(matrix.id);
    while (parent !== undefined && !survived.has(parent))
      parent = clusterParent.get(parent);
    matrix.parentId = parent;
  }

  // Stable "Series A/B/..." labels assigned before ranking so re-ranking never renames a matrix.
  matrices.forEach((matrix, i) => matrix.label = seriesLabel(i));
  // How far to trust each matrix's virtual predictions (leave-one-out over its observed cells).
  matrices.forEach((matrix) => matrix.confidence = computeMatrixConfidence(matrix));

  if (params.predictVirtual)
    await linkVirtualStructures(matrices);

  return rankMatrices(matrices, params.rankScheme, params.higherIsBetter);
}
