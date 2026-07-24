import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';

import {_package} from '../../package';
import {getRdKitService} from '../../utils/chem-common-rdkit';
import {getMmpFrags} from '../molecular-matched-pairs/mmp-analysis/mmpa-fragments';
import {SCALING_METHODS} from '../molecular-matched-pairs/mmp-viewer/mmp-constants';
import {scaleActivity} from '../molecular-matched-pairs/mmp-viewer/mmpa-utils';
import {assembleMultiPositionMatrix} from './sar-matrix-assemble';
import {computeMatrixConfidence} from './sar-matrix-confidence';
import {buildMatchedSeries, clusterRelatedCores} from './sar-matrix-clustering';
import {rankMatrices, SarRankScheme} from './sar-matrix-ranking';
import {SarMatrix, SarMatrixCell} from './sar-matrix-types';

/**
 * Build the structure of every virtual analog by linking its row core to its column
 * substituent, reusing the same `mmpLinkFragments` call the MMP Generations tab uses.
 * All matrices are linked in one batched worker call.
 */
async function linkVirtualStructures(matrices: SarMatrix[]): Promise<void> {
  const cores: string[] = [];
  const fragments: string[] = [];
  const targets: SarMatrixCell[] = [];

  for (const matrix of matrices) {
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

/** "Series A", ..., "Series Z", then "Series 27", "Series 28", ... */
function seriesLabel(i: number): string {
  return i < 26 ? `Series ${String.fromCharCode(65 + i)}` : `Series ${i + 1}`;
}

export interface SarMatrixParams {
  scaling: SCALING_METHODS;
  fragmentCutoff: number;
  predictVirtual: boolean;
  /** Core-similarity clustering threshold (lower groups more distant cores). */
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
 * cluster related cores by Murcko scaffold, assemble a potency matrix per cluster
 * (predicting virtual analogs), and return the matrices ranked by the chosen scheme.
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
  const activities = Float32Array.from(scaledCol.getRawData());
  for (let i = 0; i < activities.length; i++) {
    if (activity.isNone(i))
      activities[i] = NaN;
  }

  const [[frags], scaffolds] = await Promise.all([getMmpFrags(molList), computeMurckoScaffolds(molList)]);
  const series = buildMatchedSeries(frags, params.fragmentCutoff);
  const clusters = await clusterRelatedCores(series, scaffolds, params.threshold);

  const matrices = (await Promise.all(clusters
    .map((cluster) => assembleMultiPositionMatrix(cluster, molList, activities, params.predictVirtual, scaffolds))))
    .filter((matrix) => matrix.realCount > 0 && matrix.columns.length > 0);

  // Stable "Series A/B/..." labels assigned before ranking so re-ranking never renames a matrix.
  matrices.forEach((matrix, i) => matrix.label = seriesLabel(i));
  // How far to trust each matrix's virtual predictions (leave-one-out over its observed cells).
  matrices.forEach((matrix) => matrix.confidence = computeMatrixConfidence(matrix));

  if (params.predictVirtual)
    await linkVirtualStructures(matrices);

  return rankMatrices(matrices, params.rankScheme, params.scaling === SCALING_METHODS.MINUS_LG);
}
