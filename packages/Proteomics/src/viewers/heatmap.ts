import * as DG from 'datagrok-api/dg';
import {getGroups} from '../analysis/experiment-setup';
import {findColumn} from '../utils/column-detection';
import {SEMTYPE} from '../utils/proteomics-types';

/** Z-score normalize a value given mean and standard deviation. */
function zScore(value: number, mean: number, std: number): number {
  return std > 0 ? (value - mean) / std : 0;
}

/**
 * Creates an expression heatmap (Grid in heatmap mode) on a cloned DataFrame
 * containing only the top-N most significant proteins. The clone isolates
 * the heatmap's filter, z-score columns, and sort order from the original
 * DataFrame so that other viewers (volcano plot) are not affected.
 *
 * Rows are hierarchically clustered via the Dendrogram package service if available,
 * otherwise sorted by adj.p-value ascending (significance-based fallback).
 *
 * Z-score normalization uses statistics computed across ALL rows of the original
 * DataFrame for statistical correctness, but z-score columns are only added to the
 * cloned DataFrame.
 */
export async function createExpressionHeatmap(
  df: DG.DataFrame,
  options?: {topN?: number; labelCol?: string; title?: string},
): Promise<DG.Grid> {
  const topN = options?.topN ?? 50;
  const labelColName = options?.labelCol ?? 'Gene Name';

  // 1. Check prerequisites
  if (df.getTag('proteomics.de_complete') !== 'true')
    throw new Error('Differential expression must be run first');

  const groups = getGroups(df);
  if (!groups)
    throw new Error('Experimental groups must be annotated first');

  // 2. Select top N proteins by adj.p-value (ascending)
  const adjPCol = df.col('adj.p-value');
  if (!adjPCol)
    throw new Error('adj.p-value column not found');

  // Collect row indices with non-null adj.p-value, sorted ascending
  const indexedPValues: {idx: number; pVal: number}[] = [];
  for (let i = 0; i < df.rowCount; i++) {
    if (!adjPCol.isNone(i))
      indexedPValues.push({idx: i, pVal: adjPCol.get(i) as number});
  }
  indexedPValues.sort((a, b) => a.pVal - b.pVal);

  const topNIndices = new Set<number>();
  for (let i = 0; i < Math.min(topN, indexedPValues.length); i++)
    topNIndices.add(indexedPValues[i].idx);

  // 3. Clone the DataFrame with only top-N rows (filter isolation)
  const filter = DG.BitSet.create(df.rowCount);
  for (const idx of topNIndices)
    filter.set(idx, true);
  const heatmapDf = df.clone(filter);

  // 4. Compute z-score columns on the cloned DataFrame
  // Statistics (mean, std) are computed across ALL rows of the original df for correctness,
  // then z-scores are computed on the cloned rows only.
  const allIntensityCols = [...groups.group1.columns, ...groups.group2.columns];
  const zScoreColNames: string[] = [];

  for (const colName of allIntensityCols) {
    const zColName = `_zscore_${colName}`;
    zScoreColNames.push(zColName);

    const srcColOrig = df.col(colName);
    if (!srcColOrig) continue;

    // Compute mean and std across all non-null values in the ORIGINAL column.
    // Use direct typed-array access to skip the per-row JS->native interop.
    const origRaw = srcColOrig.getRawData() as Float32Array | Float64Array;
    let sum = 0;
    let sumSq = 0;
    let count = 0;
    for (let i = 0; i < df.rowCount; i++) {
      const v = origRaw[i];
      if (v !== DG.FLOAT_NULL) {
        sum += v;
        sumSq += v * v;
        count++;
      }
    }

    const colMean = count > 0 ? sum / count : 0;
    const variance = count > 1 ? (sumSq - sum * sum / count) / (count - 1) : 0;
    const colStd = Math.sqrt(variance);

    // Apply z-scores to the cloned DataFrame
    const srcColClone = heatmapDf.col(colName);
    if (!srcColClone) continue;

    const srcRaw = srcColClone.getRawData() as Float32Array | Float64Array;
    const zCol = heatmapDf.columns.addNewFloat(zColName);
    zCol.init((i) => srcRaw[i] === DG.FLOAT_NULL ? DG.FLOAT_NULL : zScore(srcRaw[i], colMean, colStd));
  }

  // 5. Create Grid on the cloned DataFrame
  const grid = heatmapDf.plot.grid();

  // 6. Column visibility -- hide everything, then show label + z-score columns
  for (let i = 0; i < grid.columns.length; i++) {
    const gc = grid.columns.byIndex(i);
    if (gc) gc.visible = false;
  }

  // Show label column (gene symbol preferred, protein ID as fallback)
  const labelCol = findColumn(heatmapDf, SEMTYPE.GENE_SYMBOL, [labelColName.toLowerCase(), 'gene symbol'])
    ?? findColumn(heatmapDf, SEMTYPE.PROTEIN_ID, ['protein id', 'majority protein id', 'accession']);
  if (labelCol) {
    const labelGc = grid.columns.byName(labelCol.name);
    if (labelGc) labelGc.visible = true;
  }

  // Show z-score columns in group order (control left, treatment right)
  const group1ZCols = groups.group1.columns.map((c) => `_zscore_${c}`);
  const group2ZCols = groups.group2.columns.map((c) => `_zscore_${c}`);

  for (const zColName of [...group1ZCols, ...group2ZCols]) {
    const gc = grid.columns.byName(zColName);
    if (gc) gc.visible = true;
  }

  // 7. Attempt hierarchical clustering via Dendrogram package; fall back to
  // significance-based row ordering if the package is unavailable or its
  // function signature has drifted.
  let clustered = false;
  const pi = DG.TaskBarProgressIndicator.create('Calculating distance matrix...');
  try {
    const dFunc = DG.Func.find({package: 'Dendrogram', name: 'hierarchicalClustering'})[0];
    if (!dFunc)
      throw new Error('Dendrogram package not installed (hierarchicalClustering not found)');
    if (dFunc.inputs.length !== 4)
      throw new Error(`Dendrogram.hierarchicalClustering signature changed (expected 4 inputs, got ${dFunc.inputs.length})`);
    await dFunc.apply({
      df: grid, colNameList: zScoreColNames,
      distance: 'euclidean', linkage: 'complete',
    });
    clustered = true;
  } catch (e: any) {
    console.warn(`Dendrogram clustering unavailable, falling back to significance-based row ordering: ${e?.message ?? e}`);
  } finally {
    pi.close();
  }

  if (!clustered) {
    // Sort by adjusted p-value ascending (most significant first). R scripts
    // sometimes return the column as `adj.p.value` (dots) instead of `adj.p-value`.
    const sortCol = heatmapDf.col('adj.p-value') ? 'adj.p-value' :
      (heatmapDf.col('adj.p.value') ? 'adj.p.value' : null);
    if (sortCol)
      grid.sort([sortCol]);
  }

  // 8. Configure Grid for heatmap display
  grid.props.isHeatmap = true;
  grid.props.isGrid = false;

  if (options?.title)
    grid.setOptions({title: options.title});

  return grid;
}
