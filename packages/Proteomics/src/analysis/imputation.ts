import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SEMTYPE} from '../utils/proteomics-types';
import {getGroups} from './experiment-setup';

/** Generates a random number from a normal distribution using Box-Muller transform. */
function randomNormal(mean: number, stdev: number): number {
  const u1 = Math.random();
  const u2 = Math.random();
  const z = Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  return mean + stdev * z;
}

/** Imputes missing values using the MinProb method (Perseus-style).
 * Draws random values from a downshifted normal distribution
 * simulating low-abundance proteins below detection limit.
 * Returns the total number of imputed values. */
export function imputeMinProb(
  df: DG.DataFrame,
  colNames: string[],
  downshift: number = 1.8,
  width: number = 0.3,
): number {
  let totalImputed = 0;

  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    if (col.type !== DG.COLUMN_TYPE.FLOAT) continue;

    const avg = col.stats.avg;
    const stdev = col.stats.stdev;
    if (isNaN(avg) || isNaN(stdev) || stdev === 0) continue;

    const imputeMean = avg - downshift * stdev;
    const imputeSd = width * stdev;

    const raw = col.getRawData() as Float32Array | Float64Array;
    for (let i = 0; i < df.rowCount; i++) {
      if (raw[i] === DG.FLOAT_NULL) {
        raw[i] = randomNormal(imputeMean, imputeSd);
        totalImputed++;
      }
    }
  }

  df.fireValuesChanged();
  df.setTag('proteomics.imputed', 'true');
  return totalImputed;
}

/** Imputes missing values using k-nearest neighbor averaging.
 * For each row with missing values, finds the k nearest rows by mean-squared
 * distance (RMSE) over the columns where both rows have valid values, then
 * averages those neighbors' values. Normalizing by sharedCount keeps distances
 * comparable across rows with different numbers of shared observations.
 * Falls back to the column mean when no neighbors have values for a column;
 * if a column is entirely missing, leaves the cell null.
 * Shows a progress indicator during computation.
 * Returns the total number of imputed values. */
export function imputeKnn(
  df: DG.DataFrame,
  colNames: string[],
  k: number = 10,
): number {
  const cols = colNames.map((n) => df.col(n)).filter((c) => c != null && c.type === DG.COLUMN_TYPE.FLOAT) as DG.Column[];
  if (cols.length === 0) return 0;

  // Cache the underlying typed arrays for direct in-place writes — avoids the
  // JS→native bridge per cell in the inner imputation loops below.
  const rawArrays = cols.map((c) => c.getRawData() as Float32Array | Float64Array);

  const nRows = df.rowCount;
  const nCols = cols.length;

  // Build matrix and missing mask
  const matrix: Float64Array[] = [];
  const missing: boolean[][] = [];
  for (let r = 0; r < nRows; r++) {
    const row = new Float64Array(nCols);
    const miss: boolean[] = [];
    for (let c = 0; c < nCols; c++) {
      const isNull = cols[c].isNone(r);
      miss.push(isNull);
      row[c] = isNull ? 0 : cols[c].get(r);
    }
    matrix.push(row);
    missing.push(miss);
  }

  // Compute column means (from non-missing values) as fallback.
  // Use NaN for all-missing columns so downstream consumers can detect "no signal"
  // rather than silently filling missing cells with a literal 0 (a valid log2 intensity).
  const colMeans = new Float64Array(nCols);
  for (let c = 0; c < nCols; c++) {
    let sum = 0;
    let count = 0;
    for (let r = 0; r < nRows; r++) {
      if (!missing[r][c]) {
        sum += matrix[r][c];
        count++;
      }
    }
    colMeans[c] = count > 0 ? sum / count : NaN;
  }

  // Find rows with any missing values (regardless of whether they also have valid ones —
  // the all-missing case is handled inside the loop via column means or null).
  const rowsWithMissing: number[] = [];
  for (let r = 0; r < nRows; r++) {
    if (missing[r].some((m) => m))
      rowsWithMissing.push(r);
  }

  const pi = DG.TaskBarProgressIndicator.create('kNN imputation...');
  let totalImputed = 0;

  try {
    for (let idx = 0; idx < rowsWithMissing.length; idx++) {
      const r = rowsWithMissing[idx];
      pi.update(Math.round((idx / rowsWithMissing.length) * 100), `Row ${idx + 1}/${rowsWithMissing.length}`);

      // Check if this row has any valid values for distance computation
      const hasAnyValid = missing[r].some((m) => !m);

      if (!hasAnyValid) {
        // All-missing row: use column means, but leave null if the column itself was all-missing
        for (let c = 0; c < nCols; c++) {
          if (missing[r][c] && isFinite(colMeans[c])) {
            rawArrays[c][r] = colMeans[c];
            totalImputed++;
          }
        }
        continue;
      }

      // Compute distances to all other rows on shared non-missing columns
      const distances: {row: number; dist: number}[] = [];
      for (let other = 0; other < nRows; other++) {
        if (other === r) continue;

        let sumSq = 0;
        let sharedCount = 0;
        for (let c = 0; c < nCols; c++) {
          if (!missing[r][c] && !missing[other][c]) {
            const diff = matrix[r][c] - matrix[other][c];
            sumSq += diff * diff;
            sharedCount++;
          }
        }

        if (sharedCount > 0)
          distances.push({row: other, dist: Math.sqrt(sumSq / sharedCount)});
      }

      // Sort by distance, take top k
      distances.sort((a, b) => a.dist - b.dist);
      const neighbors = distances.slice(0, k);

      // For each missing column in this row, average neighbors' values
      for (let c = 0; c < nCols; c++) {
        if (!missing[r][c]) continue;

        let sum = 0;
        let count = 0;
        for (const n of neighbors) {
          if (!missing[n.row][c]) {
            sum += matrix[n.row][c];
            count++;
          }
        }

        const val = count > 0 ? sum / count : colMeans[c];
        if (isFinite(val)) {
          rawArrays[c][r] = val;
          totalImputed++;
        }
      }
    }
  } finally {
    pi.close();
  }

  df.fireValuesChanged();
  df.setTag('proteomics.imputed', 'true');
  return totalImputed;
}

/** Imputes missing values by replacing them with zero.
 * Returns the total number of imputed values. */
export function imputeZero(df: DG.DataFrame, colNames: string[]): number {
  let totalImputed = 0;
  for (const name of colNames) {
    const col = df.col(name);
    if (!col || col.type !== DG.COLUMN_TYPE.FLOAT) continue;
    const raw = col.getRawData() as Float32Array | Float64Array;
    for (let i = 0; i < df.rowCount; i++) {
      if (raw[i] === DG.FLOAT_NULL) {
        raw[i] = 0;
        totalImputed++;
      }
    }
  }
  df.fireValuesChanged();
  df.setTag('proteomics.imputed', 'true');
  return totalImputed;
}

/** Imputes missing values by replacing them with the column mean.
 * Returns the total number of imputed values. */
export function imputeMean(df: DG.DataFrame, colNames: string[]): number {
  let totalImputed = 0;
  for (const name of colNames) {
    const col = df.col(name);
    if (!col || col.type !== DG.COLUMN_TYPE.FLOAT) continue;
    const avg = col.stats.avg;
    if (isNaN(avg)) continue;
    const raw = col.getRawData() as Float32Array | Float64Array;
    for (let i = 0; i < df.rowCount; i++) {
      if (raw[i] === DG.FLOAT_NULL) {
        raw[i] = avg;
        totalImputed++;
      }
    }
  }
  df.fireValuesChanged();
  df.setTag('proteomics.imputed', 'true');
  return totalImputed;
}

/** Imputes missing values by replacing them with the column median.
 * Returns the total number of imputed values. */
export function imputeMedian(df: DG.DataFrame, colNames: string[]): number {
  let totalImputed = 0;
  for (const name of colNames) {
    const col = df.col(name);
    if (!col || col.type !== DG.COLUMN_TYPE.FLOAT) continue;
    const med = col.stats.med;
    if (isNaN(med)) continue;
    const raw = col.getRawData() as Float32Array | Float64Array;
    for (let i = 0; i < df.rowCount; i++) {
      if (raw[i] === DG.FLOAT_NULL) {
        raw[i] = med;
        totalImputed++;
      }
    }
  }
  df.fireValuesChanged();
  df.setTag('proteomics.imputed', 'true');
  return totalImputed;
}

/** Shows a dialog for imputing missing values with method selection,
 * conditional parameters, and valid-values protein filter. */
export function showImputationDialog(df: DG.DataFrame): void {
  if (df.getTag('proteomics.imputed') === 'true') {
    grok.shell.warning('Missing values already imputed');
    return;
  }

  const log2ColNames = df.columns.toList()
    .filter((c) => c.semType === SEMTYPE.INTENSITY && c.name.startsWith('log2('))
    .map((c) => c.name);

  const log2Cols = log2ColNames.map((n) => df.col(n)!);

  // Method selector
  const methodInput = ui.input.choice('Method', {
    value: 'MinProb',
    items: ['MinProb', 'kNN', 'Zero', 'Mean', 'Median'],
    nullable: false,
  });

  // Column picker
  const colsInput = ui.input.columns('Intensity columns', {
    table: df, available: log2ColNames, value: log2Cols,
  });

  // MinProb-specific inputs
  const downshiftInput = ui.input.float('Downshift', {value: 1.8});
  downshiftInput.setTooltip('Standard deviations to shift left from mean (Perseus default: 1.8)');
  const widthInput = ui.input.float('Width', {value: 0.3});
  widthInput.setTooltip('Fraction of standard deviation for imputation distribution (Perseus default: 0.3)');
  const minProbContainer = ui.div([downshiftInput.root, widthInput.root]);

  // kNN-specific input
  const kInput = ui.input.int('k (neighbors)', {value: 10});
  const knnContainer = ui.div([kInput.root]);
  knnContainer.style.display = 'none';

  // Toggle visibility based on method
  methodInput.onChanged.subscribe(() => {
    const method = methodInput.value;
    minProbContainer.style.display = method === 'MinProb' ? '' : 'none';
    knnContainer.style.display = method === 'kNN' ? '' : 'none';
  });

  // Valid-values filter
  const minValidInput = ui.input.int('Min valid values per group', {value: 2});
  const filterCountDiv = ui.divText('');
  filterCountDiv.style.cssText = 'font-size:12px; color:#888; margin-bottom:8px;';

  const updateFilterCount = () => {
    const threshold = minValidInput.value ?? 0;
    const groups = getGroups(df);
    if (threshold === 0 || !groups) {
      filterCountDiv.textContent = `Will keep all ${df.rowCount} proteins`;
      return;
    }

    const g1Cols = groups.group1.columns.map((n) => df.col(n)).filter((c) => c != null) as DG.Column[];
    const g2Cols = groups.group2.columns.map((n) => df.col(n)).filter((c) => c != null) as DG.Column[];
    let kept = 0;

    for (let r = 0; r < df.rowCount; r++) {
      let g1Valid = 0;
      for (const col of g1Cols) {
        if (!col.isNone(r)) g1Valid++;
      }
      let g2Valid = 0;
      for (const col of g2Cols) {
        if (!col.isNone(r)) g2Valid++;
      }
      // Protein passes if ANY group has >= threshold valid values
      if (g1Valid >= threshold || g2Valid >= threshold)
        kept++;
    }

    const removed = df.rowCount - kept;
    filterCountDiv.textContent = `Will keep ${kept}/${df.rowCount} proteins (${removed} removed)`;
  };

  minValidInput.onChanged.subscribe(updateFilterCount);
  updateFilterCount();

  ui.dialog('Impute Missing Values')
    .add(methodInput)
    .add(colsInput)
    .add(minProbContainer)
    .add(knnContainer)
    .add(minValidInput)
    .add(filterCountDiv)
    .onOK(() => {
      const selected = colsInput.value.map((c: DG.Column) => c.name);
      const method = methodInput.value;

      // Apply valid-values filter if threshold > 0 and groups exist
      const threshold = minValidInput.value ?? 0;
      const groups = getGroups(df);
      if (threshold > 0 && groups) {
        const g1Cols = groups.group1.columns.map((n) => df.col(n)).filter((c) => c != null) as DG.Column[];
        const g2Cols = groups.group2.columns.map((n) => df.col(n)).filter((c) => c != null) as DG.Column[];
        const removeSet = DG.BitSet.create(df.rowCount);

        for (let r = 0; r < df.rowCount; r++) {
          let g1Valid = 0;
          for (const col of g1Cols) {
            if (!col.isNone(r)) g1Valid++;
          }
          let g2Valid = 0;
          for (const col of g2Cols) {
            if (!col.isNone(r)) g2Valid++;
          }
          if (g1Valid < threshold && g2Valid < threshold)
            removeSet.set(r, true);
        }

        const removed = removeSet.trueCount;
        if (removed > 0) {
          df.rows.removeWhere((row) => removeSet.get(row.idx));
          grok.shell.info(`Filtered: removed ${removed} proteins below threshold`);
        }
      }

      // Dispatch to selected imputation method. Imputation runs AFTER the valid-values
      // filter above, which has already mutated df.rows — if imputation throws, the
      // filter cannot be reverted, so surface the partial state clearly.
      let count = 0;
      try {
        switch (method) {
        case 'MinProb':
          count = imputeMinProb(df, selected, downshiftInput.value!, widthInput.value!);
          break;
        case 'kNN':
          count = imputeKnn(df, selected, kInput.value!);
          break;
        case 'Zero':
          count = imputeZero(df, selected);
          break;
        case 'Mean':
          count = imputeMean(df, selected);
          break;
        case 'Median':
          count = imputeMedian(df, selected);
          break;
        }
      } catch (e: any) {
        grok.shell.error(`Imputation failed: ${e?.message ?? e}. Valid-values filter was already applied; data remains filtered but un-imputed.`);
        return;
      }

      grok.shell.info(`Imputed ${count} missing values across ${selected.length} columns (${method})`);
    })
    .show();
}
