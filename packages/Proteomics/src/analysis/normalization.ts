import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';

import {SEMTYPE} from '../utils/proteomics-types';
import {unpivotIntensities} from '../viewers/qc-computations';

/** Median-centers each specified column in-place.
 * Subtracts the column median from every non-null value.
 * Uses getRawData() for bulk array access instead of per-element get/set. */
export function medianNormalize(df: DG.DataFrame, colNames: string[]): void {
  for (const name of colNames) {
    const col = df.col(name);
    if (!col) continue;
    if (col.type !== DG.COLUMN_TYPE.FLOAT) {
      console.warn(`medianNormalize: skipping non-float column "${name}" (type=${col.type})`);
      continue;
    }

    const median = col.stats.med;
    if (isNaN(median)) continue;

    const raw = col.getRawData() as Float32Array | Float64Array;
    for (let i = 0; i < raw.length; i++) {
      if (!col.isNone(i))
        raw[i] -= median;
    }
  }
  df.fireValuesChanged();
  df.setTag('proteomics.normalized', 'true');
}

/** Quantile-normalizes the specified columns in-place.
 * Aligns all column distributions to the same shape by replacing values
 * with rank-based averages. Missing values are preserved.
 * Uses getRawData() for bulk array access in tight loops. */
export function quantileNormalize(df: DG.DataFrame, colNames: string[]): void {
  if (colNames.length < 2) return;

  const cols = colNames.map((n) => df.col(n))
    .filter((c) => c != null && c.type === DG.COLUMN_TYPE.FLOAT) as DG.Column[];
  if (cols.length < 2) return;

  const nRows = df.rowCount;

  // For each column, collect indices of non-missing values and sort by value
  const colData: {indices: number[]; raw: Float32Array | Float64Array}[] = [];
  for (const col of cols) {
    const raw = col.getRawData() as Float32Array | Float64Array;
    const indices: number[] = [];
    for (let i = 0; i < nRows; i++) {
      if (!col.isNone(i))
        indices.push(i);
    }
    // Sort indices by value (ascending)
    indices.sort((a, b) => raw[a] - raw[b]);
    colData.push({indices, raw});
  }

  // Find the maximum number of non-null values across columns
  const maxValid = Math.max(...colData.map((d) => d.indices.length));
  // Need at least 2 values in the most-populated column: the rank interpolation
  // below divides by (maxValid - 1), and quantile normalization is undefined for
  // a single quantile point.
  if (maxValid < 2) return;

  // Compute rank means using scaled rank alignment
  // For each rank position r (0..maxValid-1), compute the mean of values
  // at the corresponding fractional position in each column
  const rankMeans = new Float64Array(maxValid);
  for (let r = 0; r < maxValid; r++) {
    let sum = 0;
    let count = 0;
    for (const d of colData) {
      if (d.indices.length < 2) continue;
      // Map rank r in [0..maxValid-1] to position in this column's sorted values
      const pos = (r / (maxValid - 1)) * (d.indices.length - 1);
      // Interpolate between floor and ceil positions
      const lo = Math.floor(pos);
      const hi = Math.min(Math.ceil(pos), d.indices.length - 1);
      const frac = pos - lo;
      const val = (1 - frac) * d.raw[d.indices[lo]] + frac * d.raw[d.indices[hi]];
      sum += val;
      count++;
    }
    rankMeans[r] = sum / count;
  }

  // Replace each column's sorted values with the corresponding rank mean
  for (const d of colData) {
    const n = d.indices.length;
    // n < 2: a single-value column has no rank spread; (r / (n - 1)) would be
    // 0/0 = NaN. Leave its lone value unchanged rather than corrupting it.
    if (n < 2) continue;
    for (let r = 0; r < n; r++) {
      // Map this column's rank to the global rank mean position
      const globalPos = (r / (n - 1)) * (maxValid - 1);
      const lo = Math.floor(globalPos);
      const hi = Math.min(Math.ceil(globalPos), maxValid - 1);
      const frac = globalPos - lo;
      const val = (1 - frac) * rankMeans[lo] + frac * rankMeans[hi];
      d.raw[d.indices[r]] = val;
    }
  }

  df.fireValuesChanged();
  df.setTag('proteomics.normalized', 'true');
}

/** Performs Variance-Stabilizing Normalization via server-side R script.
 * Operates on raw (non-log2) intensity columns and writes results to log2 columns.
 * Falls back to quantile normalization if R environment is unavailable. */
export async function vsnNormalize(df: DG.DataFrame, colNames: string[]): Promise<void> {
  try {
    // Find raw intensity column names by stripping log2() prefix
    const rawColNames: string[] = [];
    const log2ColNames: string[] = [];
    for (const name of colNames) {
      if (name.startsWith('log2(') && name.endsWith(')')) {
        const rawName = name.slice(5, -1);
        if (df.col(rawName)) {
          rawColNames.push(rawName);
          log2ColNames.push(name);
        }
      }
    }

    if (rawColNames.length === 0) {
      // No raw columns found, fall back to quantile
      if (colNames.length < 2) {
        grok.shell.warning('Cannot normalize — VSN needs raw intensity columns and quantile fallback needs at least 2 columns');
        return;
      }
      grok.shell.warning('No raw intensity columns found for VSN — using quantile normalization');
      quantileNormalize(df, colNames);
      return;
    }

    // Build clean DataFrame with simple column names (s1, s2, ...)
    const exprDf = DG.DataFrame.create(df.rowCount);
    for (let i = 0; i < rawColNames.length; i++) {
      const src = df.col(rawColNames[i])!;
      const srcRaw = src.getRawData() as Float32Array | Float64Array;
      const dst = exprDf.columns.addNewFloat(`s${i + 1}`);
      dst.init((r) => srcRaw[r]);
    }

    const result: DG.DataFrame = await grok.functions.call('Proteomics:VsnNormalize', {exprDf});

    // Copy VSN output (glog2 scale) into log2 intensity columns
    for (let i = 0; i < log2ColNames.length; i++) {
      const dst = df.col(log2ColNames[i])!;
      const src = result.col(`s${i + 1}`)!;
      const dstRaw = dst.getRawData() as Float32Array | Float64Array;
      const srcRaw = src.getRawData() as Float32Array | Float64Array;
      for (let r = 0; r < df.rowCount; r++) {
        if (srcRaw[r] !== DG.FLOAT_NULL)
          dstRaw[r] = srcRaw[r];
      }
    }

    df.fireValuesChanged();
    df.setTag('proteomics.normalized', 'true');
  } catch (e: any) {
    console.warn('VSN normalization failed, using quantile fallback:', e);
    if (colNames.length < 2) {
      grok.shell.warning('VSN failed and quantile fallback needs at least 2 columns — data left unnormalized');
      return;
    }
    grok.shell.warning('R environment unavailable — using quantile normalization');
    quantileNormalize(df, colNames);
  }
}

/** Shows a dialog for normalizing intensity columns.
 * Supports Median Centering, Quantile, and VSN methods with reactive box plot preview.
 * Displays an inline warning when Spectronaut pre-normalized data is detected. */
export function showNormalizationDialog(df: DG.DataFrame): void {
  if (df.getTag('proteomics.normalized') === 'true') {
    grok.shell.warning('Data already normalized');
    return;
  }

  const log2ColNames = df.columns.toList()
    .filter((c) => c.semType === SEMTYPE.INTENSITY && c.name.startsWith('log2('))
    .map((c) => c.name);

  const log2Cols = log2ColNames.map((n) => df.col(n)!);

  // Pre-normalized warning banner (NORM-04)
  const warningDiv = ui.div([], {style: {
    background: '#FFF3CD', border: '1px solid #FFEEBA', color: '#856404',
    padding: '8px', margin: '4px 0', borderRadius: '4px',
    display: df.getTag('proteomics.preNormalized') === 'true' ? '' : 'none',
  }});
  warningDiv.textContent =
    'This data may be pre-normalized (Spectronaut). Additional normalization may distort results.';

  // Method selector (NORM-03)
  const methodInput = ui.input.choice('Method', {
    value: 'Median Centering',
    items: ['Median Centering', 'Quantile', 'VSN'],
    nullable: false,
  });

  // Column picker
  const colsInput = ui.input.columns('Intensity columns', {
    table: df, available: log2ColNames, value: log2Cols,
  });

  // Box plot preview container
  const plotContainer = ui.div([], {style: {
    width: '100%', height: '250px', border: '1px solid #ddd', marginTop: '8px',
  }});

  function updatePreview(): void {
    plotContainer.innerHTML = '';
    const selected = colsInput.value.map((c: DG.Column) => c.name);
    if (selected.length === 0) return;

    // Clone DataFrame and apply selected normalization for preview
    const clone = df.clone(null, selected);
    const method = methodInput.value;
    if (method === 'Median Centering')
      medianNormalize(clone, selected);
    else if (method === 'Quantile')
      quantileNormalize(clone, selected);
    // VSN requires async R call -- show un-normalized distributions as preview

    const longDf = unpivotIntensities(clone, selected);
    const boxPlot = DG.Viewer.boxPlot(longDf, {
      valueColumnName: 'Intensity', categoryColumnName: 'Sample',
    } as any);
    boxPlot.root.style.width = '100%';
    boxPlot.root.style.height = '220px';
    plotContainer.appendChild(boxPlot.root);
  }

  methodInput.onChanged.subscribe(() => updatePreview());
  colsInput.onChanged.subscribe(() => updatePreview());
  updatePreview();

  ui.dialog('Normalize')
    .add(warningDiv)
    .add(methodInput)
    .add(colsInput)
    .add(plotContainer)
    .onOK(async () => {
      const selected = colsInput.value.map((c: DG.Column) => c.name);
      const method = methodInput.value as string;
      if (method === 'Median Centering')
        medianNormalize(df, selected);
      else if (method === 'Quantile')
        quantileNormalize(df, selected);
      else if (method === 'VSN')
        await vsnNormalize(df, selected);
      grok.shell.info(`Normalized ${selected.length} columns (${method})`);
    })
    .show();
}
