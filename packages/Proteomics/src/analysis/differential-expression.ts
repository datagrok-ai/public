import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {tTest} from '@datagrok-libraries/statistics/src/tests';
import {fdrcorrection} from '@datagrok-libraries/statistics/src/multiple-tests';
import {SEMTYPE, DEFAULT_FC_THRESHOLD, DEFAULT_P_THRESHOLD} from '../utils/proteomics-types';
import {getGroups} from './experiment-setup';

/** Guard for menu/dialog handlers that need DE results. Returns true if DE has been
 *  run on `df`; otherwise warns the user via `grok.shell.warning` and returns false.
 *  Each caller passes its own warning text so existing user-visible copy is preserved.
 *  Use as: `if (!requireDifferentialExpression(df, '...')) return;` */
export function requireDifferentialExpression(df: DG.DataFrame, message: string): boolean {
  if (df.getTag('proteomics.de_complete') === 'true')
    return true;
  grok.shell.warning(message);
  return false;
}

/** Extract non-null float values across named columns for a single row. */
function getGroupValues(df: DG.DataFrame, colNames: string[], rowIdx: number): number[] {
  const values: number[] = [];
  for (const name of colNames) {
    const col = df.columns.byName(name);
    if (col && !col.isNone(rowIdx))
      values.push(col.get(rowIdx));
  }
  return values;
}

/** Compute arithmetic mean of a number array. */
function mean(vals: number[]): number {
  let sum = 0;
  for (let i = 0; i < vals.length; i++)
    sum += vals[i];
  return sum / vals.length;
}

/** Run client-side Welch's t-test differential expression with BH FDR correction.
 *  Compares group2 (treatment) vs group1 (control). Adds log2FC, p-value, adj.p-value,
 *  and significant columns to the DataFrame in-place. */
export function runDifferentialExpression(
  df: DG.DataFrame,
  group1Cols: string[],
  group2Cols: string[],
  group1Name: string,
  group2Name: string,
  fcThreshold: number = DEFAULT_FC_THRESHOLD,
  pThreshold: number = DEFAULT_P_THRESHOLD,
): {tested: number; untestable: number} {
  // Compute per-row stats into typed arrays first, then bulk-init the columns.
  // Cheaper than column.set() per row, which hops the JS/native bridge each cell.
  const nRows = df.rowCount;
  const fcArr = new Float32Array(nRows);
  const pArr = new Float32Array(nRows);
  const adjArr = new Float32Array(nRows);
  fcArr.fill(DG.FLOAT_NULL);
  pArr.fill(DG.FLOAT_NULL);
  adjArr.fill(DG.FLOAT_NULL);

  const testableIndices: number[] = [];
  const rawPValues: number[] = [];

  for (let i = 0; i < nRows; i++) {
    const vals1 = getGroupValues(df, group1Cols, i);
    const vals2 = getGroupValues(df, group2Cols, i);
    if (vals1.length < 2 || vals2.length < 2) continue;

    fcArr[i] = mean(vals2) - mean(vals1);
    const pVal = tTest(vals1, vals2)['p-value'];
    pArr[i] = pVal;
    testableIndices.push(i);
    rawPValues.push(pVal);
  }

  // BH FDR correction on testable proteins only
  if (rawPValues.length > 0) {
    // alpha is passed for intent/consistency with the user's cutoff; it only
    // affects the (discarded) reject array — BH-adjusted p-values are
    // alpha-independent, so this does not change `corrected`.
    const [, corrected] = fdrcorrection(new Float32Array(rawPValues), pThreshold, 'i');
    for (let j = 0; j < testableIndices.length; j++)
      adjArr[testableIndices[j]] = corrected[j];
  }

  const log2fcCol = df.columns.addNewFloat('log2FC');
  const pValCol = df.columns.addNewFloat('p-value');
  const adjPValCol = df.columns.addNewFloat('adj.p-value');
  log2fcCol.init((i) => fcArr[i]);
  pValCol.init((i) => pArr[i]);
  adjPValCol.init((i) => adjArr[i]);

  const sigCol = df.columns.addNewBool('significant');
  sigCol.init((i) => {
    const adjP = adjArr[i];
    const fc = fcArr[i];
    if (adjP === DG.FLOAT_NULL || fc === DG.FLOAT_NULL) return false;
    return Math.abs(fc) >= fcThreshold && adjP <= pThreshold;
  });

  // Assign semantic types
  log2fcCol.semType = SEMTYPE.LOG2FC;
  pValCol.semType = SEMTYPE.P_VALUE;
  adjPValCol.semType = SEMTYPE.P_VALUE;

  // Mark DE as complete
  df.setTag('proteomics.de_complete', 'true');
  df.fireValuesChanged();

  const tested = testableIndices.length;
  const untestable = df.rowCount - tested;
  return {tested, untestable};
}

/** Build a clean expression dataframe with simple column names (s1, s2, ...)
 *  to avoid encoding issues with special characters in column names.
 *  Bulk-copies via getRawData/fromFloat32Array — per-row set() crosses the JS/native
 *  bridge once per cell and stalls the main thread on real-size frames (see
 *  feedback_dg_column_bulk_init memory note). */
function buildExpressionDf(
  df: DG.DataFrame, group1Cols: string[], group2Cols: string[],
): DG.DataFrame {
  const allCols = [...group1Cols, ...group2Cols];
  const cols: DG.Column[] = [];
  for (let i = 0; i < allCols.length; i++) {
    const src = df.columns.byName(allCols[i]);
    // getRawData returns the underlying Float32/Float64Array; nulls are already
    // encoded as DG.FLOAT_NULL, so a straight copy preserves missingness.
    const raw = src.getRawData();
    const buf = new Float32Array(raw.length);
    buf.set(raw);
    cols.push(DG.Column.fromFloat32Array(`s${i + 1}`, buf));
  }
  return DG.DataFrame.fromColumns(cols);
}

/** Copy log2FC / p.value / adj.p.value / significant columns from an R-script
 *  result frame into `df` as new columns, returning the count of TRUE
 *  significant rows.
 *
 *  Alignment is by the result's `row` column (1-based input index emitted by
 *  the R scripts), NOT by row position. This is load-bearing: DEqMS's
 *  `outputResult()` returns rows sorted by significance, so a positional copy
 *  silently assigns every protein another protein's statistics. When `row` is
 *  absent (defensive fallback for any R output that doesn't emit it) the copy
 *  degrades to positional.
 *
 *  R serialization sometimes flattens booleans to 0/1 or "TRUE"/"FALSE"
 *  strings — the truthiness check covers all three forms. */
export function copyDEResultsToFrame(df: DG.DataFrame, result: DG.DataFrame): number {
  const nRows = df.rowCount;
  const m = result.rowCount;

  const fcRaw  = result.columns.byName('log2FC').getRawData() as Float32Array | Float64Array;
  const pRaw   = result.columns.byName('p.value').getRawData() as Float32Array | Float64Array;
  const aRaw   = result.columns.byName('adj.p.value').getRawData() as Float32Array | Float64Array;
  const rSig   = result.columns.byName('significant');

  // result row j -> df row index. With the `row` key, realign by input index
  // (1-based); without it, identity (positional) as a defensive fallback.
  const rowCol = result.col('row');
  const target = new Int32Array(m);
  if (rowCol) {
    const rowRaw = rowCol.getRawData();
    for (let j = 0; j < m; j++) target[j] = (rowRaw[j] | 0) - 1;
  } else {
    for (let j = 0; j < m; j++) target[j] = j;
  }

  // Pre-fill with FLOAT_NULL so any df row the result doesn't cover stays null
  // rather than a spurious 0 (positional copy used to assume full coverage).
  const log2fcBuf = new Float32Array(nRows); log2fcBuf.fill(DG.FLOAT_NULL);
  const pBuf      = new Float32Array(nRows); pBuf.fill(DG.FLOAT_NULL);
  const aBuf      = new Float32Array(nRows); aBuf.fill(DG.FLOAT_NULL);
  const sigArr    = new Uint8Array(nRows);
  let sigCount = 0;

  for (let j = 0; j < m; j++) {
    const t = target[j];
    if (t < 0 || t >= nRows) continue;
    log2fcBuf[t] = fcRaw[j];
    pBuf[t] = pRaw[j];
    aBuf[t] = aRaw[j];
    const v = rSig.get(j);
    if (v === true || v === 1 || v === '1' || v === 'TRUE') {
      sigArr[t] = 1;
      sigCount++;
    }
  }

  const log2fcCol = DG.Column.fromFloat32Array('log2FC', log2fcBuf);
  const pValCol   = DG.Column.fromFloat32Array('p-value', pBuf);
  const adjPValCol = DG.Column.fromFloat32Array('adj.p-value', aBuf);
  const sigCol = df.columns.addNewBool('significant');
  sigCol.init((i) => sigArr[i] === 1);

  log2fcCol.semType = SEMTYPE.LOG2FC;
  pValCol.semType = SEMTYPE.P_VALUE;
  adjPValCol.semType = SEMTYPE.P_VALUE;

  df.columns.add(log2fcCol);
  df.columns.add(pValCol);
  df.columns.add(adjPValCol);

  return sigCount;
}

/** Run limma moderated t-test via server-side R script.
 *  Returns the number of significant proteins found. */
async function runLimmaDE(
  df: DG.DataFrame,
  group1Cols: string[],
  group2Cols: string[],
  fcThreshold: number,
  pThreshold: number,
): Promise<number> {
  const exprDf = buildExpressionDf(df, group1Cols, group2Cols);

  const result: DG.DataFrame = await grok.functions.call('Proteomics:LimmaDE', {
    exprDf: exprDf,
    nGroup1: group1Cols.length,
    fcThreshold: fcThreshold,
    pThreshold: pThreshold,
  });

  const sigCount = copyDEResultsToFrame(df, result);

  df.setTag('proteomics.de_complete', 'true');
  df.fireValuesChanged();

  return sigCount;
}

/** Run DEqMS peptide-count-weighted differential expression via server-side R script.
 *  Returns the number of significant proteins found. */
async function runDeqmsDE(
  df: DG.DataFrame,
  group1Cols: string[],
  group2Cols: string[],
  peptideCountCol: DG.Column,
  fcThreshold: number,
  pThreshold: number,
): Promise<number> {
  const exprDf = buildExpressionDf(df, group1Cols, group2Cols);

  // Build single-column peptide count DataFrame
  const peptideDf = DG.DataFrame.create(df.rowCount);
  const countCol = peptideDf.columns.addNewInt('count');
  for (let r = 0; r < df.rowCount; r++)
    countCol.set(r, peptideCountCol.isNone(r) ? 1 : Math.round(peptideCountCol.get(r) as number));

  const result: DG.DataFrame = await grok.functions.call('Proteomics:DeqmsDE', {
    exprDf: exprDf,
    nGroup1: group1Cols.length,
    peptideDf: peptideDf,
    fcThreshold: fcThreshold,
    pThreshold: pThreshold,
  });

  const sigCount = copyDEResultsToFrame(df, result);

  df.setTag('proteomics.de_complete', 'true');
  df.setTag('proteomics.de_method', 'deqms');
  df.fireValuesChanged();

  return sigCount;
}

/**
 * R3/D-09: the DE dialog's default Comparison. Returns the declared/intended
 * contrast `${g1} vs ${g2}` — group1 is the parser's first condition, which is
 * Spectronaut's declared Numerator (13-WAVE0-FINDINGS.md A3). The OK-handler's
 * `value === pairs[1]` index logic maps this string to numerator = group1.
 * Previously the default was the alphabetical `pairs[0]` (`${g2} vs ${g1}`),
 * the reverse of the declared contrast — the documented mirror defect
 * (memory project_proteomics_spectronaut_de_direction_default). The dropdown
 * still offers both orientations; only this default changed (direction-only).
 */
export function getDefaultComparison(g1Name: string, g2Name: string): string {
  return `${g1Name} vs ${g2Name}`;
}

/** Show differential expression dialog with prerequisite checks.
 *  @param onComplete Optional callback invoked after DE completes successfully. */
export function showDEDialog(df: DG.DataFrame, onComplete?: () => void): void {
  const groups = getGroups(df);
  if (!groups) {
    grok.shell.warning('Please annotate experimental groups first (Proteomics | Annotate Experiment)');
    return;
  }

  if (df.getTag('proteomics.de_complete') === 'true') {
    grok.shell.warning('Differential expression already performed');
    return;
  }

  const g1 = groups.group1;
  const g2 = groups.group2;

  const infoText = ui.divText(
    `Group 1: ${g1.name} (${g1.columns.length} samples), ` +
    `Group 2: ${g2.name} (${g2.columns.length} samples)`,
  );

  // Comparison direction picker
  const pairs = [`${g2.name} vs ${g1.name}`, `${g1.name} vs ${g2.name}`];
  const comparisonInput = ui.input.choice('Comparison', {
    value: getDefaultComparison(g1.name, g2.name), // D-09: declared contrast, not pairs[0]
    items: pairs,
    nullable: false,
  });

  // Dynamic hint text for FC interpretation. Derive direction from pairs index
  // (not string-splitting) so group names containing " vs " stay intact.
  const hintDiv = ui.div();
  hintDiv.style.cssText = 'font-style:italic; color:#888; font-size:12px; margin-bottom:8px;';
  const updateHint = () => {
    const isReversed = comparisonInput.value === pairs[1];
    const numerator = isReversed ? g1.name : g2.name;
    const denominator = isReversed ? g2.name : g1.name;
    hintDiv.textContent = `Positive log2FC = higher in ${numerator}, Negative log2FC = higher in ${denominator}`;
  };
  comparisonInput.onChanged.subscribe(updateHint);
  updateHint();

  const methodInput = ui.input.choice('Method', {
    value: 'limma',
    items: ['limma', 'DEqMS', 't-test'],
    nullable: false,
  });
  methodInput.setTooltip('limma: moderated t-test; DEqMS: peptide-count-weighted variance; t-test: client-side Welch\'s t-test');

  // Find default peptide count column
  const defaultPeptideCol = df.columns.toList().find((c) =>
    c.name === 'Unique peptides' || c.name === 'Peptides',
  ) ?? undefined;

  const peptideColInput = ui.input.column('Peptide count column', {
    table: df,
    value: defaultPeptideCol,
    filter: (col: DG.Column) => col.type === DG.COLUMN_TYPE.INT || col.type === DG.COLUMN_TYPE.FLOAT,
    nullable: false,
  });
  peptideColInput.setTooltip('Column with peptide/spectra counts per protein');

  // Show/hide peptide count picker based on method selection
  const peptideRow = peptideColInput.root;
  peptideRow.style.display = 'none';
  methodInput.onChanged.subscribe(() => {
    peptideRow.style.display = methodInput.value === 'DEqMS' ? '' : 'none';
  });

  const fcInput = ui.input.float('|log2FC| threshold', {value: DEFAULT_FC_THRESHOLD});
  fcInput.setTooltip('Minimum absolute fold change for significance');

  const pInput = ui.input.float('Adj. p-value threshold', {value: DEFAULT_P_THRESHOLD});
  pInput.setTooltip('Maximum adjusted p-value for significance');

  ui.dialog('Differential Expression')
    .add(infoText)
    .add(comparisonInput)
    .add(hintDiv)
    .add(methodInput)
    .add(peptideColInput)
    .add(fcInput)
    .add(pInput)
    .onOK(async () => {
      const fc = fcInput.value ?? DEFAULT_FC_THRESHOLD;
      const p = pInput.value ?? DEFAULT_P_THRESHOLD;
      const method = methodInput.value;

      // Determine numerator/denominator from pairs index, not string-splitting,
      // so group names containing " vs " don't break the comparison.
      const reversed = comparisonInput.value === pairs[1];
      const numeratorCols = reversed ? g1.columns : g2.columns;
      const denominatorCols = reversed ? g2.columns : g1.columns;
      const numeratorName = reversed ? g1.name : g2.name;
      const denominatorName = reversed ? g2.name : g1.name;

      const pi = DG.TaskBarProgressIndicator.create(`Running ${method} analysis...`);

      try {
        if (method === 't-test') {
          const result = runDifferentialExpression(
            df, denominatorCols, numeratorCols, denominatorName, numeratorName, fc, p,
          );
          df.setTag('proteomics.de_method', 't-test');
          grok.shell.info(
            `DE complete (t-test): ${result.tested} proteins tested, ${result.untestable} untestable`,
          );
        } else if (method === 'DEqMS') {
          const pepCol = peptideColInput.value;
          if (!pepCol) {
            grok.shell.error('Please select a peptide count column for DEqMS');
            return;
          }
          try {
            const sigCount = await runDeqmsDE(df, denominatorCols, numeratorCols, pepCol, fc, p);
            df.setTag('proteomics.de_method', 'deqms');
            grok.shell.info(`DE complete (DEqMS): ${sigCount} significant proteins`);
          } catch (deqmsErr: any) {
            console.warn('DEqMS failed, trying limma:', deqmsErr);
            pi.update(50, 'DEqMS unavailable, trying limma...');
            grok.shell.warning('DEqMS unavailable, trying limma...');
            try {
              const sigCount = await runLimmaDE(df, denominatorCols, numeratorCols, fc, p);
              df.setTag('proteomics.de_method', 'limma');
              grok.shell.info(`DE complete (limma fallback): ${sigCount} significant proteins`);
            } catch (limmaErr: any) {
              console.warn('Limma DE also failed, using client-side fallback:', limmaErr);
              pi.update(75, 'Server-side DE unavailable, using client-side t-test...');
              grok.shell.warning('R environment unavailable — using client-side t-test');
              const result = runDifferentialExpression(
                df, denominatorCols, numeratorCols, denominatorName, numeratorName, fc, p,
              );
              df.setTag('proteomics.de_method', 't-test');
              grok.shell.info(
                `DE complete (t-test): ${result.tested} proteins tested, ` +
                `${result.untestable} untestable`,
              );
            }
          }
        } else {
          try {
            const sigCount = await runLimmaDE(df, denominatorCols, numeratorCols, fc, p);
            df.setTag('proteomics.de_method', 'limma');
            grok.shell.info(`DE complete (limma): ${sigCount} significant proteins`);
          } catch (e: any) {
            pi.update(50, 'Limma unavailable, using client-side t-test...');
            grok.shell.warning('R environment unavailable — using client-side t-test');
            console.warn('Limma DE failed, using client-side fallback:', e);
            const result = runDifferentialExpression(
              df, denominatorCols, numeratorCols, denominatorName, numeratorName, fc, p,
            );
            df.setTag('proteomics.de_method', 't-test');
            grok.shell.info(
              `DE complete (t-test): ${result.tested} proteins tested, ` +
              `${result.untestable} untestable`,
            );
          }
        }
      } finally {
        pi.close();
      }

      if (onComplete)
        onComplete();
    })
    .show();
}
