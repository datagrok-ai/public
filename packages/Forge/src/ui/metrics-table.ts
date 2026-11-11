import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {textColumn} from '../columns';
import {ModelRow} from '../generated/db';
import {ClassValues, METRIC_DESCRIPTIONS, METRIC_IDS, METRIC_LABELS, MetricFinding, metricFindings, MetricId,
  MetricValues, Split} from '../metrics/metrics';
import {CLASS_PROBABILITIES, preparationOptionsOf} from '../preparation/preparation-options';
import {featuresHadGaps} from '../storage/dataset-fingerprint';
import {heldOutCount} from '../training/k-fold';
import {cautiousCaption, directionOf, limitCaption, SEARCH_CAPTIONS, SearchRecord, searchRecordOf}
  from '../training/search';
import {MetricsRecord, metricsRecordOf, Splitting, splittingOf} from '../training/train-model';
import {numberColumn, readOnlyGrid} from './data-grid';
import {reportError} from './report-error';
import {bulletList} from './rich-tooltip';

export const METRIC_HEADERS: {[column: string]: string} = {
  'Metric': 'Quality measure of the model.',
  'Train': 'Value on the training rows, from the model trained on all of them.',
  'Validation': 'Value on rows the model did not see: the pooled out-of-fold predictions of k-fold cross-validation, ' +
    'or the predictions of the held out rows by a model of the others.',
  'Test': 'Value on the rows held out before the search: the final model applied to rows neither the search nor ' +
    'the training saw.',
  'Evaluation': 'Value on the evaluated table: the saved model applied to its rows.',
};

/** A model applied to a table with known targets (Run Evaluation): its metrics there and the rows. */
export interface MetricsEvaluation {
  values: MetricValues;
  tableName: string;
  rowCount: number;
  skippedRows: number;
  /** The evaluation's per-class values. */
  classes?: ClassValues[];
}

export interface MetricsSummary {
  metrics: MetricsRecord;
  rowCount: number;
  skippedRows: number;
  /** The features had missing values and the method took them as they are. */
  keepsMissingValues?: boolean;
  /** How Validation was measured. */
  splitting?: Splitting;
  /** The split seed. */
  seed?: number;
  evaluation?: MetricsEvaluation;
  /** The search that chose the method and settings. */
  search?: SearchRecord;
  /** The classification threshold of a model saved with Predict probability (two classes). */
  threshold?: number;
  /** The metric the threshold was optimized for; absent for a threshold set by hand. */
  thresholdMetric?: MetricId;
}

const SPLIT_COLUMNS: Record<Split, string> = {train: 'Train', validation: 'Validation'};

/** The value columns of [summary] with their values: **Train**, **Validation**, after a search with a hold out
 * **Test** and, with an evaluation, **Evaluation**. */
export function valueColumns({metrics, evaluation}: Pick<MetricsSummary, 'metrics' | 'evaluation'>):
  [string, MetricValues][] {
  const columns: [string, MetricValues][] = [['Train', metrics.train], ['Validation', metrics.validation]];
  if (metrics.test !== undefined)
    columns.push(['Test', metrics.test]);
  if (evaluation !== undefined)
    columns.push(['Evaluation', evaluation.values]);
  return columns;
}

/** The metric rows of [summary]: every metric a value column other than Train has, in {@link METRIC_IDS} order. */
export function metricIdsOf({metrics, evaluation}: Pick<MetricsSummary, 'metrics' | 'evaluation'>): MetricId[] {
  return METRIC_IDS.filter((id) => metrics.validation[id] !== undefined || metrics.test?.[id] !== undefined ||
    evaluation?.values[id] !== undefined);
}

/** The [findings] that mark the value of the metric [id] in [column]. */
export function findingsOf(findings: MetricFinding[], id: MetricId, column: string): MetricFinding[] {
  return findings.filter((f) => f.id === id && f.splits.some((split) => SPLIT_COLUMNS[split] === column));
}

/** A value in full precision, what [column] measures and the texts of [findings]; '' without a value. */
export function metricTooltip(value: number | undefined, column: string, findings: MetricFinding[]): string {
  return value === undefined ? '' : [`${value}. ${METRIC_HEADERS[column]}`, ...findings.map((f) => f.text)].join(' ');
}

/** `Test: <rows> rows held out before the search`. */
export function holdOutLine(rows: number): string {
  return `Test: ${rows} rows held out before the search`;
}

/** The metrics [ids], one row each: **Metric** (the label), then the value columns of [summary]. */
function metricsFrame(summary: MetricsSummary, ids: MetricId[]): DG.DataFrame {
  const values = valueColumns(summary)
    .map(([name, metricValues]) => numberColumn(name, ids.map((id) => metricValues[id] ?? null)));
  return DG.DataFrame.fromColumns([textColumn('Metric', ids.map((id) => METRIC_LABELS[id])), ...values]);
}

/** The metrics grid (hover a header or a metric for its meaning; a value a finding marks is in the failure color,
 * the finding in its tooltip), then a bullet list of the rows, the validation, the seed, the positive class and the
 * evaluation: the **Performance** panel of the results and the model's Performance pane. {@link update} shows another
 * summary of the same metric rows in the same grid, which keeps its column widths and scroll. */
export class MetricsTable {
  readonly root: HTMLElement;
  private readonly ids: MetricId[];
  private readonly hasEvaluation: boolean;
  private readonly hasTest: boolean;
  private readonly grid: DG.Grid;
  private readonly list = ui.element('ul');
  private summary: MetricsSummary;
  private findings: MetricFinding[];

  constructor(summary: MetricsSummary) {
    this.summary = summary;
    this.ids = metricIdsOf(summary);
    this.hasEvaluation = summary.evaluation !== undefined;
    this.hasTest = summary.metrics.test !== undefined;
    this.findings = metricFindings(summary.metrics);
    // A value in full precision, with what the column measures and the findings of the cell; read from the current
    // summary.
    const valueTooltip = (column: string) => (row: number) => {
      const value = valueColumns(this.summary).find(([name]) => name === column)?.[1][this.ids[row]];
      return metricTooltip(value, column, this.findingsOf(column, row));
    };
    this.grid = readOnlyGrid(metricsFrame(summary, this.ids), {headers: METRIC_HEADERS, cells: {
      'Metric': (row) => METRIC_DESCRIPTIONS[this.ids[row]],
      'Train': valueTooltip('Train'),
      'Validation': valueTooltip('Validation'),
      'Test': valueTooltip('Test'),
      'Evaluation': valueTooltip('Evaluation'),
    }});
    this.grid.onCellPrepare((cell) => {
      const column = cell.tableColumn;
      const row = cell.tableRowIndex;
      if (cell.isTableCell && column !== null && row !== null && this.findingsOf(column.name, row).length > 0)
        cell.style.textColor = DG.Color.failure;
    });
    this.root = ui.divV([this.grid.root, this.list]);
    this.fillList();
  }

  /** Shows [summary] and returns true when it has the metric rows and the value columns of the grid; false (nothing
   * changed) when not. */
  update(summary: MetricsSummary): boolean {
    const ids = metricIdsOf(summary);
    if (ids.length !== this.ids.length || ids.some((id, i) => id !== this.ids[i]) ||
      (summary.evaluation !== undefined) !== this.hasEvaluation ||
      (summary.metrics.test !== undefined) !== this.hasTest)
      return false;
    this.summary = summary;
    this.findings = metricFindings(summary.metrics);
    for (const [name, values] of valueColumns(summary))
      this.grid.dataFrame.getCol(name).init((i) => values[ids[i]] ?? null);
    this.fillList();
    this.grid.invalidate();
    return true;
  }

  /** The findings that mark the value of [column] in [row]. */
  private findingsOf(column: string, row: number): MetricFinding[] {
    return findingsOf(this.findings, this.ids[row], column);
  }

  private fillList(): void {
    const {metrics, rowCount, skippedRows, keepsMissingValues, splitting, seed, evaluation, search, threshold,
      thresholdMetric} = this.summary;
    const items: (string | HTMLElement)[][] = [];
    if (keepsMissingValues)
      items.push(['Missing values: kept (the method handles them)']);
    if (skippedRows > 0)
      items.push([`Rows: ${rowCount} used, ${skippedRows} skipped (missing values)`]);
    const validation = splitting === undefined ? undefined : validationText(splitting, rowCount);
    if (validation !== undefined) {
      items.push([validation]);
      if (seed !== undefined && splitting?.isShuffled !== false)
        items.push(seedLine(seed));
    }
    if (search !== undefined) {
      const methods = `${search.leaders.length} ${search.leaders.length === 1 ? 'method' : 'methods'}`;
      const label = METRIC_LABELS[search.metric];
      const rule = search.mode === 'best' ? `the ${directionOf(search.metric) === 1 ? 'highest' : 'lowest'} ` +
        cautiousCaption(search.metric) : 'the one-standard-error rule';
      items.push([`Search: ${SEARCH_CAPTIONS[search.mode]} by ${label}, ` +
        `${limitCaption(search)}: ${search.budgetSeconds === undefined ? '' : `${search.trials} trials of `}` +
        `${methods}, chosen by ${rule}`]);
      items.push([search.holdOut !== undefined ? holdOutLine(search.holdOut.rows) :
        'Validation after a search is optimistic: the best of the trials was chosen on it']);
    }
    if (metrics.positiveClass !== undefined)
      items.push([`Positive class: ${metrics.positiveClass}`]);
    if (threshold !== undefined) {
      items.push([thresholdMetric === undefined ? `Classification threshold: ${threshold}` :
        `Classification threshold: ${threshold}, the best ${METRIC_LABELS[thresholdMetric]} on the validation rows: ` +
        'its validation values are optimistic']);
    }
    if (evaluation !== undefined) {
      const skipped = evaluation.skippedRows > 0 ? `, ${evaluation.skippedRows} skipped (missing values)` : '';
      items.push([`Evaluation: ${evaluation.tableName} (${evaluation.rowCount} rows${skipped})`]);
    }
    bulletList(items, this.list);
    ui.setDisplay(this.list, items.length > 0);
  }
}

/** `Validation: <k>-fold cross-validation on <n> rows` or `Validation: hold out of <p>% (<m> rows)`, then
 * `, stratified` when it was; undefined for another scheme. */
function validationText({scheme, folds, trainFraction, isStratified}: Splitting, rowCount: number):
  string | undefined {
  const stratified = isStratified ? ', stratified' : '';
  if (scheme === 'kfold' && folds !== undefined)
    return `Validation: ${folds}-fold cross-validation on ${rowCount} rows${stratified}`;
  if (scheme !== 'holdout' || trainFraction === undefined)
    return undefined;
  const testFraction = 1 - trainFraction;
  return `Validation: hold out of ${Math.round(100 * testFraction)}% ` +
    `(${heldOutCount(rowCount, testFraction)} rows)${stratified}`;
}

/** `Split seed: <seed>` with the seed selectable and a copy icon right after it. */
function seedLine(seed: number): (string | HTMLElement)[] {
  const copy = ui.icons.copy(() => {
    navigator.clipboard.writeText(`${seed}`)
      .then(() => grok.shell.info(`Split seed ${seed} copied.`))
      .catch((e) => reportError(e));
  }, 'Copy the seed');
  return ['Split seed: ', ui.span([`${seed}`], 'forge-seed'), ' ', copy];
}

/** The metrics stored with [model] and how they were measured; null without metrics. */
export function storedSummary(model: ModelRow): MetricsSummary | null {
  const metrics = metricsRecordOf(model.metrics);
  if (metrics === null)
    return null;
  const options = preparationOptionsOf(model.options);
  const {missingValues} = options;
  return {metrics, rowCount: model.row_count ?? 0, seed: model.seed, skippedRows: missingValues?.skippedRows ?? 0,
    keepsMissingValues: missingValues?.mode === 'keep' && featuresHadGaps(model.dataset_fingerprint),
    splitting: splittingOf(model.splitting), search: searchRecordOf(model.search),
    threshold: options.postprocessingInfo.includes(CLASS_PROBABILITIES) ? options.binaryClassificationThreshold :
      undefined, thresholdMetric: options.thresholdMetric};
}
