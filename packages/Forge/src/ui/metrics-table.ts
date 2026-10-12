import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelRow} from '../generated/db';
import {METRIC_DESCRIPTIONS, METRIC_IDS, METRIC_LABELS, MetricFinding, metricFindings, MetricId, MetricValues, Split}
  from '../metrics/metrics';
import {preparationOptionsOf} from '../preparation/preparation-options';
import {MetricsRecord, metricsRecordOf} from '../training/train-model';
import {numberColumn, readOnlyGrid, textColumn} from './data-grid';
import {reportError} from './report-error';
import {bulletList} from './rich-tooltip';

const METRIC_HEADERS: {[column: string]: string} = {
  'Metric': 'Quality measure of the model.',
  'Train': 'Value on the training rows, from the model trained on all of them.',
  'Validation': 'Value on rows the model did not see: the pooled out-of-fold predictions of 5-fold cross-validation.',
  'Evaluation': 'Value on the evaluated table: the saved model applied to its rows.',
};

/** A model applied to a table with known targets (Run Evaluation): its metrics there and the rows. */
export interface MetricsEvaluation {
  values: MetricValues;
  tableName: string;
  rowCount: number;
  skippedRows: number;
}

export interface MetricsSummary {
  metrics: MetricsRecord;
  rowCount: number;
  skippedRows: number;
  folds?: number;
  seed?: number;
  evaluation?: MetricsEvaluation;
}

const SPLIT_COLUMNS: Record<Split, string> = {train: 'Train', validation: 'Validation'};

/** The value columns of [summary] with their values: **Train**, **Validation** and, with an evaluation,
 * **Evaluation**. */
function valueColumns({metrics, evaluation}: MetricsSummary): [string, MetricValues][] {
  const columns: [string, MetricValues][] = [['Train', metrics.train], ['Validation', metrics.validation]];
  if (evaluation !== undefined)
    columns.push(['Evaluation', evaluation.values]);
  return columns;
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
  private readonly grid: DG.Grid;
  private readonly list = ui.element('ul');
  private summary: MetricsSummary;
  private findings: MetricFinding[];

  constructor(summary: MetricsSummary) {
    this.summary = summary;
    this.ids = MetricsTable.idsOf(summary);
    this.hasEvaluation = summary.evaluation !== undefined;
    this.findings = metricFindings(summary.metrics);
    // A value in full precision, with what the column measures and the findings of the cell; read from the current
    // summary.
    const valueTooltip = (column: string) => (row: number) => {
      const value = valueColumns(this.summary).find(([name]) => name === column)?.[1][this.ids[row]];
      return value === undefined ? '' :
        [`${value}. ${METRIC_HEADERS[column]}`, ...this.findingsOf(column, row).map((f) => f.text)].join(' ');
    };
    this.grid = readOnlyGrid(metricsFrame(summary, this.ids), {headers: METRIC_HEADERS, cells: {
      'Metric': (row) => METRIC_DESCRIPTIONS[this.ids[row]],
      'Train': valueTooltip('Train'),
      'Validation': valueTooltip('Validation'),
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
    const ids = MetricsTable.idsOf(summary);
    if (ids.length !== this.ids.length || ids.some((id, i) => id !== this.ids[i]) ||
      (summary.evaluation !== undefined) !== this.hasEvaluation)
      return false;
    this.summary = summary;
    this.findings = metricFindings(summary.metrics);
    for (const [name, values] of valueColumns(summary))
      this.grid.dataFrame.getCol(name).init((i) => values[ids[i]] ?? null);
    this.fillList();
    this.grid.invalidate();
    return true;
  }

  private static idsOf({metrics}: MetricsSummary): MetricId[] {
    return METRIC_IDS.filter((id) => metrics.validation[id] !== undefined);
  }

  /** The findings that mark the value of [column] in [row]. */
  private findingsOf(column: string, row: number): MetricFinding[] {
    const id = this.ids[row];
    return this.findings.filter((f) => f.id === id && f.splits.some((split) => SPLIT_COLUMNS[split] === column));
  }

  private fillList(): void {
    const {metrics, rowCount, skippedRows, folds, seed, evaluation} = this.summary;
    const items: (string | HTMLElement)[][] = [];
    if (skippedRows > 0)
      items.push([`Rows: ${rowCount} used, ${skippedRows} skipped (missing values)`]);
    if (folds !== undefined && seed !== undefined)
      items.push([`Validation: ${folds}-fold cross-validation on ${rowCount} rows`], seedLine(seed));
    if (metrics.positiveClass !== undefined)
      items.push([`Positive class: ${metrics.positiveClass}`]);
    if (evaluation !== undefined) {
      const skipped = evaluation.skippedRows > 0 ? `, ${evaluation.skippedRows} skipped (missing values)` : '';
      items.push([`Evaluation: ${evaluation.tableName} (${evaluation.rowCount} rows${skipped})`]);
    }
    bulletList(items, this.list);
    ui.setDisplay(this.list, items.length > 0);
  }
}

/** `Seed: <seed>` with the seed selectable and a copy icon right after it. */
function seedLine(seed: number): (string | HTMLElement)[] {
  const copy = ui.icons.copy(() => {
    navigator.clipboard.writeText(`${seed}`)
      .then(() => grok.shell.info(`Seed ${seed} copied.`))
      .catch((e) => reportError(e));
  }, 'Copy the seed');
  return ['Seed: ', ui.span([`${seed}`], 'forge-seed'), ' ', copy];
}

/** The metrics stored with [model] and how they were measured; null without metrics. */
export function storedSummary(model: ModelRow): MetricsSummary | null {
  const metrics = metricsRecordOf(model.metrics);
  if (metrics === null)
    return null;
  const folds: unknown = model.splitting?.folds;
  return {metrics, rowCount: model.row_count ?? 0, seed: model.seed,
    skippedRows: preparationOptionsOf(model.options).missingValues?.skippedRows ?? 0,
    folds: typeof folds === 'number' ? folds : undefined};
}
