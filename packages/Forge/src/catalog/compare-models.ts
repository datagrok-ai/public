import * as DG from 'datagrok-api/dg';
import {textColumn} from '../columns';
import {ForgeError} from '../forge-error';
import {ModelRow} from '../generated/db';
import {METRIC_IDS, METRIC_LABELS, Split} from '../metrics/metrics';
import {metricsRecordOf} from '../training/train-model';

/** The `model` columns a comparison reads (the system columns come with any query). */
export const COMPARE_COLUMNS = ['name', 'description', 'engine_name', 'task', 'target_name', 'row_count',
  'metrics'] as const;
export type CompareModelRow = Pick<ModelRow, typeof COMPARE_COLUMNS[number] | 'id' | 'created_on'>;

const SPLITS = ['train', 'validation'] as const;
const FORM_SUMMARY = ['Name', 'Method', 'Task', 'Target', 'Training rows', 'Created'];

/** `<metric> (<split>)` for each of [splits], for every metric at least one of the rows' stored metrics has in any of
 * them; a metric a row lacks is empty. */
export function metricColumns(rows: {metrics?: unknown}[], splits: readonly Split[]): DG.Column[] {
  const metrics = rows.map((row) => metricsRecordOf(row.metrics));
  const columns: DG.Column[] = [];
  for (const id of METRIC_IDS) {
    const values = splits.map((split) => metrics.map((m) => m?.[split][id] ?? null));
    if (values.some((list) => list.some((v) => v !== null))) {
      for (const [i, split] of splits.entries())
        columns.push(DG.Column.fromList(DG.COLUMN_TYPE.FLOAT, `${METRIC_LABELS[id]} (${split})`, values[i]));
    }
  }
  return columns;
}

/** One row per model, in the given order: the summary columns, then `<metric> (train)` and `<metric> (validation)`
 * for every metric at least one of the models has; a metric a model lacks is empty. */
export function compareModels(rows: CompareModelRow[]): DG.DataFrame {
  if (rows.length < 2)
    throw new ForgeError('Select at least two models to compare.');
  const df = DG.DataFrame.fromColumns([
    textColumn('Name', rows.map((r) => r.name)),
    textColumn('Description', rows.map((r) => r.description)),
    textColumn('Method', rows.map((r) => r.engine_name)),
    textColumn('Task', rows.map((r) => r.task)),
    textColumn('Target', rows.map((r) => r.target_name)),
    DG.Column.fromList(DG.COLUMN_TYPE.INT, 'Training rows', rows.map((r) => r.row_count ?? null)),
    DG.Column.fromList(DG.COLUMN_TYPE.DATE_TIME, 'Created', rows.map((r) => r.created_on)),
    ...metricColumns(rows, SPLITS),
  ]);
  df.name = 'Compare models';
  return df;
}

/** The columns of a {@link compareModels} table a Forms comparison shows: the summary and the validation metrics. */
export function compareFormFields(df: DG.DataFrame): string[] {
  return [...FORM_SUMMARY, ...df.columns.names().filter((name) => name.endsWith(' (validation)'))];
}
