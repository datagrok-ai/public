import * as DG from 'datagrok-api/dg';
import {apiMapping, APPLY_COLUMNS, loadedModelOf} from '../apply/apply-model';
import {EngineRegistry} from '../engines/engine-registry';
import {ForgeError} from '../forge-error';
import {ClassValues, MetricValues} from '../metrics/metrics';
import {findModelRow} from '../storage/model-store';
import {evaluateModel, EvaluationResult, storedTrainingTable} from './evaluation';
import {evaluationInsights, Insight} from './insights';
import {metricsRecordOf} from './train-model';

/** The `model` columns an evaluation reads (the system columns come with any query). */
const EVALUATE_COLUMNS = [...APPLY_COLUMNS, 'metrics', 'dataset_ref', 'dataset_table_id'] as const;

/** The inputs of `Forge:evaluateModel`; a missing optional input is null or undefined. */
export interface EvaluateModelInput {
  model: string;
  table?: DG.DataFrame | null;
  target?: string | null;
  columnNamesMap?: {[feature: string]: string} | null;
  showProgress?: boolean | null;
}

/** An insight as text: what it says, the numbers behind it and a list of them. */
export type InsightText = Omit<Insight, 'panel'>;

export interface EvaluateModelResult {
  id: string;
  name: string;
  /** The evaluated table's name. */
  table: string;
  target: string;
  /** The rows measured. */
  rowCount: number;
  skippedRows: number;
  positiveClass?: string;
  /** The stored train and validation values (empty without them) and the evaluation's. */
  metrics: {train: MetricValues; validation: MetricValues; evaluation: MetricValues};
  /** The evaluation's per-class values (a classifier with EDA's report). */
  classes?: ClassValues[];
  insights: InsightText[];
}

/** The body of `Forge:evaluateModel`: Run Evaluation without the view. The model is measured on [table], by default
 * on the data it keeps a reference to or a copy of; the features are mapped as `Forge:applyModel` maps them. Nothing
 * is added to the table and nothing is recorded. */
export async function runEvaluateModel(input: EvaluateModelInput): Promise<EvaluateModelResult> {
  const row = await findModelRow(input.model, EVALUATE_COLUMNS);
  const model = loadedModelOf(row, EngineRegistry.discover());
  const table = input.table ?? await storedTrainingTable(row);
  if (table === null)
    throw new ForgeError('The model keeps no training data. Pass a table with the actual values.');
  const targetName = input.target || row.target_name;
  const target = table.col(targetName);
  if (target === null) {
    throw new ForgeError(`The table has no column '${targetName}' with the actual values. ` +
      'Pass the column\'s name in target.');
  }
  const mapping = apiMapping(model, table, input.columnNamesMap);
  const stored = metricsRecordOf(row.metrics);
  const progress = input.showProgress === false ? undefined :
    DG.TaskBarProgressIndicator.create(`Evaluating ${row.name}`, {cancelable: true});
  let result: EvaluationResult;
  try {
    result = await evaluateModel({model, table, target, mapping, positiveClass: stored?.positiveClass}, progress);
  } catch (e) {
    throw progress?.canceled ? new ForgeError('Evaluation was cancelled.') : e;
  } finally {
    progress?.close();
  }
  const train = stored?.train ?? {};
  const validation = stored?.validation ?? {};
  const {positiveClass} = result;
  return {id: row.id, name: row.name, table: table.name, target: target.name, rowCount: result.rowCount,
    skippedRows: result.skippedRows, ...(positiveClass === undefined ? {} : {positiveClass}),
    metrics: {train, validation, evaluation: result.values},
    ...(result.classes === undefined ? {} : {classes: result.classes}),
    insights: evaluationInsights(result, {train, validation}).map(insightText)};
}

const insightText = ({panel: _, ...text}: Insight): InsightText => text;
