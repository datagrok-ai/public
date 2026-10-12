import * as DG from 'datagrok-api/dg';
import {ApplyRequest, checkApplyRequest, DEFAULT_BATCH_SIZE, LoadedModel, predictRows} from '../apply/apply-model';
import {ColumnMapping} from '../apply/column-matching';
import {LoopProgress} from '../engines/engine-calls';
import {aucOf, MetricValues, metricsOf, probabilityValues, taskOf} from '../metrics/metrics';
import {BINARY_CLASSIFICATION} from '../preparation/preparation-options';
import {columnRowCopy} from '../preparation/row-copy';

/** A saved model to measure on [table], whose [target] holds the actual values. */
export interface EvaluationRequest {
  model: LoadedModel;
  table: DG.DataFrame;
  target: DG.Column;
  mapping: ColumnMapping;
  positiveClass?: string;
}

/** The model's predictions of the kept rows of the table and their metrics against the actual values. */
export interface EvaluationResult {
  /** What the method got. */
  columns: DG.Column[];
  /** The actual values of the kept rows: the request's target, or a copy of its kept rows. */
  target: DG.Column;
  /** The labels (or values) of the kept rows. */
  prediction: DG.Column;
  /** Predict probability: the method's scores. */
  score?: DG.Column;
  keptRows: DG.BitSet | null;
  skippedRows: number;
  rowCount: number;
  values: MetricValues;
  blob: Uint8Array;
}

/** Applies the model to the table's rows (those with a missing feature value skipped) and measures the predictions
 * against the target; a row without an actual value is left out of the metrics. No column is added and no
 * application is recorded. */
export async function evaluateModel(request: EvaluationRequest, progress?: LoopProgress): Promise<EvaluationResult> {
  const {model, table, mapping} = request;
  const applyRequest: ApplyRequest = {model, table, mapping, batchSize: DEFAULT_BATCH_SIZE,
    missingValues: {mode: 'skip'}};
  checkApplyRequest(applyRequest);
  const rows = await predictRows(applyRequest, progress);
  const target = rows.keptRows === null ? request.target : columnRowCopy(request.target, rows.keptRows);
  const result = {columns: rows.columns, target, prediction: rows.labels, keptRows: rows.keptRows,
    skippedRows: rows.skippedRows, rowCount: target.length, blob: rows.blob};
  const {options} = model;
  const positiveClass = request.positiveClass ?? options.positiveClass;
  if (options.postprocessingInfo.includes(BINARY_CLASSIFICATION) && positiveClass !== undefined) {
    const values = probabilityValues(target, rows.labels, positiveClass, aucOf(target, rows.raw, positiveClass));
    return {...result, values, score: rows.raw};
  }
  return {...result, values: metricsOf(taskOf(target), target, rows.labels, positiveClass)};
}
