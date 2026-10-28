import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ApplyRequest, checkApplyRequest, DEFAULT_BATCH_SIZE, LoadedModel, missingValuesOf, predictRows}
  from '../apply/apply-model';
import {ColumnMapping} from '../apply/column-matching';
import {LoopProgress} from '../engines/engine-calls';
import {ModelRow} from '../generated/db';
import {ClassProbabilities, ClassValues, Curves, MetricValues, taskOf} from '../metrics/metrics';
import {classValuesOf, curvesOf, measure} from '../metrics/report';
import {BINARY_CLASSIFICATION} from '../preparation/preparation-options';
import {columnRowCopy} from '../preparation/row-copy';
import {openDatasetRef, storedDatasetRef} from '../storage/dataset-ref';

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
  /** The request's positive class, else the model's. */
  positiveClass?: string;
  /** The method's class probabilities of the kept rows, when it returns them. */
  probabilities?: ClassProbabilities;
  /** The curves of the positive class (two classes with a probability or score). */
  curves?: Curves;
  /** The per-class values of the evaluation. */
  classes?: ClassValues[];
}

/** Applies the model to the table's rows (missing values handled as at training, `missingValuesOf`) and measures the
 * predictions against the target; a row without an actual value is left out of the metrics. No column is added and
 * no application is recorded. */
export async function evaluateModel(request: EvaluationRequest, progress?: LoopProgress): Promise<EvaluationResult> {
  const {model, table, mapping} = request;
  const applyRequest: ApplyRequest = {model, table, mapping, batchSize: DEFAULT_BATCH_SIZE,
    missingValues: missingValuesOf(model)};
  checkApplyRequest(applyRequest);
  const rows = await predictRows(applyRequest, progress);
  const target = rows.keptRows === null ? request.target : columnRowCopy(request.target, rows.keptRows);
  const {options} = model;
  const positiveClass = request.positiveClass ?? options.positiveClass;
  const result = {columns: rows.columns, target, prediction: rows.labels, keptRows: rows.keptRows,
    skippedRows: rows.skippedRows, rowCount: target.length, blob: rows.blob, positiveClass};
  const isProbability = options.postprocessingInfo.includes(BINARY_CLASSIFICATION) && positiveClass !== undefined;
  const probabilities = isProbability ? {classes: [positiveClass], columns: [rows.raw]} : rows.probabilities;
  const {values, report} = await measure(isProbability ? 'classification' : taskOf(target), target, rows.labels,
    positiveClass, probabilities, true);
  const curves = curvesOf(report);
  return {...result, values, ...(isProbability ? {score: rows.raw} : {probabilities: rows.probabilities}),
    ...(report === undefined ? {} : {classes: classValuesOf(report)}), ...(curves === undefined ? {} : {curves})};
}

/** The data the model keeps a reference to, else its uploaded copy, opened without joining the workspace; null when
 * it keeps neither. */
export async function storedTrainingTable(model: Pick<ModelRow, 'dataset_ref' | 'dataset_table_id'>):
  Promise<DG.DataFrame | null> {
  const ref = storedDatasetRef(model.dataset_ref);
  return ref !== null ? openDatasetRef(ref) :
    model.dataset_table_id ? grok.dapi.tables.getTable(model.dataset_table_id) : null;
}
