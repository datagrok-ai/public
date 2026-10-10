import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {PREDICTION_TAG} from '../constants';
import {Engine, isComplete} from '../engines/engine';
import {apply, LoopProgress, yieldToEventLoop} from '../engines/engine-calls';
import {EngineRegistry} from '../engines/engine-registry';
import {errorMessage, ForgeError} from '../forge-error';
import {ApplicationInsert, ApplicationSource, forgeDb, ModelRow} from '../generated/db';
import {MissingValuesSettings, prepareMissingValues} from '../preparation/missing-values';
import {isRecord, PreparationOptions, preparationOptionsOf} from '../preparation/preparation-options';
import {replayPostprocessing, replayPreprocessing} from '../preparation/pipeline';
import {rowCopy} from '../preparation/row-copy';
import {ownBlob, UUID} from '../storage/model-store';
import {ColumnSchema} from '../training/train-model';
import {recordApplication} from './application-store';
import {ColumnMapping, exactMapping, MappingProblem, mappingProblems} from './column-matching';
import {featureColumns} from './feature-columns';

export const DEFAULT_BATCH_SIZE = 10000;
const YIELD_MS = 50;
// Prediction types copied through their raw data; a text column's raw data indexes its own categories.
const RAW_COPY_TYPES: string[] = [DG.COLUMN_TYPE.INT, DG.COLUMN_TYPE.FLOAT, DG.COLUMN_TYPE.QNUM];

/** The `model` columns applying reads (the system columns come with any query). */
export const APPLY_COLUMNS = ['name', 'engine_name', 'engine_namespace', 'target_name', 'features', 'options',
  'blob'] as const;
export type ApplyModelRow = Pick<ModelRow, typeof APPLY_COLUMNS[number] | 'id' | 'created_on'>;

export interface LoadedModel {
  row: ApplyModelRow;
  engine: Engine;
  /** The features a table must provide ({@link requiredFeaturesOf}). */
  features: ColumnSchema[];
  options: PreparationOptions;
  blobPath: string;
}

export interface ApplyRequest {
  model: LoadedModel;
  table: DG.DataFrame;
  mapping: ColumnMapping;
  batchSize: number;
  missingValues: MissingValuesSettings;
}

export interface ApplyResult { column: DG.Column; skippedRows: number }

const MODEL_ID = new RegExp(`^${UUID}$`, 'i');

/** The model with the id or the unique name [idOrName], ready to apply. */
export async function loadModel(idOrName: string): Promise<LoadedModel> {
  const rows: ApplyModelRow[] = await forgeDb.models.query()
    .where(MODEL_ID.test(idOrName) ? 'id' : 'name', '=', idOrName)
    .select(...APPLY_COLUMNS)
    .top(2);
  if (rows.length === 0)
    throw new ForgeError(`No Forge model '${idOrName}' is available to you.`);
  if (rows.length > 1)
    throw new ForgeError(`Several models are named '${idOrName}'. Use the model id.`);
  return loadedModelOf(rows[0], EngineRegistry.discover());
}

/** {@link loadModel} for a row already read, with the engines already discovered. */
export function loadedModelOf(row: ApplyModelRow, engines: Engine[]): LoadedModel {
  const name = row.name;
  if (!row.blob) {
    throw new ForgeError(`The model '${name}' has no model file, so it cannot be applied. ` +
      'It was not trained and saved in Forge.');
  }
  const blob = ownBlob(row.blob);
  if (blob === null)
    throw new ForgeError(`The model '${name}' points to a file outside Forge's storage and cannot be applied.`);
  const features = requiredFeaturesOf(row);
  if (features === null)
    throw new ForgeError(`The model '${name}' has no feature list.`);
  const engine = engines.find((e) => e.name === row.engine_name && isComplete(e));
  if (engine === undefined) {
    throw new ForgeError(`The method '${row.engine_name}' is not installed. ` +
      `Install the ${row.engine_namespace} package.`);
  }
  return {row, engine, features, options: preparationOptionsOf(row.options), blobPath: blob.path};
}

/** Adds the prediction column to the request's table and records the application as completed, failed or
 * cancelled; mapping problems and an empty table are refused first, without a record. */
export async function applyAndRecord(request: ApplyRequest, source: ApplicationSource,
  progress?: LoopProgress): Promise<ApplyResult> {
  checkRequest(request);
  const {model, table} = request;
  const startedOn = Date.now();
  const record = (outcome: Pick<ApplicationInsert, 'status' | 'column_name' | 'skipped_rows' | 'error'>) =>
    recordApplication({model_id: model.row.id, table_name: table.name, row_count: table.rowCount, source,
      ...outcome, duration_ms: Date.now() - startedOn});
  let result: ApplyResult;
  try {
    result = await predict(request, progress);
  } catch (e) {
    // The caller gets the application's own error, even when its record cannot be written.
    await Promise.allSettled([record(progress?.canceled ? {status: 'cancelled'} :
      {status: 'failed', error: errorMessage(e)})]);
    throw e;
  }
  await record({status: 'completed', column_name: result.column.name, skipped_rows: result.skippedRows});
  return result;
}

/** {@link applyAndRecord} under the cancellable task-bar progress "Predicting <target>". */
export async function applyWithProgress(request: ApplyRequest, source: ApplicationSource): Promise<ApplyResult> {
  const progress = DG.TaskBarProgressIndicator.create(`Predicting ${request.model.row.target_name}`,
    {cancelable: true});
  try {
    return await applyAndRecord(request, source, progress);
  } finally {
    progress.close();
  }
}

/** The body of `Forge:applyModel`: exact names plus the given pairs, rows with missing values skipped. */
export async function runApplyModel(model: string, table: DG.DataFrame,
  columnNamesMap: {[feature: string]: string} | null, showProgress: boolean): Promise<DG.DataFrame> {
  const loaded = await loadModel(model);
  const mapping = exactMapping(loaded.features, table);
  for (const [feature, column] of Object.entries(columnNamesMap ?? {}))
    mapping.set(feature, column);
  const problems = mappingProblems(loaded.features, mapping, table, loaded.engine.name);
  if (problems.length > 0)
    throw new ForgeError(apiMappingMessage(problems));
  const request: ApplyRequest = {model: loaded, table, mapping, batchSize: DEFAULT_BATCH_SIZE,
    missingValues: {mode: 'skip'}};
  await (showProgress ? applyWithProgress(request, 'api') : applyAndRecord(request, 'api'));
  return table;
}

/** `<target> (predicted)`, then `<target> (predicted 2)`, ...: the first name no column has, ignoring case. */
export function predictionName(table: DG.DataFrame, target: string): string {
  let name = `${target} (predicted)`;
  for (let i = 2; table.col(name) !== null; i++)
    name = `${target} (predicted ${i})`;
  return name;
}

function checkRequest({model, table, mapping}: ApplyRequest): void {
  const problems = mappingProblems(model.features, mapping, table, model.engine.name);
  if (problems.length > 0)
    throw new ForgeError(problems.map((p) => p.message).join(' '));
  if (table.rowCount === 0)
    throw new ForgeError('The table has no rows to predict.');
}

async function predict(request: ApplyRequest, progress?: LoopProgress): Promise<ApplyResult> {
  const {model, table, mapping} = request;
  const prepared = await prepareMissingValues(featureColumns(table, model.features, mapping), undefined,
    request.missingValues);
  const rowCount = table.rowCount - prepared.skippedRows;
  if (rowCount === 0) {
    throw new ForgeError('Every row has a missing value in the columns the model needs, so nothing can be ' +
      'predicted. Fill the missing values or choose Impute.');
  }
  const features = replayPreprocessing(prepared.features, model.options);
  const blob = await grok.dapi.files.readAsBytes(model.blobPath);
  const predictions = replayPostprocessing(
    await predictInBatches(model.engine, features, rowCount, blob, request.batchSize, progress), model.options);
  const column = prepared.keptRows === null ? predictions :
    scattered(predictions, prepared.keptRows, table.rowCount);
  column.name = predictionName(table, model.row.target_name);
  column.setTag(PREDICTION_TAG, model.row.id);
  table.columns.add(column);
  return {column, skippedRows: prepared.skippedRows};
}

/** One engine call on [features] themselves ([rowCount] rows) when they fit a batch; otherwise one call per batch of
 * rows, with a pause for the event loop at most every YIELD_MS, so a cancel and the progress repaint get through. */
async function predictInBatches(engine: Engine, features: DG.Column[], rowCount: number, blob: Uint8Array,
  batchSize: number, progress?: LoopProgress): Promise<DG.Column> {
  let yieldedAt = Date.now();
  const predictBatch = async (start: number): Promise<DG.Column> => {
    if (Date.now() - yieldedAt >= YIELD_MS) {
      await yieldToEventLoop();
      yieldedAt = Date.now();
    }
    if (progress?.canceled)
      throw new ForgeError('Application was cancelled.');
    const end = Math.min(start + batchSize, rowCount);
    const batch = rowCount <= batchSize ? features : rowCopy(features, rowRange(rowCount, start, end));
    const prediction = await apply(engine, batch, blob);
    progress?.update(100 * end / rowCount, `Rows ${start + 1}-${end} of ${rowCount}`);
    return prediction;
  };
  const first = await predictBatch(0);
  if (rowCount <= batchSize)
    return first;
  const column = DG.Column.fromType(first.type, first.name, rowCount);
  copyInto(column, first, 0);
  for (let start = batchSize; start < rowCount; start += batchSize)
    copyInto(column, await predictBatch(start), start);
  return column;
}

/** The rows [start, end) of [rowCount] rows; only the batch's bits are visited, not every row as a predicate would. */
function rowRange(rowCount: number, start: number, end: number): DG.BitSet {
  const words = new Uint32Array((rowCount + 31) >>> 5);
  for (let r = start; r < end; r++)
    words[r >>> 5] |= 1 << (r & 31);
  return DG.BitSet.fromBytes(words.buffer, rowCount);
}

function copyInto(target: DG.Column, source: DG.Column, offset: number): void {
  const targetRaw = rawNumbers(target);
  const sourceRaw = rawNumbers(source);
  if (targetRaw !== null && sourceRaw !== null && target.type === source.type)
    targetRaw.set(sourceRaw.subarray(0, source.length), offset);
  else {
    for (let r = 0; r < source.length; r++)
      target.set(offset + r, source.get(r), false);
  }
}

/** [prediction] of the kept rows spread over all [rowCount] rows; the skipped rows stay empty. */
function scattered(prediction: DG.Column, keptRows: DG.BitSet, rowCount: number): DG.Column {
  const column = DG.Column.fromType(prediction.type, prediction.name, rowCount);
  const rows = keptRows.getSelectedIndexes();
  const targetRaw = rawNumbers(column);
  const sourceRaw = rawNumbers(prediction);
  if (targetRaw !== null && sourceRaw !== null) {
    for (let i = 0; i < rows.length; i++)
      targetRaw[rows[i]] = sourceRaw[i];
  } else {
    for (let i = 0; i < rows.length; i++)
      column.set(rows[i], prediction.get(i), false);
  }
  return column;
}

/** The raw values of a numerical column, empty cells included, or null for a column of another type. */
function rawNumbers(col: DG.Column): Int32Array | Float32Array | Float64Array | null {
  if (!RAW_COPY_TYPES.includes(col.type))
    return null;
  const raw = col.getRawData();
  return raw instanceof Int32Array || raw instanceof Float32Array || raw instanceof Float64Array ? raw : null;
}

/** One sentence for the features without a column, then the other problems. */
function apiMappingMessage(problems: MappingProblem[]): string {
  const unmapped = problems.filter((p) => p.isUnmapped).map((p) => `'${p.feature}'`);
  const sentences = problems.filter((p) => !p.isUnmapped).map((p) => p.message);
  if (unmapped.length === 1)
    sentences.unshift(`The table has no column ${unmapped[0]} the model needs. Map it in columnNamesMap.`);
  else if (unmapped.length > 1) {
    sentences.unshift(`The table has no columns ${unmapped.join(', ')} the model needs. ` +
      'Map them in columnNamesMap.');
  }
  return sentences.join(' ');
}

/** The `{columns: [{name, type, semType?}]}` feature list of a model row, or null when it has none. */
export function featureSchemasOf(features: unknown): ColumnSchema[] | null {
  const columns = isRecord(features) ? features.columns : undefined;
  if (!Array.isArray(columns) || columns.length === 0)
    return null;
  const schemas = columns.map(parsedColumnSchema).filter((s): s is ColumnSchema => s !== null);
  return schemas.length === columns.length ? schemas : null;
}

/** The features a table must provide to apply the model: its feature list without the columns Skip unique categories
 * left out (`options.skippedColumns`); null when it has no feature list. */
export function requiredFeaturesOf(row: Pick<ModelRow, 'features' | 'options'>): ColumnSchema[] | null {
  const features = featureSchemasOf(row.features);
  const skipped = preparationOptionsOf(row.options).skippedColumns ?? [];
  return features === null ? null : features.filter((f) => !skipped.includes(f.name));
}

function parsedColumnSchema(value: unknown): ColumnSchema | null {
  if (!isRecord(value))
    return null;
  const {name, type, semType} = value;
  if (typeof name !== 'string' || typeof type !== 'string')
    return null;
  return typeof semType === 'string' && semType !== '' ? {name, type, semType} : {name, type};
}
