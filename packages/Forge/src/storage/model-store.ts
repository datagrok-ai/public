import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {Engine} from '../engines/engine';
import {EngineRegistry} from '../engines/engine-registry';
import {errorMessage, ForgeError} from '../forge-error';
import {forgeDb, ModelColumn, ModelInsert, ModelRow} from '../generated/db';
import {TrainingResult} from '../training/train-model';
import {DatasetFingerprint} from './dataset-fingerprint';
import {deleteTrainingCopy, uploadTrainingCopy} from './dataset-copy';
import {ModelFields, modelFieldsOf, ModelStorage, StorageChoice} from './model-fields';
import {linkTrainingRun} from './training-run-store';

export const BLOB_ROOT = 'System:DomainFiles/forge/model';
const UUID = '[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}';
const FILE_PREFIX = 'file://';
const OWN_BLOB = new RegExp(`^${FILE_PREFIX}${BLOB_ROOT}/(${UUID})/([^/]+)$`, 'i');
const MODEL_ID = new RegExp(`^${UUID}$`, 'i');
const SYNC_GAP_MS = 30000;

/** A trained model to save: what its row describes, the user's feature and target columns a `copy` uploads, and
 * the training run to link. */
export interface TrainedModel {
  name: string;
  description: string;
  tags: string[];
  engine: Engine;
  datasetName: string;
  result: TrainingResult;
  fingerprint: DatasetFingerprint;
  columns: DG.Column[];
  runId: string;
  storage: StorageChoice;
}

/** Fires after a model is saved, edited or deleted. */
export const modelsChanged = new rxjs.Subject<void>();

/** A model without a file (`engine_kind` `external`, written by {@link runRegisterModel}): its method's `apply` gets
 * the row's `options` as the model. */
export function isExternal(row: {engine_kind?: string}): boolean {
  return row.engine_kind === 'external';
}

/** The background sync: [checking] while it runs, [message] what the last one could not refresh, [changed] whether
 * it wrote a row, [finished] once one has run. */
export const syncStatus = new rxjs.BehaviorSubject({checking: false, message: '', changed: false, finished: false});
let syncRun: Promise<void> | null = null;
let syncEndedAt = 0;

/** Runs every function with `meta.mlrole: sync` (it writes the caller's models kept elsewhere through
 * `Forge:registerModel` and `Forge:deleteModel` and returns what it could not refresh), at most once in
 * {@link SYNC_GAP_MS}; callers share a running sync. Never fails. */
export function syncModels(): Promise<void> {
  if (syncRun !== null)
    return syncRun;
  const funcs = DG.Func.find({meta: {mlrole: 'sync'}});
  if (funcs.length === 0 || Date.now() - syncEndedAt < SYNC_GAP_MS)
    return Promise.resolve();
  let changed = false;
  const writes = modelsChanged.subscribe(() => changed = true);
  syncStatus.next({...syncStatus.value, checking: true, changed: false});
  syncRun = Promise.all(funcs.map((func) => func.apply({}).then((message: unknown) =>
    typeof message === 'string' ? message : '', (e: unknown) => `${func.options['mlname']}: ${errorMessage(e)}`)))
    .then((messages) => {
      writes.unsubscribe();
      syncRun = null;
      syncEndedAt = Date.now();
      syncStatus.next({checking: false, message: messages.filter((m) => m !== '').join(' '), changed,
        finished: true});
    });
  return syncRun;
}

/** The query every list of models starts from (the catalog, the Apply dialog, the Models pane, apply by name): the
 * stored rows at once, and {@link syncModels} in the background. */
export function storedModels() {
  void syncModels();
  return forgeDb.models.query();
}

/** The inputs of `Forge:registerModel`; [target], [features] and [options] are the `model` columns' json. */
export interface RegisterModelInput {
  id?: string | null;
  name: string;
  description?: string | null;
  engine: string;
  task: string;
  target: {[key: string]: unknown};
  features: {[key: string]: unknown};
  options: {[key: string]: unknown};
}

/** The body of `Forge:registerModel`: writes the caller's own row of a model without a file, a new one or, with
 * [RegisterModelInput.id], the one registered before; another user's row or a trained model is refused. Returns the
 * row id. */
export async function runRegisterModel(input: RegisterModelInput): Promise<string> {
  const engine = EngineRegistry.discover().find((e) => e.name === input.engine);
  if (engine === undefined)
    throw new ForgeError(`The method '${input.engine}' is not installed.`);
  if (input.task !== 'regression' && input.task !== 'classification')
    throw new ForgeError('The task must be regression or classification.');
  const columns: unknown = input.features?.columns;
  const fields: ModelInsert = {name: input.name, description: input.description ?? '', engine_name: engine.name,
    engine_namespace: engine.namespace, engine_kind: 'external', task: input.task,
    target_name: String(input.target?.name ?? ''), target: input.target, features: input.features,
    feature_count: Array.isArray(columns) ? columns.length : 0, options: input.options, storage_mode: 'none'};
  let id = input.id;
  if (id) {
    // A row the caller cannot see is not found.
    const row: ModelRow | null = await forgeDb.models.get(id);
    if (row?.author_id !== DG.User.current().id || !isExternal(row))
      throw new ForgeError(`The model '${id}' is not a model you registered.`);
    await forgeDb.models.update(id, fields);
  } else
    id = (await forgeDb.models.insert(fields))[0].id;
  modelsChanged.next();
  return id;
}

/** The visible model with the id or the unique name [idOrName], with its [columns] (the system columns come with any
 * query). */
export async function findModelRow<K extends ModelColumn & keyof ModelRow>(idOrName: string, columns: readonly K[]):
  Promise<Pick<ModelRow, K | DG.DomainSystemColumn>> {
  const rows = await storedModels()
    .where(MODEL_ID.test(idOrName) ? 'id' : 'name', '=', idOrName)
    .select(...columns)
    .top(2);
  if (rows.length === 0)
    throw new ForgeError(`No Forge model '${idOrName}' is available to you.`);
  if (rows.length > 1)
    throw new ForgeError(`Several models are named '${idOrName}'. Use the model id.`);
  return rows[0];
}

/** The model folder and the file path of a blob in Forge's own `<uuid>/` folder layout, or null for any other value. */
export function ownBlob(blob: string | undefined): {folder: string; path: string} | null {
  const match = OWN_BLOB.exec(blob ?? '');
  if (match === null)
    return null;
  const folder = `${BLOB_ROOT}/${match[1]}`;
  return {folder, path: `${folder}/${match[2]}`};
}

export async function saveModel(fields: ModelFields, blob: Uint8Array): Promise<string> {
  // The blob goes first, so a failed insert leaves an unreachable file, never a row without its blob.
  const path = `${BLOB_ROOT}/${DG.Utils.uuid4()}/model.bin`;
  await grok.dapi.files.write(path, blob);
  const [{id}] = await forgeDb.models.insert({...fields, blob: FILE_PREFIX + path});
  modelsChanged.next();
  return id;
}

/** Saves [model] and links its training run: `copy` uploads its columns first and deletes them again when the save
 * fails, `reference` keeps the table's origin; nothing is uploaded without `copy`. Returns the model id. */
export async function saveTrainedModel(model: TrainedModel): Promise<string> {
  const {storage: choice, columns, name, result, runId} = model;
  const tableId = choice.mode === 'copy' ? await uploadTrainingCopy(columns, name) : null;
  const storage: ModelStorage = tableId !== null ? {mode: 'copy', tableId} :
    choice.mode === 'reference' ? choice : {mode: 'none'};
  // The row's fields are built inside the chain, so their failure deletes the copy too; the caller gets the save's own
  // error, even when the copy cannot be deleted.
  const id = await Promise.resolve().then(() => saveModel(modelFieldsOf({...model, storage}), result.blob))
    .catch(async (e) => {
      if (tableId !== null)
        await deleteTrainingCopy(tableId).catch(() => undefined);
      throw e;
    });
  await linkTrainingRun(runId, id);
  return id;
}

/** The body of `Forge:deleteModel`: deletes the model with the id or the unique name [idOrName] as {@link deleteModel}
 * does, with one read of the row; returns its id. */
export async function runDeleteModel(idOrName: string): Promise<string> {
  const model = await findModelRow(idOrName, ['blob', 'dataset_table_id']);
  await deleteModelRow(model);
  return model.id;
}

/** Deletes the row, then its blob folder and its uploaded training copy, each tried; the first failure to delete
 * those is thrown after the row is gone. */
export async function deleteModel(id: string): Promise<void> {
  await deleteModelRow(await forgeDb.models.get(id));
}

async function deleteModelRow(model: Pick<ModelRow, 'id' | 'blob' | 'dataset_table_id'>): Promise<void> {
  await forgeDb.models.delete(model.id);
  // Only a blob in Forge's own per-model folder is deleted; any other file stays.
  const own = ownBlob(model.blob);
  const tableId = model.dataset_table_id;
  const outcomes = await Promise.allSettled([own === null ? undefined : grok.dapi.files.delete(own.folder),
    tableId ? deleteTrainingCopy(tableId) : undefined]);
  modelsChanged.next();
  const failed = outcomes.find((o): o is PromiseRejectedResult => o.status === 'rejected');
  if (failed !== undefined)
    throw failed.reason;
}
