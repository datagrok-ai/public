import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {Engine} from '../engines/engine';
import {errorMessage, ForgeError} from '../forge-error';
import {forgeDb, ModelColumn, ModelRow} from '../generated/db';
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

export function isExternal(row: Pick<ModelRow, 'engine_kind'>): boolean {
  return row.engine_kind === 'external';
}

export interface SyncStatus {
  checking: boolean;
  /** What the last sync could not refresh. */
  message: string;
  /** The last sync wrote a row. */
  changed: boolean;
  finished: boolean;
}

export const syncStatus = new rxjs.BehaviorSubject<SyncStatus>(
  {checking: false, message: '', changed: false, finished: false});
let syncRun: Promise<void> | null = null;
let syncEndedAt = 0;

/** Runs the `meta.mlrole: sync` functions, at most once in {@link SYNC_GAP_MS}; never fails. */
export function syncModels(): Promise<void> {
  if (syncRun !== null)
    return syncRun;
  const funcs = DG.Func.find({meta: {mlrole: 'sync'}});
  if (funcs.length === 0 || Date.now() - syncEndedAt < SYNC_GAP_MS)
    return Promise.resolve();
  return syncRun = runSync(funcs);
}

async function runSync(funcs: DG.Func[]): Promise<void> {
  let changed = false;
  const writes = modelsChanged.subscribe(() => changed = true);
  syncStatus.next({...syncStatus.value, checking: true, changed: false});
  const messages: string[] = await Promise.all(funcs.map((func) => func.apply({}).then(
    (message: unknown) => typeof message === 'string' ? message : '',
    (e: unknown) => `${func.options['mlname']}: ${errorMessage(e)}`)));
  writes.unsubscribe();
  syncRun = null;
  syncEndedAt = Date.now();
  syncStatus.next({checking: false, message: messages.filter((m) => m !== '').join(' '), changed, finished: true});
}

/** Every list of models starts here: the stored rows at once, a sync in the background. */
export function storedModels() {
  void syncModels();
  return forgeDb.models.query();
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

async function writeModelFile(bytes: Uint8Array): Promise<string> {
  const path = `${BLOB_ROOT}/${DG.Utils.uuid4()}/model.bin`;
  await grok.dapi.files.write(path, bytes);
  return FILE_PREFIX + path;
}

export async function saveModel(fields: ModelFields, blob: Uint8Array): Promise<string> {
  // The blob goes first, so a failed insert leaves an unreachable file, never a row without its blob.
  const [{id}] = await forgeDb.models.insert({...fields, blob: await writeModelFile(blob)});
  modelsChanged.next();
  return id;
}

export async function adoptLegacyFile(row: Pick<ModelRow, 'id' | 'name' | 'legacy_id'>): Promise<Uint8Array> {
  const file: DG.FileInfo = await grok.functions.call('getModelBlob', {modelId: row.legacy_id}).catch((e) => {
    throw new ForgeError(`The file of the model '${row.name}' cannot be read: ${errorMessage(e)}`);
  });
  const bytes = file.data;
  if (!bytes?.length)
    throw new ForgeError(`The file of the model '${row.name}' cannot be read: the file is empty`);
  const blob = await writeModelFile(bytes);
  try {
    await forgeDb.models.update(row.id, {blob});
  } catch (e) {
    await grok.dapi.files.delete(ownBlob(blob)!.folder);
    // A user who can only view the model applies the bytes.
    if (e instanceof DG.DomainForbiddenError)
      return bytes;
    throw e;
  }
  modelsChanged.next();
  return bytes;
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
