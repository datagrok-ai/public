import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {forgeDb} from '../generated/db';
import {deleteTrainingCopy} from './dataset-copy';
import {ModelFields} from './model-fields';

export const BLOB_ROOT = 'System:DomainFiles/forge/model';
export const UUID = '[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}';
const FILE_PREFIX = 'file://';
const OWN_BLOB = new RegExp(`^${FILE_PREFIX}${BLOB_ROOT}/(${UUID})/([^/]+)$`, 'i');

/** Fires after a model is saved, edited or deleted. */
export const modelsChanged = new rxjs.Subject<void>();

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

/** Deletes the row, then its blob folder and its uploaded training copy, each tried; the first failure to delete
 * those is thrown after the row is gone. */
export async function deleteModel(id: string): Promise<void> {
  const model = await forgeDb.models.get(id);
  await forgeDb.models.delete(id);
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
