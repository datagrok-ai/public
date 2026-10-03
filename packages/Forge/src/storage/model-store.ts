import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {forgeDb} from '../generated/db';
import {ModelFields} from './model-fields';

export const BLOB_ROOT = 'System:DomainFiles/forge/model';
const FILE_PREFIX = 'file://';
const UUID = '[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}';
const OWN_BLOB = new RegExp(`^${FILE_PREFIX}${BLOB_ROOT}/(${UUID})/[^/]+$`, 'i');

/** Fires after a model is saved or deleted. */
export const modelsChanged = new rxjs.Subject<void>();

export async function saveModel(fields: ModelFields, blob: Uint8Array): Promise<string> {
  // The blob goes first, so a failed insert leaves an unreachable file, never a row without its blob.
  const path = `${BLOB_ROOT}/${DG.Utils.uuid4()}/model.bin`;
  await grok.dapi.files.write(path, blob);
  const [{id}] = await forgeDb.models.insert({...fields, blob: FILE_PREFIX + path});
  modelsChanged.next();
  return id;
}

export async function deleteModel(id: string): Promise<void> {
  const model = await forgeDb.models.get(id);
  await forgeDb.models.delete(id);
  // Only a blob in Forge's own per-model folder is deleted; any other file stays.
  const match = OWN_BLOB.exec(model.blob ?? '');
  if (match !== null)
    await grok.dapi.files.delete(`${BLOB_ROOT}/${match[1]}`);
  modelsChanged.next();
}
