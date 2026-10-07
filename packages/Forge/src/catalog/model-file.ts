import * as grok from 'datagrok-api/grok';
import {ForgeError} from '../forge-error';
import {ModelRow} from '../generated/db';
import {ownBlob} from '../storage/model-store';

const FILE_NAME_UNSAFE = /[\\/:*?"<>|]/g;

export interface ModelFile { name: string; bytes: Uint8Array }

/** The model's own file, as `<model name>.bin` with the characters a file name cannot hold replaced by `_`. */
export async function readModelFile(row: Pick<ModelRow, 'name' | 'blob'>): Promise<ModelFile> {
  const blob = ownBlob(row.blob);
  if (blob === null) {
    throw new ForgeError(`The model '${row.name}' has no model file in Forge's storage, ` +
      'so there is nothing to download.');
  }
  return {name: `${row.name.replace(FILE_NAME_UNSAFE, '_')}.bin`, bytes: await grok.dapi.files.readAsBytes(blob.path)};
}
