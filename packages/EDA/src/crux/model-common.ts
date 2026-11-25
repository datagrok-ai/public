import * as DG from 'datagrok-api/dg';

import * as sha256 from 'fast-sha256';
import {Table, tableFromIPC} from 'apache-arrow';

import {cruxMl, CruxResource} from './crux-loader';

export const LEGACY_MODEL_ERROR =
  'This model was saved by an older EDA version and is no longer supported - retrain the model';
export const PREDICTION = 'Prediction';

/** Runs a Crux function on an Arrow IPC stream (transferred) and decodes its output */
export async function runCrux(name: string, ipc: Uint8Array, params: object, resource?: CruxResource):
  Promise<Table> {
  return tableFromIPC(await cruxMl.callFn(name, ipc, params, resource));
}

/** The columns of `df` named `names`, in that order */
export function featureColumns(df: DG.DataFrame, names: string[]): DG.Column[] {
  const columns: DG.Column[] = [];
  for (const name of names) {
    const col = df.col(name);
    if (col === null)
      throw new Error(`Crux: feature column "${name}" is missing`);
    columns.push(col);
  }
  return columns;
}

/** SHA-256 hex of a blob */
export function blobId(bytes: Uint8Array): string {
  return Array.from(sha256.hash(bytes), (b) => b.toString(16).padStart(2, '0')).join('');
}
