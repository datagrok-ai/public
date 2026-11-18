import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {onFrame} from '../preparation/shared-frame';

const COPY_SUFFIX = ' (training data)';

/** The name of the uploaded copy of a model's training data. */
export function trainingCopyName(modelName: string): string {
  return `${modelName}${COPY_SUFFIX}`;
}

/** Uploads [columns] (the features and the target, every row) as a table named after the model; returns its id. */
export async function uploadTrainingCopy(columns: DG.Column[], modelName: string): Promise<string> {
  return await onFrame(columns, (frame) => {
    frame.name = trainingCopyName(modelName);
    return grok.dapi.tables.uploadDataFrame(frame);
  });
}

/** Deletes an uploaded copy. A table that is gone already is no error; a table not named as a copy (a
 * `dataset_table_id` edited elsewhere) is left alone. The server keeps the given name as `friendlyName`. */
export async function deleteTrainingCopy(id: string): Promise<void> {
  const info: DG.TableInfo | undefined = await grok.dapi.tables.find(id);
  if (info?.friendlyName?.endsWith(COPY_SUFFIX))
    await grok.dapi.tables.delete(info);
}
