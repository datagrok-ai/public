import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ForgeError} from '../forge-error';
import {forgeDb, ModelRow} from '../generated/db';
import {tagsOf} from '../storage/model-fields';
import {modelInfoDialog, writeModelInfo} from './save-model-dialog';

const EDIT_COLUMNS = ['name', 'description', 'tags'] as const;
type EditedModel = Pick<ModelRow, typeof EDIT_COLUMNS[number] | 'id' | 'version'>;

/** Shows **Edit model** for the model's current row. */
export async function openEditModelDialog(id: string): Promise<void> {
  const model = await forgeDb.models.query().where('id', '=', id).select(...EDIT_COLUMNS).first();
  if (model === null)
    throw new ForgeError('The model is no longer available.');
  editModelDialog(model).show();
}

/** **Edit model**: Name, Description and Tags of [model]; a model changed meanwhile asks to reload or overwrite. The
 * context panel showing the model follows the write (`ForgeModelHandler.followChanges`). */
export function editModelDialog(model: EditedModel): DG.Dialog {
  return modelInfoDialog('Edit model', {name: model.name, description: model.description ?? '',
    tags: tagsOf(model.tags)}, async (info) => {
    if (await writeModelInfo(model, model.version, info, () => openEditModelDialog(model.id)) !== null)
      grok.shell.info(`Model "${info.name}" updated.`);
  });
}
