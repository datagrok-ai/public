import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ModelInfo} from '../catalog/model-edit';
import {ForgeError} from '../forge-error';
import {forgeDb, ModelRow} from '../generated/db';
import {tagsOf, tagsText} from '../storage/model-fields';
import {ForgeModelHandler, forgeModelHandler} from './model-handler';
import {modelInfoDialog, writeModelInfo} from './save-model-dialog';

// blob too: the context panel's commands of the edited model read it (Download).
const EDIT_COLUMNS = ['name', 'description', 'tags', 'blob'] as const;
type EditedModel = Pick<ModelRow, typeof EDIT_COLUMNS[number] | 'id' | 'version'> & Partial<ModelRow>;

/** Shows **Edit model** for the model's current row. */
export async function openEditModelDialog(id: string): Promise<void> {
  const model = await forgeDb.models.query().where('id', '=', id).select(...EDIT_COLUMNS).first();
  if (model === null)
    throw new ForgeError('The model is no longer available.');
  editModelDialog(model).show();
}

/** **Edit model**: Name, Description and Tags of [model]; a model changed meanwhile asks to reload or overwrite. */
export function editModelDialog(model: EditedModel): DG.Dialog {
  return modelInfoDialog('Edit model', {name: model.name, description: model.description ?? '',
    tags: tagsOf(model.tags)}, async (info) => {
    if (await writeModelInfo(model, model.version, info, () => openEditModelDialog(model.id)) === null)
      return;
    grok.shell.info(`Model "${info.name}" updated.`);
    showUpdated(model, info);
  });
}

/** The context panel showing the edited model shows it again: its title and Details have changed. */
function showUpdated(model: EditedModel, info: ModelInfo): void {
  if (forgeModelHandler.modelOf(grok.shell.o)?.id !== model.id)
    return;
  // Forced: the platform takes a row of the same id for the object it already shows.
  grok.shell.setCurrentObject(ForgeModelHandler.rowOf({...model, ...info, tags: tagsText(info.tags)}), true, true);
}
