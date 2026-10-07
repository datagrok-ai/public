import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelInfo, updateModelInfo} from '../catalog/model-edit';
import {ModelRow} from '../generated/db';
import {ButtonGate} from './button-gate';
import {reportError} from './report-error';
import {tagsInput, tagsOfInput} from './tags-input';

const NO_NAME = 'Enter a name for the model.';

export function saveModelDialog(defaultName: string, onSave: (info: ModelInfo) => Promise<void>): DG.Dialog {
  return modelInfoDialog('Save model', {name: defaultName, description: '', tags: []}, onSave);
}

/** Writes [info] over the model's [version]; a model changed meanwhile asks with the platform's conflict dialog to
 * overwrite it (written anyway) or to reload ([reload] runs). The new version, or null when nothing was written. */
export async function writeModelInfo(model: Pick<ModelRow, 'id' | 'name'>, version: number, info: Partial<ModelInfo>,
  reload: () => Promise<void>): Promise<number | null> {
  try {
    return await updateModelInfo(model.id, version, info);
  } catch (e) {
    if (!(e instanceof DG.DomainVersionConflictError))
      throw e;
    const choice = await DG.DomainObjectHandler.showConflictDialog(model.name);
    if (choice === 'overwrite')
      return updateModelInfo(model.id, undefined, info);
    if (choice === 'reload')
      await reload();
    return null;
  }
}

/** **Name**, **Description** and **Tags** prefilled from [info]; OK, unavailable while the name is blank, calls [onOK]
 * with the trimmed name. The Save model and Edit model dialogs. */
export function modelInfoDialog(title: string, info: ModelInfo, onOK: (info: ModelInfo) => Promise<void>): DG.Dialog {
  const problem = () => nameInput.value.trim() === '' ? NO_NAME : null;
  // Not nullable: the dialog's own validation then refuses an empty Name, so Enter cannot save it either.
  const nameInput = ui.input.string('Name', {value: info.name, nullable: false,
    onValueChanged: () => okGate.update()});
  nameInput.addValidator(problem);
  const descriptionInput = ui.input.textArea('Description', {value: info.description});
  const tags = tagsInput(info.tags);
  const dialog = ui.dialog(title)
    .add(nameInput)
    .add(descriptionInput)
    .add(tags)
    .onOK(async () => {
      try {
        await onOK({name: nameInput.value.trim(), description: descriptionInput.value, tags: tagsOfInput(tags)});
      } catch (e) {
        reportError(e);
      }
    });
  const okGate = new ButtonGate(dialog.getButton('OK'), problem);
  okGate.update();
  return dialog;
}
