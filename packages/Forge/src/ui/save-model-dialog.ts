import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ModelInfo, updateModelInfo} from '../catalog/model-edit';
import {ModelRow, ModelStorageMode} from '../generated/db';
import {DatasetRef, datasetRefCaption} from '../storage/dataset-ref';
import {ButtonGate} from './button-gate';
import {reportError} from './report-error';
import {tagsInput, tagsOfInput} from './tags-input';

const NO_NAME = 'Enter a name for the model.';
export const STORAGE_CAPTIONS: Record<ModelStorageMode, string> = {reference: 'Reference', none: 'None', copy: 'Copy'};

/** What the Save dialog offers to keep of the training data: [ref] the table's origin, when the platform recorded one;
 * [rowCount] the rows a copy uploads. */
export interface StorageOffer { ref: DatasetRef | null; rowCount: number }
/** The **Data storage** the user chose; `reference` carries the offered origin. */
export type StorageChoice = {mode: 'none' | 'copy'} | {mode: 'reference'; ref: DatasetRef};

/** **Save model**: [modelInfoDialog] and, with [storage], **Data storage**; OK calls [onSave] with the chosen storage
 * (`none` without [storage]). */
export function saveModelDialog(defaultName: string, onSave: (info: ModelInfo, choice: StorageChoice) => Promise<void>,
  storage?: StorageOffer): DG.Dialog {
  let choiceOf = (): StorageChoice => ({mode: 'none'});
  const dialog = modelInfoDialog('Save model', {name: defaultName, description: '', tags: []},
    (info) => onSave(info, choiceOf()));
  if (storage !== undefined)
    choiceOf = addStorage(dialog, storage);
  return dialog;
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

/** Adds **Data storage** and its lines to [dialog]; returns the chosen storage's reader. */
function addStorage(dialog: DG.Dialog, {ref, rowCount}: StorageOffer): () => StorageChoice {
  const choices: StorageChoice[] = ref === null ? [{mode: 'none'}, {mode: 'copy'}] :
    [{mode: 'reference', ref}, {mode: 'none'}, {mode: 'copy'}];
  const choiceOf = (caption: string | null) => choices.find((c) => STORAGE_CAPTIONS[c.mode] === caption) ?? choices[0];
  const lines: Record<ModelStorageMode, string> = {
    reference: ref === null ? '' : `A link to ${datasetRefCaption(ref)} is saved; the data stays where it is.`,
    none: 'Only a summary of the data is saved.',
    copy: `The training columns (${rowCount} rows) will be uploaded to the server.`,
  };
  const line = ui.divText(lines[choices[0].mode], 'forge-note');
  const input = ui.input.radio('Data storage', {items: choices.map((c) => STORAGE_CAPTIONS[c.mode]),
    value: STORAGE_CAPTIONS[choices[0].mode], tooltipText: 'What the model keeps of its training data.',
    onValueChanged: (caption) => line.textContent = lines[choiceOf(caption).mode]});
  input.root.classList.add('forge-inline-radio');
  dialog.add(input).add(line);
  return () => choiceOf(input.value);
}
