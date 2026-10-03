import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {reportError} from './report-error';

export function saveModelDialog(defaultName: string,
  onSave: (name: string, description: string) => Promise<void>): DG.Dialog {
  const nameInput = ui.input.string('Name', {value: defaultName});
  const descriptionInput = ui.input.textArea('Description', {value: ''});
  return ui.dialog('Save model')
    .add(nameInput)
    .add(descriptionInput)
    .onOK(async () => {
      try {
        await onSave(nameInput.value.trim() || defaultName, descriptionInput.value);
      } catch (e) {
        reportError(e);
      }
    });
}
