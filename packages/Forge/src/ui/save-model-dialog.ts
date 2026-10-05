import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import {ButtonGate} from './button-gate';
import {reportError} from './report-error';

const NO_NAME = 'Enter a name for the model.';

export function saveModelDialog(defaultName: string,
  onSave: (name: string, description: string) => Promise<void>): DG.Dialog {
  const problem = () => nameInput.value.trim() === '' ? NO_NAME : null;
  // Not nullable: the dialog's own validation then refuses an empty Name, so Enter cannot save it either.
  const nameInput = ui.input.string('Name', {value: defaultName, nullable: false,
    onValueChanged: () => okGate.update()});
  nameInput.addValidator(problem);
  const descriptionInput = ui.input.textArea('Description', {value: ''});
  const dialog = ui.dialog('Save model')
    .add(nameInput)
    .add(descriptionInput)
    .onOK(async () => {
      try {
        await onSave(nameInput.value.trim(), descriptionInput.value);
      } catch (e) {
        reportError(e);
      }
    });
  const okGate = new ButtonGate(dialog.getButton('OK'), problem);
  okGate.update();
  return dialog;
}
