import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import {ButtonGate} from './button-gate';
import {reportError} from './report-error';
import {TrainView} from './train-view';

/** Opens the Train view with [method] chosen (the suggested one for null): on the current table, else the only open
 * one, else on the table a dialog asks for. */
export async function launchTraining(method: string | null): Promise<void> {
  try {
    const table = grok.shell.currentTable ?? (grok.shell.tables.length === 1 ? grok.shell.tables[0] : null);
    if (table !== null) {
      await TrainView.openFor(table, method);
      return;
    }
    const tableInput = ui.input.table('Table', {value: grok.shell.tables[0] ?? null, nullable: false,
      tooltipText: 'The table to learn from.'});
    const dialog = ui.dialog(method === null ? 'Train a model' : `Train with ${method}`).add(tableInput);
    dialog.onOK(() => {
      const chosen = tableInput.value;
      if (chosen !== null)
        void TrainView.openFor(chosen, method);
    });
    const gate = new ButtonGate(dialog.getButton('OK'), () => tableInput.value === null ? 'Choose a table.' : null);
    tableInput.onChanged.subscribe(() => gate.update());
    gate.update();
    dialog.show();
  } catch (e) {
    reportError(e);
  }
}
