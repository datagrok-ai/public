import * as DG from 'datagrok-api/dg';
import {Dialog} from '@datagrok-libraries/u2';
import {funcForm} from '@datagrok-libraries/u2/dg';

export function settingsDialog(func: DG.Func, onRun: (call: DG.FuncCall) => unknown): Dialog {
  const call = func.prepare();
  const dialog = new Dialog(`${func.friendlyName ?? func.name} settings`);
  dialog.add(() => funcForm(call));
  dialog.onOK(() => onRun(call), {text: 'RUN'});
  return dialog.show();
}
