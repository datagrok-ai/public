import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ChoiceInput, Dialog, Form, Input, NumberInput, TextInput} from '@datagrok-libraries/u2';

export function settingsDialog(func: DG.Func, onRun: (call: DG.FuncCall) => unknown): Dialog {
  const dialog = new Dialog(`${func.friendlyName ?? func.name} settings`);
  const fields: [DG.Property, Input<any>][] = [];
  dialog.add(() => {
    const form = new Form();
    for (const p of func.inputs) {
      const label = p.caption ?? p.friendlyName ?? p.name;
      const input: Input<any> = p.propertyType === DG.TYPE.DATA_FRAME ?
        new ChoiceInput({label, items: grok.shell.tableNames, value: grok.shell.tableNames[0]}) :
        p.propertyType === DG.TYPE.FLOAT || p.propertyType === DG.TYPE.INT ?
          new NumberInput({label, value: p.defaultValue ?? null}) :
          new TextInput({label, value: p.defaultValue ?? ''});
      fields.push([p, input]);
      form.add(input);
    }
    return form;
  });
  dialog.onOK(() => {
    const values: Record<string, unknown> = {};
    for (const [p, input] of fields) {
      const value = input.peek();
      values[p.name] = p.propertyType === DG.TYPE.DATA_FRAME ? grok.shell.table(value) : value;
    }
    onRun(func.prepare(values));
  }, {text: 'RUN'});
  return dialog.show();
}
