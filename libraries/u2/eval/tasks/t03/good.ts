import * as DG from 'datagrok-api/dg';
import {Dialog, Form} from '@datagrok-libraries/u2';
import {FileInput, fileInput} from '@datagrok-libraries/u2/dg';

export function mergeDialog(onMerge: (file: DG.FileInfo) => unknown): Dialog {
  const dialog = new Dialog('Merge with table');
  let source!: FileInput;
  dialog.add(() => new Form().add(source = fileInput('Source file', {required: true})));
  dialog.onOK(() => onMerge(source.peek()!), {text: 'MERGE'});
  return dialog.show();
}
