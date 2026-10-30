import * as DG from 'datagrok-api/dg';
import {Dialog, Form} from '@datagrok-libraries/u2';
import {FilesInput, filesInput} from '@datagrok-libraries/u2/dg';

export function mergeDialog(onMerge: (file: DG.FileInfo) => unknown): Dialog {
  const dialog = new Dialog('Merge with table');
  let source!: FilesInput;
  dialog.add(() => new Form().add(source = filesInput('Source file', {required: true})));
  dialog.onOK(() => onMerge(source.peek()[0]), {text: 'MERGE'});
  return dialog.show();
}
