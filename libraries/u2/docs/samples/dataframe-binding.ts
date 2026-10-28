import * as DG from 'datagrok-api/dg';
import {registerAll, renderSpec} from '@datagrok-libraries/u2';
import type {Spec} from '@datagrok-libraries/u2';
import {ViewHost, platformContext} from '@datagrok-libraries/u2/dg';

const rowEditor = (table: string): Spec => ({
  $schema: 'dg-ui/1',
  components: [{tag: 'u2-table-source', name: 'table', props: {table}}],
  root: {tag: 'u2-form', children: [
    {tag: 'u2-number-input', props: {label: 'Row', mode: 'int', min: 0}, bind: {value: '$.table.currentRowIdx'}},
    {tag: 'u2-number-input', props: {label: 'Age', mode: 'int', min: 0}, bind: {value: '$.table.currentRow.AGE'}},
    {tag: 'u2-choice-input', props: {label: 'Sex', items: ['M', 'F']}, bind: {value: '$.table.currentRow.SEX'}},
    {tag: 'u2-number-input', props: {label: 'Weight', postfix: 'kg'}, bind: {value: '$.table.currentRow.WEIGHT'}},
  ]},
});

export function editCurrentRow(view: DG.TableView): void {
  registerAll();
  const host = ViewHost.of(view);
  const editor = host.build(() => renderSpec(rowEditor(view.dataFrame.name), platformContext()));
  host.dock(editor, 'right', {title: 'Current row', ratio: 0.3});
}
