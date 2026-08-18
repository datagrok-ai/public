/* Choice inputs over platform collections. The value is always the name, not the object: the
   caller resolves it (`table.col(input.value.peek())`, `grok.shell.table(...)`), which keeps the
   inputs plain `ChoiceInput`s that serialize, bind and validate like any other. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {ChoiceInput} from '../components/choice-input.js';

export interface ColumnInputOptions {
  filter?: (column: DG.Column) => boolean;
}

/** Picks a column of `table` by name; the item list follows `table.onColumnsChanged` until the
 * input is disposed. A dropped column clears the value. */
export function columnInput(label: string, table: DG.DataFrame, options: ColumnInputOptions = {}): ChoiceInput {
  const items = () => columnNames(table, options.filter);
  const input = new ChoiceInput({label, items: items()});
  const sub = table.onColumnsChanged.subscribe(() => input.setItems(items()));
  input.own(() => sub.unsubscribe());
  return input;
}

/** Picks an open table by name; the item list follows `onTableAdded`/`onTableRemoved` until the
 * input is disposed. A closed table clears the value. */
export function tableInput(label: string): ChoiceInput {
  const input = new ChoiceInput({label, items: grok.shell.tableNames});
  const subs = [grok.events.onTableAdded, grok.events.onTableRemoved]
    .map((event) => event.subscribe(() => input.setItems(grok.shell.tableNames)));
  input.own(() => {
    for (const sub of subs)
      sub.unsubscribe();
  });
  return input;
}

function columnNames(table: DG.DataFrame, filter?: (column: DG.Column) => boolean): string[] {
  return filter ? table.columns.toList().filter(filter).map((c) => c.name) : table.columns.names();
}
