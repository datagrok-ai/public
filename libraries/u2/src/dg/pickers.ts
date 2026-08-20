/* Choice inputs over platform collections. The value is always the name, not the object: the
   caller resolves it (`table.col(input.value.peek())`, `grok.shell.table(...)`), which keeps the
   inputs plain `ChoiceInput`s that serialize, bind and validate like any other. */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Input, InputOptions} from '../core/input-base.js';
import {Control} from '../core/component.js';
import {ChoiceInput, MultiChoiceInput} from '../components/choice-input.js';
import {TypeAhead} from '../components/typeahead.js';
import {iconButton} from '../components/buttons.js';
import {columnRenderer} from './column-renderer.js';

export interface ColumnInputOptions {
  filter?: (column: DG.Column) => boolean;
  /** The searchable variant: a type-ahead over the column rows (glyph, name, semantic type)
   * instead of the native select — Dart's ColumnComboBox affordances (`column_input.dart:13-17`).
   * Off by default, so the plain choice input stays what `columnInput` returns. */
  rich?: boolean;
  placeholder?: string;
}

/** A rename fires `onColumnNameChanged` and NOTHING else — `onColumnsChanged` covers additions and
 * removals only (`column.dart:68`). Every name-valued picker follows this one too and REMAPS the
 * name it holds, or it would silently point at a column that no longer answers to it. The payload
 * is the platform's `NameChangedArgs` (`data_events.dart:47-57`) behind its `EventData` wrapper. */
export function onColumnRenamed(table: DG.DataFrame,
  handler: (oldName: string, newName: string) => void): {unsubscribe(): void} {
  return table.onColumnNameChanged.subscribe((e: {args: {oldName: string, newName: string}}) =>
    handler(e.args.oldName, e.args.newName));
}

/** Picks a column of `table` by name; the item list follows `table.onColumnsChanged` until the
 * input is disposed. A dropped column clears the value, a renamed one carries it over. */
export function columnInput(label: string, table: DG.DataFrame,
  options: ColumnInputOptions & {rich: true}): ColumnPicker;
export function columnInput(label: string, table: DG.DataFrame,
  options?: ColumnInputOptions): ChoiceInput;
export function columnInput(label: string, table: DG.DataFrame,
  options: ColumnInputOptions = {}): ChoiceInput | ColumnPicker {
  if (options.rich)
    return new ColumnPicker({label, table, filter: options.filter, placeholder: options.placeholder});
  const items = () => columnNames(table, options.filter);
  const input = new ChoiceInput({label, items: items()});
  const subs = [
    table.onColumnsChanged.subscribe(() => input.setItems(items())),
    // the value moves first: setItems drops a value the new item list no longer offers
    onColumnRenamed(table, (oldName, newName) => {
      if (input.value.peek() === oldName)
        Input.system(() => input.value.value = newName);
      input.setItems(items());
    }),
  ];
  input.own(() => {
    for (const sub of subs)
      sub.unsubscribe();
  });
  return input;
}

/** Picks an open table by name; the item list follows `onTableAdded`/`onTableRemoved` until the
 * input is disposed. A closed table clears the value. The options rail carries the platform
 * input's import action (`table_input.dart:37-58`): the file it opens joins the workspace and
 * becomes the value. */
export function tableInput(label: string): ChoiceInput {
  const input = new ChoiceInput({label, items: grok.shell.tableNames});
  followOpenTables(input, (names) => input.setItems(names));
  importAction(input, (name) => {
    input.setItems(grok.shell.tableNames);
    input.value.value = name;
  });
  return input;
}

/** Picks any number of open tables — a checkbox per table, the value is their names (Dart's
 * `TablesMultiChoiceInput`, `tables_multi_choice_input.dart:48-64`). Checked tables survive an
 * open or close of another one; a closed table leaves the value with it. The import action is the
 * platform input's own, where it is a hidden-selector `TableInput` sitting under the checks
 * (`tables_multi_choice_input.dart:66-77`): the imported table joins the workspace checked. */
export function tablesInput(label: string): MultiChoiceInput {
  const input = new MultiChoiceInput({label, items: grok.shell.tableNames, emptyText: 'No open tables'});
  followOpenTables(input, (names) => input.setItems(names));
  importAction(input, (name) => {
    input.setItems(grok.shell.tableNames);
    const checked = input.value.peek();
    input.value.value = checked.includes(name) ? checked : [...checked, name];
  });
  return input;
}

export interface ColumnPickerOptions extends InputOptions<string | null> {
  table: DG.DataFrame;
  filter?: (column: DG.Column) => boolean;
  placeholder?: string;
}

/** The searchable `columnInput` variant: a {@link TypeAhead} over the live column names rendered
 * through {@link columnRenderer}. The value stays the column name, so it is interchangeable with
 * the plain choice variant everywhere. */
export class ColumnPicker extends Input<string | null, ColumnPickerOptions> {
  private _typeAhead!: TypeAhead<string>;

  constructor(options: ColumnPickerOptions) {
    super(options, null);
    this.root.dataset.u2 = 'column-picker';
  }

  protected createEditor(): HTMLElement {
    const table = this.options.table;
    const names = () => columnNames(table, this.options.filter);
    const typeAhead = new TypeAhead<string>({
      source: (query: string) => Promise.resolve(match(names(), query)),
      renderer: columnRenderer(table),
      placeholder: this.options.placeholder ?? 'Column…',
      debounceMs: 0,
    });
    this._typeAhead = typeAhead;

    // seeded before the effects, or the first one would push the type-ahead's empty selection
    // over the value the input was constructed with
    typeAhead.selected.value = this.value.peek();
    this.effect(() => {
      const picked = typeAhead.selected.value;
      if (picked !== this.value.peek())
        this.value.value = picked;
    });
    this.effect(() => {
      const value = this.value.value;
      if (value !== typeAhead.selected.peek())
        typeAhead.selected.value = value;
    });

    const subs = [
      table.onColumnsChanged.subscribe(() => {
        const value = this.value.peek();
        if (value !== null && !names().includes(value))
          Input.system(() => this.value.value = null);
      }),
      onColumnRenamed(table, (oldName, newName) => {
        if (this.value.peek() === oldName)
          Input.system(() => this.value.value = newName);
      }),
    ];
    this.own(() => {
      for (const sub of subs)
        sub.unsubscribe();
    });
    return typeAhead.root;
  }
}

function columnNames(table: DG.DataFrame, filter?: (column: DG.Column) => boolean): string[] {
  return filter ? table.columns.toList().filter(filter).map((c) => c.name) : table.columns.names();
}

function match(names: string[], query: string): string[] {
  const q = query.trim().toLowerCase();
  return q ? names.filter((name) => name.toLowerCase().includes(q)) : names;
}

function followOpenTables(input: Control, apply: (names: string[]) => void): void {
  const subs = [grok.events.onTableAdded, grok.events.onTableRemoved]
    .map((event) => event.subscribe(() => apply(grok.shell.tableNames)));
  input.own(() => {
    for (const sub of subs)
      sub.unsubscribe();
  });
}

/** What the import action opens without the platform's own `FileHandler`, which js-api does not
 * expose: text through `grok.data.parseCsv`, a d42 blob through `DG.DataFrame.fromByteArray`. */
const TEXT_FORMATS = ['csv', 'tsv', 'txt'];
const BINARY_FORMAT = 'd42';

/** The platform input's `folder-open` icon (`table_input.dart:37-58`): it opens a local file, adds
 * the table to the workspace and hands its name over — as a user edit, the way a pick from the
 * list is one. */
function importAction(input: Input<any, any>, pick: (name: string) => void): void {
  const picker = document.createElement('input');
  picker.type = 'file';
  picker.hidden = true;
  const onChange = async () => {
    const file = picker.files?.[0];
    picker.value = '';
    if (file === undefined)
      return;
    const name = await openTable(file);
    if (name !== null && !input.scope.isDisposed)
      pick(name);
  };
  picker.addEventListener('change', onChange);
  input.own(() => picker.removeEventListener('change', onChange));
  const open = input.run(() => iconButton('folder-open', () => picker.click(), {tooltip: 'Open file'}));
  optionsRail(input).append(picker, open);
}

/** The input's options rail. `Input.addOptions` materializes it with its first child and is
 * protected, so a picker built by a factory function puts the same rail there itself. */
function optionsRail(input: Input<any, any>): HTMLElement {
  const existing = input.box.querySelector('.u2-input-options') as HTMLElement | null;
  if (existing !== null)
    return existing;
  const rail = document.createElement('div');
  rail.className = 'u2-input-options';
  input.box.append(rail);
  return rail;
}

/** A local file as the workspace's newest table, and the name it ended up with. An extension
 * neither reader takes is refused with the platform input's own message
 * (`table_input.dart:41-43`) — Dart hands those to `FileHandler`, which opens xlsx, sdf and the
 * rest, and js-api has no door to it. */
async function openTable(file: File): Promise<string | null> {
  const dot = file.name.lastIndexOf('.');
  const ext = dot < 0 ? '' : file.name.substring(dot + 1).toLowerCase();
  if (ext !== BINARY_FORMAT && !TEXT_FORMATS.includes(ext)) {
    grok.shell.error(`File extension .${ext} is not supported.`);
    return null;
  }
  try {
    const table = ext === BINARY_FORMAT ?
      DG.DataFrame.fromByteArray(new Uint8Array(await file.arrayBuffer())) :
      grok.data.parseCsv(await file.text());
    table.name = dot < 0 ? file.name : file.name.substring(0, dot);
    return grok.shell.addTable(table).name;
  } catch (e) {
    grok.shell.error(`${file.name} could not be opened.`);
    return null;
  }
}
