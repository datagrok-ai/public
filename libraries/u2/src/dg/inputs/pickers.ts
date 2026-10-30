/* Choice inputs over platform collections. The value is always the name, not the object, which keeps the
   inputs plain `ChoiceInput`s that serialize, bind and validate like any other; `TableInput.table` holds the
   frame itself. What they list is what the platform's pickers list (`Platform`). */
import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Input, InputOptions, LiveOption} from '../../core/input-base.js';
import {Control} from '../../core/component.js';
import {Scope, Usage} from '../../core/scope.js';
import {ChoiceInput, MultiChoiceInput} from '../../components/inputs/choice-input.js';
import type {ChoiceInputOptions} from '../../components/inputs/choice-input.js';
import {TypeAhead} from '../../components/inputs/typeahead.js';
import {iconButton} from '../../components/actions/buttons.js';
import {columnRenderer} from '../entities/column-renderer.js';
import type {IItemsState} from '../../core/widget-like.js';
import {ColumnInput} from './column-combo.js';
import {Platform} from '../shell/platform.js';

/** The input options every u2 input takes (`value`, `nullable`, `required`, …), plus the column ones. */
export interface ColumnInputOptions extends Omit<InputOptions<string | null>, 'label'> {
  /** Narrows the listed columns; a hidden `~` column or one with no name is never listed
   * ({@link Platform.listedColumns}). */
  filter?: (column: DG.Column) => boolean;
  /** The searchable variant: a type-ahead over the column rows (glyph, name, semantic type)
   * instead of the native select — Dart's ColumnComboBox affordances (`column_input.dart:13-17`).
   * Off by default, so the plain choice input stays what `columnInput` returns. */
  rich?: boolean;
  /** The platform-parity variant: the {@link ColumnInput} combo whose dropdown is the
   * ColumnsGrid (the Dart `ColumnComboBox` popup) where the platform is present, and the
   * u2-native searchable list otherwise. */
  grid?: boolean;
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

/** Picks a column of `table` by name, among the columns the platform lists; the item list follows
 * `table.onColumnsChanged` until the input is disposed. A dropped column clears the value, a renamed one carries it
 * over. */
export function columnInput(label: string, table: DG.DataFrame,
  options: ColumnInputOptions & {grid: true}): ColumnInput;
export function columnInput(label: string, table: DG.DataFrame,
  options: ColumnInputOptions & {rich: true}): ColumnPicker;
export function columnInput(label: string, table: DG.DataFrame,
  options?: ColumnInputOptions): ChoiceInput;
export function columnInput(label: string, table: DG.DataFrame,
  options: ColumnInputOptions = {}): ChoiceInput | ColumnPicker | ColumnInput {
  const {filter, rich, grid, placeholder, ...rest} = options;
  if (grid)
    return new ColumnInput({...rest, label, table, filter, placeholder});
  if (rich)
    return new ColumnPicker({...rest, label, table, filter, placeholder});
  const items = () => columnNames(table, filter);
  const input = new ChoiceInput({...rest, label, items: items()});
  const subs = [
    table.onColumnsChanged.subscribe(() => input.setItems(items())),
    // the value moves first: setItems drops a value the new item list no longer offers
    onColumnRenamed(table, (oldName, newName) => {
      if (input.peek() === oldName)
        Input.system(() => input.value = newName);
      input.setItems(items());
    }),
  ];
  input.own(() => {
    for (const sub of subs)
      sub.unsubscribe();
  });
  return input;
}

/** One icon on a {@link TableInput}'s options rail — the platform table input's action set
 * (`table_input.dart:37-58` + the `AppEvents.onInputCreated` augmentations,
 * `inputs.dart:88-112`), per-instance manageable through {@link TableInput.actions}. */
export interface InputAction {
  /** 'open-file' | 'add-from-files' | 'query-db' | custom. */
  id: string;
  /** Font Awesome name from the platform sheet ('folder-open', 'folder-tree', 'database'). */
  icon: string;
  tooltip: string;
  run: (input: TableInput) => void | Promise<void>;
}

const pickerApi = globalThis as {
  grok_UI_PickTableFromFiles?: () => Promise<DG.DataFrame | null>,
  grok_UI_PickTableFromQuery?: () => Promise<DG.DataFrame | null>,
};

export interface TableInputOptions extends Omit<ChoiceInputOptions, 'items'> {
  /** Replaces the default action list wholesale; compose with the current one through the
   * {@link TableInput.actions} setter instead. */
  actions?: InputAction[];
}

/** Picks an open table by name ({@link Platform.tableNames}: a table with no name is not offered, a name several
 * tables share is offered once); the frame itself is {@link table}. The item list follows
 * `onTableAdded`/`onTableRemoved` until the input is disposed. A closed table clears the value. The options rail
 * carries the platform input's action icons: 'Open file' always; 'Add file from Files' and 'Query database' where
 * the platform dialogs answer (feature-detected once, at default construction — headless hosts render the
 * open-file icon alone). A picked frame joins the workspace and becomes the value as a user edit (divergence #16:
 * the Dart icons feed the input's private item list only). */
export class TableInput extends ChoiceInput {
  private _actions: InputAction[];
  private _actionsScope: Scope | undefined;
  private _held: DG.DataFrame | null = null;

  constructor(options: TableInputOptions) {
    const {actions, ...rest} = options;
    super({...rest, items: Platform.tableNames(),
      emptyText: 'No open tables'});
    followOpenTables(this, (names) => this.setItems(names));
    this._actions = actions ?? TableInput._defaults();
    this._renderActions();
    this.own(() => this._actionsScope?.dispose());
  }

  /** The picked table itself: the frame last given to this setter while it is open and its name is the value —
   * so of several open tables under one name, that one — else the open table of that name. Tracked. */
  get table(): DG.DataFrame | null {
    const name = this.value;
    if (name === null)
      return null;
    const held = this._held;
    return held !== null && held.name === name && Platform.isOpen(held) ? held : grok.shell.table(name) ?? null;
  }

  /** Picks `x`, one of the open tables the input offers; any other (a table with no name, one not open) is
   * reported and clears the pick. */
  set table(x: DG.DataFrame | null) {
    const offered = x === null || (Platform.tableNames().includes(x.name) && Platform.isOpen(x));
    if (!offered) {
      Usage.report('table-not-offered', 'a table input holds one of the open tables it offers, and this one has no ' +
        'name or is not open: name it and open it (`grok.shell.addTable`) first');
    }
    this._held = offered ? x : null;
    this.value = this._held?.name ?? null;
  }

  /** The current action list (a copy — add/remove/replace is read-modify-write through the
   * setter, which re-renders the rail). */
  get actions(): InputAction[] {
    return [...this._actions];
  }

  set actions(list: InputAction[]) {
    this._actions = [...list];
    this._renderActions();
  }

  // a per-render scope: each re-set would otherwise pile its Tooltip cleanups on the input's scope
  private _renderActions(): void {
    this._actionsScope?.dispose();
    const scope = new Scope();
    this._actionsScope = scope;
    this.box.querySelector('.u2-input-options')?.replaceChildren();
    for (const action of this._actions) {
      this.addOptions(scope.build(() =>
        iconButton(action.icon, () => {
          const run = action.run(this);
          if (run instanceof Promise)
            run.catch((e) => grok.shell.error(String(e)));
        }, {tooltip: action.tooltip})), {edits: true});
    }
  }

  private static _defaults(): InputAction[] {
    const actions = [TableInput._openFileAction()];
    const files = pickerApi.grok_UI_PickTableFromFiles;
    if (typeof files === 'function')
      actions.push(TableInput._pickAction('add-from-files', 'folder-tree', 'Add file from Files', files));
    const query = pickerApi.grok_UI_PickTableFromQuery;
    if (typeof query === 'function')
      actions.push(TableInput._pickAction('query-db', 'database', 'Query database', query));
    return actions;
  }

  /** The platform input's `folder-open` icon (`table_input.dart:37-58`) as an action: a
   * transient hidden file picker, the table into the workspace, held as a user edit. */
  private static _openFileAction(): InputAction {
    return {id: 'open-file', icon: 'folder-open', tooltip: 'Open file',
      run: async (input) => {
        const table = await TableInput._pickLocalFile(input);
        if (table !== null && !input.scope.isDisposed) {
          input.setItems(Platform.tableNames());
          input.table = table;
        }
      }};
  }

  private static _pickAction(id: string, icon: string, tooltip: string,
    pick: () => Promise<DG.DataFrame | null>): InputAction {
    return {id, icon, tooltip, run: async (input) => {
      const df = await pick();
      if (df == null || input.scope.isDisposed)
        return;
      grok.shell.addTable(df);
      input.setItems(Platform.tableNames());
      input.table = df;
    }};
  }

  /** Settles with the opened table, or null when nothing usable was picked. The picker
   * with no selection settles only where the browser reports the dialog's `cancel`. */
  private static _pickLocalFile(input: TableInput): Promise<DG.DataFrame | null> {
    return new Promise((resolve) => {
      const picker = document.createElement('input');
      picker.type = 'file';
      picker.hidden = true;
      const done = (file?: File) => {
        picker.remove();
        resolve(file === undefined ? null : openTable(file));
      };
      picker.addEventListener('change', () => done(picker.files?.[0] ?? undefined));
      picker.addEventListener('cancel', () => done());
      input.box.append(picker);
      picker.click();
    });
  }
}

export function tableInput(label: LiveOption<string>,
  options: Omit<InputOptions<string | null>, 'label'> & {actions?: InputAction[]} = {}): TableInput {
  return new TableInput({...options, label});
}

/** Picks any number of open tables — a checkbox per table, the value is their names (Dart's
 * `TablesMultiChoiceInput`, `tables_multi_choice_input.dart:48-64`). Checked tables survive an
 * open or close of another one; a closed table leaves the value with it. The import action is the
 * platform input's own, where it is a hidden-selector `TableInput` sitting under the checks
 * (`tables_multi_choice_input.dart:66-77`): the imported table joins the workspace checked. */
export function tablesInput(label: string): MultiChoiceInput {
  const input = new MultiChoiceInput({label, items: Platform.tableNames(), emptyText: 'No open tables'});
  followOpenTables(input, (names) => input.setItems(names));
  importAction(input, (name) => {
    input.setItems(Platform.tableNames());
    const checked = input.peek();
    input.value = checked.includes(name) ? checked : [...checked, name];
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
  private declare _typeAhead: TypeAhead<string>;

  constructor(options: ColumnPickerOptions) {
    super(options, null);
    this.root.dataset.u2 = 'column-picker';
  }

  /** Whether the candidates are open. */
  getState(): ReturnType<Input<string | null>['getState']> {
    return {...super.getState(), open: this._typeAhead.isOpen.peek()};
  }

  /** The columns on offer. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const names = columnNames(this.options.table, this.options.filter);
    const open = this._typeAhead.isOpen.peek();
    return Input.itemsState(names.map((name, index) => ({key: name, label: name, index, rendered: open})), opts?.all);
  }

  revealItem(keyOrLabel: string): Promise<HTMLElement | null> {
    return this._typeAhead.revealItem(keyOrLabel);
  }

  protected createEditor(): HTMLElement {
    const table = this.options.table;
    const names = () => columnNames(table, this.options.filter);
    const typeAhead = new TypeAhead<string>({
      source: (query: string) => Promise.resolve(match(names(), query)),
      renderer: columnRenderer(table),
      placeholder: this.options.placeholder ?? 'Column…',
      debounceMs: 0,
      openOnFocus: true,
    });
    this._typeAhead = typeAhead;

    // seeded before the effects, or the first one would push the type-ahead's empty selection
    // over the value the input was constructed with
    typeAhead.selected.value = this.peek();
    this.effect(() => {
      const picked = typeAhead.selected.value;
      if (picked !== this.peek())
        this.value = picked;
    });
    this.effect(() => {
      const value = this.value;
      if (value !== typeAhead.selected.peek())
        typeAhead.selected.value = value;
    });

    const subs = [
      table.onColumnsChanged.subscribe(() => {
        const value = this.peek();
        if (value !== null && !names().includes(value))
          this.systemWrite(null);
      }),
      onColumnRenamed(table, (oldName, newName) => {
        if (this.peek() === oldName)
          this.systemWrite(newName);
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
  return Platform.listedColumns(table, filter).map((c) => c.name);
}

function match(names: string[], query: string): string[] {
  const q = query.trim().toLowerCase();
  return q ? names.filter((name) => name.toLowerCase().includes(q)) : names;
}

function followOpenTables(input: Control, apply: (names: string[]) => void): void {
  const subs = [grok.events.onTableAdded, grok.events.onTableRemoved]
    .map((event) => event.subscribe(() => apply(Platform.tableNames())));
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
    const table = await openTable(file);
    if (table !== null && !input.scope.isDisposed)
      pick(table.name);
  };
  input.scope.listen(picker, 'change', onChange);
  const open = input.build(() => iconButton('folder-open', () => picker.click(), {tooltip: 'Open file'}));
  input.box.append(picker);
  input.addOptions(open, {edits: true});
}

/** A local file as the workspace's newest table, under the name it ended up with. An extension
 * neither reader takes is refused with the platform input's own message
 * (`table_input.dart:41-43`) — Dart hands those to `FileHandler`, which opens xlsx, sdf and the
 * rest, and js-api has no door to it. */
async function openTable(file: File): Promise<DG.DataFrame | null> {
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
    return grok.shell.addTable(table);
  } catch (e) {
    grok.shell.error(`${file.name} could not be opened.`);
    return null;
  }
}
