/* `domains.table(address)` — one await, then a table handle everything else hangs off (GOAL
   "What it looks like"): sources over the table, and the per-table registries an app fills once —
   actions (ribbon, row hover, context menu), the renderer (lists, pickers, chips) and validators
   (forms). Registers nothing with the platform: rendering reuses the handler `ObjectHandler.forEntity`
   already resolves for the table's rows. `domains.*` is also where every domain control is made
   (`domains.form(src)`, `domains.list(src)`, …): one object, each member constructing at call
   time, so the modules behind it may import this one back. `TRow` (an app's own row type) types
   the rows every registry and source hands out; the default reads any column as unknown. */
import * as DG from 'datagrok-api/dg';
import {Access} from '../../core/access.js';
import type {Capability} from '../../core/access.js';
import type {IProperty} from '../../core/property-like.js';
import type {ObjectRenderer} from '../../core/object-renderer.js';
import type {Action} from '../../components/actions/actions.js';
import {backends} from '../../sources/backends.js';
import type {DomainTableInfoLike, DomainTableLike} from '../../sources/domain-backend.js';
import {DomainSource} from '../../sources/domain-source.js';
import type {DomainSourceOptions} from '../../sources/domain-source.js';
import {Rows} from '../../sources/rows-like.js';
import type {ColumnOf, DomainRowLike, RowValues, RowView} from '../../sources/rows-like.js';
import {SharedSession} from '../../sources/session.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {divV, span, timestamp} from '../../core/elements.js';
import {text} from '../../core/text.js';
import {handlerRenderer} from '../entities/entity.js';
import {appView} from '../shell/app-view.js';
import {SYSTEM_COLUMNS} from './backend.js';
import {DomainErrors} from './errors.js';
// the cycles with the control modules are function-scoped only: none reads another at load time
import {DomainApp} from './app.js';
import type {DomainAppOptions} from './app.js';
import {DomainForm} from './form.js';
import type {DomainFormOptions, DomainFormTarget} from './form.js';
import {DomainList} from './list.js';
import type {DomainListOptions} from './list.js';
import {DomainPick} from './pick.js';
import type {DomainPickOptions} from './pick.js';
import {DomainSearch} from './search.js';
import type {DomainSearchOptions} from './search.js';
import {DomainFilters} from './filters.js';
import type {DomainFiltersOptions} from './filters.js';
import {DomainGrid, DomainDataTable} from './grid.js';
import type {DomainGridOptions, DomainDataTableOptions} from './grid.js';
import {DomainTree} from './tree.js';
import type {DomainTreeOptions} from './tree.js';
import {DomainHistory} from './history.js';
import type {DomainHistoryTarget} from './history.js';
import {DomainChildren} from './children.js';
import type {DomainChildrenOptions} from './children.js';
import {saveButton, discardButton, newButton} from './buttons.js';
import {route} from './routes.js';
import {mountView} from './view-sync.js';
import {bulkEdit} from './bulk.js';
import {openImport} from './import.js';
const REF_ADDRESS = /^\w+\.\w+$/;

/** An action over one row: `requires` names the capability it needs (permission ⇒ hidden),
 * `when` narrows it per row (state ⇒ absent), `run` takes the row. */
export interface DomainAction<TRow extends DomainRowLike = DomainRowLike> {
  name: string;
  icon?: string;
  requires?: Capability;
  enabled?: boolean;
  when?(row: RowView<TRow>): boolean;
  run(row: RowView<TRow>): void;
}

export type RowValidator<TRow extends DomainRowLike = DomainRowLike> =
  (value: unknown, row: RowView<TRow>) => string | null;

/** What `DomainTable.app()` takes: the app's options without the table, plus the view's name
 * (the plural name by default), its base path (`/domains/<schema>/<table>` by default) and the
 * app class — a `DomainApp` subclass with its own ribbon, presets and shortcuts. */
export interface DomainAppViewOptions extends Omit<DomainAppOptions, 'table' | 'base'> {
  name?: string;
  path?: string;
  app?: new (options: DomainAppOptions) => DomainApp;
}

export class ActionRegistry<TRow extends DomainRowLike = DomainRowLike> {
  private readonly _actions: DomainAction<TRow>[] = [];

  /** Returns the unregister function. */
  add(action: DomainAction<TRow>): () => void {
    this._actions.push(action);
    return () => {
      const at = this._actions.indexOf(action);
      if (at >= 0)
        this._actions.splice(at, 1);
    };
  }

  /** The actions that apply to `row`, bound to it — what `rowActions` and a menu take. */
  for(row: RowView<TRow>): Action[] {
    return ActionRegistry.bind(this._actions, row);
  }

  static bind<T extends DomainRowLike>(actions: readonly DomainAction<T>[], row: RowView<T>): Action[] {
    return actions.filter((a) => a.when === undefined || a.when(row)).map((a): Action =>
      ({name: a.name, icon: a.icon, requires: a.requires, enabled: a.enabled, run: () => a.run(row)}));
  }
}

export class ValidatorRegistry<TRow extends DomainRowLike = DomainRowLike> {
  private readonly _byColumn = new Map<string, RowValidator[]>();

  /** Returns the unregister function. */
  add(column: ColumnOf<TRow>, validator: RowValidator<TRow>): () => void {
    const list = this._byColumn.get(column) ?? [];
    this._byColumn.set(column, list);
    list.push(validator as RowValidator);
    return () => {
      const at = list.indexOf(validator as RowValidator);
      if (at >= 0)
        list.splice(at, 1);
    };
  }

  /** The columns something was registered for — what a whole-row check walks. */
  get columns(): string[] {
    return [...this._byColumn.keys()];
  }

  /** The first problem the column's validators report, null when the value passes. */
  check(column: string, value: unknown, row: RowView<TRow>): string | null {
    for (const validator of this._byColumn.get(column) ?? []) {
      const message = validator(value, row);
      if (message !== null)
        return message;
    }
    return null;
  }
}

export class DomainTable<TRow extends DomainRowLike = DomainRowLike> {
  /** Everything an app declares about the table's actions, once. */
  readonly actions = new ActionRegistry<TRow>();
  readonly validators = new ValidatorRegistry<TRow>();
  /** The platform handler of the table's rows — deep links, the entity view, `rowFrom`. */
  readonly handler: DG.DomainObjectHandler;
  /** How lists, pickers and chips show a row; the handler's own rendering by default. */
  renderer: ObjectRenderer<RowView<TRow>>;
  /** The renderer the handle was built with — what tells a table an app dressed from a bare one. */
  readonly defaultRenderer: ObjectRenderer<RowView<TRow>>;

  private static readonly _bySource = new WeakMap<DomainSource, DomainTable>();
  private static readonly _byRow = new WeakMap<ReadonlySignal<RowView | null>, DomainSource>();

  constructor(readonly address: string, readonly table: DomainTableLike, readonly access: Access) {
    this.handler = new DG.DomainObjectHandler(address);
    this.renderer = this.defaultRenderer = DomainTable.handlerRenderer(this);
  }

  get info(): DomainTableInfoLike {
    return this.table.info;
  }

  get properties(): IProperty[] {
    return this.table.properties;
  }

  /** The handle a source was made through, for the controls over it. */
  static of(source: DomainSource): DomainTable | undefined {
    return DomainTable._bySource.get(source);
  }

  /** The source behind a `currentRow` signal handed to a form. */
  static sourceOf(row: ReadonlySignal<RowView | null>): DomainSource | undefined {
    return DomainTable._byRow.get(row);
  }

  /** Whether the column holds another row's id rather than a value of its own — a ref column
   * (semType `<schema>.<table>`), a user or a group (the js-api editor's `isReferenceProperty`). */
  static isReference(prop: IProperty): boolean {
    const semType = prop.semType ?? '';
    return semType === 'User' || semType === 'Group' || REF_ADDRESS.test(semType);
  }

  /** A started source over this table. */
  source(options: Omit<DomainSourceOptions, 'table' | 'defaults'> & {defaults?: RowValues<TRow>} = {}):
    DomainSource<TRow> {
    if (options.deleted !== undefined && options.deleted !== 'exclude')
      DomainSource.requireRestore(this.table);
    return this._track(new DomainSource<TRow>({...options, table: this.address}));
  }

  /** Brings a soft-deleted row back at once — the Delete grant undone, with no session behind it.
   * A source stages restores into its unit of work instead ({@link DomainSource.stageRestore}). */
  restore(id: string): Promise<void> {
    DomainSource.requireRestore(this.table);
    return this.table.restore!(id);
  }

  /** A source holding one pristine draft over `values` as its current row — what a create form
   * binds to; `save()` on it inserts the row. */
  draft(values: RowValues<TRow> = {}): DomainSource<TRow> {
    return this._track(new DomainSource<TRow>({table: this.address, draft: true, defaults: values}));
  }

  /** The row as the platform sees it — the handler's `DomainRow`, built locally. */
  row(row: RowView<TRow>): DG.DomainRow {
    const values: Record<string, unknown> = {...row};
    if (Rows.isDraft(row))
      delete values.id;
    return this.handler.rowFrom(values);
  }

  /** Opens a saved row: the entity page of the app open over this table (find-or-activate), else
   * the platform's entity view. */
  open(row: RowView<TRow>): void {
    if (Rows.isDraft(row) || DomainApp.activate(DomainApp.baseOf(this.address), row.id))
      return;
    this.handler.openRow(this.row(row));
  }

  /** The table as a platform view (GOAL "What it looks like": `grok.shell.addView(issues.app())`):
   * list ⇄ entity page under one session, the ribbon, the status bar, the URL — and every gate in
   * front of dropping unsaved changes: in-app navigation, the view's ✕, the browser's unload. */
  app(options: DomainAppViewOptions = {}): DG.ViewBase {
    const {name, path, app: App = DomainApp, ...rest} = options;
    const base = path ?? `/domains/${this.address.replace('.', '/')}`;
    const session = new SharedSession();
    const app = SharedSession.runWith(session, () => new App({...rest, table: this, base}));
    const view = appView({name: name ?? DomainApp.titleOf(this.info), content: app,
      ribbon: app.ribbonGroups(), status: DomainApp.statusPanels(app), path: app.path});
    mountView({view, base, pinned: path !== undefined, session, params: ['entity', 'q', 'search', 'trash'],
      baseOf: () => app.base, rebase: (prefix) => app.rebase(prefix), open: (address) => app.open(address),
      own: (dispose) => app.own(dispose)});
    app.own(DomainApp.register(app, view));
    app.guardUnload();
    return view;
  }

  /** The handler-backed renderer: the name column as the caption where the table declares one,
   * the handler's own rendering otherwise, and the schema-driven rendering where no handler
   * claims the rows. */
  static handlerRenderer(table: DomainTable): ObjectRenderer<RowView> {
    const inner = handlerRenderer<DG.DomainRow>();
    const plain = DomainTable.schemaRenderer(() => ({properties: table.properties, info: table.info}));
    // a draft has no identity the handler could show; a handler card counts only when the handler
    // defines one — the platform's default is the generic property table
    const handled = (row: RowView): DG.DomainRow | null => {
      if (Rows.isDraft(row))
        return null;
      const x = table.row(row);
      return inner.handlerFor(x) === null ? null : x;
    };
    const by = <T>(row: RowView, on: (x: DG.DomainRow) => T, fallback: (row: RowView) => T): T => {
      const x = handled(row);
      return x === null ? fallback(row) : on(x);
    };
    // a card of the handler's own: a registered subclass that overrides renderCard — the reflective
    // default and the platform's Dart meta both paint the generic property table
    const ownCard = (x: DG.DomainRow): boolean => {
      const handler = inner.handlerFor(x);
      if (!(handler instanceof DG.DomainObjectHandler))
        return false;
      for (let p = Object.getPrototypeOf(handler); p !== null && p !== DG.DomainObjectHandler.prototype;
        p = Object.getPrototypeOf(p)) {
        if (Object.prototype.hasOwnProperty.call(p, 'renderCard'))
          return true;
      }
      return false;
    };
    return {
      caption: (row) => {
        const name = plain.caption(row);
        return name === row.id ? by(row, (x) => inner.caption(x), plain.caption) : name;
      },
      icon: (row) => by(row, (x) => inner.icon(x), (r) => span(plain.caption(r))),
      listItem: (row) => by(row, (x) => inner.listItem(x), plain.listItem!),
      markup: (row) => by(row, (x) => inner.markup(x), plain.listItem!),
      tooltip: (row) => by(row, (x) => inner.tooltip(x), plain.card!),
      card: (row) => {
        const x = handled(row);
        return x !== null && ownCard(x) ? inner.card(x) : plain.card!(row);
      },
    };
  }

  /** What a row looks like from its schema alone: the name column (else the first business-key
   * column, else the id; "New <singular>" for a draft without one), and on a card the
   * list-item-rendering recipe — the title, one muted description line, the creation time. */
  static schemaRenderer(schema: () => {properties: IProperty[], info: DomainTableInfoLike}): ObjectRenderer<RowView> {
    const label = (row: RowView): {text: string, draft: boolean} => {
      const info = schema().info;
      const column = info.nameColumn ?? info.businessKey[0];
      const name = column === undefined ? undefined : row[column];
      if (name !== null && name !== undefined && name !== '')
        return {text: String(name), draft: false};
      return Rows.isDraft(row) ? {text: `New ${info.singularName.toLowerCase() || 'row'}`, draft: true} :
        {text: row.id, draft: false};
    };
    const titleOf = (row: RowView): HTMLElement => {
      const {text: caption, draft} = label(row);
      return span(caption, draft ? 'u2-domain-list-name u2-domain-draft' : 'u2-domain-list-name');
    };
    const description = (row: RowView): string => {
      const {properties, info} = schema();
      const prose = properties.filter((p) => p.name !== info.nameColumn &&
        !SYSTEM_COLUMNS.some(([name]) => name === p.name) && !DomainTable.isReference(p) &&
        (p.propertyType ?? p.type) === 'string' && text(row[p.name!]) !== '');
      const first = prose.find((p) => p.name === 'description') ?? prose[0];
      return first === undefined ? '' : text(row[first.name!]);
    };
    return {
      caption: (row) => label(row).text,
      listItem: titleOf,
      card: (row) => {
        const created = row.created_on;
        const line = description(row);
        const title = titleOf(row);
        title.classList.add('u2-domain-card-title');
        return divV([
          title,
          ...(line === '' ? [] : [span(line, 'u2-domain-card-description')]),
          ...(created === null || created === undefined ? [] :
            [timestamp(created as Date | number | string, 'u2-domain-card-time')]),
        ], 'u2-domain-card');
      },
    };
  }

  private _track(source: DomainSource<TRow>): DomainSource<TRow> {
    DomainTable._bySource.set(source, this);
    DomainTable._byRow.set(source.currentRow, source);
    source.guard(() => this._validate(source));
    // last resort, after every guard: the writer refused over a cell and named neither
    source.nameCell = (row, column, problem) => this._nameCell(row as RowView<TRow>, column, problem);
    // every refused save is one balloon, form or not; a load failure is the list's and the hint's
    // to show
    source.effect(() => {
      const error = source.error.value;
      if (error !== undefined && source.state.peek() !== 'error' &&
          DomainErrors.codeOf(error) !== DomainSource.REFUSED)
        DomainErrors.report(error, source);
    });
    source.start();
    return source;
  }

  /** The validators over every row the batch will send — a form wires them into its inputs, but a
   * row action or a plain `row.column = …` never goes near a form, and the rule still holds. */
  private _validate(source: DomainSource<TRow>): string | null {
    const columns = this.validators.columns;
    if (columns.length === 0)
      return null;
    for (const row of source.pending()) {
      if (row[Rows.STATE] === 'deleted' || row[Rows.STATE] === 'restored')
        continue;
      for (const column of columns) {
        const problem = this.validators.check(column, row[column], row as RowView<TRow>);
        if (problem === null)
          continue;
        // the refusal names the row it is about: from a list, "Status: …" alone names nothing
        source.markProblem(row.id);
        return this._nameCell(row, column, problem);
      }
    }
    return null;
  }

  /** A refusal about one cell, named: the row's caption and the column's, as a list or a grid
   * needs them — the writer's own "Value can't be empty" says nothing about which field. */
  private _nameCell(row: RowView<TRow>, column: string, problem: string): string {
    return `${this.renderer.caption(row)}: ${DomainTable._caption(this.properties, column)}: ${problem}`;
  }

  private static _caption(properties: IProperty[], column: string): string {
    const name = properties.find((p) => p.name === column)?.friendlyName ?? column;
    return `${name.charAt(0).toUpperCase()}${name.slice(1)}`;
  }
}

export const domains = {
  /** The table's shape and the caller's access in one round of requests; a handle per call, so
   * the access is re-acquired by acquiring again. */
  async table<TRow extends DomainRowLike = DomainRowLike>(address: string): Promise<DomainTable<TRow>> {
    const backend = backends.domain;
    if (backend === undefined)
      throw new Error('no platform backend for domain tables');
    const table = await backend.table(address);
    return new DomainTable<TRow>(address, table, Access.from(await table.access()));
  },
  form: (target: DomainFormTarget, options?: DomainFormOptions): DomainForm => new DomainForm(target, options),
  list: (source: DomainSource, options?: DomainListOptions): DomainList => new DomainList(source, options),
  pick: (table: string, options?: DomainPickOptions): DomainPick => new DomainPick(table, options),
  grid: (source: DomainSource, options?: DomainGridOptions): DomainGrid => new DomainGrid(source, options),
  dataTable: (source: DomainSource, options?: DomainDataTableOptions): DomainDataTable =>
    new DomainDataTable(source, options),
  /** A tree over a table the schema declares a hierarchy; `selected` drives a collection's query
   * as `<fk> under "<id>"`. */
  tree: <TRow extends DomainRowLike = DomainRowLike>(table: DomainTable<TRow>,
    options?: DomainTreeOptions<TRow>): DomainTree<TRow> => new DomainTree(table, options),
  search: (source: DomainSource, options?: DomainSearchOptions): DomainSearch => new DomainSearch(source, options),
  filters: (source: DomainSource, options?: DomainFiltersOptions): DomainFilters => new DomainFilters(source, options),
  history: (source: DomainSource, row?: DomainHistoryTarget): DomainHistory => new DomainHistory(source, row),
  children: (parent: DomainSource, options?: DomainChildrenOptions): DomainChildren =>
    new DomainChildren(parent, options),
  app: (options: DomainAppOptions): DomainApp => new DomainApp(options),
  route,
  bulkEdit,
  /** The import wizard over a table: source → mapping → preview → the server's report. */
  import: openImport,
  saveButton,
  discardButton,
  newButton,
};
