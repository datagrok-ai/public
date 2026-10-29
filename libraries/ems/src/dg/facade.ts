import {emsBackends, capital, nounOf} from '../core/domain-backend.js';
import * as DG from 'datagrok-api/dg';
import {Access, Rows, allowedActions, div, divH, divV, span, timestamp, text, notify, plural, untracked}
  from '@datagrok-libraries/u2';
import type {Capability, IProperty, ObjectRenderer, Action, ColumnOf, DomainRowLike, RowValues, RowView,
  ReadonlySignal, NotifyHandle, Reason} from '@datagrok-libraries/u2';
import type {DomainTableInfoLike, DomainTableLike} from '../core/domain-backend.js';
import {DomainSource} from '../core/domain-source.js';
import type {DomainSourceOptions} from '../core/domain-source.js';
import {SharedSession, confirmDiscard} from '../core/session.js';
import {handlerRenderer, ViewHost, mountView, Editors} from '@datagrok-libraries/u2/dg';
import {DomainApp} from './app.js';
import {authoring} from './authoring/index.js';
import {SYSTEM_COLUMNS, DgDomainBackend} from './backend.js';
import {saveButton, discardButton, newButton} from './buttons.js';
import {DomainChildren} from './children.js';
import {DomainErrors} from './errors.js';
import {DomainFilters} from './filters.js';
import {DomainForm} from './form.js';
import {DomainGrid, DomainDataTable} from './grid.js';
import {DomainHistory} from './history.js';
import {openImport} from './import.js';
import {DomainViews} from './views.js';
import type {DomainViewsOptions} from './views.js';
import type {DomainImporter} from './import.js';
import {bulkEdit} from './bulk.js';
import {PropertyDialog} from './authoring/property-dialog.js';
import {BindingDialog, BindingWizard, createBinding} from './authoring/binding-dialog.js';
import {BindingEditDialog, editBinding} from './authoring/binding-edit-dialog.js';
// the cycles with the control modules are function-scoped only: none reads another at load time
import type {DomainAppOptions} from './app.js';
import type {DomainChildrenOptions} from './children.js';
import type {DomainFiltersOptions} from './filters.js';
import type {DomainFormOptions, DomainFormTarget} from './form.js';
import type {DomainGridOptions, DomainDataTableOptions} from './grid.js';
import type {DomainHistoryTarget, DomainHistoryOptions} from './history.js';
import {DomainList} from './list.js';
import type {DomainListOptions} from './list.js';
import {DomainPick} from './pick.js';
import type {DomainPickOptions} from './pick.js';
import {DomainSearch} from './search.js';
import type {DomainSearchOptions} from './search.js';
import {DomainTree} from './tree.js';
import type {DomainTreeOptions} from './tree.js';
import {route} from './routes.js';
export type {CreateBindingOptions, BindingResult, AccessOutcome} from './authoring/binding-dialog.js';
export type {EditBindingResult} from './authoring/binding-edit-dialog.js';
const REF_ADDRESS = /^\w+\.\w+$/;

/** An action over one row: `requires` names the capability it needs (permission ⇒ hidden),
 * `when` narrows it per row (state ⇒ absent), `disabled` says why it is off, `run` takes the row. */
export interface DomainAction<TRow extends DomainRowLike = DomainRowLike> {
  name: string;
  icon?: string;
  requires?: Capability;
  disabled?: Reason;
  when?(row: RowView<TRow>): boolean;
  run(row: RowView<TRow>): void;
}

export type RowValidator<TRow extends DomainRowLike = DomainRowLike> =
  (value: unknown, row: RowView<TRow>) => string | null;

/** An app's setup of a table — validators, renderer, actions, importer. */
export type DomainConfigure<TRow extends DomainRowLike = DomainRowLike> =
  (table: DomainTable<TRow>) => void | Promise<void>;

/** What `DomainTable.app()` takes: the app's options without the table — its own ribbon, menu,
 * shortcuts and setup among them — plus the view's name (the plural name by default) and its base
 * path (`/domains/<schema>/<table>` by default). */
export interface DomainAppViewOptions extends Omit<DomainAppOptions, 'table' | 'base'> {
  name?: string;
  path?: string;
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
      ({text: a.name, icon: a.icon, requires: a.requires, disabled: a.disabled, run: () => a.run(row)}));
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

/** The handle `domains.table(address)` answers: sources over the table, and the per-table registries an app
 * fills once (actions, the renderer, validators). Registers nothing with the platform: rendering reuses the
 * handler `ObjectHandler.forEntity` already resolves. `TRow` types the rows every registry and source hands out. */
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
  /** What every import of the table does to its rows first — derived columns, rows already there. */
  importer: DomainImporter | undefined;

  private static readonly _bySource = new WeakMap<DomainSource, DomainTable>();
  private static readonly _byRow = new WeakMap<ReadonlySignal<RowView | null>, DomainSource>();
  private static readonly _shared = new WeakMap<object, Map<string, Promise<DomainTable<any>>>>();

  constructor(readonly address: string, readonly table: DomainTableLike, readonly access: Access) {
    this.handler = new DG.DomainObjectHandler(address);
    this.renderer = this.defaultRenderer = DomainTable.handlerRenderer(this);
  }

  static forget(schema: string): void {
    const backend = emsBackends.domain;
    const handles = backend === undefined ? undefined : DomainTable._shared.get(backend);
    for (const address of [...handles?.keys() ?? []]) {
      if (address.startsWith(`${schema}.`))
        handles!.delete(address);
    }
    if (backend instanceof DgDomainBackend)
      backend.invalidate();
  }

  /** One configured handle over `address` per backend, for whoever only reads it. */
  static shared<TRow extends DomainRowLike = DomainRowLike>(address: string): Promise<DomainTable<TRow>> {
    const backend = emsBackends.domain;
    if (backend === undefined)
      return domains.table<TRow>(address);
    const handles = DomainTable._shared.get(backend) ?? new Map<string, Promise<DomainTable<any>>>();
    DomainTable._shared.set(backend, handles);
    let handle = handles.get(address);
    if (handle === undefined) {
      handles.set(address, handle = domains.table<TRow>(address));
      handle.catch(() => handles.delete(address));
    }
    return handle as Promise<DomainTable<TRow>>;
  }

  get info(): DomainTableInfoLike {
    return this.table.info;
  }

  get properties(): IProperty[] {
    return this.table.properties;
  }

  /** Registers `fn` for every handle over `address` made from now on; the same `id` again replaces
   * it, so a bundle evaluated twice does not add a second. `fn` may be async; a throw rejects the
   * handle. Returns the unregister function. */
  static configure<TRow extends DomainRowLike = DomainRowLike>(address: string, fn: DomainConfigure<TRow>,
    {id}: {id: string}): () => void {
    const entries = CONFIGURE.get(address) ?? new Map<string, DomainConfigure<any>>();
    CONFIGURE.set(address, entries);
    entries.set(id, fn);
    return () => {
      if (entries.get(id) === fn)
        entries.delete(id);
    };
  }

  /** A column as a form labels it: the declared caption, else the name made presentable. */
  static captionOf(prop: IProperty): string {
    return capital(prop.friendlyName ?? (prop.name ?? '').replace(/_/g, ' '));
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

  /** What a read-only row offers — a tree node, a data table row: Open, the table's actions and the caller's. */
  rowActions(row: RowView<TRow>, extra: readonly DomainAction<TRow>[] = [], openShortcut?: string): Action[] {
    const open: Action[] = Rows.isDraft(row) ? [] :
      [{text: 'Open', icon: 'external-link', shortcut: openShortcut, run: () => this.open(row)}];
    return allowedActions([...open, ...this.actions.for(row), ...ActionRegistry.bind(extra, row)],
      {access: this.access, row});
  }

  /** Opens a saved row: the entity page of the app open over this table (find-or-activate), else
   * the platform's entity view. Deferred: the platform brings a clicked view back to the front after the click. */
  open(row: RowView<TRow>): void {
    if (Rows.isDraft(row))
      return;
    setTimeout(() => {
      if (!DomainApp.activate(DomainApp.baseOf(this.address), row.id) && !DomainViews.openRow(this.address, row.id))
        this.handler.openRow(this.row(row));
    });
  }

  /** The table as a platform view (`grok.shell.addView(issues.app())`):
   * list ⇄ entity page under one session, the ribbon, the status bar, the URL — and every gate in
   * front of dropping unsaved changes: in-app navigation, the view's ✕, the browser's unload. */
  app(options: DomainAppViewOptions = {}): DG.ViewBase {
    const {name, path, ...rest} = options;
    const base = path ?? `/domains/${this.address.replace('.', '/')}`;
    const session = new SharedSession();
    let app!: DomainApp;
    const view = ViewHost.create({name: name ?? DomainApp.titleOf(this.info), content: () => {
      app = SharedSession.runWith(session, () => new DomainApp({...rest, table: this, base}));
      return {root: app, ribbon: app.ribbonGroups(), status: DomainApp.statusPanels(app), path: app.path};
    }}).view;
    mountView({view, base, pinned: path !== undefined, session,
      confirmLeave: () => confirmDiscard(session, {action: 'close the view'}),
      params: ['entity', 'q', 'search', 'sort', 'trash'],
      baseOf: () => app.base, rebase: (prefix) => app.rebase(prefix), open: (address) => app.open(address),
      own: (dispose) => app.own(dispose)});
    app.own(DomainApp.register(app, view));
    app.guardUnload();
    return view;
  }

  static readonly CARD_PICTURE = {width: 120, height: 80};

  /** The column a card and a table row show as a picture: the first one a display rule draws (a structure). */
  static pictureOf(properties: IProperty[]): IProperty | undefined {
    return properties.find((p) => !SYSTEM_COLUMNS.some(([name]) => name === p.name) && Editors.displays(p));
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
      // a bare number names nothing on its own: "Purchase order 9", not "9"
      if (typeof name === 'number' || typeof name === 'bigint')
        return {text: `${capital(nounOf(info, true))} ${name}`, draft: false};
      if (name !== null && name !== undefined && name !== '')
        return {text: String(name), draft: false};
      return Rows.isDraft(row) ? {text: `New ${nounOf(info, true)}`, draft: true} :
        {text: row.id, draft: false};
    };
    const titleOf = (row: RowView): HTMLElement => {
      const {text: caption, draft} = label(row);
      return span(caption, draft ? 'u2-domain-list-name u2-domain-draft' : 'u2-domain-list-name');
    };
    const descriptions = (row: RowView, count: number): string[] => {
      const {properties, info} = schema();
      const prose = properties.filter((p) => p.name !== info.nameColumn &&
        !SYSTEM_COLUMNS.some(([name]) => name === p.name) && !DomainTable.isReference(p) &&
        (p.propertyType ?? p.type) === 'string' && !Editors.displays(p) && text(row[p.name!]) !== '');
      const first = prose.find((p) => p.name === 'description');
      return [...first === undefined ? [] : [first], ...prose.filter((p) => p !== first)].slice(0, count)
        .map((p) => text(row[p.name!]));
    };
    return {
      caption: (row) => label(row).text,
      listItem: titleOf,
      card: (row) => {
        const created = row.created_on;
        const picture = DomainTable.pictureOf(schema().properties);
        const title = titleOf(row);
        title.classList.add('u2-domain-card-title');
        const body = divV([
          title,
          ...descriptions(row, picture === undefined ? 1 : 2).map((line) => span(line, 'u2-domain-card-description')),
          ...(created === null || created === undefined ? [] :
            [timestamp(created as Date | number | string, 'u2-domain-card-time')]),
        ], 'u2-domain-card');
        if (picture === undefined)
          return body;
        const value = row[picture.name!];
        const art = value === null || value === undefined || value === '' ? null :
          Editors.display(picture, value, DomainTable.CARD_PICTURE);
        body.classList.add('u2-domain-card-text');
        return divH([div([...art === null ? [] : [art]], 'u2-domain-card-picture'), body],
          'u2-domain-card u2-domain-card-pictured');
      },
    };
  }

  private _track(source: DomainSource<TRow>): DomainSource<TRow> {
    DomainTable._bySource.set(source, this);
    DomainTable._byRow.set(source.currentRow, source);
    source.guard(() => this._validate(source));
    // last resort, after every guard: the writer refused over a cell and named neither
    source.nameCell = (row, column, problem) => nameCell(this, row as RowView<TRow>, column, problem, source);
    source.caption = (row) => this.renderer.caption(row as RowView<TRow>);
    // one refusal balloon per session, closed once the refusal is settled; a load failure is the
    // list's to show
    source.effect(() => {
      const error = source.error.value;
      const shown = BALLOONS.get(source.session);
      if (shown?.source === source) {
        shown.balloon?.close();
        BALLOONS.delete(source.session);
      }
      if (error === undefined || source.state.peek() === 'error' || DomainErrors.codeOf(error) === DomainSource.REFUSED)
        return;
      BALLOONS.get(source.session)?.balloon?.close();
      const fix = untracked(() => this._fieldsToFix(source, error));
      BALLOONS.set(source.session, {source, balloon: fix.length === 0 ? DomainErrors.report(error, source) :
        notify.error(`Fix ${fix.join(', ')} to save`)});
    });
    source.start();
    return source;
  }

  private _fieldsToFix(source: DomainSource<TRow>, error: unknown): string[] {
    const row = source.currentRow.peek();
    const about = source.problemRow.peek();
    if (DomainErrors.codeOf(error) !== 'validation' || row === null || (about !== null && about !== row.id))
      return [];
    const edit = source.edit.peek();
    return this.properties.filter(({name}) => name !== undefined &&
      ((edit?.errorOf(row.id, name) ?? null) !== null || this.validators.check(name, row[name], row) !== null))
      .map((p) => DomainTable.captionOf(p));
  }

  // a row action or a plain `row.column = …` never goes near a form, and the rule still holds
  private _validate(source: DomainSource<TRow>): string | null {
    const found = firstProblem(this, source.pending().filter((row) =>
      row[Rows.STATE] !== 'deleted' && row[Rows.STATE] !== 'restored'), source);
    if (found === null)
      return null;
    source.markProblem(found.row.id);
    return found.message;
  }
}

const BALLOONS = new WeakMap<object, {source: DomainSource<any>, balloon: NotifyHandle | null}>();

// one registry per page, not per bundle: the `/domains` route runs in PowerPack's copy of this
// library; bump the version when the `DomainTable` contract changes
const CONFIGURE: Map<string, Map<string, DomainConfigure<any>>> =
  (globalThis as Record<PropertyKey, any>)[Symbol.for('ems.configure.v1')] ??= new Map();

async function configured<TRow extends DomainRowLike>(table: DomainTable<TRow>): Promise<DomainTable<TRow>> {
  for (const fn of [...CONFIGURE.get(table.address)?.values() ?? []])
    await fn(table);
  return table;
}

/** The first problem the table's validators find over `rows`, named by row and column. */
export function firstProblem<TRow extends DomainRowLike>(table: DomainTable<TRow>, rows: Iterable<RowView<TRow>>,
  source?: DomainSource<TRow>): {row: RowView<TRow>, message: string} | null {
  const columns = table.validators.columns;
  for (const row of rows) {
    for (const column of columns) {
      const problem = table.validators.check(column, row[column], row);
      // the refusal names the row it is about: from a list, "Status: …" alone names nothing
      if (problem !== null)
        return {row, message: nameCell(table, row, column, problem, source)};
    }
  }
  return null;
}

/** A refusal about one cell, named: the row's caption and the column's, as a list or a grid needs
 * them — the writer's own "Value can't be empty" says nothing about which field. */
function nameCell<TRow extends DomainRowLike>(table: DomainTable<TRow>, row: RowView<TRow>, column: string,
  problem: string, source?: DomainSource<TRow>): string {
  const caption = DomainTable.captionOf(table.properties.find((p) => p.name === column) ?? {name: column});
  // "CAS must look like 50-78-2" names its field already
  const field = problem.toLowerCase().startsWith(caption.toLowerCase()) ? '' : `${caption}: `;
  return `${rowName(table, row, source)}: ${field}${problem}`;
}

/** The row's caption — or, where that is only its id, its place among the source's rows: "Order line 2". */
function rowName<TRow extends DomainRowLike>(table: DomainTable<TRow>, row: RowView<TRow>,
  source?: DomainSource<TRow>): string {
  const caption = table.renderer.caption(row);
  const at = caption !== row.id || source === undefined ? -1 :
    source.rows.items.peek().findIndex((r) => r.id === row.id);
  return at < 0 ? caption : `${capital(nounOf(table.info, true))} ${at + 1}`;
}

/** `domains.table(address)` and every domain control (`domains.form(src)`, `domains.list(src)`, …), each made at
 * call time, so the modules behind it may import this one back. */
export const domains = {
  /** The table's shape and the caller's access in one round of requests, with everything
   * {@link configure}d for the address applied; a handle per call, so the access is re-acquired by
   * acquiring again. */
  async table<TRow extends DomainRowLike = DomainRowLike>(address: string): Promise<DomainTable<TRow>> {
    const backend = emsBackends.domain;
    if (backend === undefined)
      throw new Error('no platform backend for domain tables');
    const table = await backend.table(address);
    return configured(new DomainTable<TRow>(address, table, Access.from(await table.access())));
  },
  configure: DomainTable.configure,
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
  history: (source: DomainSource, row?: DomainHistoryTarget, options?: DomainHistoryOptions): DomainHistory =>
    new DomainHistory(source, row, options),
  children: (parent: DomainSource, options?: DomainChildrenOptions): DomainChildren =>
    new DomainChildren(parent, options),
  app: (options: DomainAppOptions): DomainApp => new DomainApp(options),
  views: (options: DomainViewsOptions): DomainViews => new DomainViews(options),
  route,
  bulkEdit,
  /** The import wizard over a table: source → mapping → preview → the server's report. */
  import: openImport,
  async editProperties(schema: string, table?: string): Promise<boolean> {
    const saved = await PropertyDialog.open(schema, table);
    if (saved)
      await DomainViews.reload(schema);
    return saved;
  },
  extensible: (schema: string): Promise<Set<string>> => PropertyDialog.extensible(schema),
  saveButton,
  discardButton,
  newButton,
  /** The manifest editor and its models — the Design step of binding authoring — and the two
   * dialogs over them: "Create domain schema" (Connection › Design › Review › Created; resolves
   * to `{name, access}`) and "Edit binding" over a registered one (Design › Review › Saved;
   * resolves to `{name, applied}`). The editor stays platform-free; the dialogs are what talk to
   * `grok.dapi`. */
  authoring: Object.assign(authoring, {createBinding, BindingDialog, editBinding, BindingEditDialog, BindingWizard}),
};
