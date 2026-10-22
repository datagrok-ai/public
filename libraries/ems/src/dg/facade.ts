/* `domains.table(address)` — one await, then a table handle everything else hangs off (GOAL
   "What it looks like"): sources over the table, and the per-table registries an app fills once —
   actions (ribbon, row hover, context menu), the renderer (lists, pickers, chips) and validators
   (forms). Registers nothing with the platform: rendering reuses the handler `ObjectHandler.forEntity`
   already resolves for the table's rows. `domains.*` is also where every domain control is made
   (`domains.form(src)`, `domains.list(src)`, …): one object, each member constructing at call
   time, so the modules behind it may import this one back. `TRow` (an app's own row type) types
   the rows every registry and source hands out; the default reads any column as unknown. */
import {emsBackends} from '../core/domain-backend.js';
import * as DG from 'datagrok-api/dg';
import {Access, Filters, KIND, Rows, divV, span, timestamp, text, isEmpty, signal} from '@datagrok-libraries/u2';
import type {Capability, IProperty, ObjectRenderer, Action, ColumnOf, DomainRowLike, RowValues, RowView,
  ReadonlySignal, Signal, FilterCondition, FilterOperatorParameter} from '@datagrok-libraries/u2';
import type {DomainTableInfoLike, DomainTableLike} from '../core/domain-backend.js';
import {DomainSource} from '../core/domain-source.js';
import type {DomainSourceOptions} from '../core/domain-source.js';
import type {EditState} from '../core/edit-state.js';
import {SharedSession, confirmDiscard} from '../core/session.js';
import {handlerRenderer, appView, mountView, ObjectForm} from '@datagrok-libraries/u2/src/dg/index.js';
import {DomainApp} from './app.js';
import {authoring} from './authoring/index.js';
import {SYSTEM_COLUMNS} from './backend.js';
import {saveButton, discardButton, newButton} from './buttons.js';
import {DomainChildren} from './children.js';
import {DomainErrors} from './errors.js';
import {DomainFilters} from './filters.js';
import {DomainForm} from './form.js';
import {DomainGrid, DomainDataTable} from './grid.js';
import {DomainHistory} from './history.js';
import {openImport, importRows} from './import.js';
import {bulkEdit} from './bulk.js';
import {BindingDialog, BindingWizard, createBinding} from './authoring/binding-dialog.js';
import {BindingEditDialog, editBinding} from './authoring/binding-edit-dialog.js';
// the cycles with the control modules are function-scoped only: none reads another at load time
import type {DomainAppOptions} from './app.js';
import type {DomainChildrenOptions} from './children.js';
import type {DomainFiltersOptions} from './filters.js';
import type {DomainFormOptions, DomainFormTarget} from './form.js';
import type {DomainGridOptions, DomainDataTableOptions} from './grid.js';
import type {DomainHistoryTarget} from './history.js';
import {DomainList} from './list.js';
import type {DomainListOptions} from './list.js';
import {DomainPick} from './pick.js';
import type {DomainPickOptions} from './pick.js';
import {DomainSearch} from './search.js';
import type {DomainSearchOptions} from './search.js';
import {DomainTree} from './tree.js';
import type {DomainTreeOptions} from './tree.js';
import {route} from './routes.js';
import {DomainViews} from './views.js';
import type {DomainViewsOptions} from './views.js';
export type {CreateBindingOptions, BindingResult, AccessOutcome} from './authoring/binding-dialog.js';
export type {EditBindingResult} from './authoring/binding-edit-dialog.js';
const REF_ADDRESS = /^\w+\.\w+$/;
const DERIVE_MS = 60000;

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

/** A row's problem with one of its cells — a derivation that failed on an imported row; `skip` where it passed the
 * row over ({@link DomainSkip}). */
export interface RowProblem {
  index: number;
  column: string;
  problem: string;
  skip?: boolean;
}

/** What a derivation throws for a row an import passes over rather than refuses — one the table holds already. A
 * form shows it as any other failure: on the field, Save refused. */
export class DomainSkip extends Error {}

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

export class Registry<T> {
  protected readonly _items: T[] = [];

  /** Returns the unregister function. */
  add(item: T): () => void {
    this._items.push(item);
    return () => {
      const at = this._items.indexOf(item);
      if (at >= 0)
        this._items.splice(at, 1);
    };
  }

  get all(): readonly T[] {
    return this._items;
  }
}

export class ActionRegistry<TRow extends DomainRowLike = DomainRowLike> extends Registry<DomainAction<TRow>> {
  /** The actions that apply to `row`, bound to it — what `rowActions` and a menu take. */
  for(row: RowView<TRow>): Action[] {
    return ActionRegistry.bind(this._items, row);
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

/** Columns the browser computes from other columns of the row: a write hook for a schema without one. The
 * server runs no hooks of its own (the EMS design leaves them out), so a derivation holds only for writes made
 * through these handles — a client writing over the REST API computes the outputs itself. */
export interface DomainDerivation<TRow extends DomainRowLike = DomainRowLike> {
  from: ColumnOf<TRow>[];
  to: ColumnOf<TRow>[];
  /** An output no two rows may share: an import passes over a row repeating an earlier row's. */
  unique?: ColumnOf<TRow>;
  derive(row: RowView<TRow>): RowValues<TRow> | Promise<RowValues<TRow>>;
}

export class DerivationRegistry<TRow extends DomainRowLike = DomainRowLike> extends Registry<DomainDerivation<TRow>> {
  get targets(): string[] {
    return this._items.flatMap((d) => d.to);
  }

  get columns(): string[] {
    return this._items.flatMap((d) => [...d.from, ...d.to]);
  }
}

/** A filter operator the app answers for a column the server cannot search that way — a structure
 * search: offered on the columns of its semType, it reaches the server as the ids of the rows it matches. */
export interface DomainHostOperator {
  id: string;
  label: string;
  semType: string;
  parameter?: FilterOperatorParameter;
  rows(condition: FilterCondition, signal: AbortSignal): Promise<readonly string[]>;
}

export class HostOperatorRegistry extends Registry<DomainHostOperator> {
  /** Registered with the page's operators, which every table of the bundle shares: not `exclusive`
   * there, or a column of another table this one's semType names would lose its core operators —
   * `DomainFilters` narrows the column of a table that answers it to this set. */
  add(op: DomainHostOperator): () => void {
    const remove = super.add(op);
    const unregister = Filters.operators.register({id: op.id, label: op.label, arity: 1, kinds: [KIND.STRING],
      semType: op.semType, editor: 'default', host: true, parameter: op.parameter});
    return () => {
      remove();
      unregister();
    };
  }
}

export class DomainTable<TRow extends DomainRowLike = DomainRowLike> {
  /** Everything an app declares about the table's actions, once. */
  readonly actions = new ActionRegistry<TRow>();
  readonly validators = new ValidatorRegistry<TRow>();
  readonly derived = new DerivationRegistry<TRow>();
  readonly operators = new HostOperatorRegistry();
  /** The platform handler of the table's rows — deep links, the entity view, `rowFrom`. */
  readonly handler: DG.DomainObjectHandler;
  /** How lists, pickers and chips show a row; the handler's own rendering by default. */
  renderer: ObjectRenderer<RowView<TRow>>;
  /** The renderer the handle was built with — what tells a table an app dressed from a bare one. */
  readonly defaultRenderer: ObjectRenderer<RowView<TRow>>;

  private static readonly _bySource = new WeakMap<DomainSource, DomainTable>();
  private static readonly _byRow = new WeakMap<ReadonlySignal<RowView | null>, DomainSource>();
  /** Per source, why each derivation last failed on a row, and a tick that moves with them. */
  private readonly _failed = new WeakMap<DomainSource,
    {failures: Map<DomainDerivation<TRow>, Map<string, string>>, changed: Signal<number>}>();

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

  /** A name as auto-mapping compares it: case, punctuation and a unit in parentheses ignored, so a file's
   * "Purity" finds "Purity (%)". */
  static matchKey(name: string): string {
    return name.toLowerCase().replace(/\(.*\)/, '').replace(/[^a-z0-9]/g, '');
  }

  /** In the column's display format: an autonumbered key reads DG-000123. */
  static captionText(properties: IProperty[], column: string, value: unknown): string {
    const prop = properties.find((p) => p.name === column);
    return prop === undefined ? text(value) : ObjectForm.displayText(prop, value);
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
   * of the app's views of it ({@link DomainViews.openRow}), else the platform's entity view. */
  open(row: RowView<TRow>): void {
    if (Rows.isDraft(row) || DomainViews.openRow(this.address, row.id))
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
    mountView({view, base, pinned: path !== undefined, session,
      confirmLeave: () => confirmDiscard(session, {action: 'close the view'}),
      params: ['entity', 'q', 'search', 'trash'],
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

  /** What a row looks like from its schema alone: the name column, else the business key (a reference
   * by its target's caption), else the id; "New <singular>" for a draft without one. A card adds one
   * muted description line and the creation time. */
  static schemaRenderer(schema: () => {properties: IProperty[], info: DomainTableInfoLike}): ObjectRenderer<RowView> {
    const label = (row: RowView): {text: string, draft: boolean} => {
      const {properties, info} = schema();
      const name = (info.nameColumn === null ? info.businessKey : [info.nameColumn]).map((column) => {
        const prop = properties.find((p) => p.name === column);
        return prop !== undefined && DomainTable.isReference(prop) ? text(row[Rows.caption(column)]) :
          isEmpty(row[column]) ? '' : DomainTable.captionText(properties, column, row[column]);
      }).filter((part) => part !== '').join(' · ');
      if (name !== '')
        return {text: name, draft: false};
      return Rows.isDraft(row) ? {text: `New ${DomainSource.noun(info, false)}`, draft: true} :
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
    source.hosted = (condition, signal) => {
      const op = this.operators.all.find((o) => o.id === condition.operator);
      if (op === undefined)
        throw new Error(`${this.address}: no operator "${condition.operator}"`);
      return op.rows(condition, signal);
    };
    this._deriving(source);
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

  /** Runs {@link derived} over plain rows in place; a row a derivation failed on keeps its outputs unset,
   * so the server refuses it by name. */
  async deriveRows(rows: Record<string, unknown>[]): Promise<RowProblem[]> {
    const failures: RowProblem[] = [];
    for (const derivation of this.derived.all) {
      await Promise.all(rows.map(async (row, index) => {
        if (DomainTable._noInputs(derivation, row))
          return;
        try {
          Object.assign(row, await DomainTable._timed(derivation.derive(row as RowView<TRow>)));
        } catch (e) {
          failures.push({index, column: derivation.from[0], problem: DomainErrors.message(e),
            skip: e instanceof DomainSkip});
        }
      }));
      if (derivation.unique === undefined)
        continue;
      const failed = new Set(failures.map((f) => f.index));
      const first = new Map<unknown, number>();
      const noun = DomainSource.noun(this.info, false);
      for (const [index, row] of rows.entries()) {
        const value = row[derivation.unique] ?? null;
        const earlier = first.get(value);
        if (value === null || failed.has(index))
          continue;
        if (earlier === undefined)
          first.set(value, index);
        else
          failures.push({index, column: derivation.from[0], problem: `Same ${noun} as row ${earlier + 1}`, skip: true});
      }
    }
    return failures.sort((a, b) => a.index - b.index);
  }

  /** Why the derivation reading `column` first failed on the row of `source`, null while none did — the
   * same words a refused Save names that column with. Read under an effect, it re-runs as a derivation
   * ends, so a form shows the failure on the field then rather than at Save. */
  derivationProblem(source: DomainSource, id: string, column: string): string | null {
    const state = this._failed.get(source);
    if (state === undefined)
      return null;
    state.changed.value;
    for (const d of this.derived.all) {
      const problem = d.from[0] === column ? state.failures.get(d)?.get(id) : undefined;
      if (problem !== undefined)
        return problem;
    }
    return null;
  }

  /** Derives a row again whenever its inputs, or the outputs last written for them, change; a result lands
   * only while those inputs and the writer still stand, and a failure refuses Save until the next good run. */
  private _deriving(source: DomainSource<TRow>): void {
    const last = new Map<DomainDerivation<TRow>, Map<string, {inputs: string, outputs: string | null}>>();
    const failures = new Map<DomainDerivation<TRow>, Map<string, string>>();
    const changed = signal(0);
    this._failed.set(source, {failures, changed});
    let writing: {d: DomainDerivation<TRow>, id: string} | null = null;
    const of = <T>(map: Map<DomainDerivation<TRow>, Map<string, T>>, d: DomainDerivation<TRow>) => {
      let rows = map.get(d);
      if (rows === undefined)
        map.set(d, rows = new Map<string, T>());
      return rows;
    };
    const fail = (d: DomainDerivation<TRow>, id: string, problem: string | null) => {
      const rows = of(failures, d);
      if ((rows.get(id) ?? null) === problem)
        return;
      if (problem === null)
        rows.delete(id);
      else
        rows.set(id, problem);
      changed.value = changed.peek() + 1;
    };
    const snapshot = (values: Record<string, unknown>, columns: string[]) =>
      JSON.stringify(columns.map((c) => values[c] ?? null));
    const write = (edit: EditState, row: RowView<TRow>, d: DomainDerivation<TRow>, inputs: string,
      values: Record<string, unknown>) => {
      const current = row as Record<string, unknown>;
      of(last, d).set(row.id, {inputs, outputs: snapshot({...current, ...values}, d.to)});
      writing = {d, id: row.id};
      try {
        for (const [column, value] of Object.entries(values)) {
          if ((current[column] ?? null) !== (value ?? null))
            edit.setValue(row.id, column, value);
        }
      } finally {
        writing = null;
      }
    };
    const run = (edit: EditState, row: RowView<TRow>, d: DomainDerivation<TRow>) => {
      if (writing?.d === d && writing.id === row.id)
        return;
      const current = row as Record<string, unknown>;
      const inputs = snapshot(current, d.from);
      const prior = of(last, d).get(row.id);
      if (prior?.inputs === inputs && (prior.outputs === null || prior.outputs === snapshot(current, d.to)))
        return;
      // a verdict is about the inputs it was given: the new ones have none until their run ends
      fail(d, row.id, null);
      if (DomainTable._noInputs(d, current)) {
        write(edit, row, d, inputs, Object.fromEntries(d.to.map((c) => [c, null])));
        return;
      }
      of(last, d).set(row.id, {inputs, outputs: null});
      const standing = () => of(last, d).get(row.id)?.inputs === inputs && source.edit.peek() === edit;
      source.hold(DomainTable._timed(Promise.resolve().then(() => d.derive(row))).then((out) => {
        if (!standing())
          return;
        write(edit, row, d, inputs, out as Record<string, unknown>);
      }, (e) => {
        if (standing())
          fail(d, row.id, DomainErrors.message(e));
      }));
    };
    source.effect(() => {
      const edit = source.edit.value;
      if (edit === undefined)
        return;
      const changes = edit.onChanged.subscribe((key) => {
        const changed = key === null ? source.pending() : [source.rows.byKey(key)];
        for (const row of changed) {
          if (row === undefined || row[Rows.STATE] === 'deleted')
            continue;
          for (const d of this.derived.all)
            run(edit, row, d);
        }
      });
      // a landed batch re-keys its drafts
      const saves = edit.onSaved.subscribe(() => {
        last.clear();
        failures.clear();
        changed.value = changed.peek() + 1;
      });
      return () => {
        changes.unsubscribe();
        saves.unsubscribe();
      };
    });
    source.guard(() => {
      for (const row of source.pending()) {
        for (const d of this.derived.all) {
          const problem = failures.get(d)?.get(row.id);
          if (problem === undefined)
            continue;
          source.markProblem(row.id);
          return this._nameCell(row, d.from[0], problem);
        }
      }
      return null;
    });
  }

  private static _noInputs(d: DomainDerivation, row: Record<string, unknown>): boolean {
    return d.from.every((c) => isEmpty(row[c]));
  }

  /** A derivation that never answers would hold every Save: it fails after {@link DERIVE_MS}. */
  private static _timed<T>(work: Promise<T> | T): Promise<T> {
    let timer: ReturnType<typeof setTimeout> | undefined;
    const expired = new Promise<never>((_, reject) => timer = setTimeout(() =>
      reject(new Error(`did not finish in ${DERIVE_MS / 1000} s`)), DERIVE_MS));
    return Promise.race([Promise.resolve(work), expired]).finally(() => clearTimeout(timer));
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
   * needs them — the writer's own "Value can't be empty" says nothing about which field. A column the
   * browser derives is named by the input it is derived from: that is the field the user can change. */
  private _nameCell(row: RowView<TRow>, column: string, problem: string): string {
    const input = this.derived.all.find((d) => d.to.includes(column as ColumnOf<TRow>))?.from[0] ?? column;
    return `${this.renderer.caption(row)}: ${DomainTable._caption(this.properties, input)}: ${problem}`;
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
    const backend = emsBackends.domain;
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
  /** Several tables as one app: a view per table under `path`, the router's answer, the tree's items, and where
   * a row of those tables opens. */
  views: (options: DomainViewsOptions): DomainViews => new DomainViews(options),
  route,
  bulkEdit,
  /** The import wizard over a table: source → mapping → preview → the server's report. */
  import: openImport,
  /** The import without the wizard: map, resolve references, derive, validate, batch; resolves to the report frame. */
  importRows,
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
