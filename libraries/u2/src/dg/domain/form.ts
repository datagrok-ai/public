/* `domains.form` — `propertyForm` over a domain row: the fields come from the table's schema,
   the access decides which are inputs and which are text, every write goes through the source's
   `EditState` (the row proxy's `set`), and validity is the union of the field rules, the table's
   validators and what the writer reports per cell. The form follows its row: a new current row
   is a new form under a fresh scope, the previous one released. Paired through the source: it
   guards the source's Save (a refusal names the field), takes the focus back when the session
   saved or discarded, and when a list bumps `activate`. A refused save is one balloon (the
   platform editor maps the server's column errors onto the cells itself). */
import * as grok from 'datagrok-api/grok';
import type * as DG from 'datagrok-api/dg';
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {Signal, ReadonlySignal, computed, signal, untracked} from '../../core/signals.js';
import {Rows} from '../../sources/rows-like.js';
import {Access} from '../../core/access.js';
import type {Input, InputOptions} from '../../core/input-base.js';
import type {IProperty} from '../../core/property-like.js';
import {Filters} from '../../core/filter/index.js';
import type {FilterGroup, FilterProperty, FilterSchema, MaskColumnLike, MaskFrameLike}
  from '../../core/filter/index.js';
import {div, span, timestamp} from '../../core/elements.js';
import {text} from '../../core/text.js';
import {loader} from '../../components/display/async-view.js';
import type {IWidgetStatus} from '../../core/widget-like.js';
import {propertyForm, ObjectForm} from '../forms/object-form.js';
import type {ObjectFormOptions} from '../forms/object-form.js';
import {Editors} from '../forms/editors.js';
import {DomainSource} from '../../sources/domain-source.js';
import type {RowView} from '../../sources/rows-like.js';
import {userInput} from '../inputs/user-input.js';
import {entityInput} from '../entities/entity.js';
import {DomainTable} from './index.js';
import {SYSTEM_COLUMNS} from './backend.js';
import {DomainErrors} from './errors.js';
import {DomainPick, PickInput} from './pick.js';
import type {DomainPickOptions} from './pick.js';

export type DomainFormTarget = RowView | ReadonlySignal<RowView | null> | DomainSource;

/** A schema constraint the form pre-validates: the parsed `expr` and the columns it names, the
 * first the field its message attaches to. */
interface Constraint {
  name: string;
  message: string;
  root: FilterGroup;
  columns: string[];
}

export interface DomainFormOptions extends Pick<ObjectFormOptions,
  'include' | 'exclude' | 'layout' | 'condensed' | 'overrides' | 'onChanged'> {
  /** What the caller may see and edit; the source's access by default. */
  access?: Access;
  /** The source a bare row or a row signal belongs to — needed unless the signal is the
   * `currentRow` of a source made through a `DomainTable`. */
  source?: DomainSource;
  /** What the form says while there is no row (default: "Select a <row> to edit."). */
  empty?: string;
  /** The system columns (id, version, created, updated, author): left out by default, or a
   * compact muted block after the fields. */
  system?: 'hidden' | 'footer';
}

const UNTOUCHED = 'u2-input-untouched';
const CHANGED = 'u2-input-changed';
const REQUIRED = 'u2-input-required';
const EMPTY = 'Value can\'t be empty';

/** The keys that save. `DomainApp` listens for them across the whole view, since `appView` puts
 * the ribbon — where Tab from the last field lands — outside the app root. */
export const SAVE_SHORTCUTS: Record<string, string> = {'Ctrl+S': 'Save', 'Ctrl+Enter': 'Save'};

/** True for a {@link SAVE_SHORTCUTS} keystroke. */
export function isSaveKey(e: KeyboardEvent): boolean {
  return (e.ctrlKey || e.metaKey) && (e.key === 's' || e.key === 'S' || e.key === 'Enter');
}

export class DomainForm extends Control {
  readonly source: DomainSource;
  readonly row: ReadonlySignal<RowView | null>;
  /** The current form's first problem, null while every field passes or there is no row. */
  readonly validity: ReadonlySignal<string | null>;
  /** The first refused field, named by its caption — what a refused save says. */
  readonly problem: ReadonlySignal<string | null>;

  private readonly _form = signal<ObjectForm | null>(null);
  private static readonly _parsed = new WeakMap<DomainTable, Constraint[]>();
  /** The schema constraints the current form checks; a violated one's message by name. */
  private _constraints: Constraint[] = [];
  private readonly _violations = new Map<string, string>();
  private readonly _violationsVersion = signal(0);
  private _attached = new Set<string>();
  private _checkGen = 0;
  private readonly _verdict: HTMLElement;
  private _shown: Scope | undefined;
  /** The system block under the fields, re-read with the row: the writer lands the server-assigned
   * columns (`number`, `version`, `created_on`, `author_id`) after the insert, not with it. */
  private _system: HTMLElement | null = null;
  /** Per readonly reference column of the current row: the id shown and the caption resolved for it. */
  private readonly _captions = new Map<string, {id: string, caption: string}>();
  private _editSub: {unsubscribe(): void} | undefined;
  private readonly _subs: {unsubscribe(): void}[];

  constructor(target: DomainFormTarget, private readonly _options: DomainFormOptions = {}) {
    super();
    this.root.classList.add('u2-domain-form');
    this.root.dataset.u2 = 'domain-form';
    const {source, row} = DomainForm.resolve(target, _options.source);
    this.source = source;
    this.row = row;
    this.validity = computed(() => this._form.value?.validity.value ?? null);
    this.problem = computed(() => {
      for (const input of this._form.value?.inputs ?? []) {
        const message = input.validity.value;
        if (message !== null) {
          const caption = `${input.label.charAt(0).toUpperCase()}${input.label.slice(1)}`;
          return message === EMPTY ? `${caption} is required` : `${caption}: ${message}`;
        }
      }
      this._violationsVersion.value;
      for (const c of this._constraints) {
        const message = this._violations.get(c.name);
        if (message !== undefined && !this._attached.has(c.name))
          return message;
      }
      return null;
    });
    // a violated constraint whose first column has no field tells under the form
    this._verdict = span('', 'u2-input-error');
    this._verdict.dataset.u2Part = 'constraints';
    this.effect(() => {
      this._violationsVersion.value;
      const lines = this._constraints.filter((c) => !this._attached.has(c.name))
        .map((c) => this._violations.get(c.name)).filter((m): m is string => m !== undefined);
      this._verdict.textContent = lines.join(' ');
      this._verdict.hidden = lines.length === 0;
    });
    // the keyboard path to Save: the ribbon is the shell's, outside the form's tab order — Ctrl+S
    // saves, and Tab past the last field lands on the session's Save button
    const onKeyDown = (e: KeyboardEvent) => {
      if (e.key === 'Tab' && !e.shiftKey) {
        const target = source.session.primaryButton as HTMLButtonElement | undefined;
        const last = this._form.peek()?.inputs.at(-1);
        if (target !== undefined && !target.disabled && last !== undefined && last.root.contains(e.target as Node)) {
          e.preventDefault();
          target.focus();
        }
        return;
      }
      if (!isSaveKey(e))
        return;
      e.preventDefault();
      if (source.session.isDirty.peek())
        void source.session.save();
    };
    this.root.addEventListener('keydown', onKeyDown);
    this.own(() => this.root.removeEventListener('keydown', onKeyDown));
    // paired through the source: Save is refused with the field named, and the focus comes back
    // after a save or a discard, and when a list hands its row over
    const unguard = source.guard(() => {
      const problem = this.problem.peek();
      if (problem !== null)
        this.validate();
      return problem;
    });
    this._subs = [source.session.onSaved.subscribe(() => this.focus()),
      source.session.onDiscarded.subscribe(() => this.focus())];
    let seen = source.activate.peek();
    this.effect(() => {
      const bumped = source.activate.value;
      if (bumped === seen)
        return;
      seen = bumped;
      this.focus();
    });
    this.own(() => {
      unguard();
      for (const s of this._subs)
        s.unsubscribe();
      this._shown?.dispose();
      this._editSub?.unsubscribe();
    });
    this.effect(() => {
      const current = this.row.value;
      const access = _options.access ?? source.access.value;
      untracked(() => this._show(current, access));
    });
    // the editing state moves under the form — a discard, a save, another control's write
    this.effect(() => {
      const edit = source.edit.value;
      this._editSub?.unsubscribe();
      this._editSub = edit?.onChanged.subscribe(() => this._sync());
    });
  }

  /** The form over the current row, null while there is none. */
  get form(): ObjectForm | null {
    return this._form.peek();
  }

  /** The input editing a column of the current row; none for a text field or without a row. */
  input(name: string): ReturnType<ObjectForm['input']> {
    return this._form.peek()?.input(name);
  }

  /** Every field's verdict shown, the first failing one focused — the explicit check a Save runs. */
  validate(): boolean {
    const form = this._form.peek();
    if (form === null)
      return false;
    for (const input of form.inputs)
      input.root.classList.remove(UNTOUCHED);
    this._checkConstraints();
    return form.validate();
  }

  /** Re-reads every field off the row. */
  refresh(): void {
    this._sync();
  }

  getWidgetStatus(): IWidgetStatus {
    return {...super.getWidgetStatus(), shortcuts: {...SAVE_SHORTCUTS}};
  }

  /** Focus to the first required field still empty, else the first — where a create form starts
   * (past what the defaults filled), and where Save returns. Deferred past the rebuild a row
   * change (a saved draft under its new key) may be running. */
  focus(): void {
    queueMicrotask(() => {
      this._form.peek()?.focusFirst(DomainForm.needsValue);
      // a Ctrl+S from inside a picker leaves its popup over the form, and the re-focus would open
      // one anyway (`openOnFocus`): the form comes back with nothing covering it
      for (const input of this._form.peek()?.inputs ?? []) {
        if (input instanceof PickInput)
          input.typeAhead.close();
      }
    });
  }

  /** A required field without a value — what the focus goes to first. */
  static needsValue(input: Input<any>): boolean {
    const value = input.value.peek();
    return !input.nullable && (value === null || value === undefined || value === '');
  }

  /** The source and the row signal behind a target: a source's current row, a signal a
   * `DomainTable` source handed out, or a bare row under `source`. */
  static resolve(target: DomainFormTarget, source?: DomainSource):
    {source: DomainSource, row: ReadonlySignal<RowView | null>} {
    if (target instanceof DomainSource)
      return {source: target, row: target.currentRow};
    if (target instanceof Signal) {
      const row = target as ReadonlySignal<RowView | null>;
      const owner = source ?? DomainTable.sourceOf(row);
      if (owner === undefined)
        throw new Error('domains.form: pass `source` — the row signal alone does not name its table');
      return {source: owner, row};
    }
    if (source === undefined)
      throw new Error('domains.form: pass `source` — the row alone does not name its table');
    return {source, row: signal(target as RowView | null)};
  }

  /** A `User` column's editor: the value is the user id, the box the user picker. */
  static userPick(options: InputOptions<string | null>): PickInput<DG.User> {
    return new PickInput<DG.User>({...options, typeAhead: () => userInput(), idOf: (user) => user.id,
      resolve: (id) => grok.dapi.users.find(id)});
  }

  static groupPick(options: InputOptions<string | null>): PickInput<DG.Group> {
    return new PickInput<DG.Group>({...options, typeAhead: () => entityInput('Group…', () => grok.dapi.groups),
      idOf: (group) => group.id, resolve: (id) => grok.dapi.groups.find(id)});
  }

  private _show(row: RowView | null, access: Access): void {
    this._shown?.dispose();
    const scope = new Scope();
    this._shown = scope;
    Scope.runWith(scope, () => {
      if (row === null) {
        this._form.value = null;
        this._system = null;
        this._constraints = [];
        this._violations.clear();
        this._violationsVersion.value = this._violationsVersion.peek() + 1;
        this.root.replaceChildren(this._hint(scope));
        return;
      }
      const form = this._build(row, access);
      this._form.value = form;
      const system = this._options.system === 'footer' ? this._footer(row) : null;
      this._system = system;
      this.root.replaceChildren(form.root, this._verdict, ...(system === null ? [] : [system]));
      this._captions.clear();
      this._resolveCaptions(form, row, scope, system === null ? [] : DomainForm.systemProps(this.source));
      this._checkConstraints();
      // a draft is where typing starts — through `focus`, which leaves the pickers' popups shut:
      // a dropdown the user did not open is in the way, not an invitation
      if (Rows.isDraft(row))
        this.focus();
    });
  }

  private _build(row: RowView, access: Access): ObjectForm {
    const source = this.source;
    const options = this._options;
    const table = DomainTable.of(source);
    const draft = Rows.isDraft(row);
    const system = new Set(SYSTEM_COLUMNS.map(([name]) => name));
    const view = access.row(row);
    // a ref column with a schema filter picks among the rows its siblings allow
    let overrides = options.overrides;
    for (const [name, filter] of Object.entries(table?.info.refFilters ?? {})) {
      const extra: Partial<DomainPickOptions> = {filter, params: () => row, siblings: source.schema};
      overrides = {...overrides, [name]: {...overrides?.[name], ...extra}};
    }
    const form = propertyForm(source.schema.properties.filter((p) => !system.has(p.name)), row, {
      include: options.include, exclude: options.exclude, layout: options.layout, condensed: options.condensed,
      overrides, access: view,
      onChanged: (name, value) => {
        options.onChanged?.(name, value);
        this._checkConstraints();
      },
    });
    form.root.dataset.u2Row = row.id;
    // the schema constraints, pre-validated here; one naming a column the caller may not see is
    // left to the server
    this._constraints = DomainForm.constraintsOf(table, source.schema)
      .filter((c) => c.columns.every((name) => view.field(name) !== 'hidden'));
    this._violations.clear();
    this._attached = new Set();
    for (const c of this._constraints) {
      const input = form.input(c.columns[0]);
      if (input === undefined)
        continue;
      this._attached.add(c.name);
      input.addValidator(() => this._violations.get(c.name) ?? null);
    }
    for (const input of form.inputs) {
      const name = input.name;
      if (name === undefined)
        continue;
      input.root.classList.toggle(REQUIRED, !input.nullable);
      input.addValidator((value) =>
        source.edit.peek()?.errorOf(row.id, name) ?? table?.validators.check(name, value, row) ?? null);
      // a pristine draft is not wrong yet: a field's verdict shows once it was edited or left
      if (draft) {
        input.root.classList.add(UNTOUCHED);
        const touched = () => input.root.classList.remove(UNTOUCHED);
        let initial = true;
        input.effect(() => {
          input.value.value;
          if (initial)
            initial = false;
          else
            touched();
        });
        input.root.addEventListener('focusout', touched);
        input.own(() => input.root.removeEventListener('focusout', touched));
      }
    }
    this._markChanged(form, row);
    return form;
  }

  /** The amber edge on every field the writer holds a change for — an edited row's cells; a
   * draft is new whole, so none of its. */
  private _markChanged(form: ObjectForm, row: RowView): void {
    const edit = this.source.edit.peek();
    const draft = Rows.isDraft(row);
    for (const input of form.inputs) {
      if (input.name !== undefined)
        input.root.classList.toggle(CHANGED, !draft && edit !== undefined && edit.isChanged(row.id, input.name));
    }
  }

  /** The system columns as a compact muted block: captions, local timestamps, the author as a
   * caption once resolved, and never a draft's key. */
  private _footer(row: RowView): HTMLElement {
    const draft = Rows.isDraft(row);
    const block = div([], 'u2-domain-form-system');
    block.dataset.u2Part = 'system';
    for (const prop of DomainForm.systemProps(this.source)) {
      const name = prop.name!;
      // a draft has none of these yet: an empty Version / Created / Author line says nothing
      if (draft && name !== 'id')
        continue;
      const caption = SYSTEM_COLUMNS.find(([column]) => column === name)?.[2] ?? name;
      const raw = row[name];
      const value = span('', 'u2-form-readonly-value');
      value.dataset.u2Part = 'readonly-value';
      if (name === 'id' && draft)
        value.textContent = 'assigned on save';
      else if ((prop.propertyType ?? prop.type) === 'datetime' && raw !== null && raw !== undefined && raw !== '')
        value.append(timestamp(raw as Date | number | string, undefined, {utcDates: true}));
      else
        value.textContent = text(raw);
      const line = div([span(caption, 'u2-input-label'), value], 'u2-form-readonly');
      line.dataset.u2 = 'readonly-field';
      line.dataset.u2Name = name;
      block.append(line);
    }
    return block;
  }

  /** The system columns the schema carries, in the platform's order. */
  static systemProps(source: DomainSource): IProperty[] {
    const props: IProperty[] = source.schema.properties;
    const found: IProperty[] = [];
    for (const [name] of SYSTEM_COLUMNS) {
      const prop = props.find((p) => p.name === name);
      if (prop !== undefined)
        found.push(prop);
    }
    return found;
  }

  /** A readonly reference shows what it points at, not the uuid it holds: resolved the way the
   * picker resolves a preset id, the id standing in until then, under the row's scope. */
  private _resolveCaptions(form: ObjectForm, row: RowView, scope: Scope, extra: IProperty[]): void {
    let live = true;
    scope.own(() => live = false);
    for (const prop of [...form.properties, ...extra]) {
      const name = prop.name!;
      const id = row[name];
      if (form.input(name) !== undefined || !DomainTable.isReference(prop) ||
          id === null || id === undefined || id === '' || this._captions.get(name)?.id === String(id))
        continue;
      void DomainForm.captionOf(prop, String(id)).then((caption) => {
        if (!live || caption === null)
          return;
        this._captions.set(name, {id: String(id), caption});
        this._applyCaptions(row);
      }, () => undefined);
    }
  }

  private _applyCaptions(row: RowView): void {
    for (const [name, {id, caption}] of this._captions) {
      if (String(row[name] ?? '') !== id)
        continue;
      const el = this.root.querySelector<HTMLElement>(`[data-u2-name="${name}"] [data-u2-part="readonly-value"]`);
      if (el !== null)
        el.textContent = caption;
    }
  }

  /** The display name behind a reference value: the target table's name column, a user's or a
   * group's friendly name; null when the platform does not answer. */
  static captionOf(prop: IProperty, id: string): Promise<string | null> {
    const semType = prop.semType ?? '';
    if (semType === 'User')
      return grok.dapi.users.find(id).then((user) => user?.friendlyName ?? null);
    if (semType === 'Group')
      return grok.dapi.groups.find(id).then((group) => group?.friendlyName ?? null);
    return DomainPick.resolve(semType, id).then((item) => item.name);
  }

  /** What stands in for the form while there is no row: the load in progress, its failure, or
   * the invitation to pick one. */
  private _hint(scope: Scope): HTMLElement {
    const source = this.source;
    const el = div([], 'u2-domain-form-empty');
    el.dataset.u2Part = 'empty';
    scope.effect(() => {
      const state = source.state.value;
      if (state === 'loading')
        el.replaceChildren(loader('Loading…'));
      else if (state === 'error')
        el.replaceChildren(span(DomainErrors.message(source.error.value), 'u2-domain-form-error'));
      else
        el.replaceChildren(span(this._options.empty ?? DomainForm.emptyHint(source), 'u2-domain-form-hint'));
    });
    return el;
  }

  private _sync(): void {
    const form = this._form.peek();
    const row = this.row.peek();
    if (form === null || row === null)
      return;
    form.refresh();
    for (const input of form.inputs)
      input.revalidate();
    this._markChanged(form, row);
    if (this._system !== null) {
      const system = this._footer(row);
      this._system.replaceWith(system);
      this._system = system;
      // the insert stamps `author_id` and the rest only now: their captions are resolved here,
      // or the footer keeps the uuid the writer landed
      if (this._shown !== undefined)
        this._resolveCaptions(form, row, this._shown, DomainForm.systemProps(this.source));
    }
    // a refresh re-reads the uuid into the text row; the caption resolved for it is put back
    this._applyCaptions(row);
    this._checkConstraints();
  }

  /** The table's grammar constraints, parsed once against the schema; one that does not parse
   * or names an unknown column is skipped — the server still enforces it. */
  static constraintsOf(table: DomainTable | undefined, schema: FilterSchema): Constraint[] {
    if (table === undefined)
      return [];
    let parsed = DomainForm._parsed.get(table);
    if (parsed !== undefined)
      return parsed;
    parsed = [];
    for (const c of table.info.constraints) {
      const {root, problems} = Filters.parse(c.expr, schema);
      if (problems.length > 0)
        continue;
      const columns: string[] = [];
      Filters.walk(root, (n) => {
        if (Filters.isGroup(n))
          return;
        const refs = (Array.isArray(n.value) ? n.value : [n.value]).filter(Filters.isColumnRef).map((v) => v.column);
        for (const name of [n.property, ...refs]) {
          if (!columns.includes(name))
            columns.push(name);
        }
      });
      parsed.push({name: c.name, message: c.message ?? `Must satisfy: ${c.expr}`, root, columns});
    }
    DomainForm._parsed.set(table, parsed);
    return parsed;
  }

  /** Evaluates every constraint over the row as the fields hold it, with the null semantics a SQL
   * CHECK has ({@link Filters.toCheckMask}): an empty column makes its own comparison unknown, not
   * the whole constraint — `a > 0 and b > 0` with a = -1 and b empty is still refused. The
   * verdicts land asynchronously on the first named column's field, or under the form. */
  private _checkConstraints(): void {
    const form = this._form.peek();
    const row = this.row.peek();
    const constraints = this._constraints;
    if (form === null || row === null || constraints.length === 0)
      return;
    const values: Record<string, unknown> = {...row, ...form.getValues()};
    const frame = DomainForm.frameOf(values, this.source.schema.properties);
    const gen = ++this._checkGen;
    void Promise.all(constraints.map(async (c): Promise<string | null> => {
      try {
        return (await Filters.toCheckMask(frame, c.root)).get(0) ? null : c.message;
      } catch {
        return null;
      }
    })).then((verdicts) => {
      if (gen !== this._checkGen || this._form.peek() !== form)
        return;
      this._violations.clear();
      constraints.forEach((c, i) => {
        if (verdicts[i] !== null)
          this._violations.set(c.name, verdicts[i]!);
      });
      for (const c of constraints)
        form.input(c.columns[0])?.revalidate();
      this._violationsVersion.value = this._violationsVersion.peek() + 1;
    });
  }

  /** One row as `Filters.toMask` reads a frame — the memory backend's column shapes, built per
   * column on demand. */
  static frameOf(values: Record<string, unknown>, properties: FilterProperty[]): MaskFrameLike {
    const columns = new Map<string, MaskColumnLike>();
    return {
      rowCount: 1,
      column: (name) => {
        const prop = properties.find((p) => p.name === name);
        if (prop === undefined)
          return null;
        let column = columns.get(name);
        if (column === undefined)
          columns.set(name, column = DomainForm._maskColumn(prop, values[name]));
        return column;
      },
    };
  }

  private static _maskColumn(prop: FilterProperty, v: unknown): MaskColumnLike {
    const type = prop.propertyType ?? prop.type ?? 'string';
    const nil = v === null || v === undefined || v === '';
    const column = (raw: ArrayLike<number>, extra: Partial<MaskColumnLike> = {}): MaskColumnLike =>
      ({name: prop.name, type, length: 1, getRawData: () => raw, ...extra});
    switch (Filters.kindOf(prop)) {
      case Filters.KIND.INT:
        return column(Int32Array.of(nil ? Filters.INT_NULL : Number(v)));
      case Filters.KIND.FLOAT:
        return column(Float64Array.of(nil ? Filters.FLOAT_NULL : Number(v)));
      case Filters.KIND.DATE_TIME:
        return column(Float64Array.of(nil ? Filters.FLOAT_NULL :
          (v instanceof Date ? v.getTime() : Date.parse(String(v))) * 1000));
      case Filters.KIND.BOOL:
        return column(Uint32Array.of(v === true ? 1 : 0));
      case Filters.KIND.BIG_INT: case Filters.KIND.STRING_LIST:
        return column(new Int32Array(0), {get: () => nil ? null : v});
      default:
        return column(Int32Array.of(0), {categories: [nil ? '' : String(v)]});
    }
  }

  static emptyHint(source: DomainSource): string {
    const name = source.schema.info.singularName.toLowerCase() || 'row';
    return `Select ${/^[aeiou]/.test(name) ? 'an' : 'a'} ${name} to edit.`;
  }
}

// the editors a schema's reference columns get everywhere a form is generated: a picker over the
// target table for a `<schema>.<table>` ref, the platform pickers for a user or a group
Editors.register({
  match: (prop) => DomainTable.isReference(prop) && prop.semType !== 'User' && prop.semType !== 'Group',
  create: (prop, options) => new DomainPick(prop.semType!, options),
});
Editors.register({
  match: (prop) => prop.semType === 'User',
  create: (prop, options) => DomainForm.userPick(options),
});
Editors.register({
  match: (prop) => prop.semType === 'Group',
  create: (prop, options) => DomainForm.groupPick(options),
});
