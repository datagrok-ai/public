/* `u2.domain.form` — `propertyForm` over a domain row: the fields come from the table's schema,
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
import type {InputOptions} from '../../core/input-base.js';
import type {IProperty} from '../../core/property-like.js';
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
import {domainPick, DomainPick, PickInput} from './pick.js';

export type DomainFormTarget = RowView | ReadonlySignal<RowView | null> | DomainSource;

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
const EMPTY = 'Value can\'t be empty';
const SHORTCUTS = {'Ctrl+S': 'Save', 'Ctrl+Enter': 'Save'};

export class DomainForm extends Control {
  readonly source: DomainSource;
  readonly row: ReadonlySignal<RowView | null>;
  /** The current form's first problem, null while every field passes or there is no row. */
  readonly validity: ReadonlySignal<string | null>;
  /** The first refused field, named by its caption — what a refused save says. */
  readonly problem: ReadonlySignal<string | null>;

  private readonly _form = signal<ObjectForm | null>(null);
  private _shown: Scope | undefined;
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
      return null;
    });
    // the keyboard path to Save: the ribbon is the shell's, outside the form's tab order
    const onKeyDown = (e: KeyboardEvent) => {
      if (!(e.ctrlKey || e.metaKey) || !(e.key === 's' || e.key === 'S' || e.key === 'Enter'))
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
    // a refused save; a load failure is the hint's, and the list's, to show
    this.effect(() => {
      const error = source.error.value;
      if (error !== undefined && source.state.peek() !== 'error')
        DomainErrors.report(error);
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
    return form.validate();
  }

  /** Re-reads every field off the row. */
  refresh(): void {
    this._sync();
  }

  getWidgetStatus(): IWidgetStatus {
    return {...super.getWidgetStatus(), shortcuts: {...SHORTCUTS}};
  }

  /** Focus to the first editable field — where a create form starts, and where Save returns.
   * Deferred past the rebuild a row change (a saved draft under its new key) may be running. */
  focus(): void {
    queueMicrotask(() => this._form.peek()?.focusFirst());
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
        throw new Error('domainForm: pass `source` — the row signal alone does not name its table');
      return {source: owner, row};
    }
    if (source === undefined)
      throw new Error('domainForm: pass `source` — the row alone does not name its table');
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
        this.root.replaceChildren(this._hint(scope));
        return;
      }
      const form = this._build(row, access);
      this._form.value = form;
      const system = this._options.system === 'footer' ? this._footer(row) : null;
      this.root.replaceChildren(form.root, ...(system === null ? [] : [system]));
      this._resolveCaptions(form, row, scope, system === null ? [] : DomainForm.systemProps(this.source));
      // a draft is where typing starts — once the form is in the document, not mid-effect
      if (Rows.isDraft(row)) {
        queueMicrotask(() => {
          if (this._form.peek() === form)
            form.focusFirst();
        });
      }
    });
  }

  private _build(row: RowView, access: Access): ObjectForm {
    const source = this.source;
    const options = this._options;
    const table = DomainTable.of(source);
    const draft = Rows.isDraft(row);
    const system = new Set(SYSTEM_COLUMNS.map(([name]) => name));
    const form = propertyForm(source.schema.properties.filter((p) => !system.has(p.name)), row, {
      include: options.include, exclude: options.exclude, layout: options.layout, condensed: options.condensed,
      overrides: options.overrides, access: access.row(row),
      onChanged: options.onChanged,
    });
    form.root.dataset.u2Row = row.id;
    for (const input of form.inputs) {
      const name = input.name;
      if (name === undefined)
        continue;
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
      const caption = SYSTEM_COLUMNS.find(([column]) => column === name)?.[2] ?? name;
      const raw = row[name];
      const value = span('', 'u2-form-readonly-value');
      value.dataset.u2Part = 'readonly-value';
      if (name === 'id' && draft)
        value.textContent = 'assigned on save';
      else if ((prop.propertyType ?? prop.type) === 'datetime' && raw !== null && raw !== undefined && raw !== '')
        value.append(timestamp(raw as Date | number | string));
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
    this._captions.clear();
    for (const prop of [...form.properties, ...extra]) {
      const name = prop.name!;
      const id = row[name];
      if (form.input(name) !== undefined || !DomainTable.isReference(prop) ||
          id === null || id === undefined || id === '')
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
    // a refresh re-reads the uuid into the text row; the caption resolved for it is put back
    this._applyCaptions(row);
  }

  static emptyHint(source: DomainSource): string {
    const name = source.schema.info.singularName.toLowerCase() || 'row';
    return `Select ${/^[aeiou]/.test(name) ? 'an' : 'a'} ${name} to edit.`;
  }
}

export function domainForm(target: DomainFormTarget, options?: DomainFormOptions): DomainForm {
  return new DomainForm(target, options);
}

// the editors a schema's reference columns get everywhere a form is generated: a picker over the
// target table for a `<schema>.<table>` ref, the platform pickers for a user or a group
Editors.register({
  match: (prop) => DomainTable.isReference(prop) && prop.semType !== 'User' && prop.semType !== 'Group',
  create: (prop, options) => domainPick(prop.semType!, options),
});
Editors.register({
  match: (prop) => prop.semType === 'User',
  create: (prop, options) => DomainForm.userPick(options),
});
Editors.register({
  match: (prop) => prop.semType === 'Group',
  create: (prop, options) => DomainForm.groupPick(options),
});
