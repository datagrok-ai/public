/* The rows of a `FilterBuilder`: condition, advanced-mode group and the read-only nested summary of simple mode.
   A row never touches the tree: edits go to the host as a patch, the host writes a new root and calls `update()`. */
import {Scope} from '../../core/scope.js';
import {Control} from '../../core/component.js';
import {div, span, button} from '../../core/elements.js';
import type {Input} from '../../core/input-base.js';
import {markSpan, resolveSpan, spanOf} from '../../core/span.js';
import {iconButton} from '../actions/buttons.js';
import {ButtonGroup} from '../actions/button-group.js';
import {icon} from '../display/icon.js';
import {ChoiceInput} from '../inputs/choice-input.js';
import {TextInput} from '../inputs/text-input.js';
import {NumberInput} from '../inputs/number-input.js';
import {BigIntInput} from '../inputs/bigint-input.js';
import {DateTimeInput} from '../inputs/date-input.js';
import {SuggestInput} from '../inputs/suggest-input.js';
import {TagsInput} from '../inputs/tags-input.js';
import {RefInput} from './ref-input.js';
import {KIND, format, isRef, isSpan, kindOf, operators, property, valueEquals} from '../../core/filter/index.js';
import type {FilterCondition, FilterGroup, FilterKind, FilterScalar, FilterValue, Lock, FilterOperator,
  FilterProperty, FilterSchema, FilterValueEditorFactory} from '../../core/filter/index.js';

/** What a row asks of the builder that owns it. */
export interface FilterRowHost {
  schema: FilterSchema;
  editors: FilterValueEditorFactory | null;
  /** The properties on offer (a template may narrow them). */
  properties(): FilterProperty[];
  operators(prop: FilterProperty): FilterOperator[];
  /** The strongest lock on the node, its ancestors and their template counterparts. */
  lockOf(id: string): Lock;
  /** Advanced mode: handles, `not`, `+ group`. */
  advanced(): boolean;
  allowAdd: boolean;
  change(id: string, patch: Partial<FilterCondition> | Partial<FilterGroup>): void;
  add(parentId: string, kind: 'condition' | 'group'): void;
  remove(id: string): void;
}

/** A mounted value editor: its element, a write of the condition's value (or a bound of it), the enabled switch. */
interface ValueEditor {
  root: HTMLElement;
  set(v: FilterValue | undefined): void;
  setEnabled(x: boolean): void;
}

function scalarText(v: FilterScalar): string {
  if (isRef(v))
    return v.name ?? v.id;
  if (isSpan(v))
    return v.span;
  return v instanceof Date ? v.toLocaleDateString() : String(v);
}

function parseScalar(text: string, kind: FilterKind): FilterScalar {
  const numeric = kind === KIND.INT || kind === KIND.FLOAT;
  return numeric && text.trim() !== '' && Number.isFinite(Number(text)) ? Number(text) : text;
}

/** An editor of one scalar: a list (another operator's shape) reads as empty. */
function scalarEditor(input: Input<any>, set: (v: FilterScalar | undefined) => void): ValueEditor {
  return {root: input.root, set: (v) => set(Array.isArray(v) ? undefined : v), setEnabled: (x) => input.enabled = x};
}

/** An editor whose input holds model values as they are; a write of what it already holds is skipped. */
function valueEditor(input: Input<any>, toInput: (v: FilterValue | undefined) => FilterValue | undefined): ValueEditor {
  return {root: input.root, setEnabled: (x) => input.enabled = x, set: (v) => {
    const next = toInput(v);
    if (!valueEquals(input.value.peek(), next))
      input.value.value = next;
  }};
}

function removeButton(onRemove: () => void, what: string = 'condition'): HTMLButtonElement {
  const el = iconButton('minus', onRemove, {tooltip: `Remove ${what}`});
  el.classList.add('u2-fb-remove');
  el.dataset.u2Part = 'remove';
  return el;
}

/** The drag grip; the builder's drag layer hears the press through delegation. */
function dragHandle(): HTMLElement {
  const el = div([icon('grip-vertical')], 'u2-fb-handle');
  el.dataset.u2Part = 'handle';
  el.title = 'Drag to move';
  el.setAttribute('role', 'button');
  el.setAttribute('aria-label', 'Drag to move');
  return el;
}

/** A hidden text part: the problem line under the controls, or what a locked picker shows instead of a select. */
function textPart(cls: string, part: string): HTMLElement {
  const el = span('', cls);
  el.dataset.u2Part = part;
  el.hidden = true;
  return el;
}

function markProblem(root: HTMLElement, line: HTMLElement, message: string | null): void {
  root.classList.toggle('u2-fb-invalid', message !== null);
  line.textContent = message ?? '';
  line.hidden = message === null;
  if (message === null)
    root.removeAttribute('title');
  else
    root.title = message;
}

export class FilterRow extends Control {
  readonly id: string;

  private readonly _host: FilterRowHost;
  private _cond: FilterCondition;
  private _handle!: HTMLElement;
  private _prop!: ChoiceInput;
  private _op!: ChoiceInput;
  private _propText!: HTMLElement;
  private _opText!: HTMLElement;
  private _remove!: HTMLButtonElement;
  private _valueHost!: HTMLElement;
  private _problem!: HTMLElement;
  private _editorKey = '';
  private _editorScope: Scope | undefined;
  private _editors: ValueEditor[] = [];

  constructor(cond: FilterCondition, host: FilterRowHost) {
    super();
    this.id = cond.id;
    this._host = host;
    this._cond = cond;
    this.root.classList.add('u2-fb-cond');
    this.root.dataset.u2 = 'filter-row';
    this.root.dataset.u2Node = cond.id;
    this.runInScope(() => this._build(cond));
    this.update(cond);
  }

  /** The node as the tree now holds it (same id). */
  update(cond: FilterCondition): void {
    this._cond = cond;
    this._prop.setItems(FilterRow._items(this._host.properties().map((p) =>
      ({value: p.name, label: p.friendlyName ?? p.name})), cond.property));
    this._prop.value.value = cond.property;
    const prop = property(this._host.schema, cond.property);
    const ops = prop ? this._host.operators(prop) : [];
    this._op.setItems(FilterRow._items(ops.map((o) => ({value: o.id, label: o.label})), cond.operator));
    this._op.value.value = cond.operator;
    const op = ops.find((o) => o.id === cond.operator);
    this._propText.textContent = prop?.friendlyName ?? cond.property;
    this._opText.textContent = op?.label ?? cond.operator;
    // editors are rebuilt on a change of shape, not of operator: `>` → `<=` keeps the number box
    const key = `${cond.property} ${op?.editor}/${op?.arity}`;
    if (key !== this._editorKey) {
      this._editorKey = key;
      this._buildEditors(cond, prop, op);
    } else
      this._syncEditors(cond);
    this._applyLock();
  }

  setProblem(message: string | null): void {
    markProblem(this.root, this._problem, message);
  }

  focus(): void {
    (this._prop.root.querySelector('select') as HTMLElement | null)?.focus();
  }

  private _build(cond: FilterCondition): void {
    this._handle = dragHandle();
    this._prop = new ChoiceInput({inline: true, nullable: false, items: [], value: cond.property,
      onChanged: (v) => this._onProperty(v)});
    this._prop.root.classList.add('u2-fb-prop');
    this._prop.root.dataset.u2Part = 'prop';
    this._op = new ChoiceInput({inline: true, nullable: false, items: [], value: cond.operator,
      onChanged: (v) => this._onOperator(v)});
    this._op.root.classList.add('u2-fb-op');
    this._op.root.dataset.u2Part = 'op';
    this._propText = textPart('u2-fb-prop-text', 'prop-text');
    this._opText = textPart('u2-fb-op-text', 'op-text');
    this._valueHost = div([], 'u2-fb-values');
    this._remove = removeButton(() => this._host.remove(this.id));
    this._problem = textPart('u2-fb-problem', 'problem');
    this.root.append(this._handle, this._prop.root, this._propText, this._op.root, this._opText, this._valueHost,
      this._remove, this._problem);
    this.own(() => this._releaseEditors());
  }

  /** The picker's items; a current value the offer lacks (an unknown name, a narrowed set) stays selectable. */
  private static _items(items: {value: string, label: string}[], current: string): {value: string, label: string}[] {
    return items.some((i) => i.value === current) ? items : [...items, {value: current, label: current}];
  }

  // onChanged fires after the batch update() wrote in, so callbacks compare with the condition, not a flag
  private _onProperty(name: string | null): void {
    if (name === null || name === this._cond.property)
      return;
    const prop = property(this._host.schema, name);
    const ops = prop ? this._host.operators(prop) : [];
    const operator = ops.some((o) => o.id === this._cond.operator) ? this._cond.operator :
      ops[0]?.id ?? this._cond.operator;
    this._host.change(this.id, {property: name, operator, value: undefined});
  }

  // the value survives an operator swap only while the editor stays the same shape (= → !=)
  private _onOperator(id: string | null): void {
    if (id === null || id === this._cond.operator)
      return;
    const [before, after] = [this._cond.operator, id].map((o) => operators.get(o));
    const keep = after !== undefined && before?.arity === after.arity && before?.editor === after.editor;
    this._host.change(this.id, keep ? {operator: id} : {operator: id, value: undefined});
  }

  private _change(value: FilterValue | undefined): void {
    if (!valueEquals(value, this._cond.value))
      this._host.change(this.id, {value});
  }

  // a locked property/operator reads as text: a disabled select says "broken", not "fixed"
  private _applyLock(): void {
    const lock = this._host.lockOf(this.id);
    const locked = lock !== 'none';
    this.root.classList.toggle('u2-fb-locked', locked);
    this._remove.hidden = locked;
    this._handle.hidden = locked || !this._host.advanced();
    this._prop.root.hidden = locked;
    this._op.root.hidden = locked;
    this._propText.hidden = !locked;
    this._opText.hidden = !locked;
    for (const editor of this._editors)
      editor.setEnabled(lock !== 'all');
  }

  private _releaseEditors(): void {
    this._editorScope?.dispose();
    this._editorScope = undefined;
    this._editors = [];
    this._valueHost.textContent = '';
  }

  private _buildEditors(cond: FilterCondition, prop: FilterProperty | null, op: FilterOperator | undefined): void {
    this._releaseEditors();
    if (!prop || !op || op.editor === 'none' || op.arity === 0)
      return;
    const scope = new Scope();
    this._editorScope = scope;
    Scope.runWith(scope, () => {
      if (op.editor === 'list')
        this._mount(this._listEditor(cond, prop), 'value');
      else if (op.editor === 'range') {
        const [lo, hi] = FilterRow._bounds(cond);
        this._mount(this._scalarEditor(prop, lo, (v) => this._writeBound(0, v)), 'value');
        this._valueHost.append(span('and', 'u2-fb-and'));
        this._mount(this._scalarEditor(prop, hi, (v) => this._writeBound(1, v)), 'value2');
      } else
        this._mount(this._singleEditor(cond, prop), 'value');
    });
  }

  private _mount(editor: ValueEditor, part: 'value' | 'value2'): void {
    editor.root.classList.add(`u2-fb-${part}`);
    editor.root.dataset.u2Part = part;
    this._valueHost.append(editor.root);
    this._editors.push(editor);
  }

  private _syncEditors(cond: FilterCondition): void {
    const values = this._editors.length === 2 ? FilterRow._bounds(cond) : [cond.value];
    for (let i = 0; i < this._editors.length; i++)
      this._editors[i].set(values[i]);
  }

  private static _bounds(c: FilterCondition): [FilterScalar | undefined, FilterScalar | undefined] {
    return Array.isArray(c.value) ? [c.value[0], c.value[1]] : [undefined, undefined];
  }

  private _writeBound(at: 0 | 1, v: FilterScalar | undefined): void {
    const pair = FilterRow._bounds(this._cond);
    pair[at] = v;
    this._change(pair.every((b) => b === undefined) ? undefined : [pair[0] ?? null, pair[1] ?? null]);
  }

  private _singleEditor(cond: FilterCondition, prop: FilterProperty): ValueEditor {
    const custom = this._host.editors?.(prop, {inline: true, value: cond.value, onChanged: (v) => this._change(v)});
    return custom ? valueEditor(custom, (v) => v) :
      this._scalarEditor(prop, Array.isArray(cond.value) ? undefined : cond.value, (v) => this._change(v));
  }

  private _listEditor(cond: FilterCondition, prop: FilterProperty): ValueEditor {
    const kind = kindOf(prop);
    const values = this._host.schema.values;
    const list = (v: FilterValue | undefined): FilterScalar[] => Array.isArray(v) ? v : [];
    const tags = new TagsInput<FilterScalar>({
      inline: true,
      value: list(cond.value),
      itemText: scalarText,
      allowNew: true,
      createItem: (text) => parseScalar(text, kind),
      placeholder: 'Type a value…',
      ...(values ?
        {source: async (q: string, signal: AbortSignal) => (await values(prop, q, signal)).map((i) =>
          isRef(i.value) && i.label ? {...i.value, name: i.label} : i.value)} :
        {items: prop.choices ?? []}),
      onChanged: (v) => this._change(v.length === 0 ? undefined : v),
    });
    return valueEditor(tags, list);
  }

  /** The core kind defaults; `onChange` gets `undefined` for a cleared editor. */
  private _scalarEditor(prop: FilterProperty, initial: FilterScalar | undefined,
    onChange: (v: FilterScalar | undefined) => void): ValueEditor {
    const kind = kindOf(prop);
    if (prop.choices) {
      const text = (v: FilterScalar | undefined) => v === undefined || v === null ? null : String(v);
      const input = new ChoiceInput({inline: true, items: prop.choices, value: text(initial),
        onChanged: (v) => onChange(v === null ? undefined : parseScalar(v, kind))});
      return scalarEditor(input, (v) => input.value.value = text(v));
    }
    const values = this._host.schema.values;
    if (kind === KIND.REF && values) {
      const input = new RefInput({inline: true, prop, schema: this._host.schema,
        value: isRef(initial) ? initial : null, onChanged: (v) => onChange(v ?? undefined)});
      return scalarEditor(input, (v) => input.value.value = isRef(v) ? v : null);
    }
    switch (kind) {
      case KIND.INT: case KIND.FLOAT: {
        const input = new NumberInput({inline: true, mode: kind, min: prop.min, max: prop.max,
          value: typeof initial === 'number' ? initial : null,
          onChanged: (v) => onChange(v === null ? undefined : v)});
        return scalarEditor(input, (v) => input.value.value = typeof v === 'number' ? v : null);
      }
      case KIND.BIG_INT: {
        const parse = (v: FilterScalar | undefined): bigint | null =>
          typeof v === 'number' || typeof v === 'string' && /^-?\d+$/.test(v) ? BigInt(v) : null;
        const input = new BigIntInput({inline: true, value: parse(initial),
          onChanged: (v) => onChange(v === null ? undefined : v.toString())});
        return scalarEditor(input, (v) => input.value.value = parse(v));
      }
      case KIND.BOOL: {
        const text = (v: FilterScalar | undefined) => v === true ? 'true' : v === false ? 'false' : null;
        const input = new ChoiceInput({inline: true, items: ['true', 'false'], value: text(initial),
          onChanged: (v) => onChange(v === null ? undefined : v === 'true')});
        return scalarEditor(input, (v) => input.value.value = text(v));
      }
      case KIND.DATE_TIME: {
        // a span is a `Date` tagged with its text; two resolutions of one span are the same value
        const toDate = (v: FilterScalar | undefined): Date | null =>
          isSpan(v) ? markSpan(resolveSpan(v.span, new Date()), v.span) : v instanceof Date ? v : null;
        const same = (a: Date | null, b: Date | null) =>
          a === b || a !== null && b !== null && (spanOf(a) ?? a.getTime()) === (spanOf(b) ?? b.getTime());
        const input = new DateTimeInput({inline: true, relative: true, value: toDate(initial),
          onChanged: (v) => onChange(v === null ? undefined : spanOf(v) === undefined ? v : {span: spanOf(v)!})});
        return scalarEditor(input, (v) => {
          if (!same(input.value.peek(), toDate(v)))
            input.value.value = toDate(v);
        });
      }
      default: {
        const text = (v: FilterScalar | undefined) => v === undefined || v === null ? '' : scalarText(v);
        const options = {inline: true, value: text(initial),
          onChanged: (v: string) => onChange(v === '' ? undefined : v)};
        const input = values && kind === KIND.STRING ?
          new SuggestInput({...options, openOnFocus: true, source: async (q, signal) =>
            (await values(prop, q, signal)).map((i) => scalarText(i.value))}) :
          new TextInput(options);
        return scalarEditor(input, (v) => input.value.value = text(v));
      }
    }
  }
}

/** A sub-group shown inside simple mode: a read-only summary of what advanced mode edits. */
export class FilterNestedRow extends Control {
  readonly id: string;

  private readonly _host: FilterRowHost;
  private readonly _text: HTMLElement;
  private readonly _remove: HTMLButtonElement;
  private readonly _problem: HTMLElement;

  constructor(group: FilterGroup, host: FilterRowHost) {
    super();
    this.id = group.id;
    this._host = host;
    this.root.classList.add('u2-fb-cond', 'u2-fb-nested');
    this.root.dataset.u2 = 'filter-row';
    this.root.dataset.u2Node = group.id;
    this._text = span('', 'u2-fb-nested-text');
    this._text.dataset.u2Part = 'value';
    this._remove = this.runInScope(() => removeButton(() => this._host.remove(this.id)));
    this._problem = textPart('u2-fb-problem', 'problem');
    this.root.append(this._text, this._remove, this._problem);
    this.update(group);
  }

  update(group: FilterGroup): void {
    this._text.textContent = `(${format(group) || '…'})`;
    const lock = this._host.lockOf(this.id);
    this.root.classList.toggle('u2-fb-locked', lock !== 'none');
    this._remove.hidden = lock !== 'none';
  }

  setProblem(message: string | null): void {
    markProblem(this.root, this._problem, message);
  }
}

/** A group in advanced mode — and, as `root`, the builder's own header and row host: the and/or
 * toggle, `not`, `+ condition`, `+ group`, `−` (not for the root) and `rows`, which the builder
 * fills with the children it reconciles. */
export class FilterGroupRow extends Control {
  readonly header: HTMLElement;
  readonly rows: HTMLElement;
  readonly connector: ButtonGroup;
  readonly not: HTMLButtonElement;
  readonly add: HTMLButtonElement;
  readonly addGroup: HTMLButtonElement;
  readonly remove: HTMLButtonElement | undefined;
  readonly handle: HTMLElement | undefined;

  private readonly _host: FilterRowHost;
  private _group: FilterGroup;

  /** Follows `update()`: the root row is told whichever root the value holds next. */
  get id(): string {
    return this._group.id;
  }

  constructor(group: FilterGroup, host: FilterRowHost, root: boolean = false) {
    super();
    this._host = host;
    this._group = group;
    this.root.classList.add(root ? 'u2-fb-root' : 'u2-fb-group');
    if (!root) {
      this.root.dataset.u2 = 'filter-group';
      this.root.dataset.u2Node = group.id;
    }
    this.connector = new ButtonGroup({toggle: 'single', density: 'toolbar',
      items: [{id: 'and', label: 'and'}, {id: 'or', label: 'or'}]});
    this.connector.root.classList.add('u2-fb-connector');
    this.connector.root.dataset.u2Part = 'connector';
    this.not = button('not', () => this._host.change(this.id, {not: this._group.not === true ? undefined : true}));
    this.not.classList.add('u2-fb-not');
    this.not.dataset.u2Part = 'not';
    this.not.title = 'Negate the group';
    this.add = this.runInScope(() => iconButton('plus', () => this._host.add(this.id, 'condition'),
      {tooltip: 'Add condition'}));
    this.add.classList.add('u2-fb-add');
    this.add.dataset.u2Part = 'add';
    this.addGroup = this.runInScope(() => iconButton('folder-plus', () => this._host.add(this.id, 'group'),
      {tooltip: 'Add group'}));
    this.addGroup.classList.add('u2-fb-add-group');
    this.addGroup.dataset.u2Part = 'add-group';
    const parts: HTMLElement[] = [this.connector.root, this.not, this.add, this.addGroup];
    if (!root) {
      this.handle = dragHandle();
      this.remove = this.runInScope(() => removeButton(() => this._host.remove(this.id), 'group'));
      parts.unshift(this.handle);
      parts.push(this.remove);
    }
    this.header = div(parts, root ? 'u2-fb-header' : 'u2-fb-group-header');
    this.rows = div([], root ? 'u2-fb-rows' : 'u2-fb-group-rows');
    this.rows.dataset.u2Part = 'rows';
    this.root.append(this.header, this.rows);
    this.effect(() => {
      const op = this.connector.selected.value[0] as 'and' | 'or' | undefined;
      if (op !== undefined && op !== this._group.op)
        this._host.change(this.id, {op});
    });
    this.update(group);
  }

  update(group: FilterGroup): void {
    this._group = group;
    if (this.connector.selected.peek()[0] !== group.op)
      this.connector.selected.value = [group.op];
    const advanced = this._host.advanced();
    const lock = this._host.lockOf(this.id);
    const locked = lock !== 'none';
    this.root.classList.toggle('u2-fb-locked', locked);
    this.root.classList.toggle('u2-fb-negated', group.not === true);
    this.not.classList.toggle('u2-fb-not-on', group.not === true);
    this.not.setAttribute('aria-pressed', String(group.not === true));
    this.not.hidden = !advanced;
    this.not.disabled = locked;
    this.connector.setEnabled('and', !locked);
    this.connector.setEnabled('or', !locked);
    this.add.hidden = locked || !this._host.allowAdd;
    this.add.disabled = this._host.properties().length === 0;
    this.addGroup.hidden = !advanced || locked || !this._host.allowAdd;
    if (this.remove)
      this.remove.hidden = locked;
    if (this.handle)
      this.handle.hidden = locked;
  }

  /** A group's problem (nesting, a lock) is the header's title: the rows below carry their own. */
  setProblem(message: string | null): void {
    this.root.classList.toggle('u2-fb-invalid', message !== null);
    if (message === null)
      this.header.removeAttribute('title');
    else
      this.header.title = message;
  }
}
