/* One-line query box over the filter grammar with completion at the caret: the text
   is the user's draft, the tree is written on commit (Enter, blur), and the bound tree re-formats
   the text unless the box is focused and dirty. Rows follow `Filters.completionContext` —
   schema properties, the property's operators, its values (or literal hints per kind) and the
   connectors — through the shared `SuggestionList` under the `u2-fq` prefix. */
import {signal, computed, Signal, ReadonlySignal} from '../../core/signals.js';
import {Input, InputOptions} from '../../core/input-base.js';
import {bindValue} from '../../core/bind.js';
import {div, span} from '../../core/elements.js';
import {SuggestionList} from '../../core/suggestion-list.js';
import type {AsyncState} from '../../core/async-source.js';
import type {IWidgetStatus} from '../../core/widget-like.js';
import {icon} from '../display/icon.js';
import {Filters} from '../../core/filter/index.js';
import type {FilterCompletion, FilterGroup, FilterKind, FilterOperator, FilterProblem, FilterProperty,
  FilterScalar, FilterSchema, FilterTarget, FilterValueItem} from '../../core/filter/index.js';

export interface FilterQueryInputOptions extends InputOptions<FilterGroup> {
  schema: FilterSchema;
  placeholder?: string;
  /** What the tree must be expressible for; validation only. */
  target?: FilterTarget;
}

export interface FilterQueryInputStatus extends IWidgetStatus {
  tree: FilterGroup;
  query: string;
  problems: FilterProblem[];
  mode: 'text';
}

/** One popup row: the grammar text it inserts and how it shows; an `info` row (a number's
 * range) inserts nothing. */
interface Option {
  text: string;
  render: () => HTMLElement;
  info?: boolean;
}

const KIND_ICON: Record<FilterKind, string> = {
  string: 'font', int: 'hashtag', float: 'hashtag', bigint: 'hashtag', datetime: 'calendar',
  bool: 'check-square', string_list: 'list', ref: 'link',
};
const KIND_HINTS: Partial<Record<FilterKind, string[]>> = {
  bool: ['true', 'false'], datetime: ['now', '-1d', '-1w', '-1m'], ref: ['@current'],
};
/** The model-only operators, spelled the way the grammar reads them. */
const SPELLED: Record<string, string> = {'is null': '= null', 'is not null': '!= null'};
/** Their value is a pattern: a picked literal is escaped the way the formatter spells it. */
const LIKE_FAMILY = new Set(['like', '!like', 'starts', 'ends']);

function row(main: HTMLElement | string, secondary?: string, glyph?: HTMLElement): HTMLElement {
  const parts: HTMLElement[] = glyph ? [glyph] : [];
  parts.push(typeof main === 'string' ? span(main, 'u2-fq-text') : main);
  if (secondary)
    parts.push(span(secondary, 'u2-fq-secondary'));
  return div(parts, 'u2-fq-row');
}

export class FilterQueryInput extends Input<FilterGroup, FilterQueryInputOptions> {
  /** The raw text; committed to `value` on Enter or blur. */
  readonly text: Signal<string>;
  /** The same text, for read-only bindings. */
  readonly query: ReadonlySignal<string>;
  readonly problems: ReadonlySignal<FilterProblem[]>;
  readonly isOpen: ReadonlySignal<boolean>;

  // every field below is built by createEditor(), which the base constructor calls before
  // subclass initializers would run — none of them may carry an inline initializer
  private _input!: HTMLInputElement;
  private _list!: SuggestionList<Option>;
  private _text!: Signal<string>;
  private _prefix!: Signal<string>;
  private _view!: Signal<AsyncState<Option>>;
  private _items!: ReadonlySignal<Option[]>;
  private _problems!: Signal<FilterProblem[]>;
  private _context: FilterCompletion | undefined;
  private _abort: AbortController | undefined;
  private _dirty!: boolean;

  constructor(options: FilterQueryInputOptions) {
    super(options, Filters.group('and'));
    this.text = this._text;
    this.query = computed(() => this._text.value);
    this.problems = this._problems;
    this.isOpen = this._list.isOpen;
    this.root.classList.add('u2-filter-query-input');
    this.root.dataset.u2 = 'filter-query-input';
  }

  /** Parses the text into the tree; false — with the first problem as validity — when it
   * cannot be, the old tree kept. A text equal to the tree re-formats without writing a new
   * root, so a bound builder keeps its rows. */
  commit(): boolean {
    const o = this.options;
    const {root, problems} = Filters.parse(this._text.peek(), o.schema, o.target);
    this._problems.value = problems;
    if (problems.length > 0)
      return false;
    this._dirty = false;
    const current = this.value.peek();
    if (Filters.equals(root, current))
      this._text.value = Filters.format(current);
    else
      this.value.value = root;
    return true;
  }

  getWidgetStatus(): FilterQueryInputStatus {
    const status = super.getWidgetStatus() as FilterQueryInputStatus;
    status.tree = this.value.peek();
    status.query = this._text.peek();
    status.problems = this._problems.peek();
    status.mode = 'text';
    return status;
  }

  protected createEditor(): HTMLElement {
    const input = document.createElement('input');
    this._input = input;
    input.type = 'text';
    input.autocomplete = 'off';
    input.spellcheck = false;
    input.placeholder = this.options.placeholder ?? '';
    input.setAttribute('role', 'combobox');
    input.setAttribute('aria-autocomplete', 'list');

    this._dirty = false;
    this._text = signal('');
    this._prefix = signal('');
    this._view = signal<AsyncState<Option>>({kind: 'idle'});
    this._problems = signal<FilterProblem[]>([]);
    this._items = computed(() => {
      const state = this._view.value;
      return state.kind === 'ready' ? state.items : [];
    });
    this._list = new SuggestionList<Option>({
      prefix: 'u2-fq',
      anchor: input,
      scope: this.scope,
      items: this._items,
      view: this._view,
      text: this._prefix,
      minChars: 0,
      render: (item) => item.render(),
      autoHighlight: true,
      onPick: (index) => this._insert(index),
      onDismiss: () => this._list.dismiss(),
      onRetry: () => this._refresh(),
    });
    this.own(() => this._abort?.abort());

    bindValue(this.scope, input, this._text);
    this._listen(input, 'input', () => {
      this._dirty = true;
      this._suggest();
    });
    this._listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    // the caret moved without typing: the context under it is another token
    const refreshIfOpen = () => {
      if (this._list.isOpen.peek())
        this._refresh();
    };
    this._listen(input, 'keyup', (e) => {
      const key = (e as KeyboardEvent).key;
      if (key === 'ArrowLeft' || key === 'ArrowRight' || key === 'Home' || key === 'End')
        refreshIfOpen();
    });
    // a click where a value is expected offers the values unasked; elsewhere it only re-aims an open list
    this._listen(input, 'click', () => {
      if (this._list.isOpen.peek())
        this._refresh();
      else if (Filters.completionContext(this._text.peek(), this._caret(), this.options.schema).expect === 'value')
        this._suggest();
    });
    this._listen(input, 'select', refreshIfOpen);
    this._listen(input, 'focus', () => {
      if (this._text.peek() === Filters.format(this.value.peek()))
        this._dirty = false;
    });
    this._listen(input, 'blur', () => {
      this._list.dismiss();
      if (this._dirty)
        this.commit();
    });
    this.effect(() => input.setAttribute('aria-expanded', String(this._list.isOpen.value)));
    this.addValidator(() => this._problems.value[0]?.message ?? null);
    this.effect(() => this._syncText());
    return input;
  }

  /** The bound tree re-formats the text — except under the user's feet. */
  private _syncText(): void {
    const root = this.value.value;
    if (this._dirty && document.activeElement === this._input)
      return;
    this._dirty = false;
    this._text.value = Filters.format(root);
    this._problems.value = Filters.validate(root, this.options.schema, this.options.target);
  }

  private _caret(): number {
    return this._input.selectionStart ?? this._text.peek().length;
  }

  /** Fresh rows for the caret, and the list open — unless a value expectation has nothing to
   * offer beyond `null`, which is not worth a popup. */
  private _suggest(): void {
    this._list.clearActive();
    if (this._refresh())
      this._list.open();
    else
      this._list.dismiss();
  }

  /** False when the rows are a lone `null` (the async values decide later, through `_ready`). */
  private _refresh(): boolean {
    const ctx = Filters.completionContext(this._text.peek(), this._caret(), this.options.schema);
    this._context = ctx;
    this._prefix.value = ctx.prefix;
    this._abort?.abort();
    this._abort = undefined;
    const q = ctx.prefix.toLowerCase();
    const starts = (s: string) => s.toLowerCase().startsWith(q);
    switch (ctx.expect) {
      case 'property':
        return this._ready(this.options.schema.properties
          .filter((p) => starts(p.name) || starts(p.friendlyName ?? ''))
          .map((p) => this._propertyOption(p)));
      case 'operator':
        return this._ready(this._operators(ctx.property).map((o) => ({spelling: SPELLED[o.id] ?? o.id, op: o}))
          .filter((x) => starts(x.spelling))
          .map((x) => ({text: x.spelling, render: () => row(span(x.spelling, 'u2-fq-code'), x.op.label)})));
      case 'connector':
        return this._ready(['and', 'or'].filter(starts)
          .map((text) => ({text, render: () => row(span(text, 'u2-fq-code'))})));
      case 'value':
        return this._values(ctx, starts);
    }
    return true;
  }

  private _ready(items: Option[]): boolean {
    this._view.value = items.length > 0 ? {kind: 'ready', items} : {kind: 'empty'};
    return !(items.length === 1 && items[0].text === 'null');
  }

  private _operators(prop: FilterProperty | undefined): FilterOperator[] {
    return prop ? Filters.operators.for(prop) : Filters.operators.all().filter((o) => o.semType === undefined);
  }

  private _propertyOption(p: FilterProperty): Option {
    const label = p.friendlyName ?? p.name;
    return {text: Filters.formatProperty(p.name), render: () => row(label, label === p.name ? undefined : p.name,
      icon(KIND_ICON[Filters.kindOf(p)], {cls: 'u2-fq-icon'}))};
  }

  /** The schema's values for the property (async, abortable) followed by the kind's literal
   * hints — a number's range first; the hints alone where the schema offers no values. */
  private _values(ctx: FilterCompletion, starts: (s: string) => boolean): boolean {
    const prop = ctx.property;
    const hints = prop ? [...(KIND_HINTS[Filters.kindOf(prop)] ?? []), 'null'] : ['null'];
    const hintOptions: Option[] = hints.filter(starts)
      .map((text) => ({text, render: () => row(span(text, 'u2-fq-code'))}));
    const range = prop ? FilterQueryInput._rangeHint(prop) : null;
    if (range !== null)
      hintOptions.unshift({text: '', info: true, render: () => row(span(range, 'u2-fq-info'))});
    const values = this.options.schema.values;
    if (!prop || !values)
      return this._ready(hintOptions);
    const abort = new AbortController();
    this._abort = abort;
    this._view.value = {kind: 'loading'};
    const renderer = this.options.schema.renderer?.(prop);
    const render = (item: FilterValueItem): HTMLElement => {
      const v = item.value;
      const custom = renderer?.listItem?.(v) ?? renderer?.markup?.(v);
      return custom ? row(custom) : row(item.label ?? (renderer ? renderer.caption(v) : String(v)),
        item.count === undefined ? undefined : String(item.count));
    };
    const failed = (e: unknown) => {
      if (!abort.signal.aborted)
        this._view.value = {kind: 'error', message: e instanceof Error ? e.message : String(e)};
    };
    values(prop, ctx.prefix, abort.signal).then((items) => {
      if (abort.signal.aborted)
        return;
      const offered = [...items.map((item) => ({text: FilterQueryInput._valueText(ctx, item.value),
        render: () => render(item)})), ...hintOptions];
      if (!this._ready(offered))
        this._list.dismiss();
    }, failed);
    return true;
  }

  private static _valueText(ctx: FilterCompletion, v: FilterScalar): string {
    const pattern = typeof v === 'string' && ctx.operator !== undefined && LIKE_FAMILY.has(ctx.operator.id);
    return Filters.formatValue(pattern ? Filters.escapeLike(v) : v);
  }

  /** `number 0–120` for a numeric property with bounds. */
  private static _rangeHint(prop: FilterProperty): string | null {
    const kind = Filters.kindOf(prop);
    if ((kind !== 'int' && kind !== 'float') || (prop.min === undefined && prop.max === undefined))
      return null;
    const range = prop.min !== undefined && prop.max !== undefined ? `${prop.min}–${prop.max}` :
      prop.min !== undefined ? `≥ ${prop.min}` : `≤ ${prop.max}`;
    return `number ${range}`;
  }

  /** Replaces the token under the caret with the row's text and one space (the one already there
   * when the pick lands mid-text), then re-opens for the next expectation. */
  private _insert(index: number): boolean {
    const item = this._items.peek()[index];
    const ctx = this._context;
    if (item === undefined || item.info || !ctx)
      return false;
    const text = this._text.peek();
    const inserted = /^\s/.test(text.slice(ctx.replace.end)) ? item.text : `${item.text} `;
    const caret = ctx.replace.start + item.text.length + 1;
    this._dirty = true;
    this._text.value = text.slice(0, ctx.replace.start) + inserted + text.slice(ctx.replace.end);
    this._input.focus();
    this._input.setSelectionRange?.(caret, caret);
    this._suggest();
    return true;
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const open = this._list.isOpen.peek();
    switch (e.key) {
      case 'ArrowDown':
        e.preventDefault();
        if (open)
          this._list.move(1);
        else
          this._suggest();
        break;
      case 'ArrowUp':
        if (open) {
          e.preventDefault();
          this._list.move(-1);
        }
        break;
      case 'Enter':
        if (e.ctrlKey || e.metaKey)
          break;
        e.preventDefault();
        if (open && this._insert(this._list.activeIndex.peek()))
          e.stopPropagation();
        else {
          this._list.dismiss();
          this.commit();
        }
        break;
      case ' ':
        if (this._text.peek() === '') {
          e.preventDefault();
          if (!open)
            this._suggest();
        }
        break;
      case 'Escape':
        if (open) {
          e.stopPropagation();
          this._list.dismiss();
        }
        break;
      case 'Tab':
        this._list.dismiss();
        break;
    }
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
