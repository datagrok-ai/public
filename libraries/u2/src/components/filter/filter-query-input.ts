/* One-line query box over the filter grammar with completion at the caret: the text
   is the user's draft, the tree is written on commit (Enter, blur), and the bound tree re-formats
   the text unless the box is focused and dirty. Rows follow `Filters.completionContext` —
   schema properties, the property's operators, its values (or literal hints per kind) and the
   connectors — through the shared `SuggestionList` under the `u2-fq` prefix. */
import {signal, computed, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {Input, InputOptions} from '../../core/input-base.js';
import {bindValue} from '../../core/bind.js';
import {div, span} from '../../core/elements.js';
import {SuggestionList} from '../../core/suggestion-list.js';
import type {AsyncState} from '../../core/async-source.js';
import type {IWidgetStatus} from '../../core/widget-like.js';
import {icon} from '../display/icon.js';
import {Filters, KIND} from '../../core/filter/index.js';
import type {FilterCompletion, FilterGroup, FilterKind, FilterOperator, FilterProblem, FilterProperty,
  FilterScalar, FilterSchema, FilterTarget, FilterValueItem} from '../../core/filter/index.js';

export interface FilterQueryInputOptions extends InputOptions<FilterGroup> {
  schema: FilterSchema;
  placeholder?: string;
  /** What the tree must be expressible for; validation only. */
  target?: FilterTarget;
  /** How the bound tree is spelled in the box — `Filters.format` when it answers null. A query a
   * control wrote bound (a preset's `$me` resolved to an id) shows the way it was written. */
  display?: (tree: FilterGroup) => string | null;
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
  [KIND.STRING]: 'font', [KIND.INT]: 'hashtag', [KIND.FLOAT]: 'hashtag', [KIND.BIG_INT]: 'hashtag',
  [KIND.DATE_TIME]: 'calendar', [KIND.BOOL]: 'check-square', [KIND.STRING_LIST]: 'list', [KIND.REF]: 'link',
};
const KIND_HINTS: Partial<Record<FilterKind, string[]>> = {
  [KIND.BOOL]: ['true', 'false'], [KIND.DATE_TIME]: ['now', '-1d', '-1w', '-1m'], [KIND.REF]: ['@current'],
};
/** The model-only operators, spelled the way the grammar reads them. */
const SPELLED: Record<string, string> = {'is null': '= null', 'is not null': '!= null'};
/** Their value is a pattern: a picked literal is escaped the way the formatter spells it. */
const LIKE_FAMILY = new Set(['like', '!like', 'starts', 'ends']);
/** `<column> under <bare word…>` at the end of the text — the shape the grammar refuses. */
const UNDER_BARE = /([A-Za-z_]\w*)\s+under\s+[^"'\s]\S*(\s+\S+)*\s*$/i;

/** A value slot whose candidates are NAMES — several words, so the whole slot is the search text
 * and the whole slot is what a pick replaces. */
function namedSlot(ctx: FilterCompletion): boolean {
  return ctx.valueSpan !== undefined && ctx.property !== undefined &&
    (ctx.operator?.id === 'under' || Filters.kindOf(ctx.property) === KIND.REF);
}

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
  /** The text {@link showText} put in the box, while it stands; null once a commit took it. */
  readonly raw: ReadonlySignal<string | null>;
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
  private _raw!: Signal<string | null>;
  private _context: FilterCompletion | undefined;
  private _abort: AbortController | undefined;
  private _dirty!: boolean;
  /** Whether the highlight is the user's — an arrow key or a click, not the auto-highlight. */
  private _moved = false;

  constructor(options: FilterQueryInputOptions) {
    super(options, Filters.group('and'));
    this.text = this._text;
    this.query = computed(() => this._text.value);
    this.problems = this._problems;
    this.raw = this._raw;
    this.isOpen = this._list.isOpen;
    this.root.classList.add('u2-filter-query-input');
    this.root.dataset.u2 = 'filter-query-input';
  }

  /** Parses the text into the tree; false — with the first problem as validity — when it
   * cannot be, the old tree kept. A text equal to the tree re-formats without writing a new
   * root, so a bound builder keeps its rows. */
  commit(): boolean {
    const o = this.options;
    // the text is still exactly what the tree in force reads as: nothing to parse, and a displayed
    // form that is not the formatter's (a preset's `$me`) is not turned into a query of its own
    if (this._raw.peek() === null && this._text.peek() === this._display(this.value.peek())) {
      // the text IS the tree in force again (the box was cleared back to it): whatever the last
      // refused commit reported is no longer about what is in the box, and holding on to it left
      // the box red and the collection marked stale over rows that were right all along
      this._problems.value = [];
      this._dirty = false;
      return true;
    }
    const text = this._text.peek();
    const {root, problems} = Filters.parse(text, o.schema, o.target);
    this._problems.value = FilterQueryInput.worded(text, o.schema, problems);
    if (problems.length > 0)
      return false;
    this._dirty = false;
    this._raw.value = null;
    const current = this.value.peek();
    if (Filters.equals(root, current))
      this._text.value = this._display(current);
    else
      this.value.value = root;
    return true;
  }

  /** `location_id under Building A` fails in the GRAMMAR (an unquoted value), long before the
   * schema check that knows what `under` takes — and "Expected a value" says nothing about the
   * one thing to do about it. The parse problems keep their places; the first one is re-worded. */
  static worded(text: string, schema: FilterSchema, problems: FilterProblem[]): FilterProblem[] {
    const match = problems.length === 0 ? null : UNDER_BARE.exec(text);
    if (match === null)
      return problems;
    const prop = Filters.property(schema, match[1]);
    const target = (prop?.ref ?? '').split('.').pop() || 'row';
    return [{...problems[0], message: `under takes a ${target} — pick one from the list`},
      ...problems.slice(1)];
  }

  /** Text the box did not produce — a `?q=` a link carried that the grammar or the schema refuses
   * — shown with its problems until a commit takes it, so the user fixes it in place instead of
   * watching it vanish; null puts the tree in force back. */
  showText(text: string | null, problems: FilterProblem[] = []): void {
    this._raw.value = text;
    if (text === null) {
      this._dirty = false;
      // the caller is an effect of its own: re-reading the tree here must not subscribe it
      untracked(() => this._syncText());
      return;
    }
    this._dirty = true;
    this._text.value = text;
    this._problems.value = problems;
    this.revalidate();
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
    this._raw = signal<string | null>(null);
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
      onPick: (index) => {
        this._moved = true;
        return this._insert(index);
      },
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
      if (this._text.peek() === this._display(this.value.peek()))
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

  /** The bound tree re-formats the text — except under the user's feet, or over text {@link
   * showText} put there that no commit has taken yet. */
  private _syncText(): void {
    const root = this.value.value;
    if (this._raw.peek() !== null || (this._dirty && document.activeElement === this._input))
      return;
    this._dirty = false;
    this._text.value = this._display(root);
    this._problems.value = Filters.validate(root, this.options.schema, this.options.target);
  }

  private _display(tree: FilterGroup): string {
    return this.options.display?.(tree) ?? Filters.format(tree);
  }

  private _caret(): number {
    return this._input.selectionStart ?? this._text.peek().length;
  }

  /** Fresh rows for the caret, and the list open — unless a value expectation has nothing to
   * offer beyond `null`, which is not worth a popup. */
  private _suggest(): void {
    this._moved = false;
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
    // `Building A` is one value: the token under the caret would search for `A` and offer every
    // row, with the wrong one highlighted
    if (namedSlot(ctx)) {
      ctx.replace = ctx.valueSpan!;
      ctx.prefix = this._text.peek().slice(ctx.valueSpan!.start, ctx.valueSpan!.end);
    }
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
    if (!prop)
      return Filters.operators.all().filter((o) => o.semType === undefined);
    const offered = Filters.operators.for(prop);
    return this.options.schema.operators?.(prop, offered) ?? offered;
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
    values(prop, ctx.prefix, abort.signal, {operator: ctx.operator?.id}).then((items) => {
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
    if ((kind !== KIND.INT && kind !== KIND.FLOAT) || (prop.min === undefined && prop.max === undefined))
      return null;
    const range = prop.min !== undefined && prop.max !== undefined ? `${prop.min}–${prop.max}` :
      prop.min !== undefined ? `≥ ${prop.min}` : `≤ ${prop.max}`;
    return `number ${range}`;
  }

  /** Replaces the token under the caret with the row's text and one space (the one already there
   * when the pick lands mid-text), APPLIES it when what it made is a filter, and re-opens for the
   * next expectation. */
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
    // an applied pick re-formats the text, so the caret goes to the end of what it became; a pick
    // that only filled a slot leaves the caret where the value ended
    const at = this._applyPick() ? this._text.peek().length : caret;
    this._input.setSelectionRange?.(at, at);
    this._suggest();
    return true;
  }

  /** A picked row is the user's answer, not a draft: when the text it made is a filter, it is
   * applied in the same action — asking for a second Enter to confirm what was just chosen is a
   * keystroke nobody expects, and the refusal on the status bar would stand over a box that is
   * already right. A text that is not a filter YET (a `between` with one bound, a property with
   * no operator) applies nothing, but the refusal that described the OLD text goes. */
  private _applyPick(): boolean {
    const o = this.options;
    if (Filters.parse(this._text.peek(), o.schema, o.target).problems.length > 0) {
      this._problems.value = [];
      return false;
    }
    return this.commit();
  }

  /** Whether Enter may take the highlighted row. It may when the user put the highlight there
   * (an arrow key or a click), and otherwise only when the row answers what is typed — an
   * auto-highlight over a search that matched nothing would otherwise swap the typed name for an
   * unrelated one on Enter, silently, and eat the keystroke that was meant to apply the filter. */
  private _acceptable(index: number): boolean {
    if (this._typed(index))
      return false;
    if (this._moved)
      return true;
    const item = this._items.peek()[index];
    const ctx = this._context;
    if (item === undefined || ctx === undefined)
      return false;
    const typed = this._text.peek().slice(ctx.replace.start, ctx.replace.end).trim().toLowerCase();
    return typed === '' || item.text.toLowerCase().includes(typed);
  }

  /** The highlighted row spells exactly the token already under the caret — auto-highlight put it
   * there, the user did not. Enter must then commit, not re-insert what is typed and do nothing. */
  private _typed(index: number): boolean {
    const item = this._items.peek()[index];
    const ctx = this._context;
    return item !== undefined && ctx !== undefined &&
      item.text === this._text.peek().slice(ctx.replace.start, ctx.replace.end);
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const open = this._list.isOpen.peek();
    switch (e.key) {
      case 'ArrowDown':
        e.preventDefault();
        if (open) {
          this._moved = true;
          this._list.move(1);
        } else
          this._suggest();
        break;
      case 'ArrowUp':
        if (open) {
          e.preventDefault();
          this._moved = true;
          this._list.move(-1);
        }
        break;
      case 'Enter':
        if (e.ctrlKey || e.metaKey)
          break;
        e.preventDefault();
        if (open && this._acceptable(this._list.activeIndex.peek()) &&
            this._insert(this._list.activeIndex.peek()))
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
