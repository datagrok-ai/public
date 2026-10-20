/* Type-ahead over items of any type: the shared `SuggestionList` popup (see its header) driven by
   a small focused/open state machine, plus the selection-to-text coupling the platform's
   look-up inputs have. The machine owns the option ROW — hover, active state and ARIA live on the
   row, never on the rendered content. */
import {signal, computed, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Input} from '../../core/input-base.js';
import {bindValue} from '../../core/bind.js';
import {AsyncSource, AsyncFetch, AsyncState} from '../../core/async-source.js';
import {SuggestionList} from '../../core/suggestion-list.js';
import {ObjectRenderer} from '../../core/object-renderer.js';
import type {IItemsState} from '../../core/widget-like.js';

export interface TypeAheadOptions<T> {
  source: T[] | AsyncFetch<T>;
  /** Input text after a pick, and the default filter for a sync source.
   * Defaults to `renderer.caption`; one of the two must be given. */
  itemText?: (item: T) => string;
  render?: (item: T) => HTMLElement;
  /** The option's key (`data-u2-item`, the items' keys); `itemText` by default. Give one where two
   * items may read the same. */
  keyOf?: (item: T) => string;
  /** Item presentation as one object; explicit `itemText`/`render` win over its members. */
  renderer?: ObjectRenderer<T>;
  placeholder?: string;
  /** Below it the popup shows a hint row and an async source is never queried. */
  minChars?: number;
  debounceMs?: number;
  /** Opens the candidates on focus, the empty query listing the first of them — the platform's
   * look-up inputs do; off by default. */
  openOnFocus?: boolean;
  /** What the empty row says, per query; "No matches" by default. */
  emptyText?: string | ((query: string) => string);
  /** A ceiling on the popup's height (a CSS length or pixels), under the room the overlay gives it. */
  maxHeight?: number | string;
  /** The accessible name (`aria-label` on the editor); a type-ahead shows no caption of its own. */
  label?: string;
  /** Enter with the popup closed: handled here, never seen by the form or dialog around. */
  onEnter?: (e: KeyboardEvent, typeAhead: TypeAhead<T>) => void;
  /** Escape with the popup closed: handled here (instead of clearing the selection). */
  onEscape?: (e: KeyboardEvent, typeAhead: TypeAhead<T>) => void;
}

export type TypeAheadState = 'idle' | 'focused' | 'open';

type TypeAheadEvent =
  {type: 'focus'} | {type: 'blur'} | {type: 'input'} | {type: 'arrowDown'} | {type: 'arrowUp'} |
  {type: 'enter'} | {type: 'escape'} | {type: 'tab'} | {type: 'dismiss'} |
  {type: 'select', index: number};

/** Picks one item by typing. `on('change', (item, typeAhead) => …)` fires on every change of
 * {@link selected}. */
export class TypeAhead<T> extends Control {
  readonly selected: Signal<T | null> = signal<T | null>(null);
  readonly text: ReadonlySignal<string>;
  readonly isOpen: ReadonlySignal<boolean>;
  readonly machineState: ReadonlySignal<TypeAheadState>;

  private readonly _input = document.createElement('input');
  private readonly _options: TypeAheadOptions<T>;
  private readonly _list: SuggestionList<T>;
  private readonly _itemText: (item: T) => string;
  private readonly _minChars: number;
  private readonly _openOnFocus: boolean;
  /** Enter arrived while the candidates were still on their way: the pick is made when they land. */
  private _pendingPick = false;
  private readonly _view: ReadonlySignal<AsyncState<T>>;
  private readonly _items: ReadonlySignal<T[]>;
  private readonly _source: AsyncSource<T> | undefined;
  private readonly _text = signal('');
  private readonly _machine = signal<TypeAheadState>('idle');
  /** The popup lists everything, whatever the box holds — until the next keystroke. */
  private readonly _all = signal(false);
  private _echo = false;

  constructor(options: TypeAheadOptions<T>) {
    super();
    this._options = options;
    this.text = this._text;
    this.machineState = this._machine;
    const renderer = options.renderer;
    const itemText = options.itemText ?? ((item: T) => renderer!.caption(item));
    this._itemText = itemText;
    const listItem = renderer?.listItem?.bind(renderer) ?? renderer?.markup?.bind(renderer);
    const render = options.render ?? listItem ??
      ((item: T) => SuggestionList.row('u2-typeahead-text', itemText(item)));
    this._minChars = options.minChars ?? 0;
    this._openOnFocus = options.openOnFocus === true;

    this.root.classList.add('u2-typeahead');
    this.root.dataset.u2 = 'typeahead';

    const input = this._input;
    input.className = 'u2-typeahead-input';
    input.type = 'text';
    input.autocomplete = 'off';
    input.spellcheck = false;
    input.placeholder = options.placeholder ?? '';
    input.setAttribute('role', 'combobox');
    input.setAttribute('aria-autocomplete', 'list');
    input.setAttribute('aria-haspopup', 'listbox');
    if (options.label)
      input.setAttribute('aria-label', options.label);
    this.root.append(input);

    if (typeof options.source === 'function') {
      const source = new AsyncSource<T>(options.source, {debounceMs: options.debounceMs});
      this._source = source;
      this._view = source.state;
      this.own(() => source.dispose());
      this.effect(() => {
        if (source.state.value.kind === 'loading')
          this.root.setAttribute('aria-busy', 'true');
        else
          this.root.removeAttribute('aria-busy');
      });
    } else {
      const items = options.source;
      this._view = computed<AsyncState<T>>(() => {
        const q = this._all.value ? '' : this._text.value.toLowerCase();
        const matched = q ? items.filter((i) => itemText(i).toLowerCase().includes(q)) : items;
        return matched.length > 0 ? {kind: 'ready', items: matched} : {kind: 'empty'};
      });
    }
    this._items = computed(() => {
      const state = this._view.value;
      return state.kind === 'ready' ? state.items : [];
    });

    const minChars = this._minChars;
    const all = this._all;
    this._list = new SuggestionList<T>({
      prefix: 'u2-typeahead',
      anchor: input,
      scope: this.scope,
      items: this._items,
      view: this._view,
      text: this._text,
      // read inside the list's computeds, so a reveal (`_all`) drops the hint row until it ends
      get minChars() { return all.value ? 0 : minChars; },
      render,
      keyOf: options.keyOf ?? itemText,
      autoHighlight: this._openOnFocus ? 'always' : true,
      onPick: (index) => this._send({type: 'select', index}),
      onDismiss: () => this._send({type: 'dismiss'}),
      onRetry: () => this._source?.retry(),
      emptyText: options.emptyText,
      maxHeight: options.maxHeight,
    });
    this.isOpen = this._list.isOpen;
    // a fast Enter: the intent was clear — the name typed, else the first candidate, once they land
    this.effect(() => {
      const view = this._view.value;
      if (!this._pendingPick || view.kind === 'loading' || view.kind === 'idle')
        return;
      this._pendingPick = false;
      if (view.kind !== 'ready')
        return;
      const typed = this._text.peek().trim().toLowerCase();
      const exact = view.items.findIndex((item) => itemText(item).toLowerCase() === typed);
      this._commit(exact >= 0 ? exact : 0, document.activeElement === input);
    });

    bindValue(this.scope, input, this._text);
    this.scope.listen(input, 'input', () => this._onInput());
    this.scope.listen(input, 'focus', () => this._send({type: 'focus'}));
    this.scope.listen(input, 'blur', () => this._send({type: 'blur'}));
    this.scope.listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));

    this.effect(() => input.setAttribute('aria-expanded', String(this._list.isOpen.value)));
    this.effect(() => this._syncText());
    let initial = true;
    this.effect(() => {
      const item = this.selected.value;
      if (initial)
        initial = false;
      else
        untracked(() => this.fireEvent('change', item));
    });
  }

  /** Focuses the editor. */
  focus(): void {
    this._input.focus();
  }

  /** The selection (JSON-projected), the typed text, whether the popup is open, whether an async
   * source is loading. */
  getState(): ReturnType<Control['getState']> {
    const placeholder = this._input.placeholder;
    return {...super.getState(), value: Control.jsonValue(this.selected.peek()) ?? null, text: this._input.value,
      open: this._list.isOpen.peek(), busy: this._view.peek().kind === 'loading', ...placeholder ? {placeholder} : {}};
  }

  /** Every candidate, not the filtered popup: the given items, or what an async source answered last. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const source = this._options.source;
    const items = typeof source === 'function' ? this._lastItems() : source;
    const rendered = new Set(this._list.isOpen.peek() ? this._items.peek() : []);
    const keyOf = this._options.keyOf ?? this._itemText;
    return Input.itemsState(items.map((item, index) =>
      ({key: keyOf(item), label: this._itemText(item), index, rendered: rendered.has(item)})), opts?.all);
  }

  /** Opens the popup over every candidate (an async source is asked for the empty query) and
   * resolves to the option's element. */
  async revealItem(keyOrLabel: string): Promise<HTMLElement | null> {
    this._all.value = true;
    if (this._source)
      this._source.query('');
    this._openPopup();
    for (let i = 0; i < 50 && this._view.peek().kind === 'loading'; i++)
      await new Promise((resolve) => setTimeout(resolve, 20));
    return super.revealItem(keyOrLabel);
  }

  protected accessibleLabel(): string | undefined {
    return super.accessibleLabel() ?? Control.labelOf(this._input) ?? (this._input.placeholder || undefined);
  }

  private _lastItems(): T[] {
    const state = this._view.peek();
    return state.kind === 'ready' ? state.items : [];
  }

  private _send(ev: TypeAheadEvent): void {
    const state = this._machine.peek();
    switch (ev.type) {
      case 'focus':
        if (state === 'idle')
          this._machine.value = 'focused';
        // the first candidates, not the held item's name as a query
        if (this._openOnFocus) {
          if (this._source && this.selected.peek() !== null)
            this._source.query('');
          else
            this._refresh();
          this._openPopup();
        }
        break;
      case 'blur':
        this._dismiss();
        this._machine.value = 'idle';
        // text that is not a pick never survives a blur; a host that restores a pick of its own
        // (its blur listener runs after this one) has had its say by then
        queueMicrotask(() => {
          if (!this.scope.isDisposed && !this._pendingPick && this.selected.peek() === null && this._text.peek() !== '')
            this.resync();
        });
        break;
      case 'input':
        this._refresh();
        this._openPopup();
        break;
      case 'arrowDown':
        if (state === 'open') {
          this._list.move(1);
          break;
        }
        this._refresh();
        this._openPopup();
        break;
      case 'arrowUp':
        if (state === 'open')
          this._list.move(-1);
        break;
      case 'enter':
        if (state !== 'open')
          break;
        if (this._view.peek().kind === 'ready')
          this._commit(this._list.activeIndex.peek());
        else
          this._pendingPick = true;
        break;
      case 'escape':
        if (state === 'open')
          this._dismiss();
        else
          this._clear();
        break;
      case 'tab':
      case 'dismiss':
        this._dismiss();
        break;
      case 'select':
        this._commit(ev.index);
        break;
    }
  }

  /** Closes the candidates without touching the selection — what a host does when something
   * outside the box (a save, a page move) must not be covered by the popup. */
  close(): void {
    this._send({type: 'dismiss'});
  }

  /** Whether an Enter is waiting for the candidates — a host does not treat the box as stray then. */
  get isPickPending(): boolean {
    return this._pendingPick;
  }

  private _onInput(): void {
    this._pendingPick = false;
    this._all.value = false;
    // the handler runs batched: the text sync sees the flag when the batch ends, and drops it
    if (this.selected.peek() !== null) {
      this._echo = true;
      this.selected.value = null;
    }
    this._send({type: 'input'});
  }

  /** Puts the text back to the selection's — what a host does with stray text on blur. */
  resync(): void {
    this._text.value = this.selected.peek() === null ? '' : this._itemText(this.selected.peek()!);
  }

  private _syncText(): void {
    const item = this.selected.value;
    if (this._echo) {
      this._echo = false;
      return;
    }
    this._text.value = item === null ? '' : this._itemText(item);
  }

  private _clear(): void {
    this.selected.value = null;
    this._text.value = '';
  }

  private _openPopup(): void {
    if (this._list.isOpen.peek())
      return;
    this._machine.value = 'open';
    this._list.open();
  }

  private _dismiss(): void {
    if (!this._list.isOpen.peek())
      return;
    this._all.value = false;
    if (this._machine.peek() === 'open')
      this._machine.value = 'focused';
    this._list.dismiss();
  }

  private _refresh(): void {
    if (this._source && this._text.peek().length >= this._minChars)
      this._source.query(this._text.peek());
  }

  private _commit(index: number, focus = true): void {
    const item = this._items.peek()[index];
    if (item === undefined)
      return;
    this.selected.value = item;
    this._dismiss();
    if (focus)
      this._input.focus();
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const own = e.key === 'Enter' ? this._options.onEnter : e.key === 'Escape' ? this._options.onEscape : undefined;
    if (own !== undefined && !this._list.isOpen.peek() && !e.isComposing) {
      e.preventDefault();
      e.stopPropagation();
      own(e, this);
      return;
    }
    switch (e.key) {
      case 'ArrowDown':
        e.preventDefault();
        this._send({type: 'arrowDown'});
        break;
      case 'ArrowUp':
        e.preventDefault();
        this._send({type: 'arrowUp'});
        break;
      case 'Enter':
        if (this._list.isOpen.peek())
          e.preventDefault();
        this._send({type: 'enter'});
        break;
      case 'Escape':
        this._send({type: 'escape'});
        break;
      case 'Tab':
        this._send({type: 'tab'});
        break;
    }
  }
}
