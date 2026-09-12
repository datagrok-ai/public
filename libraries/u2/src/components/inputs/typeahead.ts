/* Type-ahead over items of any type: the shared `SuggestionList` popup (see its header) driven by
   a small focused/open state machine, plus the selection-to-text coupling the platform's
   look-up inputs have. The machine owns the option ROW — hover, active state and ARIA live on the
   row, never on the rendered content. */
import {signal, computed, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {bindValue} from '../../core/bind.js';
import {AsyncSource, AsyncFetch, AsyncState} from '../../core/async-source.js';
import {SuggestionList} from '../../core/suggestion-list.js';
import {ObjectRenderer} from '../../core/object-renderer.js';

export interface TypeAheadOptions<T> {
  source: T[] | AsyncFetch<T>;
  /** Input text after a pick, and the default filter for a sync source.
   * Defaults to `renderer.caption`; one of the two must be given. */
  itemText?: (item: T) => string;
  render?: (item: T) => HTMLElement;
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
}

export type TypeAheadState = 'idle' | 'focused' | 'open';

type TypeAheadEvent =
  {type: 'focus'} | {type: 'blur'} | {type: 'input'} | {type: 'arrowDown'} | {type: 'arrowUp'} |
  {type: 'enter'} | {type: 'escape'} | {type: 'tab'} | {type: 'dismiss'} |
  {type: 'select', index: number};

export class TypeAhead<T> extends Control {
  readonly selected: Signal<T | null> = signal<T | null>(null);
  readonly text: ReadonlySignal<string>;
  readonly isOpen: ReadonlySignal<boolean>;
  readonly machineState: ReadonlySignal<TypeAheadState>;

  private readonly _input = document.createElement('input');
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
  private _echo = false;

  constructor(options: TypeAheadOptions<T>) {
    super();
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
    input.placeholder = options.placeholder ?? '';
    input.setAttribute('role', 'combobox');
    input.setAttribute('aria-autocomplete', 'list');
    this.root.append(input);

    if (typeof options.source === 'function') {
      const source = new AsyncSource<T>(options.source, {debounceMs: options.debounceMs});
      this._source = source;
      this._view = source.state;
      this.own(() => source.dispose());
    } else {
      const items = options.source;
      this._view = computed<AsyncState<T>>(() => {
        const q = this._text.value.toLowerCase();
        const matched = q ? items.filter((i) => itemText(i).toLowerCase().includes(q)) : items;
        return matched.length > 0 ? {kind: 'ready', items: matched} : {kind: 'empty'};
      });
    }
    this._items = computed(() => {
      const state = this._view.value;
      return state.kind === 'ready' ? state.items : [];
    });

    this._list = new SuggestionList<T>({
      prefix: 'u2-typeahead',
      anchor: input,
      scope: this.scope,
      items: this._items,
      view: this._view,
      text: this._text,
      minChars: this._minChars,
      render,
      autoHighlight: this._openOnFocus ? 'always' : true,
      onPick: (index) => this._send({type: 'select', index}),
      onDismiss: () => this._send({type: 'dismiss'}),
      onRetry: () => this._source?.retry(),
      emptyText: options.emptyText,
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
    this._listen(input, 'input', () => this._onInput());
    this._listen(input, 'focus', () => this._send({type: 'focus'}));
    this._listen(input, 'blur', () => this._send({type: 'blur'}));
    this._listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));

    this.effect(() => input.setAttribute('aria-expanded', String(this._list.isOpen.value)));
    this.effect(() => this._syncText());
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

  /** Whether an Enter is waiting for the candidates — a host does not treat the box as stray then. */
  get isPickPending(): boolean {
    return this._pendingPick;
  }

  private _onInput(): void {
    this._pendingPick = false;
    if (this.selected.peek() !== null) {
      this._echo = true;
      try {
        this.selected.value = null;
      } finally {
        this._echo = false;
      }
    }
    this._send({type: 'input'});
  }

  /** Puts the text back to the selection's — what a host does with stray text on blur. */
  resync(): void {
    this._text.value = this.selected.peek() === null ? '' : this._itemText(this.selected.peek()!);
  }

  private _syncText(): void {
    const item = this.selected.value;
    if (this._echo)
      return;
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

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
