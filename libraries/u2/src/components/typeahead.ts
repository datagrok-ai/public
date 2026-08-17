/* Type-ahead over items of any type: the same hand-written WAI-ARIA 1.2 combobox machine as
   `Combobox` (see its header for the Zag.js evaluation), generalized to `T` plus a caller-supplied
   row renderer. The machine owns the option ROW — hover, active state and ARIA live on the row,
   never on the rendered content. */
import {signal, computed, Signal, ReadonlySignal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {Scope} from '../core/scope.js';
import {bindValue} from '../core/bind.js';
import {AsyncSource, AsyncFetch, AsyncState} from '../core/async-source.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../core/overlay.js';
import {ObjectRenderer} from '../core/object-renderer.js';

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
}

export type TypeAheadState = 'idle' | 'focused' | 'open';

type TypeAheadEvent =
  {type: 'focus'} | {type: 'blur'} | {type: 'input'} | {type: 'arrowDown'} | {type: 'arrowUp'} |
  {type: 'enter'} | {type: 'escape'} | {type: 'tab'} | {type: 'dismiss'} |
  {type: 'select', index: number};

export class TypeAhead<T> extends Component {
  private static _seq = 0;

  readonly selected: Signal<T | null> = signal<T | null>(null);
  readonly text: ReadonlySignal<string>;
  readonly isOpen: ReadonlySignal<boolean>;
  readonly machineState: ReadonlySignal<TypeAheadState>;

  private readonly _id = `u2-typeahead-${++TypeAhead._seq}`;
  private readonly _input = document.createElement('input');
  private readonly _popup = document.createElement('div');
  private readonly _itemText: (item: T) => string;
  private readonly _render: (item: T) => HTMLElement;
  private readonly _minChars: number;
  private readonly _view: ReadonlySignal<AsyncState<T>>;
  private readonly _source: AsyncSource<T> | undefined;
  private readonly _text = signal('');
  private readonly _open = signal(false);
  private readonly _machine = signal<TypeAheadState>('idle');
  private readonly _activeIndex = signal(-1);
  private readonly _optionEls = signal<HTMLElement[]>([]);
  private _rows: Scope | undefined;
  private _closeOverlay: (() => void) | undefined;
  private _echo = false;

  constructor(options: TypeAheadOptions<T>) {
    super();
    this.text = this._text;
    this.isOpen = this._open;
    this.machineState = this._machine;
    const renderer = options.renderer;
    const itemText = options.itemText ?? ((item: T) => renderer!.caption(item));
    this._itemText = itemText;
    const listItem = renderer?.listItem?.bind(renderer) ?? renderer?.markup?.bind(renderer);
    this._render = options.render ?? listItem ??
      ((item) => TypeAhead._row('u2-typeahead-text', itemText(item)));
    this._minChars = options.minChars ?? 0;

    this.root.classList.add('u2-typeahead');
    this.root.dataset.u2 = 'typeahead';

    const input = this._input;
    input.className = 'u2-typeahead-input';
    input.type = 'text';
    input.autocomplete = 'off';
    input.placeholder = options.placeholder ?? '';
    input.setAttribute('role', 'combobox');
    input.setAttribute('aria-autocomplete', 'list');
    input.setAttribute('aria-controls', `${this._id}-listbox`);
    this.root.append(input);

    this._popup.className = 'u2-typeahead-popup';
    this._popup.id = `${this._id}-listbox`;
    this._popup.setAttribute('role', 'listbox');

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

    bindValue(this.scope, input, this._text);
    this._listen(input, 'input', () => this._onInput());
    this._listen(input, 'focus', () => this._send({type: 'focus'}));
    this._listen(input, 'blur', () => this._send({type: 'blur'}));
    this._listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(this._popup, 'pointerdown', (e) => this._onPopupPointerDown(e));
    this._listen(this._popup, OVERLAY_CLOSE_EVENT, () => this._send({type: 'dismiss'}));
    this.own(() => this._releaseRows());

    this.effect(() => input.setAttribute('aria-expanded', String(this._open.value)));
    this.effect(() => this._syncText());
    this.effect(() => this._renderOptions());
    this.effect(() => this._applyActive());
  }

  private _send(ev: TypeAheadEvent): void {
    const state = this._machine.peek();
    switch (ev.type) {
      case 'focus':
        if (state === 'idle')
          this._machine.value = 'focused';
        break;
      case 'blur':
        this._dismiss();
        this._machine.value = 'idle';
        break;
      case 'input':
        this._activeIndex.value = -1;
        this._refresh();
        this._openPopup();
        break;
      case 'arrowDown':
        if (state === 'open') {
          this._move(1);
          break;
        }
        this._refresh();
        this._openPopup();
        break;
      case 'arrowUp':
        if (state === 'open')
          this._move(-1);
        break;
      case 'enter':
        if (state === 'open')
          this._commit(this._activeIndex.peek());
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

  private _onInput(): void {
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
    if (this._open.peek())
      return;
    this._machine.value = 'open';
    this._open.value = true;
    this._closeOverlay = Overlay.show(this._input, this._popup, this.scope);
  }

  private _dismiss(): void {
    if (!this._open.peek())
      return;
    this._open.value = false;
    this._activeIndex.value = -1;
    if (this._machine.peek() === 'open')
      this._machine.value = 'focused';
    const close = this._closeOverlay;
    this._closeOverlay = undefined;
    if (close)
      close();
  }

  private _refresh(): void {
    if (this._source && this._text.peek().length >= this._minChars)
      this._source.query(this._text.peek());
  }

  private _move(delta: number): void {
    const count = this._optionEls.peek().length;
    if (count === 0)
      return;
    this._activeIndex.value = Math.min(Math.max(this._activeIndex.peek() + delta, 0), count - 1);
  }

  private _commit(index: number): void {
    const state = this._view.peek();
    if (state.kind !== 'ready' || index < 0 || index >= state.items.length)
      return;
    this.selected.value = state.items[index];
    this._dismiss();
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
        if (this._open.peek())
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

  // pointerdown, not click: preventDefault keeps focus (and aria-activedescendant) on the input.
  private _onPopupPointerDown(e: Event): void {
    const target = e.target as HTMLElement;
    const option = target.closest('.u2-typeahead-option') as HTMLElement | null;
    if (option) {
      e.preventDefault();
      this._send({type: 'select', index: Number(option.dataset.index)});
      return;
    }
    if (target.closest('.u2-typeahead-retry')) {
      e.preventDefault();
      if (this._source)
        this._source.retry();
    }
  }

  private _renderOptions(): void {
    const popup = this._popup;
    if (!this._open.value) {
      this._optionEls.value = [];
      return;
    }
    const short = this._text.value.length < this._minChars;
    const state = this._view.value;
    this._releaseRows();
    popup.textContent = '';
    if (!short && state.kind === 'ready') {
      const scope = new Scope();
      this._rows = scope;
      const els = Scope.runWith(scope, () => state.items.map((item, i) => this._option(item, i)));
      popup.append(...els);
      this._optionEls.value = els;
      return;
    }
    this._optionEls.value = [];
    if (short) {
      popup.append(TypeAhead._row('u2-typeahead-hint', `Type ${this._minChars} or more characters`));
      return;
    }
    if (state.kind === 'empty') {
      popup.append(TypeAhead._row('u2-typeahead-empty', 'No matches'));
      return;
    }
    if (state.kind === 'error') {
      const row = TypeAhead._row('u2-typeahead-error', state.message);
      const retry = document.createElement('button');
      retry.className = 'u2-typeahead-retry';
      retry.type = 'button';
      retry.textContent = 'Retry';
      row.append(retry);
      popup.append(row);
      return;
    }
    popup.append(TypeAhead._row('u2-typeahead-loading', 'Loading…'));
  }

  private _applyActive(): void {
    const els = this._optionEls.value;
    const active = this._activeIndex.value;
    for (let i = 0; i < els.length; i++) {
      els[i].classList.toggle('u2-typeahead-option-active', i === active);
      els[i].setAttribute('aria-selected', i === active ? 'true' : 'false');
    }
    const el = active >= 0 ? els[active] : undefined;
    if (!el) {
      this._input.removeAttribute('aria-activedescendant');
      return;
    }
    this._input.setAttribute('aria-activedescendant', el.id);
    el.scrollIntoView({block: 'nearest'});
  }

  private _option(item: T, i: number): HTMLElement {
    const el = document.createElement('div');
    el.className = 'u2-typeahead-option';
    el.id = `${this._id}-opt-${i}`;
    el.dataset.index = String(i);
    el.setAttribute('role', 'option');
    el.setAttribute('aria-selected', 'false');
    el.append(this._render(item));
    return el;
  }

  private _releaseRows(): void {
    if (this._rows)
      this._rows.dispose();
    this._rows = undefined;
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }

  private static _row(cls: string, text: string): HTMLElement {
    const el = document.createElement('div');
    el.className = cls;
    el.textContent = text;
    return el;
  }
}
