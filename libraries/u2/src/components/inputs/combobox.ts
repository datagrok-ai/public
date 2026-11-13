/* Sourcing decision (P3 spike): the ARIA/keyboard behavior below is a HAND-WRITTEN state
   machine, not vendored Zag.js.
   Evaluated @zag-js/combobox 1.43.0 as BRAINSTORM.md asks. Its vendorable closure is 12
   packages / 92 ESM modules / ~250 KB (core, popper, dismissable, dom-query, collection,
   utils, types, anatomy, live-region, focus-visible, interact-outside), and a Zag v1 machine
   cannot run on its own: `Service` requires framework-supplied `bindable` reactivity, and
   `connect()` returns prop bags that a render layer must spread onto the DOM after every
   transition. The one vanilla adapter (@zag-js/vanilla) is exactly that glue plus
   @zag-js/store — a second reactivity engine competing with our signals core. That is far
   past the work order's ~15-file / no-framework-glue bar, so the machine is ours: three
   states, typed events, WAI-ARIA 1.2 combobox pattern, and `machineState` exposed for
   getWidgetStatus()-style introspection — the property we wanted from Zag in the first place. */
import {signal, computed, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {Control} from '../../core/component.js';
import {Input} from '../../core/input-base.js';
import {bindValue} from '../../core/bind.js';
import {each} from '../../core/region.js';
import {AsyncSource, AsyncState} from '../../core/async-source.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import type {BindSource, IItemsState} from '../../core/widget-like.js';

export interface ComboboxOptions {
  items: string[] | ((query: string, abort: AbortSignal) => Promise<string[]>);
  placeholder?: string;
  /** The accessible name (`aria-label` on the editor); a combobox shows no caption of its own. */
  label?: string;
  /** Enter with the popup closed: handled here, never seen by the form or dialog around. */
  onEnter?: (e: KeyboardEvent, combobox: Combobox) => void;
  /** Escape with the popup closed: handled here (instead of clearing the text). */
  onEscape?: (e: KeyboardEvent, combobox: Combobox) => void;
}

export type ComboboxState = 'idle' | 'focused' | 'open';

type ComboboxEvent =
  {type: 'focus'} | {type: 'blur'} | {type: 'input'} | {type: 'arrowDown'} | {type: 'arrowUp'} |
  {type: 'enter'} | {type: 'escape'} | {type: 'tab'} | {type: 'dismiss'} |
  {type: 'select', index: number};

/** Free text with a list of suggestions. `on('change', (value, combobox) => …)` fires on every
 * value change. */
export class Combobox extends Control {
  private static _seq = 0;

  /** The text's signal, for wiring. */
  readonly signal: Signal<string> = signal('');
  readonly isOpen: ReadonlySignal<boolean>;
  readonly machineState: ReadonlySignal<ComboboxState>;

  private readonly _id = `u2-combobox-${++Combobox._seq}`;
  private readonly _options: ComboboxOptions;
  private readonly _input = document.createElement('input');
  private readonly _popup = document.createElement('div');
  private readonly _list = document.createElement('div');
  private readonly _status = document.createElement('div');
  private readonly _view: ReadonlySignal<AsyncState<string>>;
  private readonly _source: AsyncSource<string> | undefined;
  private readonly _open = signal(false);
  /** The popup lists everything, whatever the box holds — until the next keystroke. */
  private readonly _all = signal(false);
  private readonly _machine = signal<ComboboxState>('idle');
  private readonly _activeIndex = signal(-1);
  /** What the popup shows: the ready items while open, nothing while closed. */
  private readonly _shown: ReadonlySignal<string[]>;
  private _optionSeq = 0;
  private _closeOverlay: (() => void) | undefined;

  constructor(options: ComboboxOptions) {
    super();
    this._options = options;
    this.isOpen = this._open;
    this.machineState = this._machine;

    this.root.classList.add('u2-combobox');
    this.root.dataset.u2 = 'combobox';

    const input = this._input;
    input.className = 'u2-combobox-input';
    input.type = 'text';
    input.autocomplete = 'off';
    input.placeholder = options.placeholder ?? '';
    input.id = `${this._id}-input`;
    input.setAttribute('role', 'combobox');
    input.setAttribute('aria-autocomplete', 'list');
    input.setAttribute('aria-haspopup', 'listbox');
    input.setAttribute('aria-controls', `${this._id}-listbox`);
    if (options.label)
      input.setAttribute('aria-label', options.label);
    this.root.append(input);

    this._popup.className = 'u2-combobox-popup';
    this._popup.id = `${this._id}-listbox`;
    this._popup.setAttribute('role', 'listbox');
    this._list.className = 'u2-combobox-options';
    this._popup.append(this._list, this._status);

    if (typeof options.items === 'function') {
      const source = new AsyncSource<string>(options.items);
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
      const items = options.items;
      this._view = computed<AsyncState<string>>(() => {
        const q = this._all.value ? '' : this.value.toLowerCase();
        const matched = q ? items.filter((i) => i.toLowerCase().includes(q)) : items;
        return matched.length > 0 ? {kind: 'ready', items: matched} : {kind: 'empty'};
      });
    }
    this._shown = computed(() => {
      const state = this._view.value;
      return this._open.value && state.kind === 'ready' ? state.items : [];
    });

    bindValue(input, this.signal, {owner: this.scope});
    this.scope.listen(input, 'input', () => this._send({type: 'input'}));
    this.scope.listen(input, 'focus', () => this._send({type: 'focus'}));
    this.scope.listen(input, 'blur', () => this._send({type: 'blur'}));
    this.scope.listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this.scope.listen(this._popup, 'pointerdown', (e) => this._onPopupPointerDown(e));
    this.scope.listen(this._popup, OVERLAY_CLOSE_EVENT, () => this._send({type: 'dismiss'}));

    this.effect(() => input.setAttribute('aria-expanded', String(this._open.value)));
    this.build(() => each(this._list, this._shown, (item) => item, (item) => this._option(item)));
    this.effect(() => this._renderStatus());
    this.effect(() => this._applyActive());
    let initial = true;
    this.effect(() => {
      const value = this.value;
      if (initial)
        initial = false;
      else
        untracked(() => this.fireEvent('change', value));
    });
  }

  /** The text: followed when read inside an `effect` or a `computed`; set it to change it. */
  get value(): string {
    return this.signal.value;
  }

  set value(x: string) {
    this.signal.value = x;
  }

  /** The text, not followed. */
  peek(): string {
    return this.signal.peek();
  }

  bindStep(name: string): Signal<unknown> | BindSource | null {
    return name === '' ? this.signal as Signal<unknown> : super.bindStep(name);
  }

  /** Focuses the editor. */
  focus(): void {
    this._input.focus();
  }

  /** The value, the typed text, whether the popup is open, whether an async source is loading. */
  getState(): ReturnType<Control['getState']> {
    const placeholder = this._input.placeholder;
    return {...super.getState(), value: this.peek(), text: this._input.value, open: this._open.peek(),
      busy: this._view.peek().kind === 'loading', ...placeholder ? {placeholder} : {}};
  }

  /** Every option, not the filtered popup: the given items, or what an async source answered last. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const items = typeof this._options.items === 'function' ? this._lastItems() : this._options.items;
    const rendered = new Set(this._shown.peek());
    return Input.itemsState(items.map((item, index) => ({key: item, label: item, index,
      selected: item === this.peek(), rendered: rendered.has(item)})), opts?.all);
  }

  /** Opens the popup over the whole list (an async source is asked for the empty query) and
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

  private _lastItems(): string[] {
    const state = this._view.peek();
    return state.kind === 'ready' ? state.items : [];
  }

  private _send(ev: ComboboxEvent): void {
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
        this._all.value = false;
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
          this.value = '';
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

  private _openPopup(): void {
    if (this._open.peek())
      return;
    this._machine.value = 'open';
    this._open.value = true;
    this._closeOverlay = Overlay.show(this._input, this._popup, {owner: this.scope});
  }

  private _dismiss(): void {
    if (!this._open.peek())
      return;
    this._open.value = false;
    this._all.value = false;
    this._activeIndex.value = -1;
    if (this._machine.peek() === 'open')
      this._machine.value = 'focused';
    const close = this._closeOverlay;
    this._closeOverlay = undefined;
    if (close)
      close();
  }

  private _refresh(): void {
    if (this._source)
      this._source.query(this.peek());
  }

  private _move(delta: number): void {
    const count = this._shown.peek().length;
    if (count === 0)
      return;
    this._activeIndex.value = Math.min(Math.max(this._activeIndex.peek() + delta, 0), count - 1);
  }

  private _commit(index: number): void {
    const items = this._shown.peek();
    if (index < 0 || index >= items.length)
      return;
    this.value = items[index];
    this._dismiss();
    this._input.focus();
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const open = this._open.peek();
    const own = e.key === 'Enter' ? this._options.onEnter : e.key === 'Escape' ? this._options.onEscape : undefined;
    if (own !== undefined && !open && !e.isComposing) {
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
        if (open)
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
    const option = target.closest('.u2-combobox-option') as HTMLElement | null;
    if (option) {
      e.preventDefault();
      this._send({type: 'select', index: [...this._list.children].indexOf(option)});
      return;
    }
    if (target.closest('.u2-combobox-retry')) {
      e.preventDefault();
      if (this._source)
        this._source.retry();
    }
  }

  /** The row under the options while there are none to show: loading, empty or the error. */
  private _renderStatus(): void {
    const status = this._status;
    const state = this._view.value;
    status.textContent = '';
    // the row classes set a display, which outranks a bare `hidden`
    status.className = 'u2-combobox-status';
    status.hidden = !this._open.value || state.kind === 'ready';
    if (status.hidden)
      return;
    if (state.kind === 'error') {
      status.className = 'u2-combobox-error';
      status.textContent = state.message;
      const retry = document.createElement('button');
      retry.className = 'u2-combobox-retry';
      retry.type = 'button';
      retry.textContent = 'Retry';
      status.append(retry);
      return;
    }
    const empty = state.kind === 'empty';
    status.className = empty ? 'u2-combobox-empty' : 'u2-combobox-loading';
    status.textContent = empty ? 'No matches' : 'Loading…';
  }

  private _applyActive(): void {
    const count = this._shown.value.length;
    const active = this._activeIndex.value;
    const els = this._list.children;
    for (let i = 0; i < count; i++) {
      els[i].classList.toggle('u2-combobox-option-active', i === active);
      els[i].setAttribute('aria-selected', i === active ? 'true' : 'false');
    }
    const el = active >= 0 && active < count ? els[active] : undefined;
    if (!el) {
      this._input.removeAttribute('aria-activedescendant');
      return;
    }
    this._input.setAttribute('aria-activedescendant', el.id);
    el.scrollIntoView({block: 'nearest'});
  }

  /** Built once per item while it stays listed: the id `aria-activedescendant` names stays put. */
  private _option(item: string): HTMLElement {
    const el = document.createElement('div');
    el.className = 'u2-combobox-option';
    el.id = `${this._id}-opt-${++this._optionSeq}`;
    el.setAttribute('role', 'option');
    el.setAttribute('aria-selected', 'false');
    const text = document.createElement('span');
    text.dataset.u2Part = 'item-label';
    text.textContent = item;
    el.append(text);
    return el;
  }
}
