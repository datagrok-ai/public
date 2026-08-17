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
import {signal, computed, Signal, ReadonlySignal} from '../core/signals.js';
import {Component} from '../core/component.js';
import {bindValue} from '../core/bind.js';
import {AsyncSource, AsyncState} from '../core/async-source.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../core/overlay.js';

export interface ComboboxOptions {
  items: string[] | ((query: string, abort: AbortSignal) => Promise<string[]>);
  placeholder?: string;
}

export type ComboboxState = 'idle' | 'focused' | 'open';

type ComboboxEvent =
  {type: 'focus'} | {type: 'blur'} | {type: 'input'} | {type: 'arrowDown'} | {type: 'arrowUp'} |
  {type: 'enter'} | {type: 'escape'} | {type: 'tab'} | {type: 'dismiss'} |
  {type: 'select', index: number};

export class Combobox extends Component {
  private static _seq = 0;

  readonly value: Signal<string> = signal('');
  readonly isOpen: ReadonlySignal<boolean>;
  readonly machineState: ReadonlySignal<ComboboxState>;

  private readonly _id = `u2-combobox-${++Combobox._seq}`;
  private readonly _input = document.createElement('input');
  private readonly _popup = document.createElement('div');
  private readonly _view: ReadonlySignal<AsyncState<string>>;
  private readonly _source: AsyncSource<string> | undefined;
  private readonly _open = signal(false);
  private readonly _machine = signal<ComboboxState>('idle');
  private readonly _activeIndex = signal(-1);
  private readonly _optionEls = signal<HTMLElement[]>([]);
  private _closeOverlay: (() => void) | undefined;

  constructor(options: ComboboxOptions) {
    super();
    this.isOpen = this._open;
    this.machineState = this._machine;

    this.root.classList.add('u2-combobox');
    this.root.dataset.u2 = 'combobox';

    const input = this._input;
    input.className = 'u2-combobox-input';
    input.type = 'text';
    input.autocomplete = 'off';
    input.placeholder = options.placeholder ?? '';
    input.setAttribute('role', 'combobox');
    input.setAttribute('aria-autocomplete', 'list');
    input.setAttribute('aria-controls', `${this._id}-listbox`);
    this.root.append(input);

    this._popup.className = 'u2-combobox-popup';
    this._popup.id = `${this._id}-listbox`;
    this._popup.setAttribute('role', 'listbox');

    if (typeof options.items === 'function') {
      const source = new AsyncSource<string>(options.items);
      this._source = source;
      this._view = source.state;
      this.own(() => source.dispose());
    } else {
      const items = options.items;
      this._view = computed<AsyncState<string>>(() => {
        const q = this.value.value.toLowerCase();
        const matched = q ? items.filter((i) => i.toLowerCase().includes(q)) : items;
        return matched.length > 0 ? {kind: 'ready', items: matched} : {kind: 'empty'};
      });
    }

    bindValue(this.scope, input, this.value);
    this._listen(input, 'input', () => this._send({type: 'input'}));
    this._listen(input, 'focus', () => this._send({type: 'focus'}));
    this._listen(input, 'blur', () => this._send({type: 'blur'}));
    this._listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(this._popup, 'pointerdown', (e) => this._onPopupPointerDown(e));
    this._listen(this._popup, OVERLAY_CLOSE_EVENT, () => this._send({type: 'dismiss'}));

    this.effect(() => input.setAttribute('aria-expanded', String(this._open.value)));
    this.effect(() => this._renderOptions());
    this.effect(() => this._applyActive());
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
          this.value.value = '';
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
    if (this._source)
      this._source.query(this.value.peek());
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
    this.value.value = state.items[index];
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
    const option = target.closest('.u2-combobox-option') as HTMLElement | null;
    if (option) {
      e.preventDefault();
      this._send({type: 'select', index: Number(option.dataset.index)});
      return;
    }
    if (target.closest('.u2-combobox-retry')) {
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
    const state = this._view.value;
    popup.textContent = '';
    if (state.kind === 'ready') {
      const els = state.items.map((item, i) => this._option(item, i));
      popup.append(...els);
      this._optionEls.value = els;
      return;
    }
    this._optionEls.value = [];
    if (state.kind === 'empty') {
      popup.append(Combobox._row('u2-combobox-empty', 'No matches'));
      return;
    }
    if (state.kind === 'error') {
      const row = Combobox._row('u2-combobox-error', state.message);
      const retry = document.createElement('button');
      retry.className = 'u2-combobox-retry';
      retry.type = 'button';
      retry.textContent = 'Retry';
      row.append(retry);
      popup.append(row);
      return;
    }
    popup.append(Combobox._row('u2-combobox-loading', 'Loading…'));
  }

  private _applyActive(): void {
    const els = this._optionEls.value;
    const active = this._activeIndex.value;
    for (let i = 0; i < els.length; i++) {
      els[i].classList.toggle('u2-combobox-option-active', i === active);
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

  private _option(item: string, i: number): HTMLElement {
    const el = Combobox._row('u2-combobox-option', item);
    el.id = `${this._id}-opt-${i}`;
    el.dataset.index = String(i);
    el.setAttribute('role', 'option');
    el.setAttribute('aria-selected', 'false');
    return el;
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
