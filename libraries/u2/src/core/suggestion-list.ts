/* The suggestion popup behind `TypeAhead` and `TagsInput`: one overlay-hosted listbox, the
   active-row machine (arrows/Home/End, ARIA and the highlight on the row, never on the rendered
   content) and the standard hint/loading/empty/error-with-retry rows of an `AsyncSource`. The
   owning control keeps everything that differs — what it offers, what a pick means, and the ARIA
   on its own box — and styles the rows through its own class prefix. */
import {signal, ReadonlySignal} from './signals.js';
import {Scope} from './scope.js';
import {AsyncState} from './async-source.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from './overlay.js';

export interface SuggestionListOptions<T> {
  /** Class-name and element-id prefix of the owning control: `u2-typeahead`, `u2-tags`. */
  prefix: string;
  /** The box the popup hangs under; it also carries `aria-controls`/`aria-activedescendant`. */
  anchor: HTMLElement;
  scope: Scope;
  /** What is on offer right now — the owner filters it (`TagsInput` drops what is picked). */
  items: ReadonlySignal<T[]>;
  /** Where those items come from: its non-ready states become the loading and error rows. */
  view: ReadonlySignal<AsyncState<T>>;
  /** The typed text, against {@link minChars}, decides the hint row. */
  text: ReadonlySignal<string>;
  minChars: number;
  render: (item: T) => HTMLElement;
  /** Highlights the first row as soon as a non-empty query has matches, so plain Enter accepts it
   * (VS Code convention). Off for owners whose box holds text of its own worth. */
  autoHighlight?: boolean;
  /** A row was clicked — index into {@link items}. */
  onPick: (index: number) => void;
  /** The overlay closed on its own (outside pointerdown, Esc, anchor detached): the owner runs
   * its own dismissal path, which ends in {@link dismiss}. */
  onDismiss: () => void;
  onRetry: () => void;
}

export class SuggestionList<T> {
  private static _seq = 0;

  readonly popup = document.createElement('div');
  readonly isOpen: ReadonlySignal<boolean>;
  readonly activeIndex: ReadonlySignal<number>;

  private readonly _options: SuggestionListOptions<T>;
  private readonly _id: string;
  private readonly _open = signal(false);
  private readonly _active = signal(-1);
  private readonly _rows = signal<HTMLElement[]>([]);
  private _rowScope: Scope | undefined;
  private _closeOverlay: (() => void) | undefined;

  constructor(options: SuggestionListOptions<T>) {
    this._options = options;
    this.isOpen = this._open;
    this.activeIndex = this._active;
    this._id = `${options.prefix}-${++SuggestionList._seq}`;

    this.popup.className = `${options.prefix}-popup`;
    this.popup.id = `${this._id}-listbox`;
    this.popup.setAttribute('role', 'listbox');
    options.anchor.setAttribute('aria-controls', this.popup.id);

    this._listen(this.popup, 'pointerdown', (e) => this._onPointerDown(e));
    this._listen(this.popup, OVERLAY_CLOSE_EVENT, () => options.onDismiss());
    options.scope.own(() => this._releaseRows());
    options.scope.effect(() => this._renderRows());
    options.scope.effect(() => this._applyActive());
  }

  get rowCount(): number {
    return this._rows.peek().length;
  }

  open(): void {
    if (this._open.peek())
      return;
    this._open.value = true;
    this._closeOverlay = Overlay.show(this._options.anchor, this.popup, this._options.scope);
  }

  dismiss(): void {
    if (!this._open.peek())
      return;
    this._open.value = false;
    this._active.value = -1;
    const close = this._closeOverlay;
    this._closeOverlay = undefined;
    if (close)
      close();
  }

  /** Drops the highlight — the owner calls it when the query changes under the rows. */
  clearActive(): void {
    this._active.value = -1;
  }

  move(delta: number): void {
    this.jump(this._active.peek() + delta);
  }

  jump(index: number): void {
    const count = this._rows.peek().length;
    if (count > 0)
      this._active.value = Math.min(Math.max(index, 0), count - 1);
  }

  static row(cls: string, text: string): HTMLElement {
    const el = document.createElement('div');
    el.className = cls;
    el.textContent = text;
    return el;
  }

  // pointerdown, not click: preventDefault keeps focus (and aria-activedescendant) on the box.
  private _onPointerDown(e: Event): void {
    const target = e.target as HTMLElement;
    const option = target.closest(`.${this._options.prefix}-option`) as HTMLElement | null;
    if (option) {
      e.preventDefault();
      this._options.onPick(Number(option.dataset.index));
      return;
    }
    if (target.closest(`.${this._options.prefix}-retry`)) {
      e.preventDefault();
      this._options.onRetry();
    }
  }

  private _renderRows(): void {
    const prefix = this._options.prefix;
    const popup = this.popup;
    if (!this._open.value) {
      this._rows.value = [];
      return;
    }
    const short = this._options.text.value.length < this._options.minChars;
    const state = this._options.view.value;
    const items = this._options.items.value;
    this._releaseRows();
    popup.textContent = '';
    if (!short && state.kind === 'ready' && items.length > 0) {
      const scope = new Scope();
      this._rowScope = scope;
      const els = Scope.runWith(scope, () => items.map((item, i) => this._option(item, i)));
      popup.append(...els);
      this._rows.value = els;
      if (this._options.autoHighlight)
        this._active.value = this._options.text.value.length > 0 ? 0 : -1;
      return;
    }
    this._rows.value = [];
    if (short) {
      popup.append(SuggestionList.row(`${prefix}-hint`,
        `Type ${this._options.minChars} or more characters`));
    } else if (state.kind === 'error') {
      const row = SuggestionList.row(`${prefix}-error`, state.message);
      const retry = document.createElement('button');
      retry.className = `${prefix}-retry`;
      retry.type = 'button';
      retry.textContent = 'Retry';
      row.append(retry);
      popup.append(row);
    } else if (state.kind === 'loading' || state.kind === 'idle')
      popup.append(SuggestionList.row(`${prefix}-loading`, 'Loading…'));
    else
      popup.append(SuggestionList.row(`${prefix}-empty`, 'No matches'));
  }

  private _option(item: T, i: number): HTMLElement {
    const el = document.createElement('div');
    el.className = `${this._options.prefix}-option`;
    el.id = `${this._id}-opt-${i}`;
    el.dataset.index = String(i);
    el.setAttribute('role', 'option');
    el.setAttribute('aria-selected', 'false');
    el.append(this._options.render(item));
    return el;
  }

  private _applyActive(): void {
    const anchor = this._options.anchor;
    const els = this._rows.value;
    const active = this._active.value;
    for (let i = 0; i < els.length; i++) {
      els[i].classList.toggle(`${this._options.prefix}-option-active`, i === active);
      els[i].setAttribute('aria-selected', i === active ? 'true' : 'false');
    }
    const el = active >= 0 ? els[active] : undefined;
    if (!el) {
      anchor.removeAttribute('aria-activedescendant');
      return;
    }
    anchor.setAttribute('aria-activedescendant', el.id);
    el.scrollIntoView({block: 'nearest'});
  }

  private _releaseRows(): void {
    if (this._rowScope)
      this._rowScope.dispose();
    this._rowScope = undefined;
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this._options.scope.own(() => el.removeEventListener(type, handler));
  }
}
