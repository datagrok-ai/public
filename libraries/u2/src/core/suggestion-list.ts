/* The suggestion popup behind `TypeAhead` and `TagsInput`: one overlay-hosted listbox, the
   active-row machine (arrows/Home/End, ARIA and the highlight on the row, never on the rendered
   content) and the standard hint/loading/empty/error-with-retry rows of an `AsyncSource`. The
   owning control keeps everything that differs — what it offers, what a pick means, and the ARIA
   on its own box. Every element carries a shared `u2-suggest-*` class (the one skin, css/inputs.css)
   beside the owner's prefixed one (automation hooks and the owner's per-control extras). */
import {signal, computed, ReadonlySignal} from './signals.js';
import {Scope} from './scope.js';
import {each} from './region.js';
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
  /** The row's key (`data-u2-item`): a row is built once per key while it stays listed, so a
   * keystroke rebuilds only the rows it brings in. {@link SuggestionList.keyOf} by default. */
  keyOf?: (item: T) => string;
  /** Highlights the first row as soon as a non-empty query has matches, so plain Enter accepts it
   * (VS Code convention); `'always'` does so for the empty query too — a look-up opened on focus.
   * Off for owners whose box holds text of its own worth. */
  autoHighlight?: boolean | 'always';
  /** A row was clicked — index into {@link items}. */
  onPick: (index: number) => void;
  /** The overlay closed on its own (outside pointerdown, Esc, anchor detached): the owner runs
   * its own dismissal path, which ends in {@link dismiss}. */
  onDismiss: () => void;
  onRetry: () => void;
  /** What the empty row says, per query; "No matches" by default. */
  emptyText?: string | ((query: string) => string);
  /** A ceiling on the popup's height, under the room the overlay gives it. */
  maxHeight?: number | string;
  /** At least as wide as the anchor (default); `false` sizes the popup to its rows. */
  matchWidth?: boolean;
  /** A muted line under the rows ("10 of 20 — type to narrow"), read whenever they change; null for none. */
  footer?: () => string | null;
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
  private readonly _list = document.createElement('div');
  private readonly _status = document.createElement('div');
  private readonly _footer = document.createElement('div');
  /** The rows on show: the items while open, past the hint and ready; nothing otherwise. */
  private readonly _shown: ReadonlySignal<T[]>;
  private _rowSeq = 0;
  private _closeOverlay: (() => void) | undefined;

  constructor(options: SuggestionListOptions<T>) {
    this._options = options;
    this.isOpen = this._open;
    this.activeIndex = this._active;
    this._id = `${options.prefix}-${++SuggestionList._seq}`;

    this.popup.className = `u2-suggest-popup ${options.prefix}-popup`;
    this.popup.id = `${this._id}-listbox`;
    this.popup.setAttribute('role', 'listbox');
    options.anchor.setAttribute('aria-controls', this.popup.id);
    this._list.className = 'u2-suggest-options';
    this._footer.className = `u2-suggest-footer ${options.prefix}-footer`;
    this._footer.hidden = true;
    this.popup.append(this._list, this._status, this._footer);
    this._shown = computed(() => {
      const state = options.view.value;
      return this._open.value && options.text.value.length >= options.minChars && state.kind === 'ready' ?
        options.items.value : [];
    });

    options.scope.listen(this.popup, 'pointerdown', (e) => this._onPointerDown(e));
    options.scope.listen(this.popup, OVERLAY_CLOSE_EVENT, () => options.onDismiss());
    const keyOf = options.keyOf;
    Scope.runWith(options.scope, () => each(this._list, this._shown,
      (item, index) => keyOf ? keyOf(item) : SuggestionList.keyOf(item, index), (item) => this._option(item)));
    options.scope.effect(() => this._renderStatus());
    options.scope.effect(() => {
      const text = this._shown.value.length > 0 ? options.footer?.() ?? null : null;
      this._footer.textContent = text ?? '';
      this._footer.hidden = text === null;
    });
    options.scope.effect(() => this._highlight());
    options.scope.effect(() => this._applyActive());
  }

  get rowCount(): number {
    return this._shown.peek().length;
  }

  open(): void {
    if (this._open.peek())
      return;
    this._open.value = true;
    this._closeOverlay = Overlay.show(this._options.anchor, this.popup, this._options.scope,
      {maxHeight: this._options.maxHeight, matchWidth: this._options.matchWidth});
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
    const count = this._shown.peek().length;
    if (count > 0)
      this._active.value = Math.min(Math.max(index, 0), count - 1);
  }

  /** A string or number is its own key; an object answers its `text`, `name` or `id` string;
   * anything else its position (`#<index>`). */
  static keyOf(item: unknown, index: number): string {
    if (typeof item === 'string' || typeof item === 'number')
      return String(item);
    const fields = item as {text?: unknown, name?: unknown, id?: unknown} | null;
    for (const key of [fields?.text, fields?.name, fields?.id]) {
      if (typeof key === 'string')
        return key;
    }
    return `#${index}`;
  }

  static row(cls: string, text: string): HTMLElement {
    const el = document.createElement('div');
    el.className = cls;
    el.textContent = text;
    return el;
  }

  private _row(kind: string, text: string): HTMLElement {
    return SuggestionList.row(`u2-suggest-${kind} ${this._options.prefix}-${kind}`, text);
  }

  // pointerdown, not click: preventDefault keeps focus (and aria-activedescendant) on the box.
  private _onPointerDown(e: Event): void {
    const target = e.target as HTMLElement;
    const option = target.closest(`.${this._options.prefix}-option`) as HTMLElement | null;
    if (option) {
      e.preventDefault();
      this._options.onPick([...this._list.children].indexOf(option));
      return;
    }
    if (target.closest(`.${this._options.prefix}-retry`)) {
      e.preventDefault();
      this._options.onRetry();
    }
  }

  /** The row in place of the options: the hint, loading, the error with its retry, or empty. */
  private _renderStatus(): void {
    const status = this._status;
    const options = this._options;
    const state = options.view.value;
    const short = options.text.value.length < options.minChars;
    status.textContent = '';
    // the row classes set a display, which outranks a bare `hidden`
    status.className = 'u2-suggest-status';
    status.hidden = !this._open.value || !short && state.kind === 'ready' && options.items.value.length > 0;
    if (status.hidden)
      return;
    const kind = short ? 'hint' : state.kind === 'error' ? 'error' :
      state.kind === 'loading' || state.kind === 'idle' ? 'loading' : 'empty';
    status.className = `u2-suggest-${kind} ${options.prefix}-${kind}`;
    if (kind === 'error' && state.kind === 'error') {
      status.textContent = state.message;
      const retry = document.createElement('button');
      retry.className = `u2-suggest-retry ${options.prefix}-retry`;
      retry.type = 'button';
      retry.textContent = 'Retry';
      status.append(retry);
      return;
    }
    const empty = options.emptyText;
    status.textContent = kind === 'hint' ? `Type ${options.minChars} or more characters` :
      kind === 'loading' ? 'Loading…' :
        typeof empty === 'function' ? empty(options.text.peek()) : empty ?? 'No matches';
  }

  /** New rows on show: the first is highlighted when {@link SuggestionListOptions.autoHighlight} asks. */
  private _highlight(): void {
    const rows = this._shown.value;
    const typed = this._options.text.value.length > 0;
    const auto = this._options.autoHighlight;
    if (auto && rows.length > 0)
      this._active.value = auto === 'always' || typed ? 0 : -1;
  }

  /** Built once per key while listed: the id `aria-activedescendant` names stays put. */
  private _option(item: T): HTMLElement {
    const el = document.createElement('div');
    el.className = `u2-suggest-option ${this._options.prefix}-option`;
    el.id = `${this._id}-opt-${++this._rowSeq}`;
    el.setAttribute('role', 'option');
    el.setAttribute('aria-selected', 'false');
    const content = this._options.render(item);
    if (content.querySelector('[data-u2-part="item-label"]') === null)
      content.dataset.u2Part ??= 'item-label';
    el.append(content);
    return el;
  }

  private _applyActive(): void {
    const anchor = this._options.anchor;
    const count = this._shown.value.length;
    const els = this._list.children;
    const active = this._active.value;
    for (let i = 0; i < count; i++) {
      els[i].classList.toggle('u2-suggest-option-active', i === active);
      els[i].classList.toggle(`${this._options.prefix}-option-active`, i === active);
      els[i].setAttribute('aria-selected', i === active ? 'true' : 'false');
    }
    const el = active >= 0 && active < count ? els[active] : undefined;
    if (!el) {
      anchor.removeAttribute('aria-activedescendant');
      return;
    }
    anchor.setAttribute('aria-activedescendant', el.id);
    el.scrollIntoView({block: 'nearest'});
  }
}
