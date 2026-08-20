/* Object tags — the port of Dart's `TagsInput` (`tags_input.dart`): chips for the picked items, a
   box that suggests more as it is typed, and optional creation of items that are not on offer.
   `MultiSelect` stays the string-list control with a checkbox popup; this is its `T`-shaped twin,
   built on the machinery every other u2 async control uses (AsyncSource for suggestions,
   SuggestionList for the popup, `tag()` for the chips, ObjectRenderer for item presentation).

   Two deliberate divergences from Dart: Backspace/Delete on an empty box removes the last chip
   (Dart only closes its drop-down), and the popup never rewrites the typed text with the
   highlighted suggestion — an autocomplete-in-place that surprised people in the platform. */
import {signal, computed, Signal, ReadonlySignal} from '../core/signals.js';
import {Input, InputOptions} from '../core/input-base.js';
import {Scope} from '../core/scope.js';
import {bindValue} from '../core/bind.js';
import {AsyncSource, AsyncFetch, AsyncState} from '../core/async-source.js';
import {SuggestionList} from '../core/suggestion-list.js';
import {ObjectRenderer} from '../core/object-renderer.js';
import {tag} from './badge.js';

export interface TagsInputOptions<T> extends InputOptions<T[]> {
  /** Static suggestions, filtered by the typed text. */
  items?: T[];
  /** Async suggestions; wins over {@link items}. */
  source?: AsyncFetch<T>;
  /** Item presentation as one object: `markup` draws the chips, `listItem` the popup rows. */
  renderer?: ObjectRenderer<T>;
  /** Chip text and the default filter; defaults to `renderer.caption`, then to `String(item)`. */
  itemText?: (item: T) => string;
  /** Popup row; wins over the renderer's members. */
  renderItem?: (item: T) => HTMLElement;
  /** Enter turns the typed text into an item (Dart `allowNew`). */
  allowNew?: boolean;
  /** How that item is built; identity for `T = string`. */
  createItem?: (text: string) => T;
  /** Items it returns false for keep their chip but lose the ✕ (Dart `deleteEnabled`). */
  deleteEnabled?: (item: T) => boolean;
  placeholder?: string;
  /** Below it the popup shows a hint row and an async source is never queried. */
  minChars?: number;
  debounceMs?: number;
}

export class TagsInput<T> extends Input<T[], TagsInputOptions<T>> {
  readonly isOpen: ReadonlySignal<boolean>;

  // every field below is built by createEditor(), which the base constructor calls before
  // subclass initializers would run — none of them may carry an inline initializer
  private _chips!: HTMLElement;
  private _input!: HTMLInputElement;
  private _list!: SuggestionList<T>;
  private _text!: Signal<string>;
  private _view!: ReadonlySignal<AsyncState<T>>;
  private _available!: ReadonlySignal<T[]>;
  private _itemText!: (item: T) => string;
  private _renderItem!: (item: T) => HTMLElement;
  private _markup: ((item: T) => HTMLElement) | undefined;
  private _source: AsyncSource<T> | undefined;
  private _chipScope: Scope | undefined;

  constructor(options: TagsInputOptions<T> = {}) {
    super(options, []);
    this.isOpen = this._list.isOpen;
    this.root.dataset.u2 = 'tags-input';
  }

  protected createEditor(): HTMLElement {
    const renderer = this.options.renderer;
    this._itemText = this.options.itemText ??
      (renderer ? (item: T) => renderer.caption(item) : (item: T) => String(item));
    this._markup = renderer?.markup?.bind(renderer);
    const listItem = renderer?.listItem?.bind(renderer) ?? this._markup;
    this._renderItem = this.options.renderItem ?? listItem ??
      ((item: T) => SuggestionList.row('u2-tags-text', this._itemText(item)));

    this._text = signal('');
    this._buildView();
    this._available = computed(() => {
      const state = this._view.value;
      if (state.kind !== 'ready')
        return [];
      const taken = new Set(this.value.value.map((item) => this._itemText(item)));
      return state.items.filter((item) => !taken.has(this._itemText(item)));
    });

    const editor = document.createElement('div');
    editor.className = 'u2-tags';
    this._chips = document.createElement('div');
    this._chips.className = 'u2-tags-chips';
    const input = document.createElement('input');
    this._input = input;
    input.type = 'text';
    input.autocomplete = 'off';
    input.className = 'u2-tags-input';
    input.placeholder = this.options.placeholder ?? '';
    input.setAttribute('role', 'combobox');
    input.setAttribute('aria-autocomplete', 'list');
    editor.append(this._chips, input);

    this._list = new SuggestionList<T>({
      prefix: 'u2-tags',
      anchor: input,
      scope: this.scope,
      items: this._available,
      view: this._view,
      text: this._text,
      minChars: this.options.minChars ?? 0,
      render: (item) => this._renderItem(item),
      onPick: (index) => {
        this._accept(index);
        input.focus();
      },
      onDismiss: () => this._list.dismiss(),
      onRetry: () => this._source?.retry(),
    });

    bindValue(this.scope, input, this._text);
    this._listen(input, 'input', () => this._onInput());
    this._listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(input, 'blur', () => this._list.dismiss());
    this._listen(editor, 'click', (e) => {
      if (!(e.target as HTMLElement).closest('.u2-tag-remove'))
        input.focus();
    });
    this.own(() => this._releaseChips());

    this.effect(() => input.setAttribute('aria-expanded', String(this._list.isOpen.value)));
    this.effect(() => this._renderChips());
    return editor;
  }

  private _buildView(): void {
    const fetch = this.options.source;
    if (fetch) {
      const source = new AsyncSource<T>(fetch, {debounceMs: this.options.debounceMs});
      this._source = source;
      this._view = source.state;
      this.own(() => source.dispose());
      return;
    }
    const items = this.options.items ?? [];
    this._view = computed<AsyncState<T>>(() => {
      const query = this._text.value.trim().toLowerCase();
      const matched = query ?
        items.filter((item) => this._itemText(item).toLowerCase().includes(query)) : items;
      return matched.length > 0 ? {kind: 'ready', items: matched} : {kind: 'empty'};
    });
  }

  private _onInput(): void {
    this._list.clearActive();
    this._refresh();
    this._list.open();
  }

  private _refresh(): void {
    if (this._source && this._text.peek().length >= (this.options.minChars ?? 0))
      this._source.query(this._text.peek());
  }

  /** Adds an item the box does not already carry, then clears the text and closes the popup —
   * Dart's `displayEnteredItems` (`tags_input.dart:232-249`). */
  add(item: T): void {
    const key = this._itemText(item);
    const current = this.value.peek();
    if (!current.some((existing) => this._itemText(existing) === key))
      this.value.value = [...current, item];
    this._text.value = '';
    this._list.dismiss();
  }

  remove(item: T): void {
    const key = this._itemText(item);
    this.value.value = this.value.peek().filter((existing) => this._itemText(existing) !== key);
  }

  private _accept(index: number): boolean {
    const item = this._available.peek()[index];
    if (item === undefined)
      return false;
    this.add(item);
    return true;
  }

  private _createFromText(): boolean {
    const text = this._text.peek().trim();
    if (!this.options.allowNew || text === '')
      return false;
    const create = this.options.createItem ?? ((typed: string) => typed as unknown as T);
    this.add(create(text));
    return true;
  }

  private _removeLast(): void {
    const current = this.value.peek();
    const last = current[current.length - 1];
    if (last !== undefined && this._deletable(last))
      this.value.value = current.slice(0, current.length - 1);
  }

  private _deletable(item: T): boolean {
    const deleteEnabled = this.options.deleteEnabled;
    return deleteEnabled === undefined || deleteEnabled(item);
  }

  private _onKeyDown(e: KeyboardEvent): void {
    switch (e.key) {
      case 'ArrowDown':
        e.preventDefault();
        if (this._list.isOpen.peek())
          this._list.move(1);
        else {
          this._refresh();
          this._list.open();
        }
        break;
      case 'ArrowUp':
        if (this._list.isOpen.peek()) {
          e.preventDefault();
          this._list.move(-1);
        }
        break;
      case 'Home':
      case 'End':
        if (this._list.isOpen.peek()) {
          e.preventDefault();
          this._list.jump(e.key === 'Home' ? 0 : this._list.rowCount - 1);
        }
        break;
      // an Enter that adds a chip must not also reach the hosting dialog's OK
      // (`tags_input.dart:534-542`); an Enter that does nothing here still may, and Ctrl+Enter
      // is the dialog's own shortcut — the Dart input excludes it from the keys it handles
      // (`tags_input.dart:435-438`), so it must not add a chip either
      case 'Enter':
        if (!e.ctrlKey && !e.metaKey &&
            (this._accept(this._list.activeIndex.peek()) || this._createFromText())) {
          e.preventDefault();
          e.stopPropagation();
        }
        break;
      // the first Esc leaves the browsing state, not the dialog around it (`tags_input.dart:398`)
      case 'Escape':
        if (this._list.isOpen.peek()) {
          e.stopPropagation();
          this._list.dismiss();
        }
        break;
      case 'Tab':
        this._list.dismiss();
        break;
      case 'Backspace':
      case 'Delete':
        if (this._text.peek() === '')
          this._removeLast();
        break;
    }
  }

  private _renderChips(): void {
    const items = this.value.value;
    this._releaseChips();
    this._chips.textContent = '';
    const scope = new Scope();
    this._chipScope = scope;
    this._chips.append(...Scope.runWith(scope, () => items.map((item) => this._chip(item))));
  }

  private _chip(item: T): HTMLElement {
    const chip = tag(this._itemText(item),
      this._deletable(item) ? {onRemove: () => this.remove(item)} : undefined);
    chip.classList.add('u2-tags-chip');
    const markup = this._markup;
    if (markup)
      (chip.querySelector('.u2-tag-text') as HTMLElement).replaceWith(markup(item));
    return chip;
  }

  private _releaseChips(): void {
    if (this._chipScope)
      this._chipScope.dispose();
    this._chipScope = undefined;
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
