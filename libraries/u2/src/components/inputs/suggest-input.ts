/* Free-text input with async suggestions — the single-value reduction of TagsInput's box
   (`tags-input.ts`): every keystroke writes the value signal, a SuggestionList popup offers the
   provider's items for the typed text, and a pick replaces the text. The popup reuses the
   typeahead skin through the shared 'u2-typeahead' prefix. */
import {computed, ReadonlySignal} from '../../core/signals.js';
import {Input, InputOptions} from '../../core/input-base.js';
import {bindValue} from '../../core/bind.js';
import {AsyncSource, AsyncFetch} from '../../core/async-source.js';
import {SuggestionList} from '../../core/suggestion-list.js';

export interface SuggestInputOptions extends InputOptions<string> {
  source: AsyncFetch<string>;
  renderItem?: (item: string) => HTMLElement;
  placeholder?: string;
  /** Below it the popup shows a hint row and the source is never queried. */
  minChars?: number;
  debounceMs?: number;
}

export class SuggestInput extends Input<string, SuggestInputOptions> {
  readonly isOpen: ReadonlySignal<boolean>;

  // built by createEditor(), which the base constructor calls before subclass initializers
  private _list!: SuggestionList<string>;
  private _source!: AsyncSource<string>;
  private _items!: ReadonlySignal<string[]>;
  /** The revert point of a closed-popup Escape: the value at focus, moved by every pick and
   * every non-keystroke write — an external write must never be clobbered by Escape. */
  private _anchorValue!: string;
  private _typing = false;

  constructor(options: SuggestInputOptions) {
    super(options, '');
    this.isOpen = this._list.isOpen;
    this.root.classList.add('u2-suggest-input');
    this.root.dataset.u2 = 'suggest-input';
  }

  protected createEditor(): HTMLElement {
    const input = document.createElement('input');
    input.type = 'text';
    input.autocomplete = 'off';
    input.placeholder = this.options.placeholder ?? '';
    input.setAttribute('role', 'combobox');
    input.setAttribute('aria-autocomplete', 'list');

    const source = new AsyncSource<string>(this.options.source, {debounceMs: this.options.debounceMs});
    this._source = source;
    this.own(() => source.dispose());
    this._items = computed(() => {
      const state = source.state.value;
      return state.kind === 'ready' ? state.items : [];
    });

    this._list = new SuggestionList<string>({
      prefix: 'u2-typeahead',
      anchor: input,
      scope: this.scope,
      items: this._items,
      view: source.state,
      text: this.value,
      minChars: this.options.minChars ?? 0,
      render: (item) => this.options.renderItem?.(item) ??
        SuggestionList.row('u2-typeahead-text', item),
      onPick: (index) => {
        this._accept(index);
        input.focus();
      },
      onDismiss: () => this._list.dismiss(),
      onRetry: () => source.retry(),
    });

    // registered before bindValue so the flag is up when its handler writes the value signal
    this._listen(input, 'input', () => this._typing = true);
    bindValue(this.scope, input, this.value);
    this._listen(input, 'focus', () => this._anchorValue = this.value.peek());
    this._listen(input, 'input', () => this._onInput());
    this.effect(() => {
      const v = this.value.value;
      if (!this._typing)
        this._anchorValue = v;
    });
    this._listen(input, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));
    this._listen(input, 'blur', () => this._list.dismiss());
    this.effect(() => input.setAttribute('aria-expanded', String(this._list.isOpen.value)));
    return input;
  }

  private _onInput(): void {
    this._typing = false;
    this._list.clearActive();
    this._refresh();
    this._list.open();
  }

  private _refresh(): void {
    if (this.value.peek().length >= (this.options.minChars ?? 0))
      this._source.query(this.value.peek());
  }

  private _accept(index: number): boolean {
    const item = this._items.peek()[index];
    if (item === undefined)
      return false;
    this.value.value = item;
    this._anchorValue = item;
    this._list.dismiss();
    return true;
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
      // an Enter with the popup open must not reach the hosting dialog's OK
      // (`tags_input.dart:534-542`); Ctrl+Enter stays the dialog's own shortcut
      case 'Enter':
        if (!e.ctrlKey && !e.metaKey && this._list.isOpen.peek()) {
          e.preventDefault();
          e.stopPropagation();
          if (!this._accept(this._list.activeIndex.peek()))
            this._list.dismiss();
        }
        break;
      case 'Escape':
        if (this._list.isOpen.peek()) {
          e.stopPropagation();
          this._list.dismiss();
        }
        else if (this.value.peek() !== this._anchorValue) {
          e.stopPropagation();
          this.value.value = this._anchorValue;
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
