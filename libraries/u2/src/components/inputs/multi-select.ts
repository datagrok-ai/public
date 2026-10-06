/* Popup multi-select: the closed field carries the picks as tags, the popup is a checkbox list.
   Same hand-written machine as `Combobox` (see its header for the Zag.js evaluation), with the
   one divergence the multi-select convention demands — toggling an option never closes the popup.
   `MultiChoiceInput` (choice-input.ts) stays the inline, always-open variant. */
import {signal, computed, Signal, ReadonlySignal} from '../../core/signals.js';
import {Input, InputOptions} from '../../core/input-base.js';
import {bindValue} from '../../core/bind.js';
import {Overlay, OVERLAY_CLOSE_EVENT} from '../../core/overlay.js';
import {ChoiceItem} from './choice-input.js';
import {tag} from '../display/badge.js';

export interface MultiSelectOptions extends InputOptions<string[]> {
  items: ChoiceItem[];
  placeholder?: string;
  /** Filter box atop the popup; defaults to on above eight items, resolved at construction. */
  filterable?: boolean;
  /** Tags shown in the closed field before the '+N more' overflow. */
  maxTags?: number;
  selectAll?: boolean;
  /** Free tags (Dart `TagsInput`): Enter in the filter box adds the typed text as an item and
   * selects it, matching an existing item case-insensitively instead of duplicating it. Implies
   * a filter box, since that is where the typing happens. */
  allowCustom?: boolean;
}

export type MultiSelectState = 'idle' | 'focused' | 'open';

type MultiSelectEvent =
  {type: 'focus'} | {type: 'blur'} | {type: 'open'} | {type: 'toggleOpen'} |
  {type: 'arrowDown'} | {type: 'arrowUp'} | {type: 'space'} | {type: 'enter'} |
  {type: 'escape'} | {type: 'tab'} | {type: 'dismiss'} |
  {type: 'toggle', index: number} | {type: 'toggleAll'};

function itemValue(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.value;
}

function itemLabel(item: ChoiceItem): string {
  return typeof item === 'string' ? item : item.label;
}

export class MultiSelect extends Input<string[], MultiSelectOptions> {
  private static _seq = 0;

  readonly isOpen: ReadonlySignal<boolean>;
  readonly machineState: ReadonlySignal<MultiSelectState>;

  // every field below is built by createEditor(), which the base constructor calls before
  // subclass initializers would run — none of them may carry an inline initializer
  private _id!: string;
  private _field!: HTMLElement;
  private _tagsBox!: HTMLElement;
  private _popup!: HTMLElement;
  private _rowsBox!: HTMLElement;
  private _focusEl!: HTMLElement;
  private _filterEl: HTMLInputElement | undefined;
  private _allRow: HTMLElement | undefined;
  private _allLabel: HTMLElement | undefined;
  private _items!: Signal<ChoiceItem[]>;
  private _filter!: Signal<string>;
  private _open!: Signal<boolean>;
  private _machine!: Signal<MultiSelectState>;
  private _activeIndex!: Signal<number>;
  private _rowEls!: Signal<HTMLElement[]>;
  private _visible!: ReadonlySignal<ChoiceItem[]>;
  private _closeOverlay: (() => void) | undefined;

  constructor(options: MultiSelectOptions) {
    super(options, []);
    this.isOpen = this._open;
    this.machineState = this._machine;
    this.root.dataset.u2 = 'multi-select';
  }

  /** Replaces the item list, dropping values that no longer exist. Programmatic semantics:
   * the value is rewritten only when something actually vanished. */
  setItems(items: ChoiceItem[]): void {
    this._items.value = items;
    const current = new Set(this.value.peek());
    const pruned = items.filter((i) => current.has(itemValue(i))).map(itemValue);
    if (pruned.length !== current.size)
      this.value.value = pruned;
  }

  protected createEditor(): HTMLElement {
    this._id = `u2-multi-select-${++MultiSelect._seq}`;
    this._items = signal(this.options.items);
    this._filter = signal('');
    this._open = signal(false);
    this._machine = signal<MultiSelectState>('idle');
    this._activeIndex = signal(-1);
    this._rowEls = signal<HTMLElement[]>([]);
    this._visible = computed(() => {
      const query = this._filter.value.trim().toLowerCase();
      const items = this._items.value;
      return query ? items.filter((i) => itemLabel(i).toLowerCase().includes(query)) : items;
    });

    const field = document.createElement('div');
    this._field = field;
    field.className = 'u2-multi-select-field';
    field.tabIndex = 0;
    field.setAttribute('role', 'combobox');
    field.setAttribute('aria-haspopup', 'listbox');
    field.setAttribute('aria-controls', `${this._id}-listbox`);
    this._tagsBox = MultiSelect._row('u2-multi-select-tags', '');
    const chevron = MultiSelect._row('u2-multi-select-chevron', '');
    field.append(this._tagsBox, chevron);

    this._popup = document.createElement('div');
    this._popup.className = 'u2-multi-select-popup';
    this._popup.tabIndex = -1;
    this._focusEl = this._popup;
    if (this.options.filterable ?? (this.options.allowCustom === true || this.options.items.length > 8))
      this._popup.append(this._buildFilter());
    if (this.options.selectAll ?? true)
      this._popup.append(this._buildSelectAll());
    this._rowsBox = MultiSelect._row('u2-multi-select-options', '');
    this._rowsBox.id = `${this._id}-listbox`;
    this._rowsBox.setAttribute('role', 'listbox');
    this._rowsBox.setAttribute('aria-multiselectable', 'true');
    this._popup.append(this._rowsBox);

    this._listen(field, 'click', (e) => this._onFieldClick(e));
    this._listen(field, 'keydown', (e) => this._onFieldKeyDown(e as KeyboardEvent));
    this._listen(field, 'focus', () => this._send({type: 'focus'}));
    this._listen(field, 'blur', () => this._send({type: 'blur'}));
    this._listen(this._popup, 'keydown', (e) => this._onPopupKeyDown(e as KeyboardEvent));
    this._listen(this._popup, 'pointerdown', (e) => this._onPopupPointerDown(e));
    this._listen(this._popup, OVERLAY_CLOSE_EVENT, () => this._send({type: 'dismiss'}));

    this.effect(() => {
      const open = this._open.value;
      field.setAttribute('aria-expanded', String(open));
      field.classList.toggle('u2-multi-select-open', open);
    });
    this.effect(() => this._renderTags());
    this.effect(() => this._renderRows());
    this.effect(() => this._applySelection());
    this.effect(() => this._applyActive());
    return field;
  }

  private _buildFilter(): HTMLElement {
    const box = MultiSelect._row('u2-multi-select-filter', '');
    const filter = document.createElement('input');
    this._filterEl = filter;
    this._focusEl = filter;
    filter.className = 'u2-multi-select-filter-input';
    filter.type = 'text';
    filter.autocomplete = 'off';
    filter.placeholder = this.options.allowCustom ? 'Filter or add…' : 'Filter…';
    filter.setAttribute('aria-label', 'Filter options');
    filter.setAttribute('aria-controls', `${this._id}-listbox`);
    bindValue(this.scope, filter, this._filter);
    this._listen(filter, 'input', () => this._activeIndex.value = -1);
    box.append(filter);
    return box;
  }

  private _buildSelectAll(): HTMLElement {
    const row = MultiSelect._row('u2-multi-select-all', '');
    row.setAttribute('role', 'button');
    this._allRow = row;
    this._allLabel = MultiSelect._row('u2-multi-select-text', 'Select all');
    row.append(MultiSelect._row('u2-multi-select-check', ''), this._allLabel);
    return row;
  }

  private _send(ev: MultiSelectEvent): void {
    const state = this._machine.peek();
    switch (ev.type) {
      case 'focus':
        if (state === 'idle')
          this._machine.value = 'focused';
        break;
      case 'blur':
        if (state !== 'open')
          this._machine.value = 'idle';
        break;
      case 'open':
        this._openPopup();
        break;
      case 'toggleOpen':
        if (state === 'open')
          this._dismiss(true);
        else
          this._openPopup();
        break;
      case 'arrowDown':
        if (state === 'open')
          this._move(1);
        else
          this._openPopup();
        break;
      case 'arrowUp':
        if (state === 'open')
          this._move(-1);
        break;
      case 'space':
        this._toggleAt(this._activeIndex.peek());
        break;
      case 'enter':
        if (!this.options.allowCustom || !this._addCustom())
          this._dismiss(true);
        break;
      case 'escape':
      case 'tab':
        this._dismiss(true);
        break;
      case 'dismiss':
        this._dismiss(false);
        break;
      case 'toggle':
        this._toggleAt(ev.index);
        break;
      case 'toggleAll':
        this._toggleAll();
        break;
    }
  }

  private _openPopup(): void {
    if (this._open.peek())
      return;
    this._filter.value = '';
    this._activeIndex.value = -1;
    this._machine.value = 'open';
    this._open.value = true;
    this._closeOverlay = Overlay.show(this._field, this._popup, this.scope);
    this._focusEl.focus();
  }

  private _dismiss(refocus: boolean): void {
    if (!this._open.peek())
      return;
    this._open.value = false;
    this._activeIndex.value = -1;
    this._machine.value = refocus ? 'focused' : 'idle';
    const close = this._closeOverlay;
    this._closeOverlay = undefined;
    if (close)
      close();
    if (refocus)
      this._field.focus();
  }

  private _move(delta: number): void {
    const count = this._rowEls.peek().length;
    if (count === 0)
      return;
    this._activeIndex.value = Math.min(Math.max(this._activeIndex.peek() + delta, 0), count - 1);
  }

  private _toggleAt(index: number): void {
    const visible = this._visible.peek();
    if (index < 0 || index >= visible.length)
      return;
    const selected = new Set(this.value.peek());
    const value = itemValue(visible[index]);
    if (!selected.delete(value))
      selected.add(value);
    this._commit(selected);
  }

  private _toggleAll(): void {
    const visible = this._visible.peek();
    const selected = new Set(this.value.peek());
    const on = visible.length > 0 && visible.every((i) => selected.has(itemValue(i)));
    for (const item of visible) {
      if (on)
        selected.delete(itemValue(item));
      else
        selected.add(itemValue(item));
    }
    this._commit(selected);
  }

  /** Adds the typed text as an item and selects it; an item that already exists (by value or
   * label, ignoring case) is selected instead of being added twice. False when there is nothing
   * typed, so Enter keeps its usual "close the popup" meaning. */
  private _addCustom(): boolean {
    const text = this._filter.peek().trim();
    if (text === '')
      return false;
    const items = this._items.peek();
    const lower = text.toLowerCase();
    const existing = items.find((i) =>
      itemValue(i).toLowerCase() === lower || itemLabel(i).toLowerCase() === lower);
    if (!existing)
      this._items.value = [...items, text];
    const selected = new Set(this.value.peek());
    selected.add(existing ? itemValue(existing) : text);
    this._commit(selected);
    this._filter.value = '';
    return true;
  }

  private _commit(selected: Set<string>): void {
    this.value.value = this._items.peek().filter((i) => selected.has(itemValue(i))).map(itemValue);
  }

  // the tag's own ✕ handler already removed the value — the click must not also open the popup
  private _onFieldClick(e: Event): void {
    if (!(e.target as HTMLElement).closest('.u2-tag-remove'))
      this._send({type: 'toggleOpen'});
  }

  private _remove(value: string): void {
    const selected = new Set(this.value.peek());
    selected.delete(value);
    this._commit(selected);
  }

  private _onFieldKeyDown(e: KeyboardEvent): void {
    switch (e.key) {
      case 'ArrowDown':
      case 'Enter':
      case ' ':
        e.preventDefault();
        this._send({type: 'open'});
        break;
      case 'Escape':
        this._send({type: 'escape'});
        break;
    }
  }

  private _onPopupKeyDown(e: KeyboardEvent): void {
    switch (e.key) {
      case 'ArrowDown':
        e.preventDefault();
        this._send({type: 'arrowDown'});
        break;
      case 'ArrowUp':
        e.preventDefault();
        this._send({type: 'arrowUp'});
        break;
      // with a filter box, Space is a character until a row is active
      case ' ':
        if (this._activeIndex.peek() < 0 && this._filterEl)
          break;
        e.preventDefault();
        this._send({type: 'space'});
        break;
      case 'Enter':
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

  // pointerdown, not click: preventDefault keeps focus (and aria-activedescendant) in the popup.
  private _onPopupPointerDown(e: Event): void {
    const target = e.target as HTMLElement;
    const option = target.closest('.u2-multi-select-option') as HTMLElement | null;
    if (option) {
      e.preventDefault();
      this._send({type: 'toggle', index: Number(option.dataset.index)});
      return;
    }
    if (target.closest('.u2-multi-select-all')) {
      e.preventDefault();
      this._send({type: 'toggleAll'});
    }
  }

  private _renderTags(): void {
    const selected = this.value.value;
    const labels = new Map(this._items.value.map((i): [string, string] => [itemValue(i), itemLabel(i)]));
    const box = this._tagsBox;
    box.textContent = '';
    if (selected.length === 0) {
      const placeholder = this.options.placeholder;
      if (placeholder)
        box.append(MultiSelect._row('u2-multi-select-placeholder', placeholder));
      return;
    }
    const max = this.options.maxTags ?? 3;
    for (const value of selected.slice(0, max))
      box.append(tag(labels.get(value) ?? value, {onRemove: () => this._remove(value)}));
    if (selected.length > max)
      box.append(MultiSelect._row('u2-multi-select-more', `+${selected.length - max} more`));
  }

  private _renderRows(): void {
    const box = this._rowsBox;
    box.textContent = '';
    if (!this._open.value) {
      this._rowEls.value = [];
      return;
    }
    const visible = this._visible.value;
    if (visible.length === 0) {
      box.append(MultiSelect._row('u2-multi-select-empty', 'No matches'));
      this._rowEls.value = [];
      return;
    }
    const els = visible.map((item, i) => this._option(item, i));
    box.append(...els);
    this._rowEls.value = els;
  }

  private _option(item: ChoiceItem, i: number): HTMLElement {
    const el = MultiSelect._row('u2-multi-select-option', '');
    el.id = `${this._id}-opt-${i}`;
    el.dataset.index = String(i);
    el.dataset.value = itemValue(item);
    el.setAttribute('role', 'option');
    el.setAttribute('aria-selected', 'false');
    el.append(MultiSelect._row('u2-multi-select-check', ''),
      MultiSelect._row('u2-multi-select-text', itemLabel(item)));
    return el;
  }

  private _applySelection(): void {
    const selected = new Set(this.value.value);
    const els = this._rowEls.value;
    for (const el of els) {
      const on = selected.has(el.dataset.value as string);
      el.setAttribute('aria-selected', String(on));
      el.classList.toggle('u2-multi-select-selected', on);
    }
    if (!this._allRow || !this._allLabel)
      return;
    const all = els.length > 0 && els.every((el) => selected.has(el.dataset.value as string));
    this._allRow.classList.toggle('u2-multi-select-selected', all);
    this._allLabel.textContent = all ? 'Select none' : 'Select all';
  }

  private _applyActive(): void {
    const els = this._rowEls.value;
    const active = this._activeIndex.value;
    for (let i = 0; i < els.length; i++)
      els[i].classList.toggle('u2-multi-select-option-active', i === active);
    const el = active >= 0 ? els[active] : undefined;
    if (!el) {
      this._focusEl.removeAttribute('aria-activedescendant');
      return;
    }
    this._focusEl.setAttribute('aria-activedescendant', el.id);
    el.scrollIntoView({block: 'nearest'});
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }

  private static _row(cls: string, text: string): HTMLElement {
    const el = document.createElement('div');
    el.className = cls;
    if (text)
      el.textContent = text;
    return el;
  }
}
