/* FunctionsBrowser — the u2 counterpart of the Dart FunctionBrowser / FunctionsWidget
   (core/client/xamgle/lib/src/views/functions_view.dart): a search bar speaking the smart-search
   grammar (plain terms AND over a precomputed search string, `#tag` / `@role` terms mirrored into
   the checkbox panes and back), per-tag counts over the term-filtered set, and a `VirtualList`
   over the entries — composed, never forked (the VirtualTree precedent). Platform-free: it
   operates on `FuncItem[]`; the dg factory (src/dg/entities/functions-browser.ts) feeds it
   `DG.Func`. */
import {Control} from '../../core/component.js';
import {batch, computed, signal, Signal, ReadonlySignal} from '../../core/signals.js';
import {button, div, divH, divV, span} from '../../core/elements.js';
import type {BindProp, BindSource} from '../../core/widget-like.js';
import type {LiveOption} from '../../core/input-base.js';
import {TextInput} from '../inputs/text-input.js';
import {iconButton} from '../actions/buttons.js';
import {icon} from '../display/icon.js';
import {VirtualList} from './list.js';
import type {Action} from '../actions/actions.js';

/** One row of the browser: plain strings read once, so typing never walks a platform boundary. */
export interface FuncItem {
  /** Qualified name — the row's identity (`keyOf`). */
  name: string;
  /** Friendly name, what the row shows. */
  label: string;
  description?: string;
  /** Lowercase tags. */
  tags: string[];
  /** Lowercase roles (a platform function's tags-or-role set). */
  roles?: string[];
  /** `'(a, b) : int'` — prebuilt by the dg layer. */
  signature?: string;
  /** `${name} ${label} ${description}`.toLowerCase(), precomputed once — plain search terms
   * match descriptions too (the func-picker precedent). */
  search: string;
  /** The platform object (`DG.Func`), opaque to the core. */
  data?: unknown;
}

export interface FunctionsBrowserOptions {
  items?: FuncItem[] | ReadonlySignal<FuncItem[]>;
  /** Initial query. */
  search?: string;
  /** Pre-checked tags — merged into the query as `#tag` terms. */
  tags?: string[];
  /** Restricts the tags pane to these tags. */
  visibleTags?: string[];
  /** Tags left out of the tags pane; `['res', 'converters', 'internal']` by default. */
  ignoreTags?: string[];
  /** The roles pane; hidden when absent. */
  roles?: {role: string, description: string}[];
  showSearch?: LiveOption<boolean>;
  /** Roles + tags panes; default true. */
  showTags?: LiveOption<boolean>;
  showSignature?: LiveOption<boolean>;
  /** Needs {@link runAction}; default true. */
  showRunButton?: LiveOption<boolean>;
  itemHeight?: number;
  /** The per-row play icon. */
  runAction?: (item: FuncItem) => void;
  /** Row override; `row` is the pooled list row, for per-row attributes. */
  render?: (item: FuncItem, row: HTMLElement) => HTMLElement;
  icon?: (item: FuncItem) => HTMLElement | null;
  contextActions?: (item: FuncItem) => Action[];
  onChanged?: (item: FuncItem | null) => void;
  /** Double-click or Enter. */
  onActivate?: (item: FuncItem) => void;
}

const DEFAULT_IGNORE_TAGS = ['res', 'converters', 'internal'];

interface ParsedQuery {
  terms: string[];
  tags: string[];
  roles: string[];
}

interface PaneRow {
  row: HTMLElement;
  box: HTMLInputElement;
  count: HTMLSpanElement | null;
}

function parseQuery(query: string): ParsedQuery {
  const terms: string[] = [];
  const tags: string[] = [];
  const roles: string[] = [];
  for (const word of query.split(' ')) {
    if (word === '')
      continue;
    if (word.startsWith('#') || word.startsWith('@')) {
      const rest = word.slice(1).toLowerCase();
      if (rest !== '')
        (word.startsWith('#') ? tags : roles).push(rest);
    } else
      terms.push(word.toLowerCase());
  }
  return {terms, tags, roles};
}

/** The Dart `searchFilter` (functions_view.dart:116): every plain term over the precomputed
 * search string, with internal-tagged items always out. */
function matchesTerms(item: FuncItem, terms: string[]): boolean {
  return !item.tags.includes('internal') && terms.every((t) => item.search.includes(t));
}

function sameSet(a: ReadonlySet<string>, b: ReadonlySet<string>): boolean {
  return a.size === b.size && [...a].every((x) => b.has(x));
}

/** The Dart `filter` (functions_view.dart:119): the term filter, then any-of over the checked
 * (or `#`-typed) tags and the checked (or `@`-typed) roles. Pure — exported for tests. */
export function filterFuncItems(items: FuncItem[], query: string,
  checkedTags: ReadonlySet<string>, checkedRoles: ReadonlySet<string>): FuncItem[] {
  const parsed = parseQuery(query);
  const wantTags = new Set([...parsed.tags, ...checkedTags]);
  const wantRoles = new Set([...parsed.roles, ...checkedRoles]);
  return items.filter((item) =>
    matchesTerms(item, parsed.terms) &&
    (wantTags.size === 0 || item.tags.some((t) => wantTags.has(t))) &&
    (wantRoles.size === 0 || (item.roles ?? []).some((r) => wantRoles.has(r))));
}

/** Per-tag counts over the term-filtered set (Dart `refreshCounts`): every tag of every item is
 * a key, counting only the items the plain terms keep. Pure — exported for tests. */
export function tagCounts(items: FuncItem[], query: string): Map<string, number> {
  const {terms} = parseQuery(query);
  const counts = new Map<string, number>();
  for (const item of items) {
    const passes = matchesTerms(item, terms);
    for (const tag of item.tags)
      counts.set(tag, (counts.get(tag) ?? 0) + (passes ? 1 : 0));
  }
  return counts;
}

export class FunctionsBrowser extends Control {
  /** The search box, two-way; the `search` spec prop aliases it through {@link bindStep}. */
  readonly query: Signal<string>;
  readonly selected: ReadonlySignal<FuncItem | null>;
  /** Checked tags, lowercase; writes rewrite the query and vice versa. */
  readonly checkedTags: Signal<ReadonlySet<string>>;
  readonly checkedRoles: Signal<ReadonlySet<string>>;
  readonly showSearch: Signal<boolean>;
  readonly showTags: Signal<boolean>;
  readonly showSignature: Signal<boolean>;
  readonly showRunButton: Signal<boolean>;

  protected readonly options: FunctionsBrowserOptions;
  private readonly _list: VirtualList<FuncItem>;
  private readonly _source = signal<ReadonlySignal<FuncItem[]>>(signal<FuncItem[]>([]));
  private readonly _shown: ReadonlySignal<FuncItem[]>;
  private readonly _selected = signal<FuncItem | null>(null);
  private readonly _tagRows = new Map<string, PaneRow>();
  private readonly _roleRows = new Map<string, PaneRow>();
  private _tagList: string[] = [];
  private _lastShown: FuncItem[] = [];
  private _scrollPending = false;
  private _selectedName?: ReadonlySignal<string | null>;

  constructor(options: FunctionsBrowserOptions = {}) {
    super();
    this.options = options;
    this.root.classList.add('u2-functions-browser');
    this.root.dataset.u2 = 'functions-browser';

    let query = options.search ?? '';
    const seeded = parseQuery(query).tags;
    const preset = (options.tags ?? []).map((t) => t.toLowerCase()).filter((t) => !seeded.includes(t));
    if (preset.length > 0)
      query = [query.trim(), ...preset.map((t) => `#${t}`)].filter((s) => s !== '').join(' ');
    const parsed = parseQuery(query);
    this.query = signal(query);
    this.checkedTags = signal<ReadonlySet<string>>(new Set(parsed.tags));
    this.checkedRoles = signal<ReadonlySet<string>>(new Set(parsed.roles));
    this.showSearch = this._live(options.showSearch, true);
    this.showTags = this._live(options.showTags, true);
    this.showSignature = this._live(options.showSignature, true);
    this.showRunButton = this._live(options.showRunButton, true);
    if (options.items !== undefined)
      this.setItems(options.items);
    this._shown = computed(() => filterFuncItems(this._source.value.value, this.query.value,
      this.checkedTags.value, this.checkedRoles.value));

    const searchInput = this.run(() => new TextInput(
      {search: true, inline: true, placeholder: 'Search by name, #tag, or @role', bind: this.query}));
    // the toggle option mirrors showTags into aria-pressed and the pressed skin, however it flips
    const filterIcon = this.run(() => iconButton('filter', () => {},
      {tooltip: 'Toggle filter panel', toggle: this.showTags}));
    const searchRow = divH([searchInput, filterIcon], 'u2-fb-search');
    searchRow.dataset.u2 = 'fb-search';
    this._listen(searchInput.root, 'keydown', (e) => {
      if ((e as KeyboardEvent).key !== 'ArrowDown')
        return;
      e.preventDefault();
      this._list.root.focus();
      if (this._list.selectedIndex.peek() < 0 && this._shown.peek().length > 0)
        this._list.selectedIndex.value = 0;
    });

    const rolesPane = options.roles === undefined ? null : div([], 'u2-fb-pane');
    if (rolesPane) {
      rolesPane.dataset.u2 = 'fb-roles';
      const rows = (options.roles ?? []).map((r) => {
        const role = r.role.toLowerCase();
        const entry = this._paneRow('fb-role', this.checkedRoles, role, r.role, r.description, false);
        this._roleRows.set(role, entry);
        return entry.row;
      });
      rolesPane.replaceChildren(div(['Roles'], 'u2-fb-pane-header'), ...rows);
    }
    const tagsPane = div([], 'u2-fb-pane');
    tagsPane.dataset.u2 = 'fb-tags';
    const panes = divV(rolesPane ? [rolesPane, tagsPane] : [tagsPane], 'u2-fb-panes');
    panes.dataset.u2 = 'fb-panes';

    const contextActions = options.contextActions;
    this._list = this.run(() => new VirtualList<FuncItem>({
      itemHeight: options.itemHeight ?? 28,
      keyOf: (item) => item.name,
      render: (item, _index, row) => this._renderRow(item, row),
      contextActions: contextActions ? (item) => contextActions(item) : undefined,
    }));
    this._list.root.setAttribute('data-u2', 'fb-list');
    this._list.setItems(this._shown);

    const emptyMessage = span('', 'u2-fb-empty-message');
    const clearSearch = this.run(() => button('Clear search', () => this.query.value = ''));
    clearSearch.classList.add('u2-fb-clear');
    clearSearch.dataset.u2 = 'fb-clear';
    const empty = divV([emptyMessage, clearSearch], 'u2-async-empty u2-fb-empty');
    empty.dataset.u2 = 'fb-empty';
    const status = span('', 'u2-fb-status');
    status.dataset.u2 = 'fb-status';
    const main = divV([this._list, empty, status], 'u2-fb-main');
    this.root.append(searchRow, divH([panes, main], 'u2-fb-body'));

    this.selected = this._selected;
    // `selected` is NOT a naive (shown, selectedIndex) pairing: on an items change, the keyed
    // reselect lives in the list's own effect with no ordering guarantee against this one, so for
    // one flush the OLD index points into the NEW array — emitting whatever sits there (a function
    // the user never clicked). An items change may only keep the selection's key or clear it;
    // any other key mid-change is that stale pairing and is skipped until the reselect lands.
    this.effect(() => {
      const shown = this._shown.value;
      const index = this._list.selectedIndex.value;
      const item = shown[index] ?? null;
      const key = item?.name ?? null;
      const prevKey = this._selected.peek()?.name ?? null;
      const itemsChanged = shown !== this._lastShown;
      this._lastShown = shown;
      if (itemsChanged && key !== prevKey) {
        if (prevKey !== null && shown.some((i) => i.name === prevKey)) {
          this._scrollPending = true;
          return;
        }
        this._selected.value = null;
        this._scrollPending = false;
        return;
      }
      this._selected.value = item;
      if ((itemsChanged || this._scrollPending) && index >= 0)
        this._list.scrollToIndex(index);
      this._scrollPending = false;
    });

    // query → panes: `#tag` / `@role` terms check the boxes (Dart `_inputToCategoryChecks`)
    this.effect(() => {
      const q = parseQuery(this.query.value);
      batch(() => {
        if (!sameSet(new Set(q.tags), this.checkedTags.peek()))
          this.checkedTags.value = new Set(q.tags);
        if (!sameSet(new Set(q.roles), this.checkedRoles.peek()))
          this.checkedRoles.value = new Set(q.roles);
      });
    });
    // panes → query, only when the query does not already say the same — so typing a query that
    // checks a box is never rewritten (normalized) under the user's caret (Dart `categoryClicked`)
    this.effect(() => {
      const tags = this.checkedTags.value;
      const roles = this.checkedRoles.value;
      const q = parseQuery(this.query.peek());
      if (sameSet(new Set(q.tags), tags) && sameSet(new Set(q.roles), roles))
        return;
      this.query.value = [...q.terms, ...[...tags].map((t) => `#${t}`),
        ...[...roles].map((r) => `@${r}`)].join(' ');
    });

    // rows are (re)built only when the tag SET changes; counts and checkmarks update in place,
    // so the checkbox the user just activated keeps its focus
    this.effect(() => this._rebuildTagRows(tagsPane));
    this.effect(() => this._applyPaneState(
      tagCounts(this._source.value.value, this.query.value),
      this.checkedTags.value, this.checkedRoles.value));

    this.effect(() => searchRow.hidden = this.showSearch.value === false);
    this.effect(() => panes.hidden = this.showTags.value === false);
    this.effect(() => this.root.classList.toggle('u2-fb-no-signature', this.showSignature.value === false));
    this.effect(() => this.root.classList.toggle('u2-fb-no-run', this.showRunButton.value === false));

    this.effect(() => {
      const total = this._source.value.value.filter((i) => !i.tags.includes('internal')).length;
      const shown = this._shown.value;
      const found = shown.length > 0;
      status.textContent = `${shown.length} of ${total}`;
      emptyMessage.textContent = found ? '' : total === 0 ? 'No functions.' :
        `No functions match "${this.query.value.trim()}".`;
      clearSearch.hidden = found || total === 0;
      empty.hidden = found;
      this._list.root.hidden = !found;
    });

    let initial = true;
    this.effect(() => {
      const item = this.selected.value;
      if (initial) {
        initial = false;
        return;
      }
      this.options.onChanged?.(item);
      this.fireEvent('change', item);
    });

    this._listen(this._list.root, 'dblclick', (e) => {
      const target = e.target as Element | null;
      const row = target?.closest('.u2-list-row') as HTMLElement | null;
      if (row)
        this._activate(this._shown.peek()[Number(row.dataset.index)]);
    });
    this._listen(this._list.root, 'keydown', (e) => {
      if ((e as KeyboardEvent).key === 'Enter')
        this._activate(this.selected.peek() ?? undefined);
    });
  }

  setItems(items: FuncItem[] | ReadonlySignal<FuncItem[]>): void {
    this._source.value = Array.isArray(items) ? signal(items) : items;
  }

  /** `$.fb.search` is {@link query}; `$.fb.selected` answers the selected item's qualified name
   * (the show* meta props resolve through their same-named signal members in the base). */
  bindStep(name: string): Signal<unknown> | BindSource | null {
    if (name === 'search')
      return this.query as unknown as Signal<unknown>;
    if (name === 'selected') {
      this._selectedName ??= computed(() => this.selected.value?.name ?? null);
      return this._selectedName as unknown as Signal<unknown>;
    }
    return super.bindStep(name);
  }

  bindProps(): BindProp[] {
    const props = super.bindProps();
    props.push({name: 'search', type: 'string', writable: true});
    if (!props.some((p) => p.name === 'selected'))
      props.push({name: 'selected', type: 'string', writable: false});
    return props;
  }

  /** One-way, like `Input.liveOption` (UB-10): the member is always the control's own writable
   * signal — the filter icon writes it freely — and an option given as a signal is followed into
   * it, never adopted, so local toggles never write back into a bound context signal. */
  private _live(value: LiveOption<boolean> | undefined, fallback: boolean): Signal<boolean> {
    if (!(value instanceof Signal))
      return signal((value as boolean | undefined) ?? fallback);
    const own = signal((value.peek() as boolean | undefined) ?? fallback);
    this.effect(() => {
      const x = (value as ReadonlySignal<boolean>).value;
      if (x !== undefined)
        own.value = x; // a forward-ref proxy starts undefined; skip until linked
    });
    return own;
  }

  private _activate(item: FuncItem | undefined): void {
    if (item === undefined)
      return;
    this.options.onActivate?.(item);
    this.fireEvent('activate', item);
  }

  private _renderRow(item: FuncItem, row: HTMLElement): HTMLElement {
    row.setAttribute('data-u2-func', item.name);
    row.setAttribute('aria-label', item.signature !== undefined && this.showSignature.peek() ?
      `${item.label} — ${item.signature}` : item.label);
    if (item.description)
      row.title = item.description;
    const custom = this.options.render;
    if (custom)
      return custom(item, row);
    const parts: HTMLElement[] = [];
    const run = this.options.runAction;
    if (run) {
      // built by hand: iconButton's tooltip would own a cleanup on an ambient scope per pooled
      // render — the plain title dies with the element
      const play = document.createElement('button');
      play.type = 'button';
      play.className = 'u2-btn u2-icon-btn u2-fb-run';
      play.title = 'Run';
      play.setAttribute('aria-label', 'Run');
      play.append(icon('play', {variant: 'solid'}));
      play.addEventListener('click', () => run(item));
      parts.push(play);
    }
    const itemIcon = this.options.icon?.(item);
    if (itemIcon)
      parts.push(div([itemIcon], 'u2-fb-icon'));
    parts.push(span(item.label, 'u2-fb-label'));
    if (item.signature !== undefined)
      parts.push(span(item.signature, 'u2-fb-sig'));
    return div(parts, 'u2-fb-row');
  }

  /** Rebuilds the tags pane only when the tag SET changed (items came or went); the freshly
   * built rows are stamped with the current state, and the state effect keeps them current. */
  private _rebuildTagRows(tagsPane: HTMLElement): void {
    const items = this._source.value.value;
    const ignore = new Set((this.options.ignoreTags ?? DEFAULT_IGNORE_TAGS).map((t) => t.toLowerCase()));
    const visible = new Set((this.options.visibleTags ?? []).map((t) => t.toLowerCase()));
    const tags = [...new Set(items.flatMap((i) => i.tags))]
      .filter((t) => !ignore.has(t) && (visible.size === 0 || visible.has(t)))
      .sort();
    if (tags.length === this._tagList.length && tags.every((t, i) => t === this._tagList[i]))
      return;
    this._tagList = tags;
    this._tagRows.clear();
    const rows = tags.map((tag) => {
      const entry = this._paneRow('fb-tag', this.checkedTags, tag, tag, tag, true);
      this._tagRows.set(tag, entry);
      return entry.row;
    });
    tagsPane.replaceChildren(div(['Tags'], 'u2-fb-pane-header'), ...rows);
    this._applyPaneState(tagCounts(items, this.query.peek()),
      this.checkedTags.peek(), this.checkedRoles.peek());
  }

  private _applyPaneState(counts: Map<string, number>, checkedTags: ReadonlySet<string>,
    checkedRoles: ReadonlySet<string>): void {
    for (const [tag, entry] of this._tagRows) {
      entry.box.checked = checkedTags.has(tag);
      entry.count!.textContent = String(counts.get(tag) ?? 0);
    }
    for (const [role, entry] of this._roleRows)
      entry.box.checked = checkedRoles.has(role);
  }

  /** One pane row: the checkbox toggles membership; a click on the text exclusively checks the
   * item — unless it is the sole checked one, which unchecks it (the Dart `handleClick` toggle).
   * `tooltip` covers the M1 truncation: the full tag name, or a role's description. */
  private _paneRow(kind: string, set: Signal<ReadonlySet<string>>, value: string, text: string,
    tooltip: string, withCount: boolean): PaneRow {
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.className = 'u2-input-checkbox';
    box.setAttribute('aria-label', text);
    box.addEventListener('change', () => {
      const next = new Set(set.peek());
      if (box.checked)
        next.add(value);
      else
        next.delete(value);
      set.value = next;
    });
    const label = span(text, 'u2-fb-pane-label');
    label.title = tooltip || text;
    label.addEventListener('click', () => {
      const current = set.peek();
      set.value = current.has(value) && current.size === 1 ? new Set() : new Set([value]);
    });
    const count = withCount ? span('', 'u2-fb-count') : null;
    const row = div(count ? [box, label, count] : [box, label], 'u2-fb-pane-row');
    row.dataset.u2 = kind;
    row.dataset.value = value;
    return {row, box, count};
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void): void {
    el.addEventListener(type, handler);
    this.own(() => el.removeEventListener(type, handler));
  }
}
