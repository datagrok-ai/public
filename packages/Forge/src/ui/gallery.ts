import * as grok from 'datagrok-api/grok';
import * as ui from 'datagrok-api/ui';
import * as DG from 'datagrok-api/dg';
import * as rxjs from 'rxjs';
import {isRecord} from '../preparation/preparation-options';
import {tagsOf} from '../storage/model-fields';
import '../../css/forge.css';

export const GALLERY_MODES = ['list', 'cards', 'table'] as const;
export type GalleryMode = typeof GALLERY_MODES[number];
/** At most this many cards or list rows are rendered; the table shows every row. */
export const MAX_ITEMS = 500;
const MODE_ATTRIBUTES: Record<GalleryMode, string> = {list: 'Brief', cards: 'Card', table: 'Grid'};
const MODE_ICONS: Record<GalleryMode, [icon: string, tooltip: string]> = {
  list: ['grip-lines', 'Show as a list'],
  cards: ['grip-horizontal', 'Show as cards'],
  table: ['table', 'Show as a table'],
};
const DIRECTIONS = ['Ascending', 'Descending'];
/** Text sorts ignore case and accents: `k-Nearest Neighbors` between `Hist Gradient Boosting` and `LDA`. */
const compareText = new Intl.Collator(undefined, {sensitivity: 'base'}).compare;

export interface SortField { caption: string; column: string }
export interface GallerySort { column: string; desc: boolean }

export interface GalleryOptions {
  /** The user-settings key of the mode and the sort. */
  key: string;
  /** What the items are, in the count's tooltip: `models`. */
  noun: string;
  placeholder: string;
  /** The table mode; its frame is the state: filter, current row, selection. After the owner replaces its frame, the
   * owner calls {@link Gallery.rebind}. */
  grid: DG.Grid;
  /** Searched case-insensitively for the search text (without its `#tag` words). */
  searchColumns: string[];
  /** Each `#word` of the search must be one of this column's comma-separated tags. */
  tagsColumn?: string;
  sortFields: SortField[];
  defaultSort: GallerySort;
  card: (row: number) => HTMLElement;
  listItem: (row: number) => HTMLElement;
  tooltip: (row: number) => HTMLElement | null;
  /** Fills the item's context menu; the gallery adds nothing itself. */
  menu: (row: number, menu: DG.Menu) => void;
  dragObject?: (row: number) => unknown;
  /** Icons left of the search box. */
  actions?: HTMLElement;
  /** A line under the toolbar. */
  below?: HTMLElement;
  /** Shown instead of the items while the frame has no rows. */
  emptyHint: HTMLElement;
  extraFilter?: (row: number) => boolean;
}

/** The platform gallery's look (search, list, cards and table, sort, count) over a frame: cards and list render the
 * filtered rows in the sort order, the table is the owner's grid over the same frame. */
export class Gallery {
  /** The `grok.userSettings` storage of the mode and the sort, one key per gallery. */
  static storage = 'Forge view';
  readonly root: HTMLDivElement;
  readonly search: DG.InputBase<string>;
  private readonly options: GalleryOptions;
  private readonly body = ui.div([], 'grok-gallery-grid');
  private readonly count = ui.div([], 'grok-items-view-counts');
  private readonly modeIcons: Record<GalleryMode, HTMLElement>;
  private readonly gridSub: rxjs.Subscription;
  private frameSubs: rxjs.Subscription[] = [];
  private items = new Map<number, HTMLElement>();
  /** Every item built for the frame and the mode: a search or a sort only moves them. */
  private built = new Map<number, HTMLElement>();
  /** The sort the grid has, so a search does not sort it again. */
  private gridSortKey: string | null = null;
  private shown: number[] = [];
  private currentMode: GalleryMode;
  private currentSort: GallerySort;
  private renders = 0;

  constructor(options: GalleryOptions) {
    this.options = options;
    const stored = Gallery.stored(options.key);
    this.currentMode = GALLERY_MODES.find((m) => m === stored.mode) ?? 'cards';
    const sort = stored.sort;
    this.currentSort = isRecord(sort) && typeof sort.desc === 'boolean' &&
      options.sortFields.some((f) => f.column === sort.column) ?
      {column: String(sort.column), desc: sort.desc} : options.defaultSort;
    this.search = ui.input.search('', {placeholder: options.placeholder, onValueChanged: () => this.refilter()});
    this.modeIcons = {
      list: this.modeIcon('list'),
      cards: this.modeIcon('cards'),
      table: this.modeIcon('table'),
    };
    const sortIcon = ui.iconFA('sort-alt', (e: MouseEvent) => this.showSortMenu(e), 'Sort');
    const commands = ui.div([...GALLERY_MODES.map((m) => this.modeIcons[m]), sortIcon],
      'grok-items-view-commands grok-items-view-toggle');
    ui.tooltip.bind(this.count, () => this.countTooltip());
    const toolbar = ui.div([...(options.actions === undefined ? [] : [options.actions]), commands, this.count,
      ui.div([this.search.root], 'd4-search-bar')], 'grok-gallery-search-bar');
    options.actions?.classList.add('forge-gallery-actions');
    options.below?.classList.add('forge-gallery-below');
    this.root = ui.divV([toolbar, ...(options.below === undefined ? [] : [options.below]), this.body],
      'forge-gallery grok-card-view');
    options.grid.root.style.width = '100%';
    options.grid.root.style.height = '100%';
    this.gridSub = options.grid.onRowsSorted.subscribe(() => this.followGridSort());
    this.rebind();
  }

  get frame(): DG.DataFrame {
    return this.options.grid.dataFrame;
  }

  get mode(): GalleryMode {
    return this.currentMode;
  }

  get sort(): GallerySort {
    return {...this.currentSort};
  }

  setMode(mode: GalleryMode): void {
    if (mode === this.currentMode)
      return;
    this.currentMode = mode;
    this.built.clear();
    this.gridSortKey = null;
    this.store();
    this.render();
  }

  setSort(column: string, desc: boolean): void {
    if (column === this.currentSort.column && desc === this.currentSort.desc)
      return;
    this.currentSort = {column, desc};
    this.store();
    this.render();
  }

  /** The rows the cards or the list show, in their order. */
  shownRows(): number[] {
    return [...this.shown];
  }

  /** The card or the list row of [row]; undefined when it is not rendered. */
  itemOf(row: number): HTMLElement | undefined {
    return this.items.get(row);
  }

  /** How many times the items were rendered: once per filter, mode or sort change (tests). */
  get renderCount(): number {
    return this.renders;
  }

  /** Filters the frame by the search and the extra filter and renders the items. */
  refilter(): void {
    const {searchColumns, tagsColumn, extraFilter} = this.options;
    const frame = this.frame;
    const words = this.search.value.toLowerCase().split(/\s+/).filter((w) => w !== '');
    const isTag = (w: string) => tagsColumn !== undefined && w.length > 1 && w.startsWith('#');
    const tags = words.filter(isTag).map((w) => w.slice(1));
    const text = words.filter((w) => !isTag(w)).join(' ');
    const columns = searchColumns.map((name) => frame.col(name)).filter((c): c is DG.Column => c !== null);
    const tagsCol = tagsColumn === undefined ? null : frame.col(tagsColumn);
    const matches = (i: number) => (text === '' || columns.some((c) => c.getString(i).toLowerCase().includes(text))) &&
      (tags.length === 0 || (tagsCol !== null && tags.every((t) =>
        tagsOf(tagsCol.getString(i)).some((tag) => tag.toLowerCase() === t))));
    frame.filter.init((i) => matches(i) && (extraFilter?.(i) ?? true));
    this.render();
  }

  /** Follows the owner's new frame: its events, its text sorts, the filter and the items. */
  rebind(): void {
    this.unsubscribeFrame();
    this.built.clear();
    this.gridSortKey = null;
    const frame = this.frame;
    this.frameSubs = [
      frame.onCurrentRowChanged.subscribe(() => this.markCurrent()),
      frame.onSelectionChanged.subscribe(() => this.markCurrent()),
    ];
    for (const field of this.options.sortFields) {
      const col = frame.col(field.column);
      if (col?.type === DG.COLUMN_TYPE.STRING)
        col.valueComparer = compareText;
    }
    this.refilter();
  }

  /** Drops the subscriptions; the owner calls it from its view's `detach`. */
  detach(): void {
    this.unsubscribeFrame();
    this.gridSub.unsubscribe();
  }

  private unsubscribeFrame(): void {
    for (const sub of this.frameSubs)
      sub.unsubscribe();
    this.frameSubs = [];
  }

  private modeIcon(mode: GalleryMode): HTMLElement {
    const [icon, tooltip] = MODE_ICONS[mode];
    return ui.iconFA(icon, () => this.setMode(mode), tooltip);
  }

  private static stored(key: string): {[key: string]: unknown} {
    try {
      const parsed: unknown = JSON.parse(grok.userSettings.getValue(Gallery.storage, key) ?? '');
      return isRecord(parsed) ? parsed : {};
    } catch {
      return {};
    }
  }

  private store(): void {
    grok.userSettings.add(Gallery.storage, this.options.key,
      JSON.stringify({mode: this.currentMode, sort: this.currentSort}));
  }

  /** The filtered rows of [frame] by [column]: text case-insensitively, numbers and dates by value, an empty value
   * last; a tie keeps the frame's order. */
  // Not frame.getSortedOrder: it orders a text column by its category codes and ignores valueComparer.
  private static sortedRows(frame: DG.DataFrame, column: string, desc: boolean): number[] {
    const col = frame.getCol(column);
    const rows = Array.from(frame.filter.getSelectedIndexes());
    const isText = col.type === DG.COLUMN_TYPE.STRING;
    const raw = isText ? null : col.getRawData();
    const compare = (a: number, b: number): number => {
      const [noneA, noneB] = [col.isNone(a), col.isNone(b)];
      if (noneA || noneB)
        return noneA === noneB ? 0 : noneA ? 1 : -1;
      const order = raw === null ? compareText(col.getString(a), col.getString(b)) : raw[a] - raw[b];
      return desc ? -order : order;
    };
    return rows.sort((a, b) => compare(a, b) || a - b);
  }

  private countTooltip(): string | null {
    const frame = this.frame;
    const shown = frame.filter.trueCount;
    return shown === frame.rowCount ? null : `${shown} of ${frame.rowCount} ${this.options.noun} shown`;
  }

  private render(): void {
    const frame = this.frame;
    const {grid, emptyHint} = this.options;
    const filtered = frame.filter.trueCount;
    this.count.textContent = filtered === frame.rowCount ? `${frame.rowCount}` : `${filtered} / ${frame.rowCount}`;
    for (const mode of GALLERY_MODES)
      this.modeIcons[mode].classList.toggle('d4-current', mode === this.currentMode);
    this.body.setAttribute('mode', MODE_ATTRIBUTES[this.currentMode]);
    this.renders++;
    this.items = new Map();
    this.shown = [];
    if (frame.rowCount === 0) {
      this.body.replaceChildren(emptyHint);
      return;
    }
    const {column, desc} = this.currentSort;
    if (this.currentMode === 'table') {
      const key = Gallery.sortKey(this.currentSort);
      if (this.gridSortKey !== key) {
        grid.sort([column], [!desc]);
        this.gridSortKey = key;
      }
      if (this.body.firstChild !== grid.root)
        this.body.replaceChildren(grid.root);
      return;
    }
    if (filtered === 0) {
      this.body.replaceChildren(ui.div([ui.iconFA('empty-set')], 'grok-items-view-empty'));
      return;
    }
    const order = Gallery.sortedRows(frame, column, desc);
    this.shown = order.slice(0, MAX_ITEMS);
    const elements = this.shown.map((i) => this.itemFor(i));
    if (order.length > MAX_ITEMS) {
      elements.push(ui.divText(`Showing ${MAX_ITEMS} of ${order.length}. Search to narrow, or switch to the table.`,
        'forge-note'));
    }
    this.body.replaceChildren(...elements);
    this.markCurrent();
  }

  private itemFor(i: number): HTMLElement {
    let element = this.built.get(i);
    if (element === undefined) {
      element = this.renderItem(i);
      this.built.set(i, element);
    }
    this.items.set(i, element);
    return element;
  }

  private renderItem(i: number): HTMLElement {
    const {card, listItem, tooltip, dragObject} = this.options;
    const element = this.currentMode === 'cards' ? ui.div([card(i)], 'grok-gallery-grid-item-wrapper') : listItem(i);
    element.addEventListener('click', (e) => this.click(i, e));
    element.addEventListener('contextmenu', (e) => this.showMenu(i, e));
    ui.tooltip.bind(element, () => tooltip(i));
    if (dragObject !== undefined)
      ui.makeDraggable(element, {getDragObject: () => dragObject(i)});
    return element;
  }

  private static sortKey({column, desc}: GallerySort): string {
    return `${column}\n${desc}`;
  }

  private markCurrent(): void {
    const frame = this.frame;
    const current = frame.currentRowIdx;
    for (const [i, element] of this.items)
      element.classList.toggle('d4-current', i === current || frame.selection.get(i));
  }

  /** A plain click makes [i] the current row alone; Ctrl (Cmd) adds or removes it, the current row counting as
   * selected; Shift selects the shown rows from the current one to it. */
  private click(i: number, e: MouseEvent): void {
    const frame = this.frame;
    const current = frame.currentRowIdx;
    if (e.ctrlKey || e.metaKey) {
      if (frame.selection.trueCount === 0 && current >= 0 && current !== i)
        frame.selection.set(current, true);
      frame.selection.set(i, !frame.selection.get(i));
      // An unpicked current row leaves its place to the last picked one, so the marks are the picks.
      if (i === current && !frame.selection.get(i)) {
        const picked = frame.selection.getSelectedIndexes();
        frame.currentRowIdx = picked.length > 0 ? picked[picked.length - 1] : -1;
      }
    } else if (e.shiftKey && current >= 0) {
      const [from, to] = [this.shown.indexOf(current), this.shown.indexOf(i)];
      if (from < 0 || to < 0)
        return;
      const picked = new Set([...frame.selection.getSelectedIndexes(),
        ...this.shown.slice(Math.min(from, to), Math.max(from, to) + 1)]);
      frame.selection.init((k) => picked.has(k));
    } else {
      frame.selection.setAll(false);
      frame.currentRowIdx = i;
    }
  }

  private showMenu(i: number, e: MouseEvent): void {
    e.preventDefault();
    e.stopPropagation();
    const frame = this.frame;
    if (frame.currentRowIdx !== i && !frame.selection.get(i)) {
      frame.selection.setAll(false);
      frame.currentRowIdx = i;
    }
    const menu = DG.Menu.popup();
    this.options.menu(i, menu);
    menu.show({causedBy: e});
  }

  private showSortMenu(e: MouseEvent): void {
    const {sortFields} = this.options;
    const current = sortFields.find((f) => f.column === this.currentSort.column);
    const menu = DG.Menu.popup();
    // The platform keeps a menu of checked items open after a pick (menu.dart, shouldCloseOnClick).
    const pick = (column: string, desc: boolean) => {
      menu.hide();
      this.setSort(column, desc);
    };
    menu
      .items(sortFields.map((f) => f.caption), (caption) => {
        const field = sortFields.find((f) => f.caption === caption);
        if (field !== undefined)
          pick(field.column, this.currentSort.desc);
      }, {radioGroup: 'order by', isChecked: (caption) => caption === current?.caption})
      .separator()
      .items(DIRECTIONS, (direction) => pick(this.currentSort.column, direction === DIRECTIONS[1]),
        {radioGroup: 'order direction',
          isChecked: (direction) => (direction === DIRECTIONS[1]) === this.currentSort.desc})
      .show({causedBy: e});
  }

  /** A sort by one of the sort fields in the table is the gallery's sort; any other stays the table's own. */
  private followGridSort(): void {
    const {grid, sortFields} = this.options;
    const columns = grid.sortByColumns;
    if (columns.length !== 1 || !sortFields.some((f) => f.column === columns[0].name))
      return;
    const desc = grid.sortTypes[0] === false;
    if (columns[0].name === this.currentSort.column && desc === this.currentSort.desc)
      return;
    this.currentSort = {column: columns[0].name, desc};
    this.gridSortKey = Gallery.sortKey(this.currentSort);
    this.store();
  }
}
