/* `domains.list` — a `VirtualList` over a source's rows: rendered through the table's renderer,
   each row carrying its actions (hover block and context menu alike, permission-filtered per
   row), selection and the source's current row one thing, Enter handing the row to the form
   paired through the source, the next page loaded near the bottom, a row marked deleted kept
   struck through with Restore until the save, and the source's loading, empty and error states
   shown under the rows. */
import {Control, signal, button, iconButton, div, divH, span, timestamp, text, VirtualList, allowedActions, rowActions,
  loader, emptyState, Rows} from '@datagrok-libraries/u2';
import type {ReadonlySignal, Signal, ObjectRenderer, Action, RowView} from '@datagrok-libraries/u2';
import type {DomainSource} from '../core/domain-source.js';
import {ActionRegistry, DomainTable} from './facade.js';
import type {DomainAction} from './facade.js';
import {DomainErrors} from './errors.js';
import {collectionSelection} from '@datagrok-libraries/u2/src/dg/index.js';

export type DomainListMode = 'cards' | 'brief';

export interface DomainListOptions {
  /** `brief`: one line per row (the renderer's list item); `cards`: the renderer's card. */
  mode?: DomainListMode;
  itemHeight?: number;
  /** Replaces the renderer for the row content. */
  render?: (row: RowView) => HTMLElement;
  /** Actions on top of Open, Delete and the table's own registry. */
  actions?: DomainAction[];
  renderer?: ObjectRenderer<RowView>;
  /** What an empty result says (default: "No <rows> yet"). */
  empty?: string;
  /** What an empty, unfiltered list offers as its primary action. */
  create?: () => void;
}

/** A card is the handler's card — a title, a description line and the creation time. */
const HEIGHTS: Record<DomainListMode, number> = {brief: 28, cards: 64};
const NEAR_BOTTOM_ROWS = 5;
const PICK_KEYS = new Set(['ArrowDown', 'ArrowUp', 'Home', 'End', 'PageDown', 'PageUp', 'Enter', ' ']);

export class DomainList extends Control {
  readonly list: VirtualList<RowView>;
  readonly mode: DomainListMode;
  /** The rows are not the ones the filter box reads as — set while that filter is invalid, so
   * they are shown dimmed instead of as the answer (`DomainApp` wires it to `DomainFilters`). */
  readonly stale: Signal<boolean> = signal(false);

  private readonly _table: DomainTable | undefined;
  private _default: ObjectRenderer<RowView> | undefined;

  constructor(readonly source: DomainSource, private readonly _options: DomainListOptions = {}) {
    super();
    this.mode = _options.mode ?? 'brief';
    this._table = DomainTable.of(source);
    const itemHeight = _options.itemHeight ?? HEIGHTS[this.mode];
    this.root.classList.add('u2-domain-list', `u2-domain-list-${this.mode}`);
    this.root.dataset.u2 = 'domain-list';

    this.list = this.runInScope(() => new VirtualList<RowView>({
      itemHeight,
      keyOf: (row) => row.id,
      // a row is a live proxy over the frame: an edit changes it without a new object
      liveItems: true,
      render: (row, _index, el) => this._row(row, el),
      // with several rows picked, the selection's actions come flat, the clicked row's in a submenu
      contextActions: (row) => {
        const selected = this._selected();
        return selected.length < 2 ? this.actionsFor(row) : this.selectionActions(selected);
      },
      contextGroups: (row) => this._selected().length < 2 ? [] :
        [{label: this.renderer.caption(row), actions: this.actionsFor(row)}],
      onEnter: () => source.activate.value = source.activate.peek() + 1,
      onDelete: (row) => this.actionsFor(row).find((a) => a.name === 'Delete')?.run(),
    }));
    this.list.root.classList.add('u2-domain-list-rows');
    this.list.setItems(source.rows.items as ReadonlySignal<RowView[]>);

    collectionSelection(this, source, this.list);
    // a list over a table starts on its first row; a draft source is a form's, not a list's. The
    // seeded row is drawn current only once the user picks it, and the first Down lands on it
    let seeded = source.isDraft;
    this.effect(() => {
      const items = source.rows.items.value;
      if (seeded || source.state.value !== 'ready' || items.length === 0)
        return;
      seeded = true;
      if (source.currentRow.peek() !== null)
        return;
      this.list.selectedIndex.value = 0;
      this.root.classList.add('u2-domain-list-seeded');
    });
    const onClick = (e: Event) => {
      if ((e.target as Element).closest('.u2-list-row') !== null)
        this.root.classList.remove('u2-domain-list-seeded');
    };
    const onKey = (e: KeyboardEvent) => {
      if (!PICK_KEYS.has(e.key) || !this.root.classList.contains('u2-domain-list-seeded'))
        return;
      this.root.classList.remove('u2-domain-list-seeded');
      if (e.key === 'ArrowDown' && this.list.selectedIndex.peek() === 0) {
        e.preventDefault();
        e.stopImmediatePropagation();
      }
    };
    this.list.root.addEventListener('click', onClick);
    this.list.root.addEventListener('keydown', onKey, true);
    this.own(() => {
      this.list.root.removeEventListener('click', onClick);
      this.list.root.removeEventListener('keydown', onKey, true);
    });

    pageOnScroll(this, this.list, source);
    // a line the row cuts reads in full on hover
    const onOver = (e: Event) => {
      const el = e.target as HTMLElement;
      if (el.closest('.u2-domain-list-content') !== null && el.title === '' && el.scrollWidth > el.clientWidth)
        el.title = el.textContent ?? '';
    };
    this.list.root.addEventListener('mouseover', onOver);
    this.own(() => this.list.root.removeEventListener('mouseover', onOver));

    const status = div([], 'u2-domain-list-status');
    status.dataset.u2Part = 'status';
    this.effect(() => {
      const state = source.state.value;
      const count = source.rows.items.value.length;
      const q = source.query.value;
      const filtered = typeof q === 'string' ? q.trim() !== '' : q.nodes.length > 0;
      const clearFilter = () => button('Clear filter', () => source.query.value = '', {primary: true});
      if (state === 'loading' && count === 0)
        status.replaceChildren(loader('Loading…'));
      else if (state === 'error') {
        // a refused filter fails again on Retry: clearing it is the way out
        const failed = emptyState(DomainErrors.message(source.error.value, source), {icon: 'exclamation-circle',
          error: true, actions: [filtered ? clearFilter() : button('Retry', () => void source.refresh())]});
        failed.classList.add('u2-domain-list-error');
        status.replaceChildren(failed);
      } else if (state === 'ready' && count === 0)
        status.replaceChildren(this._empty(filtered, clearFilter));
      else
        status.replaceChildren();
      // the rows have nothing to scroll: the message takes their place
      this.root.classList.toggle('u2-domain-list-blank', count === 0 && (state === 'ready' || state === 'error'));
    });
    // a refusal names one row: marked where it is drawn, brought into view, and unmarked with the refusal
    this.effect(() => {
      const problem = source.problemRow.value;
      source.error.value;
      for (const el of Array.from(this.list.root.querySelectorAll<HTMLElement>('[data-u2-row]')))
        this._markProblem(el);
      const at = problem === null ? -1 : source.rows.items.peek().findIndex((r) => r.id === problem);
      if (at >= 0)
        this.list.scrollToIndex(at);
    });
    this.effect(() => this.root.classList.toggle('u2-domain-list-stale', this.stale.value));
    this.root.append(this.list.root, status);
  }

  /** The renderer in force: the option, else the table's (live — an app may replace it), else
   * the schema-driven default. */
  get renderer(): ObjectRenderer<RowView> {
    return this._options.renderer ?? this._table?.renderer ??
      (this._default ??= DomainTable.schemaRenderer(() => this.source.schema));
  }

  /** Every action that applies to `row` and that the caller may run on it (the row's own access):
   * Open on a saved row, Delete under the delete capability, the table's registry, the list's own.
   * A row marked deleted offers Restore alone — and so does one already in the trash, whose
   * Restore is staged into the session and lands with the Save (the Delete grant, per row). */
  actionsFor(row: RowView): Action[] {
    const source = this.source;
    const table = this._table;
    if (row[Rows.STATE] === 'restored')
      return [{name: 'Undo restore', icon: 'undo', run: () => source.unstageRestore([row.id])}];
    if (row[Rows.STATE] === 'deleted')
      return [{name: 'Restore', icon: 'undo', run: () => source.edit.peek()?.unmarkDeleted(row.id)}];
    if (Rows.isDeleted(row)) {
      return allowedActions([{name: 'Restore', icon: 'undo', requires: 'delete',
        run: () => source.stageRestore([row.id])}], {access: source.access.peek(), row});
    }
    const draft = Rows.isDraft(row);
    const actions: Action[] = [];
    if (table !== undefined && !draft)
      actions.push({name: 'Open', icon: 'external-link', shortcut: 'Enter', run: () => table.open(row)});
    if (table !== undefined)
      actions.push(...table.actions.for(row));
    actions.push(...ActionRegistry.bind(this._options.actions ?? [], row));
    actions.push({name: 'Delete', icon: 'trash-alt', requires: 'delete', shortcut: 'Del', destructive: true,
      run: () => source.edit.peek()?.markDeleted(row.id)});
    return allowedActions(actions, {access: source.access.peek(), row});
  }

  /** The context menu of a row that is part of a multi-selection: every action that applies to ALL
   * of the selected rows, marked `<N> <plural>` and run over the whole selection. An action is
   * offered there when every selected row offers it and none of them has it disabled; it runs row
   * by row in selection order, and the edits land in the session's next Save as one row's do. */
  selectionActions(rows: readonly RowView[]): Action[] {
    if (rows.length === 0)
      return [];
    // Open takes one page: everything else a row offers means the same over many
    const per = rows.map((row) => this.actionsFor(row).filter((a) => a.name !== 'Open'));
    const shared = per[0].filter((action) => per.every((actions) => {
      const same = actions.find((a) => a.name === action.name);
      return same !== undefined && same.enabled !== false;
    }));
    const detail = `${rows.length} ${DomainTable.noun(this.source.schema.info)}`;
    return shared.map((action) => ({name: action.name, icon: action.icon, detail, destructive: action.destructive,
      run: () => {
        for (const actions of per)
          actions.find((a) => a.name === action.name)!.run();
      }}));
  }

  /** An empty result says what it is empty of, and offers the way out of it. */
  private _empty(filtered: boolean, clearFilter: () => HTMLElement): HTMLElement {
    const source = this.source;
    const search = source.search.value;
    const info = source.schema.info;
    const what = DomainTable.noun(info);
    const trash = source.deleted.value === 'only';
    const create = this._options.create;
    const actions = search !== '' ? [button('Clear search', () => source.search.value = '')] :
      filtered ? [clearFilter()] :
        create !== undefined && !trash && source.access.value.can('insert') ?
          [button(`New ${DomainTable.noun(info, true)}`, () => create(), {primary: true})] : [];
    const said = this._options.empty ?? (search !== '' ? `No ${what} match "${search}"` :
      filtered ? `No ${what} match the filter` : trash ? 'Trash is empty' : `No ${what} yet`);
    const state = emptyState(said, {icon: search !== '' ? 'search' : trash ? 'trash-alt' : 'inbox', actions});
    state.classList.add('u2-domain-list-nothing');
    return state;
  }

  /** The rows the multi-selection holds, in selection order. */
  private _selected(): RowView[] {
    const items = this.source.rows.items.peek();
    return [...this.list.selectedIndices.peek()].map((i) => items[i]).filter((row) => row !== undefined);
  }

  private _row(row: RowView, el: HTMLElement): HTMLElement {
    el.dataset.u2Row = row.id;
    const state = row[Rows.STATE];
    el.classList.toggle('u2-domain-list-deleted', state === 'deleted' || Rows.isDeleted(row));
    el.classList.toggle('u2-domain-list-restored', state === 'restored');
    el.classList.toggle('u2-domain-list-pending', state === 'modified' || state === 'new' || state === 'deleted');
    this._markProblem(el);
    const render = this._options.render;
    const content = render ? render(row) : this._content(row);
    content.classList.add('u2-domain-list-content');
    const actions = this.actionsFor(row);
    const undo = state === 'deleted' ? actions.find((a) => a.name === 'Restore') : undefined;
    const tail = undo === undefined ? rowActions(actions) :
      div([iconButton('undo', () => undo.run(), {tooltip: 'Undo delete'})], 'u2-domain-list-undo');
    return div([content, tail], 'u2-domain-list-item');
  }

  /** The row a refusal names, marked while it stands — the way a pending delete is marked. */
  private _markProblem(el: HTMLElement): void {
    const problem = this.source.problemRow.peek();
    const mine = problem !== null && el.dataset.u2Row === problem;
    el.classList.toggle('u2-domain-list-invalid', mine);
    if (mine) {
      el.setAttribute('aria-invalid', 'true');
      el.title = DomainErrors.message(this.source.error.peek(), this.source);
    } else {
      el.removeAttribute('aria-invalid');
      el.removeAttribute('title');
    }
  }

  private _content(row: RowView): HTMLElement {
    const renderer = this.renderer;
    const brief = () => renderer.listItem?.(row) ?? span(renderer.caption(row), 'u2-domain-list-name');
    if (this.mode === 'cards')
      return renderer.card?.(row) ?? brief();
    const details = this._details(row);
    const line = [brief(), ...details === '' ? [] : [span(details, 'u2-domain-list-details')]];
    // the trash says WHEN, in the reader's own zone with the full moment as its title
    if (Rows.isDeleted(row) && this.source.readOnly.peek()) {
      const when = timestamp(row.updated_on as string, 'u2-domain-list-when');
      when.textContent = `Deleted ${when.textContent}`;
      line.push(when);
    }
    return line.length === 1 ? line[0] : divH(line, 'u2-domain-list-line');
  }

  /** What a bare table shows beside the name: the columns a search matches, so a hit says why it
   * is one. A table whose app brought its own renderer says what that renderer says. */
  private _details(row: RowView): string {
    const info = this.source.schema.info;
    if (this._options.renderer !== undefined || this._table === undefined ||
        this._table.renderer !== this._table.defaultRenderer)
      return '';
    return info.searchableColumns.filter((c) => c !== info.nameColumn).map((c) => text(row[c]))
      .filter((v) => v !== '').join(' · ');
  }
}

/** Loads the source's next page as `rows` scrolls near its bottom, while `owner` lives. */
export function pageOnScroll(owner: Control, rows: {root: HTMLElement, itemHeight: number},
  source: DomainSource): void {
  const scroller = rows.root;
  const onScroll = () => {
    if (scroller.scrollTop + scroller.clientHeight > scroller.scrollHeight - NEAR_BOTTOM_ROWS * rows.itemHeight)
      void source.loadMore();
  };
  scroller.addEventListener('scroll', onScroll);
  owner.own(() => scroller.removeEventListener('scroll', onScroll));
}
