/* `domains.list` — a `VirtualList` over a source's rows: rendered through the table's renderer,
   each row carrying its actions (hover block and context menu alike, permission-filtered per
   row), selection and the source's current row one thing, Enter handing the row to the form
   paired through the source, the next page loaded near the bottom, a row marked deleted kept
   struck through with Restore until the save, and the source's loading, empty and error states
   shown under the rows. */
import {Control} from '../../core/component.js';
import {signal} from '../../core/signals.js';
import type {ReadonlySignal, Signal} from '../../core/signals.js';
import type {ObjectRenderer} from '../../core/object-renderer.js';
import {button, div, divH, span} from '../../core/elements.js';
import {text} from '../../core/text.js';
import {VirtualList} from '../../components/collections/list.js';
import {allowedActions, rowActions} from '../../components/actions/actions.js';
import type {Action} from '../../components/actions/actions.js';
import {loader} from '../../components/display/async-view.js';
import type {DomainSource} from '../../sources/domain-source.js';
import {Rows} from '../../sources/rows-like.js';
import type {RowView} from '../../sources/rows-like.js';
import {ActionRegistry, DomainTable} from './index.js';
import type {DomainAction} from './index.js';
import {DomainErrors} from './errors.js';

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
  /** What an empty result says (default: "No <rows>."). */
  empty?: string;
}

/** A card is the handler's card — a title, a description line and the creation time. */
const HEIGHTS: Record<DomainListMode, number> = {brief: 28, cards: 64};
const NEAR_BOTTOM_ROWS = 5;

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
      render: (row, _index, el) => this._row(row, el),
      contextActions: (row) => this.actionsFor(row),
      onEnter: () => source.activate.value = source.activate.peek() + 1,
      onDelete: (row) => this.actionsFor(row).find((a) => a.name === 'Delete')?.run(),
    }));
    this.list.root.classList.add('u2-domain-list-rows');
    this.list.setItems(source.rows.items as ReadonlySignal<RowView[]>);

    // the list's selection and the source's current row are one thing, written in either direction —
    // seeded from the row that is current already, so the first effect does not clear it
    const current = source.currentRow.peek();
    if (current !== null)
      this.list.selectedIndex.value = source.rows.items.peek().findIndex((r) => r.id === current.id);
    this.effect(() => {
      const row = source.rows.items.peek()[this.list.selectedIndex.value] ?? null;
      if (row?.id !== source.currentRow.peek()?.id)
        source.currentRow.value = row;
    });
    this.effect(() => {
      const row = source.currentRow.value;
      const items = source.rows.items.value;
      const at = row === null ? -1 : items.findIndex((r) => r.id === row.id);
      if (at !== this.list.selectedIndex.peek())
        this.list.selectedIndex.value = at;
    });
    // a list over a table starts on its first row; a draft source is a form's, not a list's
    let seeded = source.isDraft;
    this.effect(() => {
      const items = source.rows.items.value;
      if (seeded || source.state.value !== 'ready' || items.length === 0)
        return;
      seeded = true;
      if (source.currentRow.peek() === null)
        this.list.selectedIndex.value = 0;
    });

    const onScroll = () => {
      const el = this.list.root;
      if (el.scrollTop + el.clientHeight > el.scrollHeight - NEAR_BOTTOM_ROWS * itemHeight)
        void source.loadMore();
    };
    this.list.root.addEventListener('scroll', onScroll);
    this.own(() => this.list.root.removeEventListener('scroll', onScroll));

    const status = div([], 'u2-domain-list-status');
    status.dataset.u2Part = 'status';
    this.effect(() => {
      const state = source.state.value;
      const count = source.rows.items.value.length;
      const q = source.query.value;
      const filtered = typeof q === 'string' ? q.trim() !== '' : q.nodes.length > 0;
      if (state === 'loading' && count === 0)
        status.replaceChildren(loader('Loading…'));
      else if (state === 'error') {
        // a filter the server refuses is a dead end without this: Retry runs it again
        const actions = [button('Retry', () => void source.refresh())];
        if (filtered)
          actions.push(button('Clear filter', () => source.query.value = ''));
        status.replaceChildren(divH([span(DomainErrors.message(source.error.value)), ...actions],
          'u2-domain-list-error'));
      } else if (state === 'ready' && count === 0) {
        status.replaceChildren(span(_options.empty ??
          `No ${source.schema.info.pluralName.toLowerCase() || 'rows'}.`, 'u2-domain-list-empty'));
      } else
        status.replaceChildren();
    });
    // a refusal names one row: marked where it is drawn, and unmarked with the refusal
    this.effect(() => {
      source.problemRow.value;
      source.error.value;
      for (const el of Array.from(this.list.root.querySelectorAll<HTMLElement>('[data-u2-row]')))
        this._markProblem(el);
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
   * A row marked deleted offers Restore alone. */
  actionsFor(row: RowView): Action[] {
    const source = this.source;
    const table = this._table;
    if (row[Rows.STATE] === 'deleted')
      return [{name: 'Restore', icon: 'undo', run: () => source.edit.peek()?.unmarkDeleted(row.id)}];
    const draft = Rows.isDraft(row);
    const actions: Action[] = [];
    if (table !== undefined && !draft)
      actions.push({name: 'Open', icon: 'folder-open', run: () => table.open(row)});
    actions.push({name: 'Delete', icon: 'trash-alt', requires: 'delete',
      run: () => source.edit.peek()?.markDeleted(row.id)});
    if (table !== undefined)
      actions.push(...table.actions.for(row));
    actions.push(...ActionRegistry.bind(this._options.actions ?? [], row));
    return allowedActions(actions, {access: source.access.peek(), row});
  }

  private _row(row: RowView, el: HTMLElement): HTMLElement {
    el.dataset.u2Row = row.id;
    el.classList.toggle('u2-domain-list-deleted', row[Rows.STATE] === 'deleted');
    this._markProblem(el);
    const render = this._options.render;
    const content = render ? render(row) : this._content(row);
    content.classList.add('u2-domain-list-content');
    return div([content, rowActions(this.actionsFor(row))], 'u2-domain-list-item');
  }

  /** The row a refusal names, marked while it stands — the way a pending delete is marked. */
  private _markProblem(el: HTMLElement): void {
    const problem = this.source.problemRow.peek();
    const mine = problem !== null && el.dataset.u2Row === problem;
    el.classList.toggle('u2-domain-list-invalid', mine);
    if (mine) {
      el.setAttribute('aria-invalid', 'true');
      el.title = DomainErrors.message(this.source.error.peek());
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
    return details === '' ? brief() :
      divH([brief(), span(details, 'u2-domain-list-details')], 'u2-domain-list-line');
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
