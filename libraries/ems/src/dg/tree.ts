import {Control, computed, signal, ReadonlySignal, div, span, VirtualTree, loader} from '@datagrok-libraries/u2';
import type {TreeNode, Action, DomainRowLike, RowView} from '@datagrok-libraries/u2';
import {DomainBackendError} from '../core/domain-backend.js';
import type {DomainQueryLike} from '../core/domain-backend.js';
import {domains, DomainTable} from './facade.js';
import type {DomainAction} from './facade.js';
import {DomainErrors} from './errors.js';

export interface DomainTreeOptions<TRow extends DomainRowLike = DomainRowLike> {
  /** A row to reveal once the roots are in: every ancestor opened, then the row selected. */
  expandTo?: string;
  /** How a node is labelled; the table's renderer caption by default. */
  render?: (row: RowView<TRow>) => HTMLElement;
  /** Actions on top of Open and the table's own registry, per node. */
  actions?: DomainAction<TRow>[];
  /** A first row standing for the whole table ('All locations'): selecting it is `selected`
   * null, which is what a collection bound to the tree reads as "no subtree filter"; it is selected
   * whenever no node is. Without it the only way back to everything is clicking the selected node again. */
  allNode?: string;
  /** How many children one request reads (default 500); a larger level is read in pages. */
  pageSize?: number;
  /** A live source the tree rides: its probe follows the tree's table too ({@link DomainSource.watch}),
   * and a peer's rename, insert, move or delete there — or a save in this browser — reads the open
   * tree again. Without it the tree is read once, and again only on {@link DomainTree.refresh}. While the
   * source or its session holds pending changes, a selected node a peer deleted stays `selected` (marked
   * deleted) instead of handing the source to its parent; it moves on once they are saved or discarded. */
  follow?: {watch(address: string, moved: () => void): () => void, isDirty?: ReadonlySignal<boolean>,
    session?: {isDirty: ReadonlySignal<boolean>}};
}

const PAGE_SIZE = 500;
/** The id of the synthetic {@link DomainTreeOptions.allNode} row — no row of the table carries it. */
const ALL = '~all';
/** An id no node carries: `expandPath` to it reads every branch on the way and selects nothing. */
const NONE = '~none';

/** A `VirtualTree` over a hierarchy table, a branch loaded when it is opened. Nodes are plain rows: no frame, no
 * writer and no editor per node, which is what makes a tree over a large table affordable. `selected` feeds a
 * source's query as `<fk> under "<id>"`. */
export class DomainTree<TRow extends DomainRowLike = DomainRowLike> extends Control {
  readonly tree: VirtualTree<RowView<TRow>>;
  /** The row the tree has selected, null while it has none. */
  readonly selected: ReadonlySignal<RowView<TRow> | null>;
  /** Why the roots or an `expandTo` did not arrive — a table that is not a hierarchy included,
   * when the tree was made from an address. A branch reports its own failure on its row. */
  readonly error: ReadonlySignal<string | null>;
  /** What a refresh did to the selection — "Lab 101 was deleted" — until another node is selected. */
  readonly notice: ReadonlySignal<string | null>;

  private readonly _error = signal<string | null>(null);
  private readonly _notice = signal<string | null>(null);
  // the selection a refresh is re-finding: `selected` stays on it while the rows below shift
  private readonly _held = signal<{row: RowView<TRow> | null} | null>(null);
  private _read = false;
  private _refreshing: Promise<void> | undefined;
  private _again = false;
  private _noticeAt: string | null = null;
  // the selected row a peer deleted while `follow` had pending changes: still `selected`, so nothing rebinds
  private readonly _kept = signal<RowView<TRow> | null>(null);
  private readonly _loading = signal(true);
  private readonly _status = div([], 'u2-domain-tree-status');
  private readonly _ready: Promise<DomainTable<TRow>>;
  private _table: DomainTable<TRow> | undefined;

  /** `target` is the table handle, or its `'<schema>.<table>'` address for a spec-built tree —
   * a handle is checked here and a bad one refused at once, an address once it has resolved. */
  constructor(target: DomainTable<TRow> | string, private readonly _options: DomainTreeOptions<TRow> = {}) {
    super();
    if (typeof target === 'string')
      this._ready = domains.table<TRow>(target);
    else {
      DomainTree.requireHierarchy(target);
      this._table = target;
      this._ready = Promise.resolve(target);
    }
    this.error = this._error;
    this.notice = this._notice;
    this.root.classList.add('u2-domain-tree');
    this.root.dataset.u2 = 'domain-tree';
    this._status.dataset.u2Part = 'status';

    const render = _options.render;
    this.tree = this.build(() => new VirtualTree<RowView<TRow>>({
      contextActions: (node) => node.data === undefined ? [] : this.actionsFor(node.data),
      render: (node) => node.id === ALL ? span(node.label, 'u2-domain-tree-all') :
        render === undefined ? span(node.label) : render(node.data!),
    }));
    // a node with no row behind it (the "all" row) reads as no selection, which is what makes
    // selecting it clear a collection's subtree filter
    this.selected = computed(() => {
      const held = this._held.value;
      if (held !== null)
        return held.row;
      return this._kept.value ?? this.tree.selectedNode.value?.data ?? null;
    });
    this.root.append(this.tree.root, this._status);
    // clicking the selected node again, and Escape, clear the selection — a tree is a navigator,
    // and "nothing selected" (the whole table) is one of its states. Capture, so the list's own
    // click handler does not re-select the row on the way down.
    this.scope.listen(this.root, 'click', (e) => {
      const target = e.target as Element | null;
      if (target === null || target.closest('.u2-tree-twistie') !== null)
        return;
      const node = this.tree.nodeForRow(target);
      if (node === null || node.id !== this.tree.selectedNode.peek()?.id)
        return;
      e.stopPropagation();
      void this._clear();
    }, {capture: true});
    this.scope.listen(this.root, 'keydown', (e) => {
      if ((e as KeyboardEvent).key !== 'Escape' || this.tree.selectedNode.peek() === null)
        return;
      e.stopPropagation();
      void this._clear();
    }, {capture: true});
    this.effect(() => {
      const problem = this._error.value;
      const notice = this._notice.value;
      if (this._loading.value)
        this._status.replaceChildren(loader('Loading…'));
      else if (problem !== null)
        this._status.replaceChildren(span(problem, 'u2-domain-tree-error'));
      else if (notice !== null)
        this._status.replaceChildren(span(notice, 'u2-domain-tree-notice'));
      else
        this._status.replaceChildren();
    });
    // the notice is about the node a refresh picked instead; picking another one retires it
    this.effect(() => {
      const id = this.selected.value?.id ?? null;
      if (this._held.peek() === null && this._notice.peek() !== null && id !== this._noticeAt)
        this._notice.value = null;
    });
    const follow = _options.follow;
    // a kept row leaves once the pending changes are saved or discarded — or at once when the user picks a node
    this.effect(() => {
      const kept = this._kept.value;
      if (kept === null)
        return;
      if (this.tree.selectedNode.value !== null && this._held.peek() === null)
        this._kept.value = null;
      else if (!this._pending(true))
        void this._leave(kept);
    });
    if (follow !== undefined) {
      void this._ready.then((table) => {
        if (!this.scope.isDisposed)
          this.own(follow.watch(table.address, () => void this.refresh()));
      }, () => {});
    }
    void this.refresh();
  }

  /** The table handle, once resolved. */
  get table(): DomainTable<TRow> | undefined {
    return this._table;
  }

  /** The self-referencing column the tree walks, once the table is known. */
  get parentColumn(): string | null {
    return this._table?.info.parentColumn ?? null;
  }

  /** Reads the roots again, and every branch that is open, keeping the selection on its row: a renamed
   * or moved row stays selected where it now is, a deleted one gives way to its parent (or to
   * nothing) with a {@link notice}. Calls made while one runs are one more read after it. */
  refresh(): Promise<void> {
    if (this._refreshing !== undefined) {
      this._again = true;
      return this._refreshing;
    }
    const run = async (): Promise<void> => {
      do {
        this._again = false;
        await this._refresh();
      } while (this._again && !this.scope.isDisposed);
    };
    return this._refreshing = run().finally(() => this._refreshing = undefined);
  }

  private async _refresh(): Promise<void> {
    const first = !this._read;
    const was = this.tree.selectedNode.peek();
    const row = was?.data ?? null;
    // the selection is a list index: rows that land above it would move it to another node; a kept row
    // has no node, and is what `selected` is
    if (!first)
      this._held.value = {row: this._kept.peek() ?? row};
    if (first)
      this._loading.value = true;
    this._error.value = null;
    try {
      const table = await this._ready;
      DomainTree.requireHierarchy(table);
      this._table = table;
      const roots = await this._children(null);
      if (this.scope.isDisposed)
        return;
      this.tree.setRoots(roots);
      this._read = true;
      const reveal = this._options.expandTo;
      if (first && reveal !== undefined)
        await this.expandTo(reveal);
      else if (this._kept.peek() !== null)
        await this._reopen();
      else if (first || was === null)
        await this._clear();
      else {
        await this._reopen();
        // the index still points where the row was: placed afresh, so a row that left is not found "there"
        this.tree.clearSelection();
        if (row === null)
          await this.tree.expandPath([was.id]);
        else
          await this._reselect(row);
      }
    } catch (e) {
      this._error.value = DomainErrors.message(e);
    } finally {
      this._held.value = null;
      this._loading.value = false;
    }
  }

  // every open branch read again under the new roots before the selection is placed, parents first, so no
  // answer lands above it afterwards: `expandPath` to an id no node has awaits the tree's own reads. One read
  // per open branch, one after another, per change the probe sees; changes during a refresh make one more
  private async _reopen(): Promise<void> {
    const done = new Set<string>();
    for (;;) {
      const items = this.tree.getItems({all: true}).list;
      const chain: string[] = [];
      let next: string[] | undefined;
      for (const item of items) {
        chain.length = (item.level ?? 1) - 1;
        chain.push(item.key);
        if (item.expanded === true && item.key !== ALL && !done.has(item.key)) {
          next = [...chain];
          break;
        }
      }
      if (next === undefined || this.scope.isDisposed)
        return;
      done.add(next[next.length - 1]);
      await this.tree.expandPath([...next, NONE]);
    }
  }

  // the row selected before the refresh, found again: where it now is, or its parent when a peer deleted it
  // (two to four reads, on top of `_reopen`'s one per open branch: a tree's refresh costs what it shows)
  private async _reselect(row: RowView<TRow>): Promise<void> {
    const table = this._table!;
    const [now] = await table.table.query({filter: {property: 'id', operator: '=', value: row.id}, limit: 1});
    if (now !== undefined) {
      const ancestors = table.table.ancestors;
      const path = ancestors === undefined ? [] : (await ancestors.call(table.table, row.id)).map((a) => a.id);
      await this.tree.expandPath([...path, row.id]);
      if (this.tree.selectedNode.peek()?.id === row.id)
        return;
    }
    const caption = table.renderer.caption(row);
    if (now === undefined && this._pending(false)) {
      // as a list keeps its rows over pending edits and says "Data changed": what is bound to the row stays
      this.tree.clearSelection();
      this._kept.value = row;
      this._noticeAt = row.id;
      this._notice.value = `${caption} was deleted — save or discard your changes to leave it`;
      return;
    }
    await this._toParent(row,
      now === undefined ? `${caption} was deleted` : `${caption} is no longer reachable from here`);
  }

  private _pending(track: boolean): boolean {
    const follow = this._options.follow;
    const read = (s: ReadonlySignal<boolean> | undefined): boolean =>
      s === undefined ? false : track ? s.value : s.peek();
    return read(follow?.isDirty) || read(follow?.session?.isDirty);
  }

  private async _leave(row: RowView<TRow>): Promise<void> {
    this._held.value = {row};
    this._kept.value = null;
    try {
      await this._toParent(row, `${this._table!.renderer.caption(row)} was deleted`);
    } finally {
      this._held.value = null;
    }
  }

  private async _toParent(row: RowView<TRow>, gone: string): Promise<void> {
    const table = this._table!;
    const parent = (row as Record<string, unknown>)[table.info.parentColumn!];
    if (typeof parent === 'string' && parent !== '') {
      const [up] = await table.table.query({filter: {property: 'id', operator: '=', value: parent}, limit: 1});
      if (up !== undefined) {
        const path = table.table.ancestors === undefined ? [] :
          (await table.table.ancestors(parent)).map((a) => a.id);
        await this.tree.expandPath([...path, parent]);
      }
    }
    if (this.tree.selectedNode.peek()?.id !== parent)
      await this._clear();
    this._noticeAt = this.tree.selectedNode.peek()?.data?.id ?? null;
    this._notice.value = gone;
  }

  /** Opens every ancestor of the row and selects it. The chain is the backend's answer, so one
   * running through a node the caller cannot see stops there, as the server's path does. */
  async expandTo(id: string): Promise<void> {
    const handle = await this._ready;
    const table = handle.table;
    if (table.ancestors === undefined) {
      throw new DomainBackendError('unsupported',
        `${handle.address}: the backend cannot answer a row's ancestors`);
    }
    const path = await table.ancestors(id);
    await this.tree.expandPath([...path.map((a) => a.id), id]);
    // `expandPath` walks down from the roots and stops where a step is missing — which is exactly
    // what a chain truncated at an ancestor the caller cannot see looks like from here
    if (this.tree.selectedNode.peek()?.data?.id !== id)
      this._error.value = `${id} is not reachable from here — an ancestor is not visible to you`;
  }

  /** Every action that applies to the node's row and that the caller may run on it: Open, the
   * table's registry, the tree's own. */
  actionsFor(row: RowView<TRow>): Action[] {
    return this._table?.rowActions(row, this._options.actions) ?? [];
  }

  /** The named refusal for a table the registry does not declare a hierarchy: there is no column
   * for a tree to walk, whichever backend holds the rows. */
  static requireHierarchy<T extends DomainRowLike>(table: DomainTable<T>): void {
    const parent = table.info.parentColumn;
    if (table.info.hierarchy !== true || parent === null || parent === undefined) {
      throw new DomainBackendError('filter',
        `${table.address} is not a hierarchy table — declare "hierarchy": true and one ref column to itself`);
    }
  }

  private async _clear(): Promise<void> {
    if (this._options.allNode === undefined)
      this.tree.clearSelection();
    else
      await this.tree.expandPath([ALL]);
  }

  // `id` null is the roots: what the `is empty` operator compiles to
  private async _children(id: string | null): Promise<TreeNode<RowView<TRow>>[]> {
    const table = this._table!;
    const limit = this._options.pageSize ?? PAGE_SIZE;
    const spec: DomainQueryLike = {filter: {property: table.info.parentColumn!, operator: '=', value: id},
      limit, withAccess: true};
    const name = table.info.nameColumn;
    if (name !== null)
      spec.sort = name;
    const rows: Record<string, unknown>[] = [];
    for (;;) {
      const page = await table.table.query({...spec, offset: rows.length});
      rows.push(...page);
      if (page.length < limit)
        break;
    }
    const nodes = rows.map((row) => this._node(row as unknown as RowView<TRow>));
    const all = this._options.allNode;
    return id === null && all !== undefined ? [{id: ALL, label: all, tooltip: all}, ...nodes] : nodes;
  }

  private _node(row: RowView<TRow>): TreeNode<RowView<TRow>> {
    const label = this._table!.renderer.caption(row);
    return {id: row.id, label, tooltip: label, data: row, children: () => this._children(row.id)};
  }
}
