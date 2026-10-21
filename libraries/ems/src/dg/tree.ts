/* `domains.tree` — a `VirtualTree` over a table the schema declares a hierarchy (`hierarchy:
   true`, one self-referencing ref column): the roots are the rows whose parent column is empty,
   a branch loads its children when it is opened, and both go through `DomainTableLike.query` as
   plain rows — no frame, no writer and no editor per node, which is what makes a tree over a
   large table affordable. `expandTo` walks the row's ancestors (the platform's `GET …/{id}/path`,
   which EXCLUDES the row itself), opens every one of them and then selects the row. What the
   tree is FOR is driving a collection: `selected` feeds a source's query as `<fk> under "<id>"`. */
import {Control, computed, signal, ReadonlySignal, div, span, VirtualTree, allowedActions,
  loader} from '@datagrok-libraries/u2';
import type {TreeNode, Action, DomainRowLike, RowView} from '@datagrok-libraries/u2';
import {DomainBackendError} from '../core/domain-backend.js';
import type {DomainQueryLike} from '../core/domain-backend.js';
import {domains, ActionRegistry, DomainTable} from './facade.js';
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
   * null, which is what a collection bound to the tree reads as "no subtree filter". Without it
   * the only way back to everything is clicking the selected node again. */
  allNode?: string;
  /** How many children one level reads at a time (default 500); a full page ends in a
   * "Load more…" row. */
  pageSize?: number;
}

const PAGE_SIZE = 500;
/** The id of the synthetic {@link DomainTreeOptions.allNode} row — no row of the table carries it. */
const ALL = '~all';
/** The id prefix of a level's "Load more…" row; the parent's id follows, none for the roots. */
const MORE = '~more:';

/** What one level has read: its nodes, and the read of its next page while the last one was full. */
interface Level<T> {
  nodes: TreeNode<T>[];
  next?: () => Promise<void>;
}

export class DomainTree<TRow extends DomainRowLike = DomainRowLike> extends Control {
  readonly tree: VirtualTree<RowView<TRow>>;
  /** The row the tree has selected, null while it has none. */
  readonly selected: ReadonlySignal<RowView<TRow> | null>;
  /** Why the roots, an `expandTo` or a "Load more…" did not arrive — a table that is not a
   * hierarchy included, when the tree was made from an address. A branch reports its own failure
   * on its row. */
  readonly error: ReadonlySignal<string | null>;

  private readonly _error = signal<string | null>(null);
  private readonly _loading = signal(true);
  private readonly _status = div([], 'u2-domain-tree-status');
  private readonly _ready: Promise<DomainTable<TRow>>;
  private _table: DomainTable<TRow> | undefined;
  private _roots: TreeNode<RowView<TRow>>[] = [];
  /** Every level read since the roots, by its parent's id ('' for the roots). */
  private _levels = new Map<string, Level<RowView<TRow>>>();

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
    this.root.classList.add('u2-domain-tree');
    this.root.dataset.u2 = 'domain-tree';
    this._status.dataset.u2Part = 'status';

    const render = _options.render;
    this.tree = this.runInScope(() => new VirtualTree<RowView<TRow>>({
      contextActions: (node) => node.data === undefined ? [] : this.actionsFor(node.data),
      render: (node) => node.data === undefined || render === undefined ?
        span(node.label, node.id.startsWith(MORE) ? 'u2-link' : undefined) : render(node.data),
    }));
    // a node with no row behind it (the "all" row) reads as no selection, which is what makes
    // selecting it clear a collection's subtree filter
    this.selected = computed(() => this.tree.selectedNode.value?.data ?? null);
    this.root.append(this.tree.root, this._status);
    // clicking the selected node again, and Escape, clear the selection — a tree is a navigator,
    // and "nothing selected" (the whole table) is one of its states. Capture, so the list's own
    // click handler does not re-select the row on the way down.
    this.scope.listen(this.root, 'click', (e) => {
      const more = this.tree.nodeForRow(e.target)?.id;
      if (more?.startsWith(MORE)) {
        e.stopPropagation();
        void this._loadMore(more);
        return;
      }
      const target = e.target as Element | null;
      if (target === null || target.closest('.u2-tree-twistie') !== null)
        return;
      const node = this.tree.nodeForRow(target);
      if (node === null || node.id !== this.tree.selectedNode.peek()?.id)
        return;
      e.stopPropagation();
      this.tree.clearSelection();
    }, {capture: true});
    this.scope.listen(this.root, 'keydown', (e) => {
      const key = (e as KeyboardEvent).key;
      const node = this.tree.selectedNode.peek();
      if (key === 'Enter' && node?.id.startsWith(MORE)) {
        e.stopPropagation();
        void this._loadMore(node.id);
        return;
      }
      if (key !== 'Escape' || node === null)
        return;
      e.stopPropagation();
      this.tree.clearSelection();
    }, {capture: true});
    this.effect(() => {
      const problem = this._error.value;
      if (this._loading.value)
        this._status.replaceChildren(loader('Loading…'));
      else if (problem !== null)
        this._status.replaceChildren(span(problem, 'u2-domain-tree-error'));
      else
        this._status.replaceChildren();
    });
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

  /** Reads the roots again; the children of a branch are re-read when it is opened again. */
  async refresh(): Promise<void> {
    this._loading.value = true;
    this._error.value = null;
    try {
      const table = await this._ready;
      DomainTree.requireHierarchy(table);
      this._table = table;
      this._levels = new Map();
      this._roots = await this._children(null);
      this.tree.setRoots(this._roots);
      const reveal = this._options.expandTo;
      if (reveal !== undefined)
        await this.expandTo(reveal);
    } catch (e) {
      this._error.value = DomainErrors.message(e);
    } finally {
      this._loading.value = false;
    }
  }

  /** Opens every ancestor of the row and selects it, reading the further pages of a level the
   * next step is not in yet. The chain is the backend's answer, so one running through a node the
   * caller cannot see stops there, as the server's path does. */
  async expandTo(id: string): Promise<void> {
    const handle = await this._ready;
    const table = handle.table;
    if (table.ancestors === undefined) {
      throw new DomainBackendError('unsupported',
        `${handle.address}: the backend cannot answer a row's ancestors`);
    }
    const path = await table.ancestors(id);
    const ids = [...path.map((a) => a.id), id];
    for (;;) {
      await this.tree.expandPath(ids);
      if (this.tree.selectedNode.peek()?.data?.id === id)
        return;
      const next = this._nextOnPath(ids);
      if (next === undefined)
        break;
      await next();
    }
    // every page of the level the walk stops in is read, so what is left is exactly what a chain
    // truncated at an ancestor the caller cannot see looks like from here
    this._error.value = `${id} is not reachable from here — an ancestor is not visible to you`;
  }

  /** Every action that applies to the node's row and that the caller may run on it: Open, the
   * table's registry, the tree's own. */
  actionsFor(row: RowView<TRow>): Action[] {
    const table = this._table;
    if (table === undefined)
      return [];
    const actions: Action[] = [{name: 'Open', icon: 'folder-open', run: () => table.open(row)}];
    actions.push(...table.actions.for(row));
    actions.push(...ActionRegistry.bind(this._options.actions ?? [], row));
    return allowedActions(actions, {access: table.access, row});
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

  /** One level: the rows whose parent column is `parent` — null for the roots, which is what the
   * `is empty` operator compiles to — read one page past `loaded`; a full page ends in a
   * "Load more…" row. */
  private async _children(parent: TreeNode<RowView<TRow>> | null,
    loaded: TreeNode<RowView<TRow>>[] = []): Promise<TreeNode<RowView<TRow>>[]> {
    const table = this._table!;
    const levels = this._levels;
    const limit = this._options.pageSize ?? PAGE_SIZE;
    const spec: DomainQueryLike = {filter: {property: table.info.parentColumn!, operator: '=',
      value: parent?.id ?? null}, limit, offset: loaded.length, withAccess: true};
    const name = table.info.nameColumn;
    // the id breaks ties between equal names, so consecutive pages neither repeat nor skip a row
    if (name !== null)
      spec.sort = `${name},id`;
    const rows = await table.table.query(spec);
    const nodes = [...loaded, ...rows.map((row) => this._node(row as unknown as RowView<TRow>))];
    const key = parent?.id ?? '';
    const full = rows.length === limit;
    levels.set(key, {nodes, next: full ? () => this._next(levels, parent, nodes) : undefined});
    const level = full ? [...nodes, {id: MORE + key, label: 'Load more…'}] : nodes;
    const all = this._options.allNode;
    return parent === null && all !== undefined ? [{id: ALL, label: all, tooltip: all}, ...level] : level;
  }

  /** Reads the level's next page and shows it in place of the level's "Load more…" row. */
  private async _next(levels: Map<string, Level<RowView<TRow>>>,
    parent: TreeNode<RowView<TRow>> | null, nodes: TreeNode<RowView<TRow>>[]): Promise<void> {
    const next = await this._children(parent, nodes);
    if (levels !== this._levels)
      return;
    if (parent === null)
      this._roots = next;
    else
      parent.children = () => Promise.resolve(next);
    this.tree.setRoots(this._roots);
  }

  /** What a "Load more…" click or Enter runs, with the status line saying it is busy or failed. */
  private async _loadMore(id: string): Promise<void> {
    const next = this._levels.get(id.slice(MORE.length))?.next;
    if (next === undefined || this._loading.peek())
      return;
    this._loading.value = true;
    this._error.value = null;
    try {
      await next();
    } catch (e) {
      this._error.value = DomainErrors.message(e);
    } finally {
      this._loading.value = false;
    }
  }

  /** The next-page read of the level a walk down `ids` stops in, while that level has one. */
  private _nextOnPath(ids: string[]): (() => Promise<void>) | undefined {
    let key = '';
    for (const step of ids) {
      const level = this._levels.get(key);
      if (level === undefined || !level.nodes.some((n) => n.id === step))
        return level?.next;
      key = step;
    }
    return undefined;
  }

  private _node(row: RowView<TRow>): TreeNode<RowView<TRow>> {
    const label = this._table!.renderer.caption(row);
    const node: TreeNode<RowView<TRow>> = {id: row.id, label, tooltip: label, data: row};
    node.children = () => this._children(node);
    return node;
  }
}
