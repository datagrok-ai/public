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
   * "Load more…" row. The server clamps a page to 10 000 rows, so a larger size reads as the
   * level's last page, and an external table cannot page past its first 10 000 rows. */
  pageSize?: number;
}

const PAGE_SIZE = 500;
/** The id of the synthetic {@link DomainTreeOptions.allNode} row — no row of the table carries it. */
const ALL = '~all';
/** The id prefix of a level's "Load more…" row; the parent's id follows, none for the roots. */
const MORE = '~more:';

/** What one level has read: its parent and the ids down to it, its nodes, and whether its last
 * page was full and brought something new. */
interface Level<T> {
  parent: TreeNode<T> | null;
  path: string[];
  nodes: TreeNode<T>[];
  full: boolean;
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
  private readonly _selected = signal<RowView<TRow> | null>(null);
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
    // selecting it clear a collection's subtree filter; the "Load more…" row is not a place, so the
    // row selected before it stands, also while a read replaces it
    this.selected = computed(() => this._selected.value);
    this.effect(() => {
      const node = this.tree.selectedNode.value;
      if (node?.id.startsWith(MORE) || node === null && this._loading.peek())
        return;
      this._selected.value = node?.data ?? null;
    });
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
    this._levels = new Map();
    try {
      const table = await this._ready;
      DomainTree.requireHierarchy(table);
      this._table = table;
      this._roots = await this._children(null, []);
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
    const [top] = await table.query({filter: {property: 'id', operator: '=', value: ids[0]}, limit: 1});
    // a chain that does not start at a root stops at an ancestor the caller cannot see: no page holds it
    if (top !== undefined && top[handle.info.parentColumn!] == null) {
      for (;;) {
        await this.tree.expandPath(ids);
        if (this.tree.selectedNode.peek()?.data?.id === id)
          return;
        const level = this._unread(ids);
        if (level === undefined)
          break;
        await this._next(level);
      }
    }
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
   * "Load more…" row. `path` is the ids from a root down to `parent`. */
  private async _children(parent: TreeNode<RowView<TRow>> | null, path: string[],
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
    const rows = await table.table.query(spec) as unknown as RowView<TRow>[];
    // a row inserted above the page boundary comes back; a page with nothing new ends the level
    const known = new Set(loaded.map((n) => n.id));
    const fresh = rows.filter((row) => !known.has(row.id)).map((row) => this._node(row, path));
    const nodes = [...loaded, ...fresh];
    const full = rows.length === limit && fresh.length > 0;
    const key = parent?.id ?? '';
    levels.set(key, {parent, path, nodes, full});
    const level = full ? [...nodes, {id: MORE + key, label: 'Load more…'}] : nodes;
    const all = this._options.allNode;
    return parent === null && all !== undefined ? [{id: ALL, label: all, tooltip: all}, ...level] : level;
  }

  private async _next(level: Level<RowView<TRow>>): Promise<TreeNode<RowView<TRow>>[] | undefined> {
    const levels = this._levels;
    const next = await this._children(level.parent, level.path, level.nodes);
    if (levels !== this._levels)
      return undefined;
    if (level.parent === null)
      this._roots = next;
    else
      level.parent.children = () => Promise.resolve(next);
    this.tree.setRoots(this._roots);
    return levels.get(level.parent?.id ?? '')!.nodes;
  }

  private async _loadMore(id: string): Promise<void> {
    const levels = this._levels;
    const level = levels.get(id.slice(MORE.length));
    if (level === undefined || !level.full || this._loading.peek())
      return;
    const onRow = this.tree.selectedNode.peek()?.id === id;
    this._loading.value = true;
    this._error.value = null;
    try {
      const nodes = await this._next(level);
      if (nodes === undefined)
        return;
      // the first row read takes the selection; with none, a cursor that was on the row stays on the level
      const target = nodes[level.nodes.length] ?? (onRow ? nodes[nodes.length - 1] : undefined);
      if (target !== undefined)
        await this.tree.expandPath([...level.path, target.id]);
    } catch (e) {
      if (levels === this._levels)
        this._error.value = DomainErrors.message(e);
    } finally {
      if (levels === this._levels)
        this._loading.value = false;
    }
  }

  private _unread(ids: string[]): Level<RowView<TRow>> | undefined {
    let key = '';
    for (const step of ids) {
      const level = this._levels.get(key);
      if (level === undefined || !level.nodes.some((n) => n.id === step))
        return level?.full ? level : undefined;
      key = step;
    }
    return undefined;
  }

  private _node(row: RowView<TRow>, path: string[]): TreeNode<RowView<TRow>> {
    const label = this._table!.renderer.caption(row);
    const node: TreeNode<RowView<TRow>> = {id: row.id, label, tooltip: label, data: row};
    node.children = () => this._children(node, [...path, node.id]);
    return node;
  }
}
