/* Virtualized tree — a visible-node projection over `VirtualList` (src/components/collections/list.ts):
   roots + expanded descendants flatten into a row array fed to the list as a signal, so the list
   keeps owning virtualization, recycling, selection and vertical keyboard motion. The tree adds
   expansion, lazy branches, inline rename and the tree ARIA layer on top. A flattened row is kept
   while nothing it shows changed, so the list renders again only the nodes that did. */
import {Control} from '../../core/component.js';
import {Scope} from '../../core/scope.js';
import {signal, computed, untracked, Signal, ReadonlySignal} from '../../core/signals.js';
import {VirtualList} from './list.js';
import type {VirtualListOptions} from './list.js';
import {icon} from '../display/icon.js';
import type {Action} from '../actions/actions.js';
import {ITEMS_CAP} from '../../core/widget-like.js';
import type {IItemState, IItemsState} from '../../core/widget-like.js';

export interface TreeNode<T = unknown> {
  id: string;
  label: string;
  /** The row's `title`; the label itself when omitted, so a truncated row always says its name. */
  tooltip?: string;
  data?: T;
  children?: TreeNode<T>[] | (() => Promise<TreeNode<T>[]>);
  /** A checkbox before the label with this state; absent, the row carries none. */
  checked?: boolean;
  /** The checkbox cannot be toggled — a state the node's owner decides (a key column). */
  locked?: boolean;
  /** The row is greyed out (`u2-tree-row-disabled`, `aria-disabled`): still selectable, so the
   * reason in its tooltip and its detail can be read. */
  disabled?: boolean;
}

export interface TreeOptions<T> {
  itemHeight?: number;
  /** The node's label element; the label text in a `u2-tree-label` span by default. Whatever it
   * answers carries that class, so the row keeps its ellipsis and its rename handle. Runs in the
   * row's scope, again only when the node (or what its row shows) changes — `refresh(id)` for a
   * change it cannot see. */
  render?: (node: TreeNode<T>) => HTMLElement;
  onRename?: (node: TreeNode<T>, newLabel: string) => void;
  /** The node's full action list: right-click selects the row and opens it as a menu at the cursor. */
  contextActions?: (node: TreeNode<T>) => Action[];
  /** A checkbox toggled — by click or by Space on the selected row; the click never moves the
   * selection. The box never flips on its own: it shows `node.checked`, so the owner answers by
   * new roots (or `refresh(id)` after changing the node in place). Without it the boxes are disabled. */
  onCheck?: (node: TreeNode<T>, checked: boolean) => void;
  /** The checkbox's accessible name; the node's label by default. */
  checkLabel?: (node: TreeNode<T>) => string;
}

interface FlatRow<T> {
  node: TreeNode<T>;
  depth: number;
  kind: 'node' | 'loading' | 'error';
  message?: string;
  /** A branch that may still have children — one that answered none is a leaf. */
  branch: boolean;
  open: boolean;
  /** What the row shows of the node, so a node changed in place is still rendered again. */
  label: string;
  tooltip?: string;
  checked?: boolean;
  locked?: boolean;
  disabled?: boolean;
}

/** A lazy branch's read: the loader it ran, the scope it ran in, and its answer once there. */
interface Load<T> {
  loader: () => Promise<TreeNode<T>[]>;
  scope: Scope;
  promise: Promise<TreeNode<T>[]>;
  /** The answer — or, while a new loader runs, the previous one, shown until it is replaced. */
  children?: TreeNode<T>[];
  error?: string;
}

const INDENT = 12;

/** The tree's scroller: its `[role=tree]` root answers for the tree, so whoever lands on it reads
 * the nodes, not the flattened rows. */
class TreeRows<T> extends VirtualList<FlatRow<T>> {
  constructor(options: VirtualListOptions<FlatRow<T>>, private readonly _tree: VirtualTree<T>) {
    super(options);
  }

  getItems(opts?: {all?: boolean}): IItemsState {
    return this._tree.getItems(opts);
  }

  revealItem(keyOrLabel: string): Promise<HTMLElement | null> {
    return this._tree.revealItem(keyOrLabel);
  }

  /** The list's own reveal of a row by its key. */
  revealRow(key: string): Promise<HTMLElement | null> {
    return super.revealItem(key);
  }

  getState(): ReturnType<Control['getState']> {
    return this._tree.getState();
  }
}

export class VirtualTree<T = unknown> extends Control {
  readonly expanded: Signal<Set<string>> = signal(new Set<string>());
  readonly selectedNode: ReadonlySignal<TreeNode<T> | null>;
  /** The full multi-selection in selection order (F5); {@link selectedNode} is its lead. */
  readonly selectedNodes: ReadonlySignal<TreeNode<T>[]>;

  private readonly _list: TreeRows<T>;
  private readonly _roots = signal<TreeNode<T>[]>([]);
  private readonly _flat = signal<FlatRow<T>[]>([]);
  private readonly _loads = new Map<string, Load<T>>();
  /** Bumped when a load settles: what the flattening reads instead of the map itself. */
  private readonly _loaded = signal(0);
  private readonly _rename: ((node: TreeNode<T>, newLabel: string) => void) | undefined;
  private readonly _render: ((node: TreeNode<T>) => HTMLElement) | undefined;
  private readonly _check: ((node: TreeNode<T>, checked: boolean) => void) | undefined;
  private readonly _checkLabel: ((node: TreeNode<T>) => string) | undefined;
  private _shown = new Map<string, FlatRow<T>>();
  private _endRename: (() => void) | undefined;

  constructor(options: TreeOptions<T> = {}) {
    super();
    this._rename = options.onRename;
    this._render = options.render;
    this._check = options.onCheck;
    this._checkLabel = options.checkLabel;
    const contextActions = options.contextActions;
    this._list = new TreeRows<T>({
      itemHeight: options.itemHeight ?? 22,
      rowRole: 'treeitem',
      keyOf: VirtualTree._key,
      labelOf: (row) => row.node.label,
      render: (row, index, el) => this._renderRow(row, el),
      contextActions: contextActions ?
        (row) => row.kind === 'node' ? contextActions(row.node) : [] :
        undefined,
    }, this);
    this.own(() => this._list.dispose());
    this.own(() => this._endRename?.());

    this.root.classList.add('u2-tree');
    this.root.dataset.u2 = 'tree';
    this._list.root.setAttribute('role', 'tree');
    this._list.setItems(this._flat);
    this.root.append(this._list.root);

    this.selectedNode = computed(() => {
      const row = this._flat.value[this._list.selectedIndex.value];
      return row && row.kind === 'node' ? row.node : null;
    });
    this.selectedNodes = computed(() => {
      const rows = this._flat.value;
      const nodes: TreeNode<T>[] = [];
      for (const index of this._list.selectedIndices.value) {
        const row = rows[index];
        if (row && row.kind === 'node')
          nodes.push(row.node);
      }
      return nodes;
    });

    this.scope.listen(this._list.root, 'click', (e) => this._onClick(e as MouseEvent), {capture: true});
    // the whole row is the rename handle — short auto-names leave most of the row bare, and a
    // dblclick aimed at "the row" must not silently miss the text; the twistie keeps toggling
    this.scope.listen(this._list.root, 'dblclick', (e) => {
      const target = e.target as Element | null;
      if (target && !target.closest('.u2-tree-twistie'))
        this._startRename(VirtualTree._rowIndex(target));
    });
    this.scope.listen(this._list.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));

    this.effect(() => this._rebuild());
  }

  /** The tree shows `roots`; a node changed in place shows its change once the roots are set again. */
  setRoots(roots: TreeNode<T>[]): void {
    this._roots.value = [...roots];
    this._prune();
  }

  /** Renders the row of the node `id` again — or every row on screen: what a change the node's
   * fields do not show needs (whatever a custom `render` reads from `data`). */
  refresh(id?: string): void {
    this._list.refresh(id);
  }

  clearSelection(): void {
    this._list.selectedIndex.value = -1;
  }

  /** The node behind whatever element an event hit — null off the rows and on loading/error rows. */
  nodeForRow(target: EventTarget | null): TreeNode<T> | null {
    const row = this._flat.peek()[VirtualTree._rowIndex(target)];
    return row && row.kind === 'node' ? row.node : null;
  }

  /** Expands every ancestor in `ids` (awaiting lazy loads), then selects and reveals the last one. */
  async expandPath(ids: string[]): Promise<void> {
    if (ids.length === 0)
      return;
    let nodes = this._roots.peek();
    for (let i = 0; i < ids.length - 1; i++) {
      const node = nodes.find((n) => n.id === ids[i]);
      if (!node)
        return;
      nodes = await this._childrenOf(node);
    }
    // read after the loads, so the clicks on other branches made meanwhile are kept
    const current = this.expanded.peek();
    const expanded = new Set([...current, ...ids.slice(0, -1)]);
    // an unchanged set is never written: the rebuild would recycle every row between the two
    // clicks of a double-click, and the browser only fires dblclick on a stable target
    if (expanded.size !== current.size)
      this.expanded.value = expanded;
    const target = ids[ids.length - 1];
    const index = this._flat.peek().findIndex((r) => r.kind === 'node' && r.node.id === target);
    if (index >= 0)
      this._select(index);
  }

  /** The visible nodes in order (a collapsed branch's descendants are not items until it opens),
   * capped at 500 unless `all`: id, label, depth, expansion, checkbox and selection state. */
  getItems(opts?: {all?: boolean}): IItemsState {
    const rows = this._flat.peek();
    const selected = this._list.selectedIndices.peek();
    const list: IItemState[] = [];
    let count = 0;
    for (const [at, row] of rows.entries()) {
      if (row.kind !== 'node')
        continue;
      if (opts?.all || count < ITEMS_CAP) {
        list.push({key: row.node.id, label: row.label, index: count, level: row.depth + 1,
          expanded: row.branch ? row.open : undefined, checked: row.checked,
          selected: selected.has(at), disabled: row.disabled === true, rendered: this._list.rowAt(at) !== null});
      }
      count++;
    }
    return {count, list, truncated: list.length < count};
  }

  /** The base state, with the selection as the ids of the selected nodes, in selection order. */
  getState(): ReturnType<Control['getState']> {
    return {...super.getState(), selection: this.selectedNodes.peek().map((node) => node.id)};
  }

  /** Brings a node — by id, else by label among the nodes read so far — into view: opens its
   * ancestors, scrolls to it and resolves to its row once rendered; null for no match, a rejection
   * naming the candidates when several labels match. */
  async revealItem(keyOrLabel: string): Promise<HTMLElement | null> {
    const path = this._pathTo(keyOrLabel);
    if (path === null)
      return null;
    const current = this.expanded.peek();
    const expanded = new Set([...current, ...path.slice(0, -1)]);
    if (expanded.size !== current.size)
      this.expanded.value = expanded;
    const id = path[path.length - 1];
    return this._flat.peek().some((r) => r.kind === 'node' && r.node.id === id) ? this._list.revealRow(id) : null;
  }

  /** The ids from a root down to the node `keyOrLabel` names — its id first, else its label —
   * among the nodes known now (arrays and answered lazy reads, so every ancestor's children are
   * there to show once it opens); null for none. */
  private _pathTo(keyOrLabel: string): string[] | null {
    const wanted = keyOrLabel.trim().toLowerCase();
    const byId: string[][] = [];
    const byLabel: string[][] = [];
    const walk = (nodes: TreeNode<T>[], path: string[]): void => {
      for (const node of nodes) {
        const here = [...path, node.id];
        if (node.id === keyOrLabel)
          byId.push(here);
        else if (node.label.trim().toLowerCase() === wanted)
          byLabel.push(here);
        const children = Array.isArray(node.children) ? node.children : this._loads.get(node.id)?.children;
        if (children)
          walk(children, here);
      }
    };
    walk(this._roots.peek(), []);
    if (byId.length > 0)
      return byId[0];
    if (byLabel.length > 1) {
      throw new Error(`u2: ${byLabel.length} items match "${keyOrLabel}": ` +
        byLabel.map((path) => path[path.length - 1]).join(', '));
    }
    return byLabel[0] ?? null;
  }

  private _rebuild(): void {
    const rows: FlatRow<T>[] = [];
    const load: TreeNode<T>[] = [];
    this._loaded.value;
    this._flatten(this._roots.value, 0, rows, load);
    const shown = new Map<string, FlatRow<T>>();
    const kept = rows.map((row) => {
      const key = VirtualTree._key(row);
      const old = this._shown.get(key);
      const same = old !== undefined && VirtualTree._same(old, row) ? old : row;
      shown.set(key, same);
      return same;
    });
    this._shown = shown;
    this._flat.value = kept;
    // read outside the flattening: a loader neither subscribes it nor lives in its run
    untracked(() => {
      for (const node of load)
        void this._load(node);
    });
  }

  private _flatten(nodes: TreeNode<T>[], depth: number, rows: FlatRow<T>[], load: TreeNode<T>[]): void {
    const expanded = this.expanded.value;
    for (const node of nodes) {
      const lazy = typeof node.children === 'function' ? this._loads.get(node.id) : undefined;
      const current = lazy !== undefined && lazy.loader === node.children;
      // a branch that has answered with no children is a leaf: its twistie would open nothing
      const branch = VirtualTree._isBranch(node) && !(current && lazy!.children?.length === 0);
      const open = branch && expanded.has(node.id);
      rows.push({node, depth, kind: 'node', branch, open, label: node.label, tooltip: node.tooltip,
        checked: node.checked, locked: node.locked, disabled: node.disabled});
      if (!expanded.has(node.id))
        continue;
      if (Array.isArray(node.children)) {
        this._flatten(node.children, depth + 1, rows, load);
        continue;
      }
      if (typeof node.children !== 'function')
        continue;
      if (!current)
        load.push(node);
      if (lazy?.children !== undefined)
        this._flatten(lazy.children, depth + 1, rows, load);
      else if (current && lazy!.error !== undefined)
        rows.push(VirtualTree._status(node, depth + 1, 'error', lazy!.error));
      else
        rows.push(VirtualTree._status(node, depth + 1, 'loading'));
    }
  }

  private _childrenOf(node: TreeNode<T>): Promise<TreeNode<T>[]> {
    if (Array.isArray(node.children))
      return Promise.resolve(node.children);
    if (typeof node.children !== 'function')
      return Promise.resolve([]);
    return this._load(node);
  }

  /** Reads a lazy branch once per loader: the read runs in its own scope, owned by the tree — an
   * abort signal or a view it takes lives until the tree drops the answer — and a failed read is
   * tried again when asked again (expanding the branch after collapsing it). New roots with new
   * loaders read again, showing the previous answer until the new one is in. */
  private _load(node: TreeNode<T>): Promise<TreeNode<T>[]> {
    const loader = node.children as () => Promise<TreeNode<T>[]>;
    const previous = this._loads.get(node.id);
    if (previous !== undefined && previous.loader === loader && previous.error === undefined)
      return previous.promise;
    previous?.scope.dispose();
    // not lazy: an abort signal alone does not link a lazy scope, and the tree must still abort it
    const scope = this.scope.child('tree-load');
    const load: Load<T> = {loader, scope, promise: Promise.resolve([]), children: previous?.children};
    this._loads.set(node.id, load);
    let started: Promise<TreeNode<T>[]>;
    try {
      started = Scope.runWith(scope, () => untracked(loader));
    } catch (e) {
      started = Promise.reject(e);
    }
    const settle = (children: TreeNode<T>[] | undefined, error?: string): void => {
      if (this._loads.get(node.id) !== load)
        return;
      const replaced = load.children !== undefined;
      load.children = children;
      load.error = error;
      if (replaced)
        this._prune();
      this._loaded.value++;
    };
    load.promise = started.then((children) => {
      settle(children);
      return children;
    }, (e: unknown) => {
      settle(undefined, e instanceof Error ? e.message : String(e));
      return [] as TreeNode<T>[];
    });
    return load.promise;
  }

  /** Drops the reads of nodes no longer in the tree: under the roots, through the arrays and the
   * answered reads, open or not — a collapsed branch keeps what it read. */
  private _prune(): void {
    if (this._loads.size === 0)
      return;
    const known = new Set<string>();
    const walk = (nodes: TreeNode<T>[]): void => {
      for (const node of nodes) {
        if (known.has(node.id))
          continue;
        known.add(node.id);
        const children = Array.isArray(node.children) ? node.children : this._loads.get(node.id)?.children;
        if (children)
          walk(children);
      }
    };
    walk(this._roots.peek());
    for (const [id, load] of this._loads) {
      if (!known.has(id)) {
        load.scope.dispose();
        this._loads.delete(id);
      }
    }
  }

  private _renderRow(row: FlatRow<T>, rowEl: HTMLElement): HTMLElement {
    rowEl.setAttribute('aria-level', String(row.depth + 1));
    const el = document.createElement('div');
    el.className = 'u2-tree-row';
    el.style.setProperty('--u2-tree-indent', `${row.depth * INDENT}px`);
    if (row.kind !== 'node') {
      // a status row is not an item: nothing addresses "Loading…"
      rowEl.removeAttribute('data-u2-item');
      el.append(VirtualTree._span('u2-tree-twistie u2-tree-twistie-hidden', '▸'),
        VirtualTree._span(`u2-tree-${row.kind}`, row.kind === 'loading' ? 'Loading…' : row.message ?? ''));
      return el;
    }
    if (row.branch)
      rowEl.setAttribute('aria-expanded', String(row.open));
    if (row.checked !== undefined)
      rowEl.setAttribute('aria-checked', String(row.checked));
    const twistie = row.branch ? 'u2-tree-twistie' : 'u2-tree-twistie u2-tree-twistie-hidden';
    const label = this._render === undefined ? VirtualTree._span('u2-tree-label', row.label) :
      this._render(row.node);
    label.classList.add('u2-tree-label');
    label.dataset.u2Part = 'item-label';
    el.append(VirtualTree._span(twistie, row.open ? '▾' : '▸'));
    if (row.checked !== undefined)
      el.append(this._checkbox(row));
    el.append(label);
    if (row.disabled === true) {
      el.classList.add('u2-tree-row-disabled');
      rowEl.setAttribute('aria-disabled', 'true');
    }
    el.title = row.tooltip ?? row.label;
    return el;
  }

  private _onClick(e: MouseEvent): void {
    const target = e.target as Element | null;
    if (!target)
      return;
    const row = this._flat.peek()[VirtualTree._rowIndex(target)];
    if (!row || row.kind !== 'node')
      return;
    if (target.closest('.u2-tree-check')) {
      // the checkbox is not the row: toggling it leaves the selection where it is; and the box
      // never flips itself — it shows the node, and only the owner's answer changes the node
      e.stopPropagation();
      e.preventDefault();
      this._toggleCheck(row.node);
    } else if (target.closest('.u2-tree-twistie')) {
      e.stopPropagation();
      this._setExpanded(row.node.id, !this.expanded.peek().has(row.node.id));
    }
  }

  private _toggleCheck(node: TreeNode<T>): void {
    if (this._check !== undefined && node.checked !== undefined && node.locked !== true)
      this._check(node, !node.checked);
  }

  private _onKeyDown(e: KeyboardEvent): void {
    const rows = this._flat.peek();
    const index = this._list.selectedIndex.peek();
    const row = rows[index];
    if (!row)
      return;
    switch (e.key) {
    case 'ArrowRight':
      e.preventDefault();
      if (row.kind === 'node' && VirtualTree._isBranch(row.node) && !this.expanded.peek().has(row.node.id))
        this._setExpanded(row.node.id, true);
      else if (index + 1 < rows.length && rows[index + 1].depth > row.depth)
        this._select(index + 1);
      break;
    case 'ArrowLeft':
      e.preventDefault();
      if (row.kind === 'node' && this.expanded.peek().has(row.node.id)) {
        this._setExpanded(row.node.id, false);
        break;
      }
      for (let i = index - 1; i >= 0; i--) {
        if (rows[i].depth < row.depth) {
          this._select(i);
          break;
        }
      }
      break;
    case 'Enter':
      e.preventDefault();
      this._select(index);
      break;
    case 'F2':
      e.preventDefault();
      this._startRename(index);
      break;
    case ' ':
      if (row.kind !== 'node' || row.node.checked === undefined)
        break;
      e.preventDefault();
      this._toggleCheck(row.node);
      break;
    }
  }

  private _startRename(index: number): void {
    const rename = this._rename;
    const row = this._flat.peek()[index];
    if (!rename || !row || row.kind !== 'node')
      return;
    this._endRename?.();
    const label = this._list.rowAt(index)?.querySelector('.u2-tree-label') as HTMLElement | null;
    if (!label)
      return;

    const input = document.createElement('input');
    input.className = 'u2-tree-rename';
    input.value = row.node.label;
    label.replaceWith(input);
    input.focus();
    input.select();

    let active = true;
    const end = (newLabel?: string) => {
      if (!active)
        return;
      active = false;
      if (newLabel && newLabel !== row.node.label) {
        label.textContent = newLabel;
        rename(row.node, newLabel);
      }
      if (input.isConnected)
        input.replaceWith(label);
    };
    this._endRename = end;
    input.addEventListener('blur', () => end());
    input.addEventListener('keydown', (ev) => {
      ev.stopPropagation();
      if (ev.key !== 'Enter' && ev.key !== 'Escape')
        return;
      end(ev.key === 'Enter' ? input.value.trim() : undefined);
      this._list.root.focus();
    });
  }

  /** Collapsing a branch whose read failed forgets the failure, so opening it again reads again. */
  private _setExpanded(id: string, expand: boolean): void {
    const next = new Set(this.expanded.peek());
    if (expand)
      next.add(id);
    else {
      next.delete(id);
      const load = this._loads.get(id);
      if (load?.error !== undefined) {
        load.scope.dispose();
        this._loads.delete(id);
      }
    }
    this.expanded.value = next;
  }

  private _select(index: number): void {
    this._list.selectedIndex.value = index;
    this._list.scrollToIndex(index);
  }

  /** The disabled box is the no-owner case: nothing would answer its click. */
  private _checkbox(row: FlatRow<T>): HTMLInputElement {
    const box = document.createElement('input');
    box.type = 'checkbox';
    box.className = 'u2-tree-check';
    box.checked = row.checked === true;
    // a natively disabled box paints grey even when checked: a locked checked one only refuses the click
    box.disabled = this._check === undefined || row.locked === true && row.checked !== true;
    if (row.locked === true && row.checked === true)
      box.setAttribute('aria-disabled', 'true');
    box.tabIndex = -1;
    box.setAttribute('aria-label', this._checkLabel?.(row.node) ?? row.label);
    return box;
  }

  /** The list key of a row: a node row is its node's id, the `data-u2-item` tests address it by. */
  private static _key<T>(row: FlatRow<T>): string {
    return row.kind === 'node' ? row.node.id : `${row.kind}:${row.node.id}`;
  }

  private static _same<T>(a: FlatRow<T>, b: FlatRow<T>): boolean {
    return a.node === b.node && a.depth === b.depth && a.kind === b.kind && a.message === b.message &&
      a.branch === b.branch && a.open === b.open && a.label === b.label && a.tooltip === b.tooltip &&
      a.checked === b.checked && a.locked === b.locked && a.disabled === b.disabled;
  }

  private static _status<T>(node: TreeNode<T>, depth: number, kind: 'loading' | 'error', message?: string): FlatRow<T> {
    return {node, depth, kind, message, branch: false, open: false, label: node.label};
  }

  private static _isBranch(node: TreeNode<unknown>): boolean {
    return node.children !== undefined;
  }

  private static _rowIndex(target: EventTarget | null): number {
    const row = target instanceof Element ? target.closest('.u2-list-row') as HTMLElement | null : null;
    return row ? Number(row.dataset.index) : -1;
  }

  /** The platform's caret, not a text glyph: it inherits the icon sheet every other twistie in
   * the app is drawn from, and keeps its box when it is hidden so the rows stay aligned. */
  private static _twistie(branch: boolean, open: boolean): HTMLElement {
    const el = icon(open ? 'caret-down' : 'caret-right', {variant: 'solid', cls: 'u2-tree-twistie'});
    if (!branch)
      el.classList.add('u2-tree-twistie-hidden');
    return el;
  }

  private static _span(cls: string, text: string): HTMLElement {
    const el = document.createElement('span');
    el.className = cls;
    el.textContent = text;
    return el;
  }
}
