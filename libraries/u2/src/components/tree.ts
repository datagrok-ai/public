/* Virtualized tree — a visible-node projection over `VirtualList` (src/components/list.ts):
   roots + expanded descendants flatten into a row array fed to the list as a signal, so the list
   keeps owning virtualization, recycling, selection and vertical keyboard motion. The tree adds
   expansion, lazy branches, inline rename and the tree ARIA layer on top. */
import {Component} from '../core/component.js';
import {signal, computed, Signal, ReadonlySignal} from '../core/signals.js';
import {VirtualList} from './list.js';

export interface TreeNode<T = unknown> {
  id: string;
  label: string;
  data?: T;
  children?: TreeNode<T>[] | (() => Promise<TreeNode<T>[]>);
}

export interface TreeOptions<T> {
  itemHeight?: number;
  onRename?: (node: TreeNode<T>, newLabel: string) => void;
}

interface FlatRow<T> {
  node: TreeNode<T>;
  depth: number;
  kind: 'node' | 'loading' | 'error';
  message?: string;
}

const INDENT = 12;

export class VirtualTree<T = unknown> extends Component {
  readonly expanded: Signal<Set<string>> = signal(new Set<string>());
  readonly selectedNode: ReadonlySignal<TreeNode<T> | null>;

  private readonly _list: VirtualList<FlatRow<T>>;
  private readonly _roots = signal<TreeNode<T>[]>([]);
  private readonly _flat = signal<FlatRow<T>[]>([]);
  private readonly _loaded = signal(new Map<string, TreeNode<T>[]>());
  private readonly _errors = signal(new Map<string, string>());
  private readonly _pending = new Map<string, Promise<TreeNode<T>[]>>();
  private readonly _rename: ((node: TreeNode<T>, newLabel: string) => void) | undefined;
  private _endRename: (() => void) | undefined;

  constructor(options: TreeOptions<T> = {}) {
    super();
    this._rename = options.onRename;
    this._list = new VirtualList<FlatRow<T>>({
      itemHeight: options.itemHeight ?? 22,
      rowRole: 'treeitem',
      keyOf: (row) => `${row.kind}:${row.node.id}`,
      render: (row, index, el) => this._renderRow(row, el),
    });
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

    this._listen(this._list.root, 'click', (e) => this._onClick(e as MouseEvent), true);
    this._listen(this._list.root, 'dblclick', (e) => this._startRename(VirtualTree._rowIndex(e.target)));
    this._listen(this._list.root, 'keydown', (e) => this._onKeyDown(e as KeyboardEvent));

    this.effect(() => this._rebuild());
  }

  setRoots(roots: TreeNode<T>[]): void {
    this._roots.value = roots;
  }

  /** Expands every ancestor in `ids` (awaiting lazy loads), then selects and reveals the last one. */
  async expandPath(ids: string[]): Promise<void> {
    if (ids.length === 0)
      return;
    const expanded = new Set(this.expanded.peek());
    let nodes = this._roots.peek();
    for (let i = 0; i < ids.length - 1; i++) {
      const node = nodes.find((n) => n.id === ids[i]);
      if (!node)
        return;
      expanded.add(node.id);
      nodes = await this._childrenOf(node);
    }
    this.expanded.value = expanded;
    const target = ids[ids.length - 1];
    const index = this._flat.peek().findIndex((r) => r.kind === 'node' && r.node.id === target);
    if (index >= 0)
      this._select(index);
  }

  private _rebuild(): void {
    const rows: FlatRow<T>[] = [];
    this._flatten(this._roots.value, 0, rows);
    this._flat.value = rows;
  }

  private _flatten(nodes: TreeNode<T>[], depth: number, rows: FlatRow<T>[]): void {
    const expanded = this.expanded.value;
    for (const node of nodes) {
      rows.push({node, depth, kind: 'node'});
      if (!expanded.has(node.id))
        continue;
      if (Array.isArray(node.children)) {
        this._flatten(node.children, depth + 1, rows);
        continue;
      }
      if (typeof node.children !== 'function')
        continue;
      const message = this._errors.value.get(node.id);
      if (message !== undefined) {
        rows.push({node, depth: depth + 1, kind: 'error', message});
        continue;
      }
      const children = this._loaded.value.get(node.id);
      if (children) {
        this._flatten(children, depth + 1, rows);
        continue;
      }
      rows.push({node, depth: depth + 1, kind: 'loading'});
      this._load(node);
    }
  }

  private _childrenOf(node: TreeNode<T>): Promise<TreeNode<T>[]> {
    if (Array.isArray(node.children))
      return Promise.resolve(node.children);
    if (typeof node.children !== 'function')
      return Promise.resolve([]);
    return this._load(node);
  }

  private _load(node: TreeNode<T>): Promise<TreeNode<T>[]> {
    const cached = this._loaded.peek().get(node.id);
    if (cached)
      return Promise.resolve(cached);
    const started = this._pending.get(node.id);
    if (started)
      return started;
    const loader = node.children as () => Promise<TreeNode<T>[]>;
    const promise = loader().then((children) => {
      this._pending.delete(node.id);
      this._loaded.value = new Map(this._loaded.peek()).set(node.id, children);
      return children;
    }, (e: unknown) => {
      this._pending.delete(node.id);
      const message = e instanceof Error ? e.message : String(e);
      this._errors.value = new Map(this._errors.peek()).set(node.id, message);
      return [] as TreeNode<T>[];
    });
    this._pending.set(node.id, promise);
    return promise;
  }

  private _renderRow(row: FlatRow<T>, rowEl: HTMLElement): HTMLElement {
    rowEl.setAttribute('aria-level', String(row.depth + 1));
    const el = document.createElement('div');
    el.className = 'u2-tree-row';
    el.style.setProperty('--u2-tree-indent', `${row.depth * INDENT}px`);
    if (row.kind !== 'node') {
      el.append(VirtualTree._span('u2-tree-twistie u2-tree-twistie-hidden', '▸'),
        VirtualTree._span(`u2-tree-${row.kind}`, row.kind === 'loading' ? 'Loading…' : row.message ?? ''));
      return el;
    }
    const branch = VirtualTree._isBranch(row.node);
    const open = branch && this.expanded.peek().has(row.node.id);
    if (branch)
      rowEl.setAttribute('aria-expanded', String(open));
    const twistie = branch ? 'u2-tree-twistie' : 'u2-tree-twistie u2-tree-twistie-hidden';
    el.append(VirtualTree._span(twistie, open ? '▾' : '▸'), VirtualTree._span('u2-tree-label', row.node.label));
    return el;
  }

  private _onClick(e: MouseEvent): void {
    const target = e.target as Element | null;
    if (!target || !target.closest('.u2-tree-twistie'))
      return;
    const row = this._flat.peek()[VirtualTree._rowIndex(target)];
    if (!row || row.kind !== 'node')
      return;
    e.stopPropagation();
    this._setExpanded(row.node.id, !this.expanded.peek().has(row.node.id));
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
    }
  }

  private _startRename(index: number): void {
    const rename = this._rename;
    const row = this._flat.peek()[index];
    if (!rename || !row || row.kind !== 'node')
      return;
    this._endRename?.();
    const rowEl = this._list.root.querySelector(`.u2-list-row[data-index="${index}"]`);
    const label = rowEl?.querySelector('.u2-tree-label') as HTMLElement | null;
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

  private _setExpanded(id: string, expand: boolean): void {
    const next = new Set(this.expanded.peek());
    if (expand)
      next.add(id);
    else
      next.delete(id);
    this.expanded.value = next;
  }

  private _select(index: number): void {
    this._list.selectedIndex.value = index;
    this._list.scrollToIndex(index);
  }

  private _listen(el: EventTarget, type: string, handler: (e: Event) => void, capture = false): void {
    el.addEventListener(type, handler, capture);
    this.own(() => el.removeEventListener(type, handler, capture));
  }

  private static _isBranch(node: TreeNode<unknown>): boolean {
    return node.children !== undefined;
  }

  private static _rowIndex(target: EventTarget | null): number {
    const row = target instanceof Element ? target.closest('.u2-list-row') as HTMLElement | null : null;
    return row ? Number(row.dataset.index) : -1;
  }

  private static _span(cls: string, text: string): HTMLElement {
    const el = document.createElement('span');
    el.className = cls;
    el.textContent = text;
    return el;
  }
}
