import {Control, computed, ReadonlySignal, div, link, span, icon, plural, VirtualTree, badge} from '@datagrok-libraries/u2';
import type {TreeNode} from '@datagrok-libraries/u2';
import {ManifestModel} from './manifest-model.js';
import type {ColumnView, DriftView, ManifestDiagnostic, ManifestSelection, TableView} from './manifest-model.js';

export interface ManifestNode {
  selection: ManifestSelection;
  table?: TableView;
  column?: ColumnView;
  /** The diagnostics addressed to this node. */
  problems: ManifestDiagnostic[];
  add?: boolean;
}

export interface ManifestTreeOptions {
  /** Whether the checkboxes may be toggled; a view-mode editor locks them. */
  editable?: boolean;
  /** Diagnostics addressed by manifest path; the row that owns one is marked. */
  diagnostics?: ReadonlySignal<ManifestDiagnostic[]>;
}

const SCHEMA_ID = 'schema';

/** The Design step's left pane: schema › tables › columns. A view, a keyless table or an unsupported column
 * stays in the tree greyed with its reason, never hidden. */
export class ManifestTree extends Control {
  readonly tree: VirtualTree<ManifestNode>;
  /** What is selected — the schema until a row is chosen. */
  readonly selected: ReadonlySignal<ManifestSelection>;

  private readonly _editable: boolean;
  private readonly _diagnostics: ReadonlySignal<ManifestDiagnostic[]> | undefined;

  constructor(readonly model: ManifestModel, options: ManifestTreeOptions = {}) {
    super();
    this._editable = options.editable !== false;
    this._diagnostics = options.diagnostics;
    this.root.classList.add('u2-manifest-tree');
    this.root.dataset.u2 = 'manifest-tree';
    this.tree = this.build(() => new VirtualTree<ManifestNode>({
      render: (node) => this._label(node),
      onCheck: (node, checked) => this._check(node.data!, checked),
    }));
    this.selected = computed(() => this.tree.selectedNode.value?.data?.selection ?? {kind: 'schema'});
    this.root.append(this.build(() => this._header()), this.tree.root);
    this.effect(() => this.tree.setRoots([this._schemaNode()]));
    const first = model.tables.peek().find((t) => t.included && (!model.extending || t.extensible));
    this.tree.expanded.value = new Set(first === undefined ? [SCHEMA_ID] :
      [SCHEMA_ID, ManifestTree.tableId(first.remote)]);
    void this.tree.expandPath([SCHEMA_ID]);
  }

  static tableId(remote: string): string {
    return `t:${remote}`;
  }

  static columnId(table: string, column: string): string {
    return `c:${table}:${column}`;
  }

  /** The node ids from the root down to the selection's node. */
  static pathOf(selection: ManifestSelection): string[] {
    switch (selection.kind) {
    case 'schema': return [SCHEMA_ID];
    case 'table': return [SCHEMA_ID, ManifestTree.tableId(selection.table)];
    case 'column': return [SCHEMA_ID, ManifestTree.tableId(selection.table),
      ManifestTree.columnId(selection.table, selection.column)];
    }
  }

  /** Opens the path to the node and selects it. */
  select(selection: ManifestSelection): Promise<void> {
    return this.tree.expandPath(ManifestTree.pathOf(selection));
  }

  reveal(table: string): Promise<void> {
    this.tree.expanded.value = new Set([SCHEMA_ID, ManifestTree.tableId(table)]);
    return this.select({kind: 'table', table});
  }

  static ownership(column: ColumnView): HTMLElement {
    return !column.owned ? badge('plugin') : !column.registered ? badge('new', {variant: 'accent'}) :
      column.included ? badge('custom', {variant: 'success'}) : badge('removed', {variant: 'error'});
  }

  private _header(): HTMLElement {
    const model = this.model;
    const count = computed(() => {
      const tables = model.tables.value;
      if (model.extending)
        return plural(tables.filter((t) => t.extensible).length, 'table accepts', 'tables accept') + ' custom properties';
      const bindable = tables.filter((t) => t.bindable).length;
      return `${tables.filter((t) => t.included).length} of ${bindable} bindable tables`;
    });
    const header = div([span(count, 'u2-manifest-tree-count')], 'u2-manifest-tree-header');
    if (this._editable && !model.extending) {
      header.append(link('Check all', () => model.includeTables(true)),
        link('Clear', () => model.includeTables(false)));
    }
    return header;
  }

  private _schemaNode(): TreeNode<ManifestNode> {
    const model = this.model;
    const problems = this._problemsOf({kind: 'schema'});
    const name = model.name.value;
    return {
      id: SCHEMA_ID, label: name, tooltip: ManifestTree._tooltip(name, problems),
      data: {selection: {kind: 'schema'}, problems},
      children: model.tables.value.filter((t) => !model.extending || t.extensible).map((t) => this._tableNode(t)),
    };
  }

  private _tableNode(t: TableView): TreeNode<ManifestNode> {
    const selection: ManifestSelection = {kind: 'table', table: t.remote};
    const problems = this._problemsOf(selection);
    const children = this.model.columns(t.remote).value.map((c) => this._columnNode(t, c));
    const node: TreeNode<ManifestNode> = {
      id: ManifestTree.tableId(t.remote), label: t.remote,
      tooltip: ManifestTree._tooltip(t.reason ?? t.remote, problems),
      data: {selection, table: t, problems}, children,
    };
    const disabled = t.drafted ? null : t.reason ?? 'Not in this draft';
    if (!this.model.extending)
      return {...node, checked: t.included, locked: !t.drafted || !this._editable, disabled};
    if (this._editable)
      children.push({id: `a:${t.remote}`, label: 'Add property', data: {selection, table: t, problems: [], add: true}});
    return node;
  }

  private _columnNode(t: TableView, c: ColumnView): TreeNode<ManifestNode> {
    const selection: ManifestSelection = {kind: 'column', table: t.remote, column: c.remote};
    const problems = this._problemsOf(selection);
    const node: TreeNode<ManifestNode> = {
      id: ManifestTree.columnId(t.remote, c.remote), label: c.remote,
      tooltip: ManifestTree._tooltip(c.reason ?? c.relations[0]?.reason ?? c.remote, problems),
      data: {selection, table: t, column: c, problems},
    };
    return this.model.extending ? {...node, label: c.logical} : {...node,
      checked: c.supported && t.included && (c.isKey || c.included),
      locked: c.isKey && this._editable ? 'The primary key is always bound: the row id encodes it' :
        c.isKey || !c.supported || !t.included || !this._editable,
      disabled: c.supported ? null : c.reason ?? 'Its warehouse type cannot be bound'};
  }

  private _label(node: TreeNode<ManifestNode>): HTMLElement {
    const data = node.data!;
    if (data.add === true)
      return div([link('+ Add property', () => void this._add(data.table!.remote))], 'u2-manifest-node');
    const el = div([span(node.label, 'u2-manifest-node-name')], 'u2-manifest-node');
    if (data.problems.length > 0)
      el.classList.add('u2-manifest-node-problem');
    const column = data.column;
    const table = data.table;
    if (table !== undefined && (!table.included || column?.included === false))
      el.classList.add('u2-manifest-node-excluded');
    const editing = this.model.editing;
    if (column !== undefined && this.model.extending) {
      if (!column.owned)
        el.prepend(icon('lock', {cls: 'u2-manifest-node-hint'}));
      el.classList.toggle('u2-manifest-node-excluded', !column.owned || !column.included);
      el.append(span(column.type, 'u2-manifest-node-hint'));
      if (column.isKey)
        el.append(badge('key'));
      if (column.owned)
        el.append(ManifestTree.ownership(column));
    } else if (column !== undefined) {
      el.append(span(column.supported ? column.type : `${column.dbType ?? ''} · not bindable`,
        'u2-manifest-node-hint'));
      if (column.isKey)
        el.append(badge('key'));
      for (const relation of column.relations)
        el.append(badge(`→ ${relation.targetTable}`, {variant: relation.ref ? 'accent' : 'warning'}));
      if (column.isName)
        el.append(badge('name', {variant: 'success'}));
      if (editing && column.supported && !column.registered)
        el.append(badge('new', {variant: 'accent'}));
      ManifestTree._drift(el, column.drift);
    } else if (table !== undefined) {
      if (table.included && table.logical !== table.remote)
        el.append(span(`· ${table.logical}`, 'u2-manifest-node-hint'));
      if (!table.bindable)
        el.append(badge(ManifestTree.shortReason(table.code), {variant: 'error'}));
      else if (!table.drafted)
        el.append(badge('not in this draft', {variant: 'warning'}));
      else if (table.included && table.key.length > 0)
        el.append(badge(table.key.join(', ')));
      if (editing && table.drafted && !table.registered)
        el.append(badge('new', {variant: 'accent'}));
      ManifestTree._drift(el, table.drift);
    } else {
      if (this.model.storage?.kind === 'external')
        el.append(span(`over ${this.model.storage.schema ?? ''}`, 'u2-manifest-node-hint'));
      if (editing && this.model.catalog === 'unknown')
        el.append(badge('catalog not read', {variant: 'warning'}));
    }
    return el;
  }

  /** What the catalog says about a registered item, as a badge label; null for "unknown", which is
   * said once, on the schema. */
  static driftLabel(drift: DriftView | undefined): string | null {
    if (drift === undefined || drift.kind === 'unknown')
      return null;
    return drift.kind === 'missing' ? 'missing remotely' : drift.kind === 'unbindable' ? 'not bindable' :
      drift.kind === 'key' ? 'key changed' : 'type changed';
  }

  private static _drift(el: HTMLElement, drift: DriftView | undefined): void {
    const label = ManifestTree.driftLabel(drift);
    if (label !== null)
      el.append(badge(label, {variant: drift!.blocks ? 'error' : 'warning'}));
  }

  private async _add(table: string): Promise<void> {
    const column = this.model.addColumn(table);
    if (column !== null)
      await this.select({kind: 'column', table, column});
  }

  private _check(node: ManifestNode, checked: boolean): void {
    const selection = node.selection;
    if (selection.kind === 'table')
      this.model.includeTable(selection.table, checked);
    else if (selection.kind === 'column')
      this.model.includeColumn(selection.table, selection.column, checked);
  }

  private _problemsOf(selection: ManifestSelection): ManifestDiagnostic[] {
    const all = this._diagnostics?.value ?? [];
    return all.filter((d) => ManifestModel.sameSelection(this.model.resolvePath(d.path), selection));
  }

  /** The badge an unbindable table wears, from the draft's code. */
  static shortReason(code: string | undefined): string {
    return code === 'external-table-view' ? 'view' : code === 'external-key-missing' ? 'no key' : 'not bindable';
  }

  private static _tooltip(text: string, problems: ManifestDiagnostic[]): string {
    return problems.length === 0 ? text : `${text}\n${problems.map((p) => p.message).join('\n')}`;
  }
}
