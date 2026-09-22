/* The Design step's left pane: schema › tables › columns as a `VirtualTree` over the manifest
   model. A checkbox per table and column includes it (keys and unsupported columns locked), the
   badges say what the manifest will carry (key, `→ target` as a ref or a plain value, name), and
   a view, a keyless table or an unsupported column stays in the tree greyed with its reason —
   never hidden. The roots are rebuilt from the model on every change; the tree keeps the
   selection and the expansion by node id. */
import {Control} from '../../../core/component.js';
import {computed, ReadonlySignal} from '../../../core/signals.js';
import {div, link, span} from '../../../core/elements.js';
import {VirtualTree} from '../../../components/collections/tree.js';
import type {TreeNode} from '../../../components/collections/tree.js';
import {badge} from '../../../components/display/badge.js';
import {ManifestModel} from './manifest-model.js';
import type {ColumnView, ManifestDiagnostic, ManifestSelection, TableView} from './manifest-model.js';

export interface ManifestNode {
  selection: ManifestSelection;
  table?: TableView;
  column?: ColumnView;
  /** The diagnostics addressed to this node. */
  problems: ManifestDiagnostic[];
}

export interface ManifestTreeOptions {
  /** Whether the checkboxes may be toggled; a view-mode editor locks them. */
  editable?: boolean;
  /** Diagnostics addressed by manifest path; the row that owns one is marked. */
  diagnostics?: ReadonlySignal<ManifestDiagnostic[]>;
}

const SCHEMA_ID = 'schema';

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
    this.tree = this.runInScope(() => new VirtualTree<ManifestNode>({
      render: (node) => this._label(node),
      onCheck: (node, checked) => this._check(node.data!, checked),
    }));
    this.selected = computed(() => this.tree.selectedNode.value?.data?.selection ?? {kind: 'schema'});
    this.root.append(this.runInScope(() => this._header()), this.tree.root);
    this.effect(() => this.tree.setRoots([this._schemaNode()]));
    const first = model.tables.peek().find((t) => t.included);
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

  private _header(): HTMLElement {
    const model = this.model;
    const count = computed(() => {
      const tables = model.tables.value;
      const bindable = tables.filter((t) => t.bindable).length;
      return `${tables.filter((t) => t.included).length} of ${bindable} bindable tables`;
    });
    const header = div([span(count, 'u2-manifest-tree-count')], 'u2-manifest-tree-header');
    if (this._editable) {
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
      children: model.tables.value.map((t) => this._tableNode(t)),
    };
  }

  private _tableNode(t: TableView): TreeNode<ManifestNode> {
    const selection: ManifestSelection = {kind: 'table', table: t.remote};
    const problems = this._problemsOf(selection);
    const columns = this.model.columns(t.remote).value;
    return {
      id: ManifestTree.tableId(t.remote), label: t.remote,
      tooltip: ManifestTree._tooltip(t.reason ?? t.remote, problems),
      data: {selection, table: t, problems},
      checked: t.bindable ? t.included : false,
      locked: !t.bindable || !this._editable,
      disabled: !t.bindable,
      children: columns.map((c) => this._columnNode(t, c)),
    };
  }

  private _columnNode(t: TableView, c: ColumnView): TreeNode<ManifestNode> {
    const selection: ManifestSelection = {kind: 'column', table: t.remote, column: c.remote};
    const problems = this._problemsOf(selection);
    return {
      id: ManifestTree.columnId(t.remote, c.remote), label: c.remote,
      tooltip: ManifestTree._tooltip(c.reason ?? c.relation?.reason ?? c.remote, problems),
      data: {selection, table: t, column: c, problems},
      checked: c.supported && c.included && t.included,
      locked: c.isKey || !c.supported || !t.included || !this._editable,
      disabled: !c.supported,
    };
  }

  private _label(node: TreeNode<ManifestNode>): HTMLElement {
    const data = node.data!;
    const el = div([span(node.label, 'u2-manifest-node-name')], 'u2-manifest-node');
    if (data.problems.length > 0)
      el.classList.add('u2-manifest-node-problem');
    const column = data.column;
    const table = data.table;
    if (table !== undefined && (!table.included || column?.included === false))
      el.classList.add('u2-manifest-node-excluded');
    if (column !== undefined) {
      el.append(span(column.supported ? column.type : `${column.dbType ?? ''} · not bindable`,
        'u2-manifest-node-hint'));
      if (column.isKey)
        el.append(badge('key'));
      const relation = column.relation;
      if (relation !== undefined)
        el.append(badge(`→ ${relation.targetTable}`, {variant: relation.ref ? 'accent' : 'warning'}));

      if (column.isName)
        el.append(badge('name', {variant: 'success'}));
    } else if (table !== undefined) {
      if (table.included && table.logical !== table.remote)
        el.append(span(`· ${table.logical}`, 'u2-manifest-node-hint'));
      if (!table.bindable)
        el.append(badge(ManifestTree.shortReason(table.code), {variant: 'error'}));
      else if (table.included)
        el.append(badge(table.key.join(', ')));
    } else
      el.append(span(`over ${this.model.storage?.schema ?? ''}`, 'u2-manifest-node-hint'));
    return el;
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
