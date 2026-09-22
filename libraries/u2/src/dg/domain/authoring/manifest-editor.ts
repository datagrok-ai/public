/* The reusable manifest editor: the tree and the context panel side by side over one
   `ManifestModel` and one `AccessModel`, the field offer decided by the editor context.
   Platform-free — the dialog lane (PowerPack) drafts, dry-runs, creates and applies the access;
   this control only edits and answers `plan()`. Diagnostics the dry run returns are set on
   `diagnostics` and land on the row and the panel that own their manifest path. */
import {Control} from '../../../core/component.js';
import {signal, Signal, ReadonlySignal} from '../../../core/signals.js';
import {Splitter} from '../../../components/containers/splitter.js';
import type {EditorContext, FieldOffer} from './editor-context.js';
import {AccessModel, ManifestModel} from './manifest-model.js';
import type {AccessJson, DraftEnvelope, ManifestDiagnostic, ManifestJson, ManifestSelection} from './manifest-model.js';
import {ManifestTree} from './manifest-tree.js';
import {ManifestContextPanel} from './manifest-panel.js';

export interface ManifestEditorOptions {
  context: EditorContext;
  /** The groups and users the access pickers offer. */
  groups?: string[];
  /** The schema's friendly name, carried beside the manifest. */
  friendlyName?: string;
  /** Access rows to start from (an editor reopened on what it produced). */
  access?: Partial<AccessJson>;
  /** The left pane's share of the width (0.42 by default). */
  treeShare?: number;
}

/** A grant to apply after Create, by the names the API takes. */
export interface PlannedGrant {
  /** The table's logical name, or null for the schema (every table). */
  table: string | null;
  group: string;
  view: boolean;
  edit: boolean;
  delete: boolean;
}

/** A column restriction to apply after Create, by logical names. */
export interface PlannedRestriction {
  table: string;
  column: string;
  groups: string[];
}

/** Everything the dialog submits: the create envelope and what to apply afterwards. Grants and
 * restrictions of tables and columns the manifest no longer carries are left out. */
export interface ManifestPlan {
  name: string;
  friendlyName: string;
  manifest: ManifestJson;
  grants: PlannedGrant[];
  restrictions: PlannedRestriction[];
}

export class ManifestEditor extends Control {
  readonly model: ManifestModel;
  readonly access: AccessModel;
  readonly tree: ManifestTree;
  readonly panel: ManifestContextPanel;
  readonly selected: ReadonlySignal<ManifestSelection>;
  /** The dry run's findings, addressed by manifest path; set them and the rows light up. */
  readonly diagnostics: Signal<ManifestDiagnostic[]>;
  readonly context: EditorContext;

  constructor(draft: DraftEnvelope, options: ManifestEditorOptions) {
    super();
    this.context = options.context;
    this.model = new ManifestModel(draft, {friendlyName: options.friendlyName});
    this.access = new AccessModel(options.access);
    this.diagnostics = signal<ManifestDiagnostic[]>(draft.diagnostics ?? []);
    this.root.classList.add('u2-manifest-editor');
    this.root.dataset.u2 = 'manifest-editor';
    const editable = options.context.mode !== 'view';
    this.tree = this.runInScope(() => new ManifestTree(this.model, {editable, diagnostics: this.diagnostics}));
    this.selected = this.tree.selected;
    this.panel = this.runInScope(() => new ManifestContextPanel(this.model, {selected: this.selected,
      access: this.access, context: options.context, groups: options.groups, diagnostics: this.diagnostics}));
    const share = options.treeShare ?? 0.42;
    const splitter = this.runInScope(() => new Splitter([this.tree, this.panel],
      {direction: 'horizontal', sizes: [share, 1 - share], minSize: 200}));
    this.root.append(splitter.root);
  }

  get offer(): FieldOffer {
    return this.panel.offer;
  }

  /** Opens the path to the node and selects it — where a diagnostic's row is taken to. */
  select(selection: ManifestSelection): Promise<void> {
    return this.tree.select(selection);
  }

  /** The create envelope plus the access to apply, resolved to the manifest's logical names. */
  plan(): ManifestPlan {
    const model = this.model;
    const manifest = model.toJSON();
    const logicalOf = (remote: string): string | null => {
      const table = model.table(remote);
      return table !== undefined && table.included ? table.logical : null;
    };
    const grants: PlannedGrant[] = [];
    for (const g of this.access.grants.peek()) {
      const table = g.scope === 'schema' ? null : logicalOf(g.scope);
      if (g.scope === 'schema' || table !== null)
        grants.push({table, group: g.group, view: g.view, edit: g.edit, delete: g.delete});
    }
    const restrictions: PlannedRestriction[] = [];
    for (const v of this.access.visibility.peek()) {
      const table = logicalOf(v.table);
      const column = model.column(v.table, v.column);
      if (v.groups !== null && table !== null && column !== undefined && column.included)
        restrictions.push({table, column: column.logical, groups: [...v.groups]});
    }
    return {name: model.name.peek(), friendlyName: model.friendlyName.peek(), manifest, grants, restrictions};
  }
}
