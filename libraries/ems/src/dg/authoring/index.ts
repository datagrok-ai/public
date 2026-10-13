/* `domains.authoring` — the reusable manifest editor (tree + context panel) and its models, built
   against fixture drafts headlessly. Platform-free: nothing here imports `grok`; the dialog lane
   adds `createBinding` (draft → dry run → create → grants) to this namespace. */
import {ManifestEditor} from './manifest-editor.js';
import {ManifestModel, AccessModel} from './manifest-model.js';
import {ManifestTree} from './manifest-tree.js';
import {ManifestContextPanel} from './manifest-panel.js';
import {ManifestRules} from './manifest-rules.js';
import {fieldOffer} from './editor-context.js';

export {ManifestEditor, ManifestModel, AccessModel, ManifestTree, ManifestContextPanel, ManifestRules, fieldOffer};
export type {ManifestEditorOptions, ManifestPlan, PlannedGrant, PlannedRestriction, ApplyPayload, ApplyAccess,
  EditPlan} from './manifest-editor.js';
export type {DraftEnvelope, DraftInventory, InventoryTable, InventoryColumn, InventoryRelation, ManifestJson,
  ManifestTableJson, ManifestColumnJson, ManifestStorageJson, ManifestDiagnostic, ManifestSelection, TableView,
  ColumnView, RelationView, DriftView, ManifestChange, ManifestPatch, RebaseOp, RebaseReport, ManifestModelOptions,
  AccessGrant, AccessCapability, ColumnVisibility, AccessJson, AccessSnapshot, AccessSnapshotTable,
  AccessSnapshotColumn, AccessSnapshotGrant, AccessSnapshotGroup, AccessDelta, AccessTriple, AccessRestriction,
  AccessEditOptions}
  from './manifest-model.js';
export type {ManifestTreeOptions, ManifestNode} from './manifest-tree.js';
export type {ManifestContextPanelOptions, PrincipalPicker} from './manifest-panel.js';
export type {EditorContext, EditorMode, EditorStorage, FieldOffer} from './editor-context.js';

/** The namespace `domains.authoring` exposes; a plain object, so the dialog lane extends it. */
export const authoring = {
  ManifestEditor, ManifestModel, AccessModel, ManifestTree, ManifestContextPanel, ManifestRules, fieldOffer,
};
