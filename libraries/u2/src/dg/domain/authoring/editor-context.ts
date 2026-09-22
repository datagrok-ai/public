/* The two-part seam the editor's field offer follows (Astra 3): the manifest STORAGE decides the
   vocabulary — what the parser accepts for a table and a column — and the editor CONTEXT decides
   what may be changed here: create, edit, view, extend. Nothing in the editor branches on who
   made the schema. Only the arms this feature ships are built; the others refuse by name so a
   later feature fills them in rather than inheriting a wrong offer silently. */

export type EditorMode = 'create' | 'edit' | 'view' | 'extend';
export type EditorStorage = 'external' | 'domain';

export interface EditorContext {
  mode: EditorMode;
  storage: EditorStorage;
}

export interface FieldOffer {
  /** Whether anything may be changed at all — `view` shows every field as text. */
  editable: boolean;
  /** The mapped type is read-only (external: the draft owns the mapping). */
  type: 'readonly' | 'editable';
  /** Refs come from the warehouse's foreign keys, never from an edit. */
  refs: 'relations' | 'editable';
  /** The remote names and the binding step exist. */
  remoteNames: boolean;
  /** Schema-level `storage.writable` and the per-table read-only opt-out. */
  writable: boolean;
  required: boolean;
  nameColumn: boolean;
  searchable: boolean;
  /** Write-side vocabulary an external storage refuses. */
  defaultValue: boolean;
  autoNumber: boolean;
  immutable: boolean;
  unique: boolean;
  choices: boolean;
}

export function fieldOffer(context: EditorContext): FieldOffer {
  let vocabulary: Omit<FieldOffer, 'editable'>;
  switch (context.storage) {
  case 'external':
    vocabulary = {type: 'readonly', refs: 'relations', remoteNames: true, writable: true, required: true,
      nameColumn: true, searchable: true, defaultValue: false, autoNumber: false, immutable: false,
      unique: false, choices: false};
    break;
  case 'domain':
    throw new Error('u2: the manifest editor over a platform-stored schema is not built yet (ui2/ems ruling 6)');
  default:
    throw new Error(`u2: unknown manifest storage "${String(context.storage)}"`);
  }
  switch (context.mode) {
  case 'create':
    return {editable: true, ...vocabulary};
  case 'view':
    return {editable: false, ...vocabulary};
  case 'edit':
  case 'extend':
    throw new Error(`u2: the manifest editor's "${context.mode}" mode is not built yet`);
  default:
    throw new Error(`u2: unknown manifest editor mode "${String(context.mode)}"`);
  }
}
