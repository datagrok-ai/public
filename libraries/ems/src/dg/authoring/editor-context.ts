export type EditorMode = 'create' | 'edit' | 'view' | 'extend';
export type EditorStorage = 'external' | 'domain';

export interface EditorContext {
  mode: EditorMode;
  storage: EditorStorage;
}

export interface FieldOffer {
  /** Whether anything may be changed at all — `view` shows every field as text. */
  editable: boolean;
  /** The schema identifier may be set (a registered one is locked). */
  identifier: boolean;
  /** The schema description travels with the change (the entity's, beside the friendly name). */
  description: boolean;
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

/** The manifest storage decides the vocabulary, the editor context what may be changed here; the per-item edit
 * locks are the model's. An arm not built yet refuses by name, so a later feature fills it in rather than
 * inheriting a wrong offer silently. */
export function fieldOffer(context: EditorContext): FieldOffer {
  let vocabulary: Omit<FieldOffer, 'editable' | 'identifier' | 'description'>;
  switch (context.storage) {
  case 'external':
    vocabulary = {type: 'readonly', refs: 'relations', remoteNames: true, writable: true, required: true,
      nameColumn: true, searchable: true, defaultValue: false, autoNumber: false, immutable: false,
      unique: false, choices: false};
    break;
  case 'domain':
    if (context.mode !== 'extend')
      throw new Error('u2: the manifest editor over a platform-stored schema is not built yet (ui2/ems ruling 6)');
    vocabulary = {type: 'editable', refs: 'relations', remoteNames: false, writable: false, required: false,
      nameColumn: false, searchable: false, defaultValue: true, autoNumber: false, immutable: false,
      unique: false, choices: true};
    break;
  default:
    throw new Error(`u2: unknown manifest storage "${String(context.storage)}"`);
  }
  switch (context.mode) {
  case 'create':
    return {editable: true, identifier: true, description: false, ...vocabulary};
  case 'edit':
    return {editable: true, identifier: false, description: true, ...vocabulary};
  case 'view':
    return {editable: false, identifier: false, description: false, ...vocabulary};
  case 'extend':
    if (context.storage !== 'domain')
      throw new Error('u2: only a platform-stored schema is extended');
    return {editable: true, identifier: false, description: false, ...vocabulary};
  default:
    throw new Error(`u2: unknown manifest editor mode "${String(context.mode)}"`);
  }
}
