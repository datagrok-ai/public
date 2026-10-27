/* The two-part seam the editor's field offer follows (Astra 3): the manifest STORAGE decides the
   vocabulary — what the parser accepts for a table and a column — and the editor CONTEXT decides
   what may be changed here: create, edit, view, extend. Nothing in the editor branches on who
   made the schema. Only the arms this feature ships are built; the others refuse by name so a
   later feature fills them in rather than inheriting a wrong offer silently. What edit mode locks
   per item — the logical names and remote mappings of registered tables and columns — the model
   decides by baseline membership; the offer holds only the schema-wide locks. Extend mode locks
   per item by ownership: the plugin's tables and columns are read, the user's own are edited. */

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
  /** A column's description, range, format and semantic type. */
  attributes: boolean;
  /** Columns are declared here — added and removed — rather than picked from a catalog. */
  columns: boolean;
  /** The access grids and the column visibility. */
  access: boolean;
  /** A user's column on a plugin table stays optional, non-unique and never the name column; the
   * panel shows those switches off, with the reason. */
  optionalOnly: boolean;
}

export function fieldOffer(context: EditorContext): FieldOffer {
  let vocabulary: Omit<FieldOffer, 'editable' | 'identifier' | 'description' | 'optionalOnly'>;
  switch (context.storage) {
  case 'external':
    vocabulary = {type: 'readonly', refs: 'relations', remoteNames: true, writable: true, required: true,
      nameColumn: true, searchable: true, defaultValue: false, autoNumber: false, immutable: false,
      unique: false, choices: false, attributes: false, columns: false, access: true};
    break;
  case 'domain':
    vocabulary = {type: 'editable', refs: 'editable', remoteNames: false, writable: false, required: true,
      nameColumn: true, searchable: true, defaultValue: true, autoNumber: true, immutable: true,
      unique: true, choices: true, attributes: true, columns: true, access: true};
    break;
  default:
    throw new Error(`u2: unknown manifest storage "${String(context.storage)}"`);
  }
  switch (context.mode) {
  case 'create':
  case 'edit':
    if (context.storage === 'domain')
      throw new Error(`u2: the manifest editor's "${context.mode}" mode over a platform-stored schema is not built yet`);
    return {editable: true, identifier: context.mode === 'create', description: context.mode === 'edit',
      optionalOnly: false, ...vocabulary};
  case 'view':
    return {editable: false, identifier: false, description: false, optionalOnly: false, ...vocabulary};
  case 'extend':
    if (context.storage !== 'domain')
      throw new Error('u2: only a platform-stored schema is extended — an external binding is edited');
    // the server refuses a required, unique or naming user column, a ref is not offered here, and the
    // apply's extend section carries no access
    return {...vocabulary, editable: true, identifier: false, description: false, refs: 'relations',
      required: false, nameColumn: false, searchable: false, autoNumber: false, immutable: false, unique: false,
      access: false, optionalOnly: true};
  default:
    throw new Error(`u2: unknown manifest editor mode "${String(context.mode)}"`);
  }
}
