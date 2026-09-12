/* The permission protocol every control consults (GOAL "Access"): what the caller may do with a
   table and with each of its fields, refined per row where the rows carry the per-row access
   columns. Platform-free — the EMS server answers this exact shape, dg passes it through, and
   every non-EMS caller uses one of the two constants. */
import {Rows} from '../sources/rows-like.js';

export type FieldAccess = 'hidden' | 'readonly' | 'editable';
/** The five table capabilities, plus whatever custom `permissions` a schema declares. */
export type Capability = 'view' | 'insert' | 'edit' | 'delete' | 'share' | (string & {});

export interface AccessData {
  can: Record<string, boolean>;
  /** A field absent here is hidden — the server lists every column the caller may see. */
  fields: Record<string, FieldAccess>;
}

/** The two rules every control degrades by: permission ⇒ hidden, state ⇒ disabled; readonly is
 * text, never a dead input. {@link full} and {@link readOnly} answer every name; an access built
 * by {@link from} answers exactly what its data lists — an unlisted capability is denied, an
 * unlisted field hidden.
 *
 * The split with the server: `fields[col]` is column security only, `can.*` the table-level
 * row rights (false negatives on a row-mode table), and a row's `~can_<name>` columns the per-row
 * truth. A listed `editable` field means the column is not restricted for the caller, not that
 * the caller may write the row: {@link field} folds the write capability in — `edit` for an
 * existing row (narrowed per row by {@link row}), `insert` for a draft ({@link forDraft}) — and
 * answers `readonly` where it is denied. */
export class Access {
  /** The per-row service columns a `withAccess` query adds (js-api `DOMAIN_ACCESS_COLUMNS`), each
   * with the capability it carries — the documented server set; {@link row} reads any
   * `~can_<name>` beyond it. */
  static readonly ROW_COLUMNS: readonly (readonly [column: string, capability: Capability])[] =
    [['~can_edit', 'edit'], ['~can_delete', 'delete'], ['~can_share', 'share']];
  static readonly ROW_PREFIX = '~can_';

  static readonly full: Access = Object.freeze(new Access({}, {}, true, 'editable')) as Access;
  static readonly readOnly: Access = Object.freeze(new Access({view: true}, {}, false, 'readonly')) as Access;

  private constructor(private readonly _can: Record<string, boolean>,
    private readonly _fields: Record<string, FieldAccess>, private readonly _unlistedCan: boolean,
    private readonly _unlistedField: FieldAccess, private readonly _writes: Capability = 'edit') {}

  static from(data: AccessData): Access {
    return new Access({...data.can}, {...data.fields}, false, 'hidden');
  }

  can(capability: Capability): boolean {
    return this._can[capability] ?? this._unlistedCan;
  }

  /** `hidden` where the field is unlisted, `editable` where it is listed so AND the row may be
   * written, `readonly` otherwise. */
  field(name: string): FieldAccess {
    const listed = this._fields[name] ?? this._unlistedField;
    return listed === 'editable' && !this.can(this._writes) ? 'readonly' : listed;
  }

  /** Whether {@link field} gates on `insert` — the view of a draft. */
  get isDraft(): boolean {
    return this._writes === 'insert';
  }

  /** The same access seen from a draft: a field is editable under the `insert` capability, since
   * a row that does not exist yet has nothing to `edit`. */
  forDraft(): Access {
    return this.isDraft ? this :
      new Access(this._can, this._fields, this._unlistedCan, this._unlistedField, 'insert');
  }

  /** Row-level security: a `~can_<name>` column the row carries with a boolean is the server's
   * per-row truth for that capability and REPLACES the table's answer — a row-mode table's
   * table-level `edit`/`delete`/`share` are false negatives, so a row may say yes where the
   * table said no, and no where it said yes. A column that is absent or not a boolean (the
   * server sends `~can_share` as null off row mode, and omits it from a frame) is not carried
   * and answers the table's access; `insert` is never per row. A draft (the editor's `~state`
   * reads `'new'`) has no row truth yet — its `~can_*` cells are the frame's defaults — and
   * answers {@link forDraft}. */
  row(row: unknown): Access {
    const r = row as Record<string, unknown> | null | undefined;
    if (r === null || r === undefined)
      return this;
    if (r[Rows.STATE] === 'new')
      return this.forDraft();
    // rows are read by name (a row proxy never enumerates its service columns), so the
    // capabilities looked for are the table's own plus the documented server set
    const names = new Set([...Object.keys(this._can), ...Access.ROW_COLUMNS.map(([, c]) => c)]);
    let can: Record<string, boolean> | undefined;
    for (const capability of names) {
      const value = r[`${Access.ROW_PREFIX}${capability}`];
      if (typeof value !== 'boolean')
        continue;
      can ??= {...this._can};
      can[capability] = value;
    }
    return can === undefined ? this :
      new Access(can, this._fields, this._unlistedCan, this._unlistedField, this._writes);
  }
}
