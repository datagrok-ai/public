/**
 * {@link DomainSession} — the unit of work over several {@link DomainFrameEditor}s:
 * ONE `/transaction` for all of their pending batches, so a draft row of one
 * editor can be referenced by a row of another before either exists, and a
 * parent form plus its child grids save (or fail) together.
 *
 * The editors stay per table; the session is an aggregator that concatenates
 * their ops, posts them once, and hands each editor its slice of the results
 * through the save participant protocol (`prepareSave` / `applyResults` /
 * `resolveConflict` / `mapValidationError` / `setSaving`). An
 * editor's own `save()` is a session of one.
 *
 * @module domains-session
 */

import * as rxjs from 'rxjs';

import {IDartApi} from '../../api/grok_api.g';
import {DomainsDataSource} from '../../dapi';
import {DomainTransactionOp, DomainValidationError, DomainVersionConflictError} from '../../domains';
import {Balloon} from '../../widgets/menu';
import {DomainFrameEditor} from './domains-editor';
import type {DomainPendingOp, DomainSaveResult, IEditorHost} from './domains-editor';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;
const balloon = new Balloon();

function domains(): DomainsDataSource {
  return new DomainsDataSource(api.grok_Dapi_Domains());
}

/** How many conflict round trips one {@link DomainSession.save} tolerates. Each
 * outcome resolves exactly one row, so the chain is finite; the cap only guards
 * against a server that keeps reporting the same op. */
const CONFLICT_RETRY_LIMIT = 32;

export interface DomainSessionOptions {
  /** Suppresses the 'Saved N rows' balloon; by default the session is quiet
   * when every editor is (see `DomainFrameEditorOptions.quiet`). */
  quiet?: boolean;
  /** URL schema of the transaction; the first editor's by default. Ops of an
   * editor on another schema are addressed as `<schema>.<table>`. */
  schema?: string;
}

/** One editor's share of a {@link DomainSession.save}. */
export interface DomainSessionPart {
  editor: DomainFrameEditor;
  pending: DomainPendingOp[];
}

/** Why a {@link DomainSession.save} did not land — see {@link DomainSession.lastRefusal}. */
export interface DomainSessionRefusal {
  /** The sentence the save reported; the same one the cell says where the refusal reaches
   * cells. Bare — a surface that leads with its own "Cannot save:" (the balloon does) prefixes
   * it itself. */
  message: string;
  /** The editor that owns the failing row; null for a refusal the batch as a whole
   * answers for (nothing located it, or no single participant caused it). */
  editor: DomainFrameEditor | null;
  /** The transaction-wide index of the failing op, where the server named one. */
  opIndex?: number;
}

/**
 * A set of editors that save as ONE transaction.
 *
 * ```ts
 * const session = new DG.DomainSession([parent, children]);
 * const draft = parent.addRow({name: 'Study 7'});         // id = '~new:…'
 * children.addRow({study_id: parent.dataFrame.get('id', draft)});
 * await session.save();                                     // one /transaction, the ref resolved server-side
 * ```
 *
 * `isDirty` / `changeCount` / `isSaving` aggregate over the editors; the
 * signals mirror theirs. {@link editors} is live — {@link add} / {@link remove}
 * resubscribe — so a page can grow its session as widgets appear.
 */
export class DomainSession implements IEditorHost {
  /** What {@link lastRefusal} says when an editor refused its own batch: the blocking cells
   * name themselves, and `prepareSave` has already ballooned the first of them. */
  static readonly CELL_ERRORS = 'Fix the cell errors first';

  readonly editors: DomainFrameEditor[] = [];

  private readonly _subs = new Map<DomainFrameEditor, rxjs.Subscription[]>();
  private readonly _quiet?: boolean;
  private readonly _schema?: string;
  private _dirty = false;
  private _saving = false;
  private _lastRefusal: DomainSessionRefusal | null = null;

  private readonly _onChanged = new rxjs.Subject<DomainSession>();
  private readonly _onDirtyChanged = new rxjs.Subject<boolean>();
  private readonly _onSavingChanged = new rxjs.Subject<boolean>();
  private readonly _onSaved = new rxjs.Subject<DomainSaveResult>();

  constructor(editors?: DomainFrameEditor[], options?: DomainSessionOptions) {
    this._quiet = options?.quiet;
    this._schema = options?.schema;
    for (const editor of editors ?? [])
      this.add(editor);
  }

  /** The transaction's URL schema (see {@link DomainSessionOptions.schema}). */
  get schema(): string | undefined { return this._schema ?? this.editors[0]?.client.schema; }

  get quiet(): boolean { return this._quiet ?? this.editors.every((e) => e.quiet); }

  /** Whether any editor has a pending change. */
  get isDirty(): boolean { return this.editors.some((e) => e.isDirty); }

  /** Sum of the editors' pending cell changes. */
  get changeCount(): number { return this.editors.reduce((n, e) => n + e.changeCount, 0); }

  /** Whether a transaction is in flight (every editor is closed while it is). */
  get isSaving(): boolean { return this.editors.some((e) => e.isSaving); }

  /** Why the last {@link save} did not land; null after a successful save or before the first
   * one. A host that shows the refusal in its own words (a status line, a form footer) reads the
   * sentence here instead of intercepting the balloon. */
  get lastRefusal(): DomainSessionRefusal | null { return this._lastRefusal; }

  /** Fires on every service-state write of any editor. */
  get onChanged(): rxjs.Observable<DomainSession> { return this._onChanged; }
  /** Fires when the aggregate {@link isDirty} flips. */
  get onDirtyChanged(): rxjs.Observable<boolean> { return this._onDirtyChanged; }
  /** Fires when the aggregate {@link isSaving} flips. */
  get onSavingChanged(): rxjs.Observable<boolean> { return this._onSavingChanged; }
  /** Fires after a successful {@link save} with the counters summed over the
   * editors and every draft id the server assigned. */
  get onSaved(): rxjs.Observable<DomainSaveResult> { return this._onSaved; }

  add(editor: DomainFrameEditor): void {
    if (this.editors.includes(editor))
      return;
    this.editors.push(editor);
    this._subs.set(editor, [
      editor.onChanged.subscribe(() => this._fire()),
      editor.onSavingChanged.subscribe(() => this._fireSaving()),
    ]);
    this._fire();
    this._fireSaving();
  }

  remove(editor: DomainFrameEditor): void {
    const index = this.editors.indexOf(editor);
    if (index < 0)
      return;
    for (const sub of this._subs.get(editor) ?? [])
      sub.unsubscribe();
    this._subs.delete(editor);
    this.editors.splice(index, 1);
    this._fire();
    this._fireSaving();
  }

  /** Every editor's batch (`prepareSave`), in editor order — what {@link save}
   * concatenates; null when any editor refused (a blocking cell error, already
   * reported). */
  buildOps(): DomainSessionPart[] | null {
    const parts: DomainSessionPart[] = [];
    for (const editor of this.editors) {
      const pending = editor.prepareSave();
      if (pending == null)
        return null;
      parts.push({editor: editor, pending: pending});
    }
    return parts;
  }

  /**
   * Writes every editor's pending batch as ONE `/transaction`: audit rows share
   * a `tx_id`, and any failure rolls all of it back. Resolves to whether the
   * batch landed.
   *
   * Blocking cell errors refuse the save (naming the first one). A version
   * conflict goes through the platform's standard reload/overwrite dialog on
   * the editor that owns the failing row, and the WHOLE batch is rebuilt and
   * retried. A server validation error — a business-key clash included, which the
   * `onDuplicate: 'error'` of every insert makes an atomic refusal — lands on the
   * offending editor's cells and everything stays pending. Every editor is
   * CLOSED while this runs (see
   * `DomainFrameEditor.isSaving`). Every unsuccessful return leaves its reason in
   * {@link lastRefusal}.
   */
  async save(): Promise<boolean> {
    if (this.isSaving) {
      const busy = 'The batch is already being saved.';
      balloon.warning(busy);
      return this._refuse(busy);
    }
    let parts = this.buildOps();
    if (parts == null)
      return this._refuse(DomainSession.CELL_ERRORS);
    const overlapping = DomainSession._overlapping(parts);
    if (overlapping != null) {
      balloon.error(overlapping);
      return this._refuse(overlapping);
    }
    if (parts.every((p) => p.pending.length === 0)) {
      this._lastRefusal = null;
      return true;
    }
    const schema = this.schema!;
    for (const editor of this.editors)
      editor.setSaving(true);
    try {
      for (let attempt = 0; ; attempt++) {
        if (attempt > CONFLICT_RETRY_LIMIT) {
          balloon.error('Cannot save: too many version conflicts in a row');
          return this._refuse('Too many version conflicts in a row');
        }
        let results: any[];
        try {
          results = await domains().transaction(schema, this._ops(parts, schema));
        } catch (e: any) {
          const at = this._locate(parts, e?.opIndex);
          if (e instanceof DomainVersionConflictError && at != null) {
            if (!(await at.editor.resolveConflict(e, at.failing)))
              return this._refuse(e?.message ?? `${e}`, at.editor, e?.opIndex);
            parts = this.buildOps();
            if (parts == null)
              return this._refuse(DomainSession.CELL_ERRORS);
            continue;
          }
          if (at == null) {
            balloon.error(e?.message ?? `${e}`);
            return this._refuse(e?.message ?? `${e}`, null, e?.opIndex);
          }
          // ONE sentence for all four surfaces: the cells, the status line that reads them,
          // the balloon, and `lastRefusal` for a host that words it itself
          const message = await at.editor.refusalFor(e, at.failing);
          if (e instanceof DomainValidationError)
            at.editor.mapValidationError(e, at.failing, message);
          balloon.error(message);
          return this._refuse(message, at.editor, e?.opIndex);
        }
        const total: DomainSaveResult = {inserted: 0, updated: 0, deleted: 0, assigned: Object.create(null)};
        const slices: any[][] = [];
        let offset = 0;
        for (const part of parts) {
          slices.push(results.slice(offset, offset + part.pending.length));
          offset += part.pending.length;
        }
        // Read from the WHOLE transaction before any editor applies its slice: a
        // frame holding another editor's draft id resolves it from here, so the
        // post-save re-read only enriches and may fail without leaving a dangling ref.
        const assigned: {[draftId: string]: string} = Object.create(null);
        for (let i = 0; i < parts.length; i++)
          Object.assign(assigned, DomainFrameEditor.assignedOf(parts[i].pending, slices[i]));
        for (let i = 0; i < parts.length; i++) {
          const r = await parts[i].editor.applyResults(parts[i].pending, slices[i], assigned);
          total.inserted += r.inserted;
          total.updated += r.updated;
          total.deleted += r.deleted;
          Object.assign(total.assigned, r.assigned);
        }
        this._lastRefusal = null;
        this._onSaved.next(total);
        const n = total.inserted + total.updated + total.deleted;
        if (!this.quiet)
          balloon.info(`Saved ${n} row${n === 1 ? '' : 's'}`);
        return true;
      }
    } finally {
      for (const editor of this.editors)
        editor.setSaving(false);
    }
  }

  /** Releases every editor: the session's subscriptions go, the editors stay. An ad-hoc session —
   * {@link DomainFrameEditor.save}, one `saveAll` — MUST be disposed when its transaction is done,
   * or its subscriptions accumulate on the editors for the life of the frame. */
  dispose(): void {
    for (const editor of [...this.editors])
      this.remove(editor);
  }

  /** Drops every editor's pending batch. */
  discard(): void {
    for (const editor of this.editors)
      editor.discard();
  }

  private _refuse(message: string, editor?: DomainFrameEditor | null, opIndex?: unknown): false {
    this._lastRefusal = {message: message, editor: editor ?? null,
      opIndex: typeof opIndex === 'number' ? opIndex : undefined};
    return false;
  }

  /** The same persisted row addressed by two participants: the transaction would carry two ops
   * for it, the second losing to the first's version check and rolling everything back. Named
   * and refused before anything is sent — a host with its own vocabulary (u2's `SharedSession`)
   * says it in the row's own words above this. Null when the batch is disjoint. */
  private static _overlapping(parts: DomainSessionPart[]): string | null {
    const seen = new Set<string>();
    for (const {editor, pending} of parts)
      for (const {op} of pending) {
        if (op.id == null)
          continue;
        const address = `${editor.client.schema}.${op.table}`;
        const key = `${address}/${op.id}`;
        if (seen.has(key))
          return `${address} ${op.id} is edited in two places`
            + ' — discard one of the two edits and save again';
        seen.add(key);
      }
    return null;
  }

  /** The concatenated ops, an editor on another schema addressed as `<schema>.<table>`. */
  private _ops(parts: DomainSessionPart[], schema: string): DomainTransactionOp[] {
    const ops: DomainTransactionOp[] = [];
    for (const {editor, pending} of parts)
      for (const p of pending)
        ops.push(editor.client.schema === schema ? p.op
          : {...p.op, table: `${editor.client.schema}.${p.op.table}`});
    return ops;
  }

  /** The editor and pending op a transaction-wide `opIndex` points at. */
  private _locate(parts: DomainSessionPart[], index: unknown):
      {editor: DomainFrameEditor; failing: DomainPendingOp} | null {
    if (typeof index !== 'number' || index < 0)
      return null;
    let offset = 0;
    for (const part of parts) {
      if (index < offset + part.pending.length)
        return {editor: part.editor, failing: part.pending[index - offset]};
      offset += part.pending.length;
    }
    return null;
  }

  private _fire(): void {
    this._onChanged.next(this);
    const dirty = this.isDirty;
    if (dirty !== this._dirty) {
      this._dirty = dirty;
      this._onDirtyChanged.next(dirty);
    }
  }

  private _fireSaving(): void {
    const saving = this.isSaving;
    if (saving !== this._saving) {
      this._saving = saving;
      this._onSavingChanged.next(saving);
    }
  }
}
