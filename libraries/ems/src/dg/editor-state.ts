import type * as DG from 'datagrok-api/dg';
import {signal, computed, ReadonlySignal, Emitter, FrameRows} from '@datagrok-libraries/u2';
import type {EditSaved, EditState} from '../core/edit-state.js';
import type {DataFrameLike} from '@datagrok-libraries/u2';

/** `EditState` over the js-api `DomainFrameEditor`: the editor stays the single writer of the frame, and its
 * `onChanged` is the one signal the mirrored state follows. Keys are what `FrameRows` uses. */
export class EditorEditState implements EditState {
  readonly isDirty: ReadonlySignal<boolean>;
  readonly changeCount: ReadonlySignal<number>;
  readonly validity: ReadonlySignal<string | null>;
  readonly blocking: ReadonlySignal<string | null>;
  readonly isSaving: ReadonlySignal<boolean>;
  readonly onChanged = new Emitter<string | null>();
  readonly onSaved = new Emitter<EditSaved>();
  /** The sentence the platform session worked out for the last refusal — see `DgDomainBackend`,
   * which is where it is captured, and `SharedSession._refused`, which says it. */
  problem: string | null = null;
  problemCode: string | null = null;

  private readonly _subs: {unsubscribe(): void}[];
  // bumped on every editor change: the validity walk runs when it is read, not per keystroke
  private readonly _version = signal(0);
  private _index: Map<string, number> | undefined;
  // the cells the server's last refusal marked (row key → column → what it said): shown as any
  // cell error, yet not the client's verdict — the platform marks them `kind: 'error'` as it does
  // its own validation, and only an edit of that cell would take them back
  private readonly _refused = new Map<string, Map<string, string>>();

  constructor(readonly editor: DG.DomainFrameEditor) {
    const dirty = signal(editor.isDirty);
    const count = signal(editor.changeCount);
    const saving = signal(editor.isSaving);
    this.isDirty = dirty;
    this.changeCount = count;
    this.validity = computed(() => {
      this._version.value;
      return this._validity(false);
    });
    this.blocking = computed(() => {
      this._version.value;
      return this._validity(true);
    });
    this.isSaving = saving;
    const df = editor.dataFrame;
    const drop = () => this._index = undefined;
    this._subs = [
      editor.onChanged.subscribe(() => {
        // a saved delete leaves the frame without a row event (`removeAt(row, 1, false)`)
        drop();
        this._prune();
        count.value = editor.changeCount;
        dirty.value = editor.isDirty;
        this._version.value = this._version.peek() + 1;
        this.onChanged.fire(null);
      }),
      editor.onDirtyChanged.subscribe((d) => dirty.value = d),
      editor.onSavingChanged.subscribe((s) => saving.value = s),
      editor.onSaved.subscribe((r) => {
        // an editor that reports no counters (an older js-api) wrote everything it had pending
        const written = r.inserted + r.updated + r.deleted;
        this.onSaved.fire(Number.isNaN(written) ? {assigned: r.assigned} : {assigned: r.assigned, written});
      }),
      df.onValuesChanged.subscribe(drop),
      df.onRowsAdded.subscribe(drop),
      df.onRowsRemoved.subscribe(drop),
    ];
  }

  get df(): DG.DataFrame {
    return this.editor.dataFrame;
  }

  /** The frame row behind a key, -1 when the frame no longer holds it. */
  indexOf(key: string): number {
    const df = this.df;
    if (this._index === undefined) {
      this._index = new Map();
      for (let i = 0; i < df.rowCount; i++)
        this._index.set(FrameRows.keyOf(df as unknown as DataFrameLike, i), i);
    }
    return this._index.get(key) ?? -1;
  }

  isChanged(key: string, column: string): boolean {
    const at = this.indexOf(key);
    return at >= 0 && this.editor.isChanged(at, column);
  }

  errorOf(key: string, column: string): string | null {
    const at = this.indexOf(key);
    return at < 0 ? null : this.editor.errorOf(at, column)?.message ?? null;
  }

  refused(key: string, column: string): boolean {
    this._version.value;
    const message = this._refused.get(key)?.get(column);
    return message !== undefined && message === this.errorOf(key, column);
  }

  /** Every blocking cell error standing now, by row key and column: what {@link markRefused}
   * compares against. */
  errorCells(): Map<string, Map<string, string>> {
    const found = new Map<string, Map<string, string>>();
    const df = this.df;
    for (let row = 0; row < df.rowCount; row++) {
      const errors = this.editor.errorsOf(row);
      for (const column of Object.keys(errors)) {
        if (errors[column].kind !== 'error')
          continue;
        const key = FrameRows.keyOf(df as unknown as DataFrameLike, row);
        let cells = found.get(key);
        if (cells === undefined)
          found.set(key, cells = new Map());
        cells.set(column, errors[column].message);
      }
    }
    return found;
  }

  /** After a save the server refused: the cell errors that were not there when the save began
   * ([before], from {@link errorCells}) are the server's, mapped onto the cells it named, and are
   * kept apart from the client's verdict. An error the client had already found is never one. */
  markRefused(before: Map<string, Map<string, string>>): void {
    this._refused.clear();
    for (const [key, cells] of this.errorCells()) {
      for (const [column, message] of cells) {
        if (before.get(key)?.get(column) === message)
          continue;
        let marked = this._refused.get(key);
        if (marked === undefined)
          this._refused.set(key, marked = new Map());
        marked.set(column, message);
      }
    }
    this._version.value = this._version.peek() + 1;
  }

  /** Before a save: the server's marks still standing are dropped, so the editor sends the batch
   * and the server decides again — re-marking what it still refuses. `validate()` re-runs the
   * client's own checks over the pending rows and takes back every error they do not find. */
  retry(): void {
    if (this._refused.size === 0)
      return;
    this._refused.clear();
    this.editor.validate();
  }

  setValue(key: string, column: string, value: unknown): void {
    const at = this.indexOf(key);
    if (at >= 0)
      this.editor.setValue(at, column, value);
  }

  newRow(values: Record<string, unknown> = {}, options?: {pristine?: boolean}): string {
    const at = this.editor.addRow(values, options);
    if (at < 0)
      throw new Error(`${this.editor.table}: cannot add a row while the batch is being saved`);
    return FrameRows.keyOf(this.df as unknown as DataFrameLike, at);
  }

  markDeleted(key: string): void {
    const at = this.indexOf(key);
    if (at >= 0)
      this.editor.markDeleted(at);
  }

  unmarkDeleted(key: string): void {
    const at = this.indexOf(key);
    if (at >= 0)
      this.editor.unmarkDeleted(at);
  }

  markRestored(key: string): void {
    const at = this.indexOf(key);
    if (at >= 0)
      this.editor.markRestored(at);
  }

  unmarkRestored(key: string): void {
    const at = this.indexOf(key);
    if (at >= 0)
      this.editor.unmarkRestored(at);
  }

  discard(): void {
    this.editor.discard();
  }

  setSaving(saving: boolean): void {
    this.editor.setSaving(saving);
  }

  rebind(assigned: Record<string, string>): void {
    this.editor.rebind(assigned);
  }

  save(): Promise<boolean> {
    return this.editor.save();
  }

  dispose(): void {
    for (const s of this._subs)
      s.unsubscribe();
    this.onChanged.clear();
    this.onSaved.clear();
    this.editor.detach();
  }

  // what the editor's own save refuses on; `client`: the server's marks left out
  private _validity(client: boolean): string | null {
    const editor = this.editor;
    const df = this.df;
    for (let row = 0; row < df.rowCount; row++) {
      const state = editor.stateOf(row) as string;
      if (state === 'deleted' || state === 'restored')
        continue;
      const errors = editor.errorsOf(row);
      const refused = client ? this._refused.get(FrameRows.keyOf(df as unknown as DataFrameLike, row)) : undefined;
      for (const column of Object.keys(errors)) {
        if (errors[column].kind === 'error' && refused?.get(column) !== errors[column].message)
          return errors[column].message;
      }
    }
    return null;
  }

  // a mark the editor has since taken back or replaced is no longer the server's: an equal
  // sentence the client finds later blocks
  private _prune(): void {
    for (const [key, cells] of this._refused) {
      for (const [column, message] of cells) {
        if (this.errorOf(key, column) !== message)
          cells.delete(column);
      }
      if (cells.size === 0)
        this._refused.delete(key);
    }
  }
}
