/* `EditState` over the js-api `DomainFrameEditor` — the third host of the editing state
   (STATE-CONTRACT "u2 DomainSource"): the editor stays the single writer of the frame (H3), every
   write here goes through it by row index, and its `onChanged` is the one signal the mirrored
   state follows (H10). Keys are what `FrameRows` uses — the `id` cell, the draft key for a row
   that has no id yet — so the source, the form and the list address a row the same way. */
import type * as DG from 'datagrok-api/dg';
import {signal, computed, ReadonlySignal} from '../../core/signals.js';
import {Emitter} from '../../core/emitter.js';
import type {EditState} from '../../sources/edit-state.js';
import {FrameRows} from '../../sources/df-rows.js';
import {Rows} from '../../sources/rows-like.js';
import type {DataFrameLike} from '../../sources/df-bindings.js';

export class EditorEditState implements EditState {
  readonly isDirty: ReadonlySignal<boolean>;
  readonly changeCount: ReadonlySignal<number>;
  readonly validity: ReadonlySignal<string | null>;
  readonly isSaving: ReadonlySignal<boolean>;
  readonly onChanged = new Emitter<string | null>();

  private readonly _subs: {unsubscribe(): void}[];
  /** Bumped on every editor change: the validity walk runs when it is read, not per keystroke. */
  private readonly _version = signal(0);
  /** id → row index, built on demand and dropped whenever the frame's cells or rows move. */
  private _index: Map<string, number> | undefined;

  constructor(readonly editor: DG.DomainFrameEditor) {
    const dirty = signal(editor.isDirty);
    const count = signal(editor.changeCount);
    const saving = signal(editor.isSaving);
    this.isDirty = dirty;
    this.changeCount = count;
    this.validity = computed(() => {
      this._version.value;
      return this._validity();
    });
    this.isSaving = saving;
    const df = editor.dataFrame;
    const drop = () => this._index = undefined;
    this._subs = [
      editor.onChanged.subscribe(() => {
        count.value = editor.changeCount;
        dirty.value = editor.isDirty;
        this._version.value = this._version.peek() + 1;
        this.onChanged.fire(null);
      }),
      editor.onDirtyChanged.subscribe((d) => dirty.value = d),
      editor.onSavingChanged.subscribe((s) => saving.value = s),
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
    const draft = Rows.draftRow(key, df.rowCount);
    if (draft !== null)
      return draft;
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

  discard(): void {
    this.editor.discard();
  }

  save(): Promise<boolean> {
    return this.editor.save();
  }

  dispose(): void {
    for (const s of this._subs)
      s.unsubscribe();
    this.onChanged.clear();
    this.editor.detach();
  }

  /** The first blocking problem on a row that is not deleted — what the editor's own save refuses on. */
  private _validity(): string | null {
    const editor = this.editor;
    for (let row = 0; row < this.df.rowCount; row++) {
      if (editor.stateOf(row) === 'deleted')
        continue;
      const errors = editor.errorsOf(row);
      for (const column of Object.keys(errors)) {
        if (errors[column].kind === 'error')
          return errors[column].message;
      }
    }
    return null;
  }
}
