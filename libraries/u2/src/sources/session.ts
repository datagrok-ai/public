/* The unit of work a Save button drives (GOAL "Unit of work"): what is pending across the sources
   it covers, one `save` for all of it — every writer's batch as ONE transaction through
   `DomainBackend.saveAll` — one `discard`, and the events a paired form follows. A source alone
   is a session of one; the sources built for one spec instance or one app share the ambient
   session (`SharedSession.runWith`), and a hand-built page passes `{session}` explicitly. */
import {signal, computed, ReadonlySignal} from '../core/signals.js';
import {Emitter} from '../core/emitter.js';
import type {ObservableLike} from '../core/widget-like.js';
import {notify} from '../components/display/notify.js';
import {Dialog} from '../components/containers/dialog.js';
import {backends} from './backends.js';
import {Rows} from './rows-like.js';
import type {RowView} from './rows-like.js';
import type {DomainSource} from './domain-source.js';

export interface DomainSession {
  readonly isDirty: ReadonlySignal<boolean>;
  readonly changeCount: ReadonlySignal<number>;
  /** The first blocking problem, null when the batch may be saved. */
  readonly validity: ReadonlySignal<string | null>;
  readonly isSaving: ReadonlySignal<boolean>;
  /** A batch landed — the paired form returns the focus to its fields. */
  readonly onSaved: ObservableLike<void>;
  readonly onDiscarded: ObservableLike<void>;
  /** Every pending change as one transaction; resolves to whether it landed. A refusal — the
   * form's own gate, the writer's validity, the server — is the source's `error`. */
  save(): Promise<boolean>;
  discard(): void;
  /** A source joins at construction; the unregister runs at its disposal. */
  add?(source: DomainSource): () => void;
  /** The Save button of the page, when one registered — where Tab from a form's last field goes. */
  primaryButton?: HTMLElement;
}

export class SharedSession implements DomainSession {
  /** The session every source built while {@link runWith} runs joins — set per spec instance
   * and per app; a source built outside any gets a session of its own. */
  static ambient: SharedSession | undefined;

  readonly sources: ReadonlySignal<readonly DomainSource[]>;
  readonly isDirty: ReadonlySignal<boolean>;
  readonly changeCount: ReadonlySignal<number>;
  readonly validity: ReadonlySignal<string | null>;
  readonly isSaving: ReadonlySignal<boolean>;
  /** `''` when clean, `'N unsaved changes'`, `' in M tables'` appended when several are dirty. */
  readonly summary: ReadonlySignal<string>;
  readonly onSaved = new Emitter<void>();
  readonly onDiscarded = new Emitter<void>();
  primaryButton: HTMLElement | undefined;

  private readonly _sources = signal<readonly DomainSource[]>([]);
  private readonly _saving = signal(false);

  constructor() {
    this.sources = this._sources;
    this.isDirty = computed(() => this._sources.value.some((s) => s.isDirty.value));
    this.changeCount = computed(() => this._sources.value.reduce((n, s) => n + s.changeCount.value, 0));
    this.validity = computed(() => this._sources.value.map((s) => s.validity.value).find((v) => v !== null) ?? null);
    this.isSaving = computed(() => this._saving.value || this._sources.value.some((s) => s.isSaving.value));
    this.summary = computed(() => {
      const changes = this.changeCount.value;
      if (changes === 0)
        return '';
      const tables = new Set(this._sources.value.filter((s) => s.isDirty.value).map((s) => s.table)).size;
      return `${changes} unsaved change${changes === 1 ? '' : 's'}${tables > 1 ? ` in ${tables} tables` : ''}`;
    });
  }

  static runWith<T>(session: SharedSession, fn: () => T): T {
    const outer = SharedSession.ambient;
    SharedSession.ambient = session;
    try {
      return fn();
    } finally {
      SharedSession.ambient = outer;
    }
  }

  add(source: DomainSource): () => void {
    this._sources.value = [...this._sources.peek(), source];
    return () => this._sources.value = this._sources.peek().filter((s) => s !== source);
  }

  /** Every dirty source checked first (a problem is that source's `error`), then one `saveAll`;
   * the one "saved" balloon: "<Singular> saved", "N changes saved", "N changes saved in M tables". */
  async save(): Promise<boolean> {
    if (this.isSaving.peek())
      return false;
    const dirty = this._sources.peek().filter((s) => s.isDirty.peek());
    if (dirty.length === 0)
      return false;
    this._saving.value = true;
    try {
      const changes = this.changeCount.peek();
      const batch = this._withReferenced(dirty);
      for (const source of batch) {
        if (source.check() !== null)
          return false;
      }
      if (SharedSession._overlapping(batch))
        return false;
      // every draft id the transaction assigned, over the whole batch: a source's query or
      // defaults may name one (a child collection under a draft parent), and what it reads next
      // must be about the id that parent was given
      const assigned: Record<string, string> = {};
      const heard = batch.map((s) => s.edit.peek()!.onSaved.subscribe((r) => Object.assign(assigned, r.assigned)));
      try {
        if (!await backends.domain!.saveAll(batch.map((s) => s.edit.peek()!))) {
          SharedSession._refused(batch);
          return false;
        }
      } catch (e) {
        for (const source of batch)
          source.fail(e);
        return false;
      } finally {
        for (const sub of heard)
          sub.unsubscribe();
      }
      const rebased = new Set(batch);
      for (const source of this._sources.peek()) {
        if (source.rebind(assigned) && !rebased.has(source))
          source.refresh().catch((e) => source.fail(e));
      }
      for (const source of batch)
        await source.afterSave();
      const what = batch[0].schema.info.singularName || 'Row';
      const tables = new Set(batch.map((s) => s.table)).size;
      notify.info(tables > 1 ? `${changes} change${changes === 1 ? '' : 's'} saved in ${tables} tables` :
        changes > 1 ? `${changes} changes saved` : `${what.charAt(0).toUpperCase()}${what.slice(1)} saved`);
      this.onSaved.fire();
      return true;
    } finally {
      this._saving.value = false;
    }
  }

  discard(): void {
    for (const source of this._sources.peek())
      source.revert();
    this.onDiscarded.fire();
  }

  /** A backend that refused without throwing (the platform editor reports the server's error
   * itself and answers false): the summary still has to say the changes are stuck, and the
   * column errors the editor mapped onto the cells name why. */
  private static _refused(batch: readonly DomainSource[]): void {
    const problem = batch.map((s) => s.validity.peek()).find((v) => v !== null) ??
      'the changes were refused';
    for (const source of batch)
      source.refuse(problem);
  }

  /** The save set: the dirty sources, plus every source holding a draft one of them refers to — a
   * pristine parent created with defaults only is not dirty, yet its child's `$~new:` needs it.
   * Over the pending rows, not the visible ones: what the writer sends is the frame, and a parent
   * the frame's filter hides is still the row the child points at. */
  private _withReferenced(dirty: readonly DomainSource[]): DomainSource[] {
    const batch = [...dirty];
    const referenced = new Set<string>();
    for (let i = 0; i < batch.length; i++) {
      for (const row of batch[i].pending()) {
        for (const value of Object.values(row)) {
          for (const v of Array.isArray(value) ? value : [value]) {
            if (typeof v === 'string' && Rows.isDraft(v))
              referenced.add(v);
          }
        }
      }
      for (const source of this._sources.peek()) {
        if (!batch.includes(source) && source.pending().some((r) => Rows.isDraft(r) && referenced.has(r.id)))
          batch.push(source);
      }
    }
    return batch;
  }

  /** The same persisted row pending in two sources of the batch: they would send two ops for it,
   * the second losing to the first's version check, and the whole transaction would roll back.
   * Refused before anything is sent, naming the row on both sources — never merged, however
   * disjoint the columns are. */
  private static _overlapping(batch: readonly DomainSource[]): boolean {
    const seen = new Map<string, {source: DomainSource, row: RowView}>();
    for (const source of batch) {
      for (const row of source.pending()) {
        if (Rows.isDraft(row))
          continue;
        const key = `${source.table}/${row.id}`;
        const first = seen.get(key);
        if (first === undefined) {
          seen.set(key, {source, row});
          continue;
        }
        const problem = `${SharedSession._nameOf(first.source, first.row)} is edited in two places` +
          ' — discard one of the two edits and save again';
        first.source.refuseRow(first.row.id, problem);
        source.refuseRow(row.id, problem);
        return true;
      }
    }
    return false;
  }

  /** How a refusal names a row without a renderer: the table's singular name and the row's own
   * name column (its business key, else its id). */
  private static _nameOf(source: DomainSource, row: RowView): string {
    const info = source.schema.info;
    const column = info.nameColumn ?? info.businessKey[0];
    const name = column === undefined ? undefined : row[column];
    const what = info.singularName || source.table;
    return `${what.charAt(0).toUpperCase()}${what.slice(1)} "${name ?? row.id}"`;
  }
}

/** THE gate in front of whatever would drop pending changes — a navigation, a filter change, a
 * view close: resolves to whether the caller may proceed. Clean → true without a dialog; a save
 * in flight → a warning and false; else the user saves (false when the save is refused),
 * discards, or cancels. */
export function confirmDiscard(session: DomainSession, options: {action?: string, subject?: string} = {}):
  Promise<boolean> {
  // saving first: a batch in flight is clean for a while before its write-back lands, and the
  // caller must not be let through that window either
  if (session.isSaving.peek()) {
    notify.warning('Wait for the batch being saved to finish.');
    return Promise.resolve(false);
  }
  if (!session.isDirty.peek())
    return Promise.resolve(true);
  const changes = session.changeCount.peek();
  const action = options.action ?? 'continue';
  const subject = options.subject ?? 'this view';
  return new Promise((resolve) => {
    const text = document.createElement('div');
    const them = changes === 1 ? 'it' : 'them';
    for (const line of [`${changes} unsaved change${changes === 1 ? '' : 's'} in ${subject}.`,
      `Save ${them}, discard ${them}, or cancel and do not ${action}.`]) {
      const p = document.createElement('p');
      p.textContent = line;
      text.append(p);
    }
    let decided = false;
    const dialog = Dialog.create('Unsaved changes', {name: 'unsaved-changes'});
    const decide = (answer: boolean | Promise<boolean>) => {
      if (decided)
        return;
      decided = true;
      dialog.dispose();
      resolve(answer);
    };
    dialog.add(text)
      .addButton('SAVE', () => decide(session.save()), {primary: true})
      .addButton('DISCARD', () => {
        session.discard();
        decide(true);
      })
      .onCancel(() => decide(false))
      .show({modal: true});
  });
}
