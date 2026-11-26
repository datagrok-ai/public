import {emsBackends, capital, nounOf, errorCode, hierarchyCycle, DomainBackendError} from './domain-backend.js';
import {signal, computed, ReadonlySignal, Emitter, plural, notify, Dialog, Rows, div} from '@datagrok-libraries/u2';
import type {ObservableLike, SpecInstance, RowView} from '@datagrok-libraries/u2';
import type {DomainSource} from './domain-source.js';

export interface DomainSession {
  readonly isDirty: ReadonlySignal<boolean>;
  readonly changeCount: ReadonlySignal<number>;
  /** Why Save is off: the first source's {@link DomainSource.validity}, null when the batch may be sent. */
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
  add(source: DomainSource): () => void;
  /** The Save button of the page, when one registered — where Tab from a form's last field goes. */
  primaryButton?: HTMLElement;
}

const specSessions = new WeakMap<SpecInstance, SharedSession>();

/** The unit of work a Save button drives: every writer's batch as ONE transaction through `DomainBackend.saveAll`.
 * A source alone is a session of one; the sources of one spec instance or one app share the ambient session
 * (`SharedSession.runWith`), and a hand-built page passes `{session}` explicitly. */
export class SharedSession implements DomainSession {
  /** The session every source built while {@link runWith} runs joins — set per spec instance
   * and per app; a source built outside any gets a session of its own. */
  static ambient: SharedSession | undefined;

  /** The session the sources of one rendered spec share — none unless the spec's registry has ems
   * installed ({@link sessionAmbient}). */
  static of(instance: SpecInstance): SharedSession | undefined {
    return specSessions.get(instance);
  }

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
    // a source with nothing pending is not in the batch: what it shows holds no Save
    this.validity = computed(() => this._sources.value.filter((s) => s.isDirty.value).map((s) => s.validity.value)
      .find((v) => v !== null) ?? null);
    this.isSaving = computed(() => this._saving.value || this._sources.value.some((s) => s.isSaving.value));
    this.summary = computed(() => {
      const changes = this.changeCount.value;
      if (changes === 0)
        return '';
      const tables = new Set(this._sources.value.filter((s) => s.isDirty.value).map((s) => s.table)).size;
      return `${plural(changes, 'unsaved change', 'unsaved changes')}${tables > 1 ? ` in ${tables} tables` : ''}`;
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
    if (!this.isDirty.peek())
      return false;
    this._saving.value = true;
    try {
      let dirty: DomainSource[];
      let batch: DomainSource[];
      // every source is checked (a check takes back an earlier refusal), and an edit made meanwhile again
      for (;;) {
        dirty = this._sources.peek().filter((s) => s.isDirty.peek());
        batch = this._withReferenced(dirty);
        const seen = batch.map((s) => s.revision.peek());
        if (batch.map((s) => s.check()).some((problem) => problem !== null) || this._overlapping(batch) ||
            await this._cycle(batch).catch((e) => (batch.forEach((s) => s.fail(e)), true)))
          return false;
        const now = this._sources.peek().filter((s) => s.isDirty.peek());
        if (now.length === dirty.length && now.every((s) => dirty.includes(s)) &&
            batch.every((s, i) => s.revision.peek() === seen[i]))
          break;
      }
      // what the batch DOES, read before it lands and empties the pending rows: a batch of deletes
      // alone is not "saved", it is deleted, and one of restores is restored — the balloon has to
      // say so
      const states = new Set(batch.flatMap((s) => s.pending().map((r) => r[Rows.STATE])));
      const verb = states.size !== 1 ? 'saved' : states.has('deleted') ? 'deleted' :
        states.has('restored') ? 'restored' : 'saved';

      // dirty, yet no row carries a change: refused, never announced as saved
      if (batch.every((s) => s.pending().length === 0)) {
        SharedSession._unwritten(dirty);
        return false;
      }
      // every draft id the transaction assigned, over the whole batch: a source's query, defaults
      // or cells may name one (a child collection under a draft parent), and what it reads next
      // must be about the id that parent was given. No prototype: what is read against this map is
      // any string a cell or a literal may hold, `constructor` among them
      const assigned: Record<string, string> = Object.create(null);
      const before = new Map(batch.map((s) => [s, {changes: s.changeCount.peek(), rows: s.pending().length}]));
      const written = new Map<DomainSource, number>();
      const heard = batch.map((s) => s.edit.peek()!.onSaved.subscribe((r) => {
        Object.assign(assigned, r.assigned);
        written.set(s, (written.get(s) ?? 0) + r.written);
      }));
      try {
        await emsBackends.domain!.saveAll(batch.map((s) => s.edit.peek()!));
      } catch (e) {
        for (const source of batch)
          source.fail(e);
        // every control over those tables follows the re-read access, not only the sources that sent the batch
        if (errorCode(e) === 'forbidden') {
          emsBackends.domain!.dropAccess();
          await Promise.all(this._sources.peek().filter((s) => batch.some((b) => b.table === s.table))
            .map((s) => s.reloadAccess()));
        }
        return false;
      } finally {
        for (const sub of heard)
          sub.unsubscribe();
      }
      // the backend re-opened its writers when the transaction landed, and the re-base below
      // still has to swap the frame: an edit made in that window would be swapped away with it,
      // so the writers stay closed — refused, as an edit during the transaction is
      SharedSession._hold(batch, true);
      try {
        const rebased = new Set(batch);
        for (const source of this._sources.peek()) {
          if (source.rebind(assigned) && !rebased.has(source))
            source.refresh().catch((e) => source.fail(e));
        }
        // the batch has landed: a re-read that fails afterwards is that source's news to carry,
        // never a failed save
        for (const source of batch) {
          try {
            await source.afterSave();
          } catch (e) {
            source.markStale(e);
          }
        }
      } finally {
        SharedSession._hold(batch, false);
      }
      // a row the conflict dialog reloaded is not written; a writer that built no op for its changes still holds them
      const unwritten = batch.filter((s) => (written.get(s) ?? 0) === 0);
      SharedSession._unwritten(unwritten.filter((s) => s.isDirty.peek()));
      const landing = batch.filter((s) => !unwritten.includes(s)).map((s) => {
        const was = before.get(s)!;
        const rows = written.get(s)!;
        return {source: s, changes: rows >= was.rows ? was.changes : Math.min(rows, was.changes)};
      });
      if (landing.length === 0) {
        if (!this.isDirty.peek())
          notify.info('Reloaded; your changes were discarded');
        return false;
      }
      const landed = landing.reduce((n, l) => n + l.changes, 0);
      const tables = new Set(landing.map((l) => l.source.table)).size;
      notify.info(tables > 1 ? `${plural(landed, 'change', 'changes')} saved in ${tables} tables` :
        SharedSession._landed(landing[0].source, landed, verb));
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

  // the writer each source holds now: after a re-base, the frame the batch landed in
  private static _hold(batch: readonly DomainSource[], saving: boolean): void {
    for (const source of batch)
      source.edit.peek()?.setSaving(saving);
  }

  private static _unwritten(sources: readonly DomainSource[]): void {
    for (const source of sources) {
      source.fail(new DomainBackendError('validation', 'Cannot save: the unsaved changes could not be ' +
        'written to the rows; reload the page to see the table as it is now'));
    }
  }

  private static _landed(source: DomainSource, changes: number, verb: string): string {
    const info = source.schema.info;
    const one = changes === 1;
    const what = nounOf(info, one);
    if (verb === 'saved' && !one)
      return `${changes} changes saved`;
    const subject = one ? capital(what) : `${changes} ${what}`;
    return `${subject} ${verb}`;
  }

  // a pristine parent draft is not dirty, yet its child's `$~new:` needs it; over the pending rows, not the
  // visible ones: a parent the frame's filter hides is still the row the child points at
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

  // two ops for one row: the second loses to the first's version check and the whole transaction rolls back
  private _overlapping(batch: readonly DomainSource[]): boolean {
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
        const problem = `${SharedSession._subject(first.source, first.row)} is edited in two places` +
          ' — discard one of the two edits and save again';
        first.source.refuseRow(first.row.id, problem);
        source.refuseRow(row.id, problem);
        return true;
      }
    }
    return false;
  }

  // the server keeps no hierarchy acyclic: a row moved under its own descendant drops out of every tree
  private async _cycle(batch: readonly DomainSource[]): Promise<boolean> {
    for (const address of new Set(batch.map((s) => s.table))) {
      const sources = batch.filter((s) => s.table === address);
      const column = sources[0].schema.info.parentColumn;
      if (column == null || sources[0].schema.info.hierarchy !== true)
        continue;
      const rows = sources.flatMap((source) => source.pending().map((row) => ({source, row})));
      const moves = rows.filter(({source, row}) => row[Rows.STATE] === 'modified' &&
        source.edit.peek()!.isChanged(row.id, column));
      if (moves.length === 0)
        continue;
      const moved = new Map<unknown, unknown>(rows.filter(({row}) => row[Rows.STATE] === 'new' ||
        moves.some((m) => m.row === row)).map(({row}) => [row.id, row[column]]));
      const table = await emsBackends.domain!.table(address);
      for (const {source, row} of moves) {
        const problem = await hierarchyCycle(table, row[column], (ids) => ids.includes(row.id), moved,
          source.nameOf(row));
        if (problem !== null) {
          source.refuseRow(row.id, problem);
          return true;
        }
      }
    }
    return false;
  }

  /** The unsaved-changes dialog's title and body: one record by its caption and edited fields, a new
   * one by its noun, several by the summary. */
  static atStake(session: DomainSession, action: string): {title: string, body: string} {
    const pending = session instanceof SharedSession ?
      session._sources.peek().flatMap((source) => source.pending().map((row) => ({source, row}))) : [];
    const lost = `will be lost if you ${action}.`;
    if (pending.length !== 1) {
      const summary = session instanceof SharedSession ? session.summary.peek() : '';
      return {title: 'Save changes?', body: `${summary === '' ? 'The changes' : capital(summary)} ${lost}`};
    }
    const [{source, row}] = pending;
    const info = source.schema.info;
    if (Rows.isDraft(row)) {
      const noun = nounOf(info, true);
      return {title: `Save the new ${noun}?`, body: `The new ${noun} ${lost}`};
    }
    const caption = source.nameOf(row);
    const edit = source.edit.peek();
    const fields = source.schema.properties.filter((p) => edit?.isChanged(row.id, p.name!))
      .map((p) => p.friendlyName ?? p.name!);
    const named = fields.length > 3 ? [...fields.slice(0, 2), `${fields.length - 2} more fields`] : fields;
    const list = named.length < 2 ? named.join('') : `${named.slice(0, -1).join(', ')} and ${named.at(-1)}`;
    return {title: caption === null ? `Save changes to the ${nounOf(info, true)}?` : `Save changes to "${caption}"?`,
      body: fields.length === 0 ? `Your changes ${lost}` :
        `Your ${fields.length === 1 ? 'change' : 'changes'} to ${list} ${lost}`};
  }

  private static _subject(source: DomainSource, row: RowView): string {
    const name = source.nameOf(row);
    return `${capital(nounOf(source.schema.info, true))}${name === null ? '' : ` "${name}"`}`;
  }
}

/** THE gate in front of whatever would drop pending changes — a navigation, a filter change, a
 * view close: resolves to whether the caller may proceed. Clean → true without a dialog; a save
 * in flight → a warning and false; else the user saves (false when the save is refused; SAVE is
 * off while the session's validity says why), discards, or cancels. */
export function confirmDiscard(session: DomainSession, options: {action?: string} = {}): Promise<boolean> {
  // saving first: a batch in flight is clean for a while before its write-back lands, and the
  // caller must not be let through that window either
  if (session.isSaving.peek()) {
    notify.warning('Wait for the batch being saved to finish.');
    return Promise.resolve(false);
  }
  if (!session.isDirty.peek())
    return Promise.resolve(true);
  const {title, body} = SharedSession.atStake(session, options.action ?? 'continue');
  return new Promise((resolve) => {
    let decided = false;
    const dialog = Dialog.create(title, {name: 'unsaved-changes'});
    const decide = (answer: boolean | Promise<boolean>) => {
      if (decided)
        return;
      decided = true;
      dialog.close();
      resolve(answer);
    };
    dialog.add(div([body]))
      .addButton('DISCARD', () => {
        session.discard();
        decide(true);
      })
      .addButton('SAVE', () => decide(session.save().then((saved) => saved || !session.isDirty.peek())), {primary: true,
        disabled: session.validity})
      .onCancel(() => decide(false))
      .show({modal: true});
  });
}

/** The spec ambient the ems plugins install: every source a spec instance declares joins the
 * instance's one session — one Save, one transaction. */
export function sessionAmbient<T>(instance: SpecInstance, build: () => T): T {
  let session = specSessions.get(instance);
  if (session === undefined)
    specSessions.set(instance, session = new SharedSession());
  return SharedSession.runWith(session, build);
}
