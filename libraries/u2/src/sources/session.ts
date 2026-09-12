/* The unit of work a Save button drives (GOAL "Unit of work"): what is pending across the sources
   it covers, one `save` for all of it, one `discard`, and the events a paired form follows. Phase
   1 ships `SingleSession` — a source is a session of one; the phase-2 `DomainSession` over several
   tables implements the same interface, so nothing built on `src.session` moves. */
import {signal, computed, ReadonlySignal} from '../core/signals.js';
import {Emitter} from '../core/emitter.js';
import type {ObservableLike} from '../core/widget-like.js';
import {notify} from '../components/display/notify.js';
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
}

export class SingleSession implements DomainSession {
  readonly isDirty: ReadonlySignal<boolean>;
  readonly changeCount: ReadonlySignal<number>;
  readonly validity: ReadonlySignal<string | null>;
  readonly isSaving: ReadonlySignal<boolean>;
  readonly onSaved = new Emitter<void>();
  readonly onDiscarded = new Emitter<void>();

  private readonly _saving = signal(false);

  constructor(private readonly _source: DomainSource) {
    this.isDirty = _source.isDirty;
    this.changeCount = _source.changeCount;
    this.validity = _source.validity;
    this.isSaving = computed(() => this._saving.value || _source.isSaving.value);
  }

  /** The one "saved" balloon: "<Singular> saved" for one change, "<N> changes saved" for more. */
  async save(): Promise<boolean> {
    if (this._saving.peek())
      return false;
    this._saving.value = true;
    try {
      const changes = this.changeCount.peek();
      const saved = await this._source.commit();
      if (saved) {
        const what = this._source.schema.info.singularName || 'Row';
        notify.info(changes > 1 ? `${changes} changes saved` : `${what.charAt(0).toUpperCase()}${what.slice(1)} saved`);
        this.onSaved.fire();
      }
      return saved;
    } finally {
      this._saving.value = false;
    }
  }

  discard(): void {
    this._source.revert();
    this.onDiscarded.fire();
  }
}
