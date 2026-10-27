/* `domains.history` — a row's audit trail, newest first: who did what when, and for an update
   one line per changed column, `<caption>  <before> → <after>`. Reads `table.audit(id)` through
   the seam a page at a time, so the memory backend's history shows the same way; "Show older"
   under the list loads the next page, and the pane is refreshed when the row's session saves. A
   draft has no history yet, and says so. */
import {emsBackends, nounOf} from '../core/domain-backend.js';
import * as grok from 'datagrok-api/grok';
import {Component, Control, signal, AsyncSource, div, divH, divV, span, timestamp, text, badge, button, AsyncView,
  Accordion, Rows, Scope, isEmpty} from '@datagrok-libraries/u2';
import type {ReadonlySignal, BadgeVariant, IProperty, RowView} from '@datagrok-libraries/u2';
import type {AuditEntryLike} from '../core/domain-backend.js';
import type {DomainSource} from '../core/domain-source.js';
import {SYSTEM_COLUMNS} from './backend.js';
import {DomainTable} from './facade.js';
import {DomainForm} from './form.js';
import {DomainDataTable} from './grid.js';
import {DomainErrors} from './errors.js';

/** The current row of the source (default), a row signal, or one row. */
export type DomainHistoryTarget = RowView | ReadonlySignal<RowView | null>;

export interface DomainHistoryOptions {
  /** The columns the user edits (the form's `include`); changes to derived columns are left out. */
  include?: string[];
  exclude?: string[];
}

const OPS: Record<string, [label: string, variant: BadgeVariant]> = {
  insert: ['created', 'success'], update: ['updated', 'default'], delete: ['deleted', 'error'],
  promote: ['shared', 'accent'], undelete: ['restored', 'success'],
};
const PAGE_SIZE = 50;
const SHOWN_CHANGES = 5;
const NO_VALUE = '(empty)';

export class DomainHistory extends Control {
  readonly row: ReadonlySignal<RowView | null>;
  readonly view: AsyncView<AuditEntryLike>;

  private readonly _hint: HTMLElement;
  /** Whether the table keeps no history at all — its handle declares it, so the pane says so
   * instead of showing an empty list nothing will ever fill. */
  private readonly _noHistory = signal(false);
  private readonly _summary = signal('');
  private _failure = '';
  /** id → the name behind it, resolved once per pane: the actors and the ref cells alike. */
  private readonly _names = new Map<string, Promise<string>>();

  constructor(readonly source: DomainSource, row?: DomainHistoryTarget,
    private readonly _options: DomainHistoryOptions = {}) {
    super();
    this.root.classList.add('u2-domain-history', 'u2-async-quiet');
    this.root.dataset.u2 = 'domain-history';
    this.row = DomainForm.resolve(row ?? source, source).row;
    this._hint = div([], 'u2-domain-history-hint');
    this._hint.dataset.u2Part = 'hint';
    const entries = new AsyncSource<AuditEntryLike>((id) => this._load(id), {debounceMs: 0});
    this.own(() => entries.dispose());
    this.view = this.runInScope(() => new AsyncView<AuditEntryLike>(entries, (items) => this._lines(items),
      {empty: 'No history yet'}));
    this.effect(() => {
      const failed = entries.state.value.kind === 'error';
      this.view.root.title = failed ? this._failure : '';
      if (failed)
        this._summary.value = 'unavailable';
    });
    const accordion = this.runInScope(() => new Accordion());
    accordion.addPane('History', div([this._hint, this.view.root]),
      {hint: this._summary, expanded: false, persistKey: 'u2-domain-history'});
    this.root.append(accordion.root);
    this.effect(() => {
      const row = this.row.value;
      // a draft id, and any other key the frame stood in for a row it has no id for, names
      // nothing the server could be asked about
      const hint = this._noHistory.value ? 'This table keeps no history.' :
        row === null ? `Select a ${nounOf(source.schema.info, true)} to see its history.` :
          Rows.isService(row.id) ? 'Not saved yet' : null;
      this._hint.textContent = hint ?? '';
      this._hint.hidden = hint === null;
      this.view.root.hidden = hint !== null;
      this._summary.value = hint === 'Not saved yet' ? 'not saved yet' : '';
      if (hint === null)
        this.view.refresh(row!.id);
    });
    const saved = source.session.onSaved.subscribe(() => {
      if (!this.view.root.hidden)
        this.view.refresh();
    });
    this.own(() => saved.unsubscribe());
  }

  /** The changed columns of an update, captions from the schema; the system columns a write
   * always touches (version, updated) are not changes the user made. */
  changesOf(entry: AuditEntryLike): string {
    return this.changedColumns(entry)
      .map(({prop, before, after}) =>
        `${prop.friendlyName ?? prop.name}: ${DomainHistory._display(before)} → ${DomainHistory._display(after)}`)
      .join(', ');
  }

  /** What an update changed, as the schema describes each column — the ids a reference cell holds
   * are resolved to names when the line is drawn. */
  changedColumns(entry: AuditEntryLike): {prop: IProperty, before: unknown, after: unknown}[] {
    if (entry.op !== 'update')
      return [];
    const before = entry.before ?? {};
    const after = entry.after ?? {};
    return this.source.schema.properties
      .filter((p) => p.name! in after && !SYSTEM_COLUMNS.some(([name]) => name === p.name) &&
        !Component.sameValue(before[p.name!], after[p.name!]))
      .map((prop) => ({prop, before: before[prop.name!], after: after[prop.name!]}));
  }

  /** A value at midnight UTC is a date, never shifted a day by the local zone. */
  private static _display(raw: unknown, prop?: IProperty): string {
    if (isEmpty(raw))
      return NO_VALUE;
    if (raw !== null && typeof raw === 'object' && raw.constructor === Object)
      return JSON.stringify(raw);
    switch (prop?.propertyType ?? prop?.type) {
    case 'datetime': return timestamp(raw as string, undefined, {utcDates: true}).textContent || text(raw);
    case 'double': case 'float': case 'num': return DomainDataTable.number(raw as number, prop!.format);
    default: return text(raw);
    }
  }

  private async _load(id: string): Promise<AuditEntryLike[]> {
    const table = await emsBackends.domain!.table(this.source.table);
    if (table.audit === undefined) {
      this._noHistory.value = true;
      return [];
    }
    try {
      const entries = await table.audit(id, {limit: PAGE_SIZE});
      return entries.slice().reverse();
    } catch (e) {
      this._failure = DomainErrors.message(e, this.source);
      throw new Error('History is unavailable');
    }
  }

  private _lines(items: AuditEntryLike[]): HTMLElement {
    const list = divV(items.map((entry) => this._line(entry)), 'u2-domain-history-lines');
    list.setAttribute('role', 'list');
    this._summarize(items[0]);
    const scope = Scope.ambient!;
    let oldest = items.at(-1)!;
    // a button, not a scroll trigger: on the entity page the page scrolls, not the list
    const more: HTMLButtonElement = button('Show older', async () => {
      more.disabled = true;
      try {
        const table = await emsBackends.domain!.table(this.source.table);
        const page = await table.audit!(oldest.id, {limit: PAGE_SIZE, before: oldest.seq});
        if (scope.isDisposed)
          return;
        const older = page.slice().reverse();
        Scope.runWith(scope, () => list.append(...older.map((entry) => this._line(entry))));
        oldest = older.at(-1) ?? oldest;
        more.hidden = page.length < PAGE_SIZE;
      } catch (e) {
        if (!scope.isDisposed)
          DomainErrors.report(e, this.source);
      } finally {
        more.disabled = false;
      }
    });
    more.classList.add('u2-domain-history-more');
    more.dataset.u2Part = 'more';
    more.hidden = items.length < PAGE_SIZE;
    return divV([list, more], 'u2-domain-history-trail');
  }

  private _summarize(entry: AuditEntryLike): void {
    const [label] = OPS[entry.op] ?? [entry.op];
    const when = `${label} ${timestamp(entry.ts).textContent ?? ''}`;
    this._summary.value = when;
    if (entry.actor_id !== null)
      void this._actor(entry.actor_id).then((name) => this._summary.value = `${when} by ${name}`);
  }

  private _actor(id: string): Promise<string> {
    return this._resolve(id, () => grok.dapi.users.find(id).then((user) => user?.friendlyName ?? null));
  }

  private _line(entry: AuditEntryLike): HTMLElement {
    const [label, variant] = OPS[entry.op] ?? [entry.op, 'default'];
    const actorId = entry.actor_id;
    const actor = span(actorId ?? 'system', 'u2-domain-history-actor');
    if (actorId !== null)
      void this._actor(actorId).then((name) => actor.textContent = name);
    const head = divH([timestamp(entry.ts, 'u2-domain-history-time'), actor, badge(label, {variant})],
      'u2-domain-history-head');
    const changes = this._changes(entry);
    const line = divV(changes === null ? [head] : [head, changes], 'u2-domain-history-line');
    line.setAttribute('role', 'listitem');
    return line;
  }

  private _changes(entry: AuditEntryLike): HTMLElement | null {
    const {include, exclude} = this._options;
    const own = this.changedColumns(entry).filter(({prop}) => (include === undefined || include.includes(prop.name!)) &&
      !exclude?.includes(prop.name!));
    if (own.length === 0)
      return null;
    const lines = own.map(({prop, before, after}) => divH([span(DomainTable.captionOf(prop), 'u2-domain-history-field'),
      this._value(prop, before), span('→', 'u2-domain-history-arrow'), this._value(prop, after)],
    'u2-domain-history-change'));
    const box = divV(lines, 'u2-domain-history-changes');
    const folded = lines.slice(SHOWN_CHANGES);
    if (folded.length > 0) {
      for (const line of folded)
        line.hidden = true;
      const more: HTMLButtonElement = button(`${folded.length} more`, () => {
        for (const line of folded)
          line.hidden = false;
        more.remove();
      });
      more.classList.add('u2-domain-history-expand');
      box.append(more);
    }
    return box;
  }

  /** A cell as the line shows it: a reference is the name it points at, an id standing in until
   * the platform answers. This one keeps its own lookup where the rest of the domain path reads
   * `~caption_<column>` off the frame: an audit line's `before`/`after` values are ids out of a
   * jsonb snapshot, not frame cells, and no caption column will ever carry them. */
  private _value(prop: IProperty, raw: unknown): HTMLElement {
    const el = span(DomainHistory._display(raw, prop), 'u2-domain-history-value');
    if (isEmpty(raw)) {
      el.classList.add('u2-domain-history-empty');
      return el;
    }
    if (DomainTable.isReference(prop)) {
      const id = text(raw);
      void this._resolve(id, () => DomainForm.captionOf(prop, id)).then((name) => el.textContent = name);
    }
    return el;
  }

  private _resolve(id: string, lookUp: () => Promise<string | null>): Promise<string> {
    let name = this._names.get(id);
    if (name === undefined) {
      name = lookUp().then((n) => n ?? id, () => id);
      this._names.set(id, name);
    }
    return name;
  }
}
