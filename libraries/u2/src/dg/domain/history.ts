/* `domains.history` — a row's audit trail, newest first: who did what when, and for an update
   the columns it changed as `<caption>: <before> → <after>`. Reads `table.audit(id)` through the
   seam, so the memory backend's history shows the same way; refreshed when the row's session
   saves. A draft has no history yet, and says so. */
import * as grok from 'datagrok-api/grok';
import {Control} from '../../core/component.js';
import {signal} from '../../core/signals.js';
import type {ReadonlySignal} from '../../core/signals.js';
import {AsyncSource} from '../../core/async-source.js';
import {div, divH, span, timestamp} from '../../core/elements.js';
import {text} from '../../core/text.js';
import {badge} from '../../components/display/badge.js';
import type {BadgeVariant} from '../../components/display/badge.js';
import {AsyncView} from '../../components/display/async-view.js';
import {Section} from '../../components/containers/section.js';
import {VirtualList} from '../../components/collections/list.js';
import {backends} from '../../sources/backends.js';
import type {AuditEntryLike} from '../../sources/domain-backend.js';
import type {DomainSource} from '../../sources/domain-source.js';
import type {IProperty} from '../../core/property-like.js';
import {Rows} from '../../sources/rows-like.js';
import type {RowView} from '../../sources/rows-like.js';
import {SYSTEM_COLUMNS} from './backend.js';
import {DomainTable} from './index.js';
import {DomainForm} from './form.js';

/** The current row of the source (default), a row signal, or one row. */
export type DomainHistoryTarget = RowView | ReadonlySignal<RowView | null>;

const OPS: Record<string, [label: string, variant: BadgeVariant]> = {
  insert: ['created', 'success'], update: ['updated', 'default'], delete: ['deleted', 'error'],
  promote: ['shared', 'accent'],
};
const ITEM_HEIGHT = 40;

export class DomainHistory extends Control {
  readonly row: ReadonlySignal<RowView | null>;
  readonly view: AsyncView<AuditEntryLike>;

  private readonly _hint: HTMLElement;
  /** Whether the table keeps no history at all — its handle declares it, so the pane says so
   * instead of showing an empty list nothing will ever fill. */
  private readonly _noHistory = signal(false);
  /** id → the name behind it, resolved once per pane: the actors and the ref cells alike. */
  private readonly _names = new Map<string, Promise<string>>();

  constructor(readonly source: DomainSource, row?: DomainHistoryTarget) {
    super();
    this.root.classList.add('u2-domain-history');
    this.root.dataset.u2 = 'domain-history';
    this.row = DomainForm.resolve(row ?? source, source).row;
    this._hint = div([], 'u2-domain-history-hint');
    this._hint.dataset.u2Part = 'hint';
    const entries = new AsyncSource<AuditEntryLike>((id) => this._load(id), {debounceMs: 0});
    this.own(() => entries.dispose());
    this.view = this.runInScope(() => new AsyncView<AuditEntryLike>(entries, (items) => this._lines(items),
      {empty: 'No history yet'}));
    const section = this.runInScope(() => new Section({title: 'History'}));
    section.add(this._hint, this.view.root);
    this.root.append(section.root);
    this.effect(() => {
      const row = this.row.value;
      const singular = source.schema.info.singularName.toLowerCase() || 'row';
      // a draft id, and any other key the frame stood in for a row it has no id for, names
      // nothing the server could be asked about
      const hint = this._noHistory.value ? 'This table keeps no history.' :
        row === null ? `Select a ${singular} to see its history.` :
          Rows.isService(row.id) ? 'Not saved yet' : null;
      this._hint.textContent = hint ?? '';
      this._hint.hidden = hint === null;
      this.view.root.hidden = hint !== null;
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
      .map(({prop, before, after}) => `${prop.friendlyName ?? prop.name}: ${text(before)} → ${text(after)}`)
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
        before[p.name!] !== after[p.name!])
      .map((prop) => ({prop, before: before[prop.name!], after: after[prop.name!]}));
  }

  private async _load(id: string): Promise<AuditEntryLike[]> {
    const table = await backends.domain!.table(this.source.table);
    if (table.audit === undefined) {
      this._noHistory.value = true;
      return [];
    }
    const entries = await table.audit(id);
    return entries.slice().reverse();
  }

  private _lines(items: AuditEntryLike[]): HTMLElement {
    const list = new VirtualList<AuditEntryLike>({
      itemHeight: ITEM_HEIGHT, rowRole: 'listitem',
      render: (entry) => this._line(entry),
    });
    list.root.classList.add('u2-domain-history-lines');
    list.setItems(items);
    return list.root;
  }

  private _line(entry: AuditEntryLike): HTMLElement {
    const [label, variant] = OPS[entry.op] ?? [entry.op, 'default'];
    const actorId = entry.actor_id;
    const actor = span(actorId ?? 'system', 'u2-domain-history-actor');
    if (actorId !== null) {
      void this._resolve(actorId, () => grok.dapi.users.find(actorId).then((user) => user?.friendlyName ?? null))
        .then((name) => actor.textContent = name);
    }
    const changes = div([], 'u2-domain-history-changes');
    const retitle = () => changes.title = changes.textContent ?? '';
    let first = true;
    for (const {prop, before, after} of this.changedColumns(entry)) {
      changes.append(span(`${first ? '' : ', '}${prop.friendlyName ?? prop.name}: `),
        this._value(prop, before, retitle), span(' → '), this._value(prop, after, retitle));
      first = false;
    }
    retitle();
    const line = divH([timestamp(entry.ts, 'u2-domain-history-time'), actor, badge(label, {variant}), changes],
      'u2-domain-history-line');
    return line;
  }

  /** A cell as the line shows it: a reference is the name it points at, an id standing in until
   * the platform answers. This one keeps its own lookup where the rest of the domain path reads
   * `~caption_<column>` off the frame: an audit line's `before`/`after` values are ids out of a
   * jsonb snapshot, not frame cells, and no caption column will ever carry them. */
  private _value(prop: IProperty, raw: unknown, resolved: () => void): HTMLElement {
    const el = span(text(raw), 'u2-domain-history-value');
    const id = text(raw);
    if (id !== '' && DomainTable.isReference(prop)) {
      void this._resolve(id, () => DomainForm.captionOf(prop, id)).then((name) => {
        el.textContent = name;
        resolved();
      });
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
