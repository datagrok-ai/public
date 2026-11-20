import {emsBackends} from '../core/domain-backend.js';
import * as grok from 'datagrok-api/grok';
import {Rows, notify} from '@datagrok-libraries/u2';
import type {NotifyHandle, RowView} from '@datagrok-libraries/u2';
import {DomainSource} from '../core/domain-source.js';
import {DgDomainBackend} from './backend.js';

/** How both backends word the veto on a row that would be live under a deleted parent: a restore, an insert or
 * a reference moved onto it. Matched by its shape, whatever the code says: over the platform it arrives as the
 * session's refusal text (code `'refused'`), and the sentence is the only thing both paths share. */
const RESTRICT_PARENT = /Column "(.+?)" references a deleted row in "(.+?)"/;

/** What a domain control says when the source reports a failure. Over the platform the js-api editor never
 * throws from `save()` (it balloons, maps column errors onto the cells, resolves a 409 itself), so nothing here
 * duplicates that. */
export class DomainErrors {
  /** The seam's one code, from either backend: the memory one throws `DomainBackendError`, the
   * platform's `DomainError.code` is the server's `body.error` — so `'unsupported'` reads the
   * same whether the table refused the op here or the engine refused it there. */
  static codeOf(e: unknown): string {
    const code = (e as {code?: unknown} | null)?.code;
    return typeof code === 'string' ? code : '';
  }

  static isNotFound(e: unknown): boolean {
    return (e as {status?: unknown} | null)?.status === 404 || DomainErrors.codeOf(e) === 'not-found';
  }

  /** `source` is what the sentence is worded from where the refusal is about one of its rows. */
  static message(e: unknown, source?: DomainSource): string {
    const raw = e instanceof Error ? e.message : String((e as {message?: unknown} | null)?.message ?? e);
    const parts = RESTRICT_PARENT.exec(raw);
    if (parts === null)
      return raw;
    // the backend names a column and a table; the user asked to bring a ROW back, and the way
    // forward is the parent, not the column
    const column = parts[1];
    const staged = (source?.pending() ?? []).filter((r) => r[Rows.STATE] === 'restored');
    const writer = staged.length > 0 ? undefined : DomainErrors._writer(source, column);
    const owner = writer?.source ?? source;
    const prop = owner?.schema.properties.find((p) => p.name === column);
    const target = (prop?.friendlyName ?? parts[2]).replace(/_/g, ' ').toLowerCase();
    if (writer === undefined) {
      const row = staged.length !== 1 ? null : DomainErrors._nameOf(source!, staged[0]);
      return `Cannot restore${row === null ? '' : ` "${row}"`}: its ${target} was deleted too. ` +
        `Restore the ${target} first.`;
    }
    // a new row, or a reference moved: the parent is the way forward, so it is named as the user knows it
    const child = DomainErrors._nameOf(writer.source, writer.row);
    const subject = `Cannot save${child === null ? '' : ` "${child}"`}: its ${target}`;
    const parentId = String(writer.row[column]);
    for (const live of DomainSource.live) {
      const parent = live.rows.byKey(parentId);
      if (parent === undefined)
        continue;
      const named = DomainErrors._nameOf(live, parent);
      return parent[Rows.STATE] === 'deleted' ?
        `${subject} "${named}" is being deleted. Undo the deletion, or remove the new reference to it.` :
        `${subject} "${named}" was deleted. Restore the ${target} first.`;
    }
    return `${subject} was deleted. Restore the ${target} first.`;
  }

  // the pending row of the session that points at the deleted parent through `column`
  private static _writer(source: DomainSource | undefined, column: string):
    {source: DomainSource, row: RowView} | undefined {
    if (source === undefined)
      return undefined;
    for (const s of [source, ...[...DomainSource.live].filter((x) => x !== source && x.session === source.session)]) {
      if (!s.schema.properties.some((p) => p.name === column))
        continue;
      const row = s.pending().find((r) =>
        (r[Rows.STATE] === 'new' || r[Rows.STATE] === 'modified') && r[column] != null);
      if (row !== undefined)
        return {source: s, row};
    }
    return undefined;
  }

  private static _nameOf(source: DomainSource, row: RowView): string {
    const name = source.schema.info.nameColumn;
    const value = name == null ? null : row[name];
    return String(value ?? source.caption?.(row) ?? row.id);
  }

  /** One balloon; after a 403 the affordances were built from an access snapshot the server just
   * contradicted, so every cache of it goes. */
  static report(e: unknown, source?: DomainSource): NotifyHandle | null {
    if (DomainErrors.codeOf(e) === 'forbidden') {
      grok.dapi.domains.invalidateUiCaches();
      if (emsBackends.domain instanceof DgDomainBackend)
        emsBackends.domain.invalidate();
    }
    return notify.error(DomainErrors.message(e, source));
  }
}
