import {emsBackends, nounOf} from '../core/domain-backend.js';
import * as grok from 'datagrok-api/grok';
import {Rows, notify, plural} from '@datagrok-libraries/u2';
import type {NotifyHandle, RowView} from '@datagrok-libraries/u2';
import {DomainSource} from '../core/domain-source.js';
import {DgDomainBackend} from './backend.js';

/** How both backends word the veto on a row that would be live under a deleted parent: a restore, an insert or
 * a reference moved onto it. Matched by its shape, whatever the code says: over the platform it arrives as the
 * session's refusal text (code `'refused'`), and the sentence is the only thing both paths share. */
// the parent's id follows where the backend gives it (an insert's or an update's veto); an older server's restore
// veto ends at the table
const RESTRICT_PARENT = /Column "(.+?)" references a deleted row in "(.+?)"(?: \("(.+?)"\))?/;

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
    const writer = staged.length > 0 ? undefined : DomainErrors._writer(source, column, parts[3]);
    const owner = writer?.source ?? source;
    const prop = owner?.schema.properties.find((p) => p.name === column);
    const target = (prop?.friendlyName ?? parts[2]).replace(/_/g, ' ').toLowerCase();
    if (writer === null)
      return `Cannot save: a reference points at a deleted ${target}. Restore the ${target} first.`;
    if (writer === undefined) {
      const row = staged.length !== 1 ? null : DomainErrors._nameOf(source!, staged[0]);
      return `Cannot restore${row === null ? '' : ` "${row}"`}: its ${target} was deleted too. ` +
        `Restore the ${target} first.`;
    }
    // a new row, or a reference moved: the parent is the way forward, so it is named as the user knows it
    const child = DomainErrors._nameOf(writer.source, writer.row);
    const subject = `Cannot save${child === null ? '' : ` "${child}"`}: its ${target}`;
    const parent = DomainErrors._held(String(writer.row[column]));
    if (parent === undefined)
      return `${subject} was deleted. Restore the ${target} first.`;
    const named = DomainErrors._nameOf(parent.source, parent.row);
    const quoted = named === null ? '' : ` "${named}"`;
    return parent.row[Rows.STATE] === 'deleted' ?
      `${subject}${quoted} is being deleted. Undo the deletion, or remove the new reference to it.` :
      `${subject}${quoted} was deleted. Restore the ${target} first.`;
  }

  // the row naming `parentId`; without it (an older server) the one whose parent is known deleted, null if ambiguous
  private static _writer(source: DomainSource | undefined, column: string, parentId: string | undefined):
    {source: DomainSource, row: RowView} | null | undefined {
    if (source === undefined)
      return undefined;
    const candidates: {source: DomainSource, row: RowView}[] = [];
    for (const s of [source, ...[...DomainSource.live].filter((x) => x !== source && x.session === source.session)]) {
      if (!s.schema.properties.some((p) => p.name === column))
        continue;
      for (const row of s.pending()) {
        if ((row[Rows.STATE] === 'new' || row[Rows.STATE] === 'modified') && row[column] != null &&
          (parentId === undefined || String(row[column]) === parentId))
          candidates.push({source: s, row});
      }
    }
    if (candidates.length <= 1 || parentId !== undefined)
      return candidates[0];
    const dead = candidates.filter((c) => {
      const parent = DomainErrors._held(String(c.row[column]))?.row;
      return parent !== undefined && (parent[Rows.STATE] === 'deleted' || Rows.isDeleted(parent));
    });
    return dead.length === 1 ? dead[0] : null;
  }

  // the row a live source of this page holds under `id`
  private static _held(id: string): {source: DomainSource, row: RowView} | undefined {
    for (const source of DomainSource.live) {
      const row = source.rows.byKey(id);
      if (row !== undefined)
        return {source, row};
    }
    return undefined;
  }

  // the renderer's caption: "New <noun>" for a draft with no name; a saved row's bare id names nothing, so null
  private static _nameOf(source: DomainSource, row: RowView): string | null {
    const name = source.schema.info.nameColumn;
    const caption = source.caption?.(row) ?? String((name == null ? null : row[name]) ?? '');
    return caption === '' || caption === row.id ? null : caption;
  }

  /** What a page over `source` says in place of its empty or not-found state when the table denies the caller
   * View: the server answers such a reader's reads with no rows, by design, so only the access tells "no access"
   * from "nothing there". Null while the caller may view, or holds rows: View is table-wide, and a row-mode
   * table still answers the rows granted one by one. Null too once a search, or a user's filter (`filtered`: a
   * child tab's link or a row page's id is the page's own), narrowed the read: granted rows may just not match. */
  static noView(source: DomainSource, filtered = false): string | null {
    if (source.access.value.can('view') || source.rows.items.value.length > 0 || filtered ||
      source.search.value !== '')
      return null;
    return DomainErrors.noAccess(source.schema.info);
  }

  /** The no-View wording for a table whose access is known without a source (a picker's target). */
  static noAccess(info: Parameters<typeof nounOf>[0]): string {
    return `You do not have access to ${nounOf(info)}. Ask the table's owner to share it with you.`;
  }

  /** The wording for a child tab whose FK to the parent is a column the caller cannot read: the server refuses
   * every read and write through that column, so the tab cannot show or add the children. */
  static noLink(child: Parameters<typeof nounOf>[0], parent: Parameters<typeof nounOf>[0]): string {
    return `You do not have access to ${nounOf(child)}: you cannot see how ${nounOf(child)} link to this ` +
      `${nounOf(parent, true)}. Ask the table's owner to share the link column with you.`;
  }

  /** The wording for creating on a table with required columns the caller cannot set: every Save would be refused.
   * A hidden column is asked to be shared, a View-only one to be made editable; a column the caller may not
   * read (`unseen` of them) is counted, never named. */
  static noRequired(info: Parameters<typeof nounOf>[0], hidden: readonly string[],
    viewOnly: readonly string[], unseen = 0): string {
    const names = [...hidden, ...viewOnly,
      ...(unseen > 0 ? [`${plural(unseen, 'column', 'columns')} you cannot see`] : [])].join(', ');
    const asks = [...(hidden.length > 0 ? [`share ${hidden.join(', ')} with you`] : []),
      ...(unseen > 0 ? [`share ${unseen === 1 ? 'it' : 'them'} with you`] : []),
      ...(viewOnly.length > 0 ? [`let you edit ${viewOnly.join(', ')}`] : [])];
    return `You cannot create ${nounOf(info)}: they need a value for ${names}, which you cannot set. ` +
      `Ask the table's owner to ${asks.join(' and ')}.`;
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
