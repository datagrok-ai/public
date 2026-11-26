import {emsBackends, errorCode, nounOf} from '../core/domain-backend.js';
import {Rows, notify, plural} from '@datagrok-libraries/u2';
import type {NotifyHandle} from '@datagrok-libraries/u2';
import {DomainSource} from '../core/domain-source.js';

/** The deleted-parent veto, matched by its sentence: the two backends share nothing else. */
const RESTRICT_PARENT = /Column "(.+?)" references a deleted row in "(.+?)" \("(.+?)"\)/;

/** What a domain control says when the source reports a failure. Over the platform the js-api editor never
 * throws from `save()` (it balloons, maps column errors onto the cells, resolves a 409 itself), so nothing here
 * duplicates that. */
export class DomainErrors {
  /** The seam's one code, from either backend — {@link errorCode}. */
  static readonly codeOf = errorCode;

  static isNotFound(e: unknown): boolean {
    return (e as {status?: unknown} | null)?.status === 404 || errorCode(e) === 'not-found';
  }

  /** `source` is what the sentence is worded from where the refusal is about one of its rows. */
  static message(e: unknown, source?: DomainSource): string {
    const raw = e instanceof Error ? e.message : String((e as {message?: unknown} | null)?.message ?? e);
    const parts = RESTRICT_PARENT.exec(raw);
    if (parts === null)
      return raw;
    // the backend names a column and a table; the way forward is the parent row, not the column
    const [, column, table, parentId] = parts;
    // the refusal is about the batch: the row may be pending in any source of the session
    const rows = source === undefined ? [] : [...DomainSource.live].filter((s) => s.session === source.session &&
      s.schema.properties.some((p) => p.name === column)).flatMap((s) => s.pending()
      .filter((r) => String(r[column]) === parentId).map((r) => ({source: s, row: r})));
    const found = rows.length === 1 ? rows[0] : undefined;
    const name = found === undefined ? null : found.source.nameOf(found.row);
    const prop = (found?.source ?? source)?.schema.properties.find((p) => p.name === column);
    const target = (prop?.friendlyName ?? table).replace(/_/g, ' ').toLowerCase();
    const quoted = name === null ? '' : ` "${name}"`;
    if (found?.row[Rows.STATE] === 'restored')
      return `Cannot restore${quoted}: its ${target} was deleted too. Restore the ${target} first.`;
    const parent = DomainSource.draftOf(parentId);
    const named = parent === undefined ? null : parent.source.nameOf(parent.row);
    const fix = parent?.row[Rows.STATE] !== 'deleted' ? `was deleted. Restore the ${target} first.` :
      'is being deleted. Undo the deletion, or remove the new reference to it.';
    return `Cannot save${quoted}: its ${target}${named === null ? '' : ` "${named}"`} ${fix}`;
  }

  /** What a page says in place of its empty state when the table denies the caller View; null otherwise. */
  static noView(source: DomainSource, filtered = false): string | null {
    if (source.access.value.can('view') || source.rows.items.value.length > 0 || filtered ||
      source.search.value !== '')
      return null;
    return DomainErrors.noAccess(source.schema.info);
  }

  static noAccess(info: Parameters<typeof nounOf>[0]): string {
    return `You do not have access to ${nounOf(info)}. Ask the table's owner to share it with you.`;
  }

  static noLink(child: Parameters<typeof nounOf>[0], parent: Parameters<typeof nounOf>[0]): string {
    return `You do not have access to ${nounOf(child)}: you cannot see how ${nounOf(child)} link to this ` +
      `${nounOf(parent, true)}. Ask the table's owner to share the link column with you.`;
  }

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
    if (errorCode(e) === 'forbidden')
      emsBackends.domain?.dropAccess();
    return notify.error(DomainErrors.message(e, source));
  }
}
