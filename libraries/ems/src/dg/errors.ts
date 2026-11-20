import {emsBackends, nounOf} from '../core/domain-backend.js';
import * as grok from 'datagrok-api/grok';
import {Rows, notify} from '@datagrok-libraries/u2';
import type {NotifyHandle} from '@datagrok-libraries/u2';
import type {DomainSource} from '../core/domain-source.js';
import {DgDomainBackend} from './backend.js';

/** How both backends word the veto on restoring a child whose parent is still deleted. Matched by its shape,
 * whatever the code says: over the platform it arrives as the session's refusal text (code `'refused'`), and
 * the sentence is the only thing both paths share. */
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
    const prop = source?.schema.properties.find((p) => p.name === parts[1]);
    const target = (prop?.friendlyName ?? parts[2]).replace(/_/g, ' ').toLowerCase();
    const staged = (source?.pending() ?? []).filter((r) => r[Rows.STATE] === 'restored');
    const name = source?.schema.info.nameColumn;
    const row = staged.length !== 1 ? null : (name == null ? null : staged[0][name]) ?? staged[0].id;
    return `Cannot restore${row === null ? '' : ` "${String(row)}"`}: its ${target} was deleted too. ` +
      `Restore the ${target} first.`;
  }

  /** What a page over `source` says in place of its empty or not-found state when the table denies the caller
   * View: the server answers such a reader's reads with no rows, by design, so only the access tells "no access"
   * from "nothing there". Null while the caller may view. */
  static noView(source: DomainSource): string | null {
    if (source.access.value.can('view'))
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
   * A hidden column is asked to be shared, a View-only one to be made editable. */
  static noRequired(info: Parameters<typeof nounOf>[0], hidden: readonly string[],
    viewOnly: readonly string[]): string {
    const names = [...hidden, ...viewOnly].join(', ');
    const asks = [...(hidden.length > 0 ? [`share ${hidden.join(', ')} with you`] : []),
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
