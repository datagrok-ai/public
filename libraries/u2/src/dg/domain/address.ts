/* Everything about how a domain row is spelled in an address and read back out of one (GOAL R2):
   `/domains/<schema>/<table>[/<keyOrId>]` addresses a row by a path segment, an app mounted at
   `/apps/…` by `?entity=<id>`. One module, so `app.ts` and `routes.ts` cannot disagree about a link,
   and an external binding that addresses rows differently has ONE place to say so. */
import {Filters} from '../../core/filter/index.js';
import type {FilterGroup, FilterSchema} from '../../core/filter/index.js';
import type {RowView} from '../../sources/rows-like.js';

export class DomainAddress {
  /** A row id as the platform spells one — what tells an id from a business key in an address. */
  static readonly ID = /^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$/i;

  /** `/domains/<schema>/<table>[/<segment>]`, the platform's own row route; the query is
   * {@link DomainApp.open}'s to read. */
  static readonly ROUTE = /^\/domains\/(\w+)\/(\w+)(\/[^/?#]*)?$/i;

  /** Whether `base` addresses rows by a path segment (`${base}/${key}`) rather than by
   * `?entity=`: the platform's own `/domains/…` routes are, an app mounted at `/apps/…` is not
   * (ruling R2). */
  static entityPath(base: string): boolean {
    return base.startsWith('/domains/');
  }

  /** How a row is spelled in a `/domains` path: the business key when it is unambiguous, the id
   * otherwise — no key declared, a null component, or a composite key whose values carry a '-'
   * (the TS half of `domain_row_meta.dart` `deepLink`). */
  static keyOf(row: RowView, businessKey: readonly string[]): string {
    const parts: string[] = [];
    for (const column of businessKey) {
      const value = row[column];
      if (value === null || value === undefined)
        return row.id;
      parts.push(String(value));
    }
    if (parts.length === 0 || (parts.length > 1 && parts.some((part) => part.includes('-'))))
      return row.id;
    return parts.join('-');
  }

  /** Whether `path` is `base` itself or continues it at a segment boundary — a view at
   * `/domains/grit/issue` does not answer for `/domains/grit/issue_label`. Both are compared as
   * given: a caller that does not know the case of either lower-cases them first. */
  static under(path: string, base: string): boolean {
    if (!path.startsWith(base))
      return false;
    const rest = path.slice(base.length);
    return rest === '' || rest.startsWith('/') || rest.startsWith('?');
  }

  /** What `path` carries under any of `bases`: `''` for a base itself, `'/<row>…'` below one,
   * null when the address is another view's. Both of an app's routes count — a `/domains/…` link
   * still reaches an app that has rebased onto `/apps/…`. */
  static restOf(path: string, bases: readonly string[]): string | null {
    const here = path.toLowerCase();
    for (const base of bases) {
      if (!here.startsWith(base.toLowerCase()))
        continue;
      const rest = path.slice(base.length);
      if (rest === '' || rest.startsWith('/'))
        return rest;
    }
    return null;
  }

  /** The trailing segment of an address under `bases` — the row it names; null when the address
   * is a base itself, or another view's. */
  static segmentOf(path: string, bases: readonly string[]): string | null {
    const rest = DomainAddress.restOf(path, bases);
    if (rest === null || rest === '')
      return null;
    const segment = rest.slice(1).split('/')[0];
    return segment === '' ? null : decodeURIComponent(segment);
  }

  /** The query finding the row an address segment names, null where the segment IS the id: a uuid
   * is one, and so is anything under a table with no business key; one key column takes the whole
   * segment, a composite key splits on '-' and must match its arity (`domain_entity_view.dart`
   * `_resolve`). A miss shows the form's not-found, as the Dart view does. */
  static keyQuery(segment: string, businessKey: readonly string[], schema: FilterSchema): FilterGroup | null {
    if (businessKey.length === 0 || DomainAddress.ID.test(segment))
      return null;
    const parts = businessKey.length === 1 ? [segment] : segment.split('-');
    if (parts.length !== businessKey.length)
      return null;
    return Filters.group('and', businessKey.map((column, at) => {
      const property = Filters.property(schema, column);
      const kind = property === null ? Filters.KIND.STRING : Filters.kindOf(property);
      const numeric = kind === Filters.KIND.INT || kind === Filters.KIND.FLOAT;
      return Filters.cond(column, '=', numeric ? Number(parts[at]) : parts[at]);
    }));
  }
}
