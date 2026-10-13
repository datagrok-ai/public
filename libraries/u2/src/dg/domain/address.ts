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

  /** Every literal the platform's row routes register at the key position (`GET …/<table>/version`
   * answers the table version, never a row), the dot segments a browser folds away, and the app's
   * draft sentinel — the one list the Dart `DomainRowKey` and the Node dapi keep too. */
  static readonly RESERVED = ['access', 'aggregate', 'audit', 'batch', 'delete', 'facets', 'filters', 'query',
    'update', 'version', 'watch', 'new', '.', '..'];

  /** A key as a path segment spells it, before percent-encoding: a reserved one, or one led by `~`,
   * gets a `~` in front — `watch` ↔ `~watch`, `~x` ↔ `~~x`, a uuid never changes — so a segment
   * names one row and never a route. The key is the row's canonical id (an external table's is its
   * business key, components percent-encoded and comma-joined) or its business-key spelling. */
  static escape(key: string): string {
    return DomainAddress.RESERVED.includes(key) || key.startsWith('~') ? `~${key}` : key;
  }

  /** {@link escape}'s inverse: one leading `~` off a percent-decoded segment. */
  static unescape(segment: string): string {
    return segment.startsWith('~') ? segment.slice(1) : segment;
  }

  /** The path segment of a key: escaped, then percent-encoded — the write twin of
   * {@link segmentOf} + {@link unescape}. */
  static segment(key: string): string {
    return encodeURIComponent(DomainAddress.escape(key));
  }

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

  /** The trailing segment of an address under `bases`, percent-decoded but still escaped — the
   * draft sentinel is read off it as is, {@link unescape} gives the row's key; null when the
   * address is a base itself, or another view's. */
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
