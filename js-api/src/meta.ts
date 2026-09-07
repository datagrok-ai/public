/**
 * Discovery of the properties a platform entity type or a domain table exposes —
 * the same curated catalogs the platform's own grids, filter panels and server-side
 * facets are built from, so what is discoverable here is exactly what is filterable
 * and renderable there.
 */
import {DomainRegistryClient} from './dapi';
import {DomainValidationError} from './domains';
import {Property} from './entities/property';
import {RELATION_KIND, TYPE} from './const';
import {toJs} from './wrappers';
import {IDartApi} from './api/grok_api.g';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/** One discovered property. `refType` and `relationKind` are answered by platform
 * types only — a domain-table column carries neither. */
export interface EntityPropertyInfo {
  name: string;
  type: TYPE;
  semType: string | null;
  friendlyName: string;
  description: string | null;
  /** The type this property points at (`'User'` for `UserSession.user`), else null. */
  refType: string | null;
  relationKind: RELATION_KIND | null;
}

/** Address of the read-only domain table that serves a platform type's rows. */
export interface CoreLocation {
  schema: string;
  table: string;
}

/** Reflection over what the platform knows about its own types (`grok.meta`). */
export class Meta {
  /** Properties of [type] — a platform entity type name (`'User'`, `'UserSession'`) or a
   * domain table addressed as `'<schema>.<table>'` (`'Core.users'`, `'plates.plate'`).
   * Null when the platform curates no catalog for the type, and null for a domain table
   * that does not exist — a caller cannot tell an uncurated type from an unknown one,
   * and neither can the platform.
   *
   * With `filterable`, narrows a platform type to its database-backed subset — the
   * properties a filter string may name. A domain table answers the same full column
   * list either way: every domain column is filterable.
   *
   * These are the descriptors the platform renders and filters by. For the richer
   * per-column metadata of a domain table (nullable, choices, min/max, defaultValue,
   * get/set-bound to row values), read {@link DomainRegistryClient.rowProperties}
   * through `grok.dapi.domains.registry`.
   * @example
   * const props = await grok.meta.propertiesOf('User', {filterable: true});
   * const users = await grok.dapi.users.filter(`${props[0].name} = "admin"`).list();
   */
  async propertiesOf(type: string, options?: {filterable?: boolean}): Promise<EntityPropertyInfo[] | null> {
    let props: Property[] | null;
    if (type.includes('.')) {
      try {
        props = await new DomainRegistryClient().rowProperties(type);
      } catch (e) {
        if (!(e instanceof DomainValidationError))
          throw e;
        return null;
      }
    } else
      props = toJs(options?.filterable ? api.grok_Meta_FilterPropertiesOf(type) : api.grok_Meta_GridPropertiesOf(type));
    return props == null ? null : props.map((p) => ({
      name: p.name,
      type: p.propertyType,
      semType: p.semType ?? null,
      friendlyName: p.friendlyName,
      description: p.description ?? null,
      refType: p.refType,
      relationKind: p.relationKind,
    }));
  }

  /** Where [type]'s rows are served as a read-only domain table, or null when the
   * platform does not serve them — pass the result to `grok.dapi.domains.table()`
   * to query, facet or filter them with the names {@link propertiesOf} returned. */
  async coreLocationOf(type: string): Promise<CoreLocation | null> {
    return (await api.grok_Meta_CoreLocationOf(type)) ?? null;
  }
}
