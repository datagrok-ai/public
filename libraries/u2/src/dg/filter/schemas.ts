import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {propertyFields} from '../../core/property-like.js';
import {Filters, FilterError} from '../../core/filter/index.js';
import type {FilterRef} from '../../core/filter/model.js';
import type {FilterProperty, FilterSchema, FilterValueItem} from '../../core/filter/schema.js';
import type {ObjectRenderer} from '../../core/object-renderer.js';
import {moleculeRenderer} from '../inputs/molecule.js';
import {dapiSource} from '../entities/dapi-source.js';
import type {DapiSourceLike} from '../entities/dapi-source.js';
import {HandlerRenderer} from '../entities/entity.js';

const VALUES_LIMIT = 50;
const REF_ADDRESS = /^\w+\.\w+$/;

/** The entity types whose rows a dapi collection lists — the ref values of {@link FilterSchemas.forEntityType}. */
const ENTITY_SOURCES: Record<string, () => DapiSourceLike<DG.Entity>> = {
  User: () => grok.dapi.users,
  Group: () => grok.dapi.groups,
  UserGroup: () => grok.dapi.groups,
  Project: () => grok.dapi.projects,
  DataConnection: () => grok.dapi.connections,
  DataQuery: () => grok.dapi.queries,
  Script: () => grok.dapi.scripts,
  Package: () => grok.dapi.packages,
};

const REF_RENDERER: ObjectRenderer<unknown> = {
  caption: (x) => Filters.isRef(x) ? x.name ?? x.id : String(x),
};

/** Entity values render through their ObjectHandler (icon, markup, tooltip); a `FilterRef` — the
 * id and display name the schema's values carry — is captioned, never handed to a handler. */
class EntityRefRenderer extends HandlerRenderer<unknown> {
  handlerFor(x: unknown): DG.ObjectHandler | null {
    return Filters.isRef(x) ? null : super.handlerFor(x);
  }

  caption(x: unknown): string {
    return Filters.isRef(x) ? x.name ?? x.id : super.caption(x);
  }
}

function rendererFor(prop: FilterProperty, refRenderer: ObjectRenderer<unknown> = REF_RENDERER):
  ObjectRenderer<any> | undefined {
  if (prop.semType === 'Molecule')
    return moleculeRenderer();
  return prop.ref ? refRenderer : undefined;
}

function ref(type: string, id: unknown, name?: string | null): FilterRef {
  const r: FilterRef = {type, id: String(id)};
  if (name != null)
    r.name = name;
  return r;
}

export interface DataFrameFilterSchema extends FilterSchema {
  /** Re-reads the columns — a semantic type detected after the snapshot, a rename, an added
   * column; `properties` is replaced, so a holder compares it by identity. */
  refresh(): void;
}

export class FilterSchemas {
  /** One property per column; string columns offer their categories as values, molecules render
   * through the sketcher depiction. The frame's min/max and categories never become
   * `min`/`max`/`choices` — those validate as hard bounds, and a filter may name what the frame lacks. */
  static forDataFrame(df: DG.DataFrame): DataFrameFilterSchema {
    const categories = new Map<string, string[]>();
    const categoriesOf = (col: DG.Column) => {
      let list = categories.get(col.name);
      if (list === undefined)
        categories.set(col.name, list = col.categories.filter((c) => c));
      return list;
    };
    const properties = () => df.columns.toList().map((col) => {
      const prop: FilterProperty = {name: col.name, type: col.type};
      if (col.semType)
        prop.semType = col.semType;
      return prop;
    });
    const schema: DataFrameFilterSchema = {
      properties: properties(),
      refresh: () => {
        categories.clear();
        schema.properties = properties();
      },
      values: (prop, query) => {
        const col = df.columns.byName(prop.name);
        const q = query.toLowerCase();
        const items = col?.type === DG.COLUMN_TYPE.STRING ?
          categoriesOf(col).filter((c) => c.toLowerCase().includes(q)) : [];
        return Promise.resolve(items.map((value) => ({value, label: value})));
      },
      renderer: rendererFor,
    };
    return schema;
  }

  /** The registry's row properties of `'<schema>.<table>'`; a ref column's semType is its
   * target address. Values come from the table's `categories` facet under the typed text. */
  static async forDomainTable(address: string): Promise<FilterSchema> {
    const registry = grok.dapi.domains.registry;
    const properties = (await registry.rowProperties(address)).map((p) => {
      const prop: FilterProperty = {...propertyFields(p), name: p.name};
      if (REF_ADDRESS.test(p.semType ?? ''))
        prop.ref = p.semType;
      return prop;
    });
    return {
      properties,
      values: async (prop, query) => {
        const {facets} = await grok.dapi.domains.table(address).facets({facets: [{
          id: 'v', kind: 'categories', column: prop.name, search: query, limit: VALUES_LIMIT}]});
        return facets.v.categories.map((c): FilterValueItem => ({
          value: prop.ref ? ref(prop.ref, c.value, c.display) : c.value,
          label: c.display ?? String(c.value),
          count: c.total,
        }));
      },
      resolveRef: (prop) => FilterSchemas.forDomainTable(prop.ref!),
      renderer: rendererFor,
    };
  }

  /** The properties a filter string may name on a platform type (`'User'`, `'Project'`); a ref
   * property points at its `refType`, and its values are that type's dapi collection. */
  static async forEntityType(type: string): Promise<FilterSchema> {
    const infos = await grok.meta.propertiesOf(type, {filterable: true});
    if (infos == null)
      throw new FilterError(`No filterable properties for type "${type}"`);
    const properties = infos.map((p) => {
      const prop: FilterProperty = {name: p.name, type: p.type, friendlyName: p.friendlyName};
      if (p.semType)
        prop.semType = p.semType;
      if (p.description)
        prop.description = p.description;
      if (p.refType)
        prop.ref = p.refType;
      return prop;
    });
    const refRenderer = new EntityRefRenderer();
    return {
      properties,
      values: async (prop, query, signal) => {
        const source = prop.ref ? ENTITY_SOURCES[prop.ref] : undefined;
        if (!source)
          return [];
        const entities = await dapiSource(source)(query, signal);
        return entities.map((e) => ({value: ref(prop.ref!, e.id, e.friendlyName), label: e.friendlyName}));
      },
      resolveRef: (prop) => FilterSchemas.forEntityType(prop.ref!),
      renderer: (prop) => rendererFor(prop, refRenderer),
    };
  }
}
