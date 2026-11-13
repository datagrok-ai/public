import * as grok from 'datagrok-api/grok';
import * as DG from 'datagrok-api/dg';
import {Filters, FilterError} from '../../core/filter/index.js';
import type {FilterRef} from '../../core/filter/model.js';
import type {FilterProperty, FilterSchema} from '../../core/filter/schema.js';
import type {ObjectRenderer} from '../../core/object-renderer.js';
import {moleculeRenderer} from '../inputs/molecule.js';
import {dapiSource} from '../entities/dapi-source.js';
import type {DapiSourceLike} from '../entities/dapi-source.js';
import {HandlerRenderer} from '../entities/entity.js';
import {Platform} from '../shell/platform.js';

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

export function filterRenderer(prop: FilterProperty, refRenderer: ObjectRenderer<unknown> = REF_RENDERER):
  ObjectRenderer<any> | undefined {
  if (prop.semType === DG.SEMTYPE.MOLECULE)
    return moleculeRenderer();
  return prop.ref ? refRenderer : undefined;
}

export function filterRef(type: string, id: unknown, name?: string | null): FilterRef {
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
  /** One property per listed column ({@link Platform.listedColumns}); string columns offer their categories as
   * values, molecules render through the sketcher depiction. The frame's min/max and categories never become
   * `min`/`max`/`choices` — those validate as hard bounds, and a filter may name what the frame lacks. */
  static forDataFrame(df: DG.DataFrame): DataFrameFilterSchema {
    const categories = new Map<string, string[]>();
    const categoriesOf = (col: DG.Column) => {
      let list = categories.get(col.name);
      if (list === undefined)
        categories.set(col.name, list = col.categories.filter((c) => c));
      return list;
    };
    const properties = () => Platform.listedColumns(df).map((col) => {
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
      renderer: filterRenderer,
    };
    return schema;
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
        return entities.map((e) => ({value: filterRef(prop.ref!, e.id, e.friendlyName), label: e.friendlyName}));
      },
      resolveRef: (prop) => FilterSchemas.forEntityType(prop.ref!),
      renderer: (prop) => filterRenderer(prop, refRenderer),
    };
  }
}
