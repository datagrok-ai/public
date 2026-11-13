// was `FilterSchemas.forDomainTable` in u2; moved so u2 never calls `grok.dapi.domains`
import * as grok from 'datagrok-api/grok';
import {propertyFields} from 'datagrok-api/u2core';
import {Filters} from '@datagrok-libraries/u2';
import type {FilterProperty, FilterSchema, FilterValueItem} from '@datagrok-libraries/u2';
import {filterRenderer, filterRef} from '@datagrok-libraries/u2/dg';

const VALUES_LIMIT = 50;
const REF_ADDRESS = /^\w+\.\w+$/;

/** The registry's row properties of `'<schema>.<table>'`; a ref column's semType is its
 * target address. Values come from the table's `categories` facet under the typed text — on a
 * `'basic'` backend (`support.filters`) for string and ref columns only, the ones it searches. */
export async function domainFilterSchema(address: string): Promise<FilterSchema> {
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
      const table = grok.dapi.domains.table(address);
      const kind = Filters.kindOf(prop);
      if (kind !== Filters.KIND.STRING && kind !== Filters.KIND.REF &&
          (await table.access()).support.filters === 'basic')
        return [];
      const {facets} = await table.facets({facets: [{
        id: 'v', kind: 'categories', column: prop.name, search: query, limit: VALUES_LIMIT}]});
      return facets.v.categories.map((c): FilterValueItem => ({
        value: prop.ref ? filterRef(prop.ref, c.value, c.display) : c.value,
        label: c.display ?? String(c.value),
        ...(typeof c.total === 'number' ? {count: c.total} : {}),
      }));
    },
    resolveRef: (prop) => domainFilterSchema(prop.ref!),
    renderer: filterRenderer,
  };
}
