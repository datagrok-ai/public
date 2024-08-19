import {toDart, toJs} from "./wrappers";
import {Schema, EntityType, EntityProperty} from "./entities";
import {IDartApi} from "./api/grok_api.g";
import { Column, DataFrame } from "./dataframe";

const api: IDartApi = <any>window;

/**
 * API for {@link https://datagrok.ai/help/govern/catalog/sticky-meta | Sticky meta}. Allows attaching arbitrary metadata on custom entities.
 * Can be applied to a column of entities based on semtype or column tags.
 * 
 * See samples: {@link https://public.datagrok.ai/js/samples/dapi/sticky-meta-tags}
 * 
 * @example
 * var schema = await dapi.stickyMeta.createSchema('test-schema', {name: 'molecule', matchBy: 'semtype=molecule'}, {name: 'prop', type: 'string'});
 * var metadata = await dapi.stickyMeta.getAllValues(schema, DG.Column.fromList('string', 'key', ['h2o']));
 */
export class StickyMeta {

  /**
   * Fetch available sticky meta configurations
   * @returns list of schemas
   */
  async getSchemas(): Promise<Schema[]> {
    return toJs(await api.grok_Sticky_GetSchemas());
  }

  /**
   * Saves sticky meta values for list of entities.
   * @param {Schema} schema - a configuration of sticky meta to be applied
   * @param {Column} keys - a column with identifiers of entities.
   * This column should be configured with tags so type of entity would be recognizable.
   * @param {DataFrame} values - a dataframe with metadata values to set.
   * Names of columns in this dataframe should match property names of schema.
   */
  async setAllValues(schema: Schema, keys: Column, values: DataFrame): Promise<void> {
    return api.grok_Sticky_SetAllValues(toDart(schema), toDart(keys), toDart(values));
  }

  /**
   * Fetches sticky meta values for list of entities.
   * @param {Schema} schema - configuration of sticky meta to fetch data for.
   * @param {Column} keys - column with identifiers of entities.
   * This column should be configured with tags so type of entity would be recognizable.
   * @returns {DataFrame}, where properties for keys are stored in columns with respective name
   */
  async getAllValues(schema: Schema, keys: Column): Promise<DataFrame> {
    return toJs(api.grok_Sticky_GetAllValues(toDart(schema), toDart(keys)));
  }

  /**
   * Creates a new schema, which is an instance of sticky meta configuration
   * @param {string} name - specifies name of new schema
   * @param {{name: string, matchBy: string[]}} types - list of types that schema is applied to.
   * Every type has a name and a matching expression it applies to.
   * @param {{name: string, type: string[]}} properties - list of typed properties that can be stored within the schema 
   * @returns {Schema}.
   */
  async createSchema(name: string, types: {name: string, matchBy: string}[], properties: {name: string, type: string}[]): Promise<Schema> {
    var schema = Schema.create(name);
    let entityTypes: Array<EntityType> = [];
    for (var typeIdx in types) {
      var typeDesc = types[typeIdx];
      entityTypes.push(EntityType.create(typeDesc.name, typeDesc.matchBy));
    }
    let entityProperties: Array<EntityProperty> = [];
    for (var propIdx in properties) {
      var propDesc = properties[propIdx];
      entityProperties.push(EntityProperty.create(propDesc.name, propDesc.type));
    }
    schema.entityTypes = entityTypes;
    schema.properties = entityProperties;
    return this.saveSchema(schema).then(() => schema);
  }

  /**
   * Saves modified schema.
   * @param {Schema} schema - modified schema
   */
  async saveSchema(schema: Schema): Promise<void> {
    return api.grok_Sticky_SaveSchema(toDart(schema));
  }

  /**
   * Deletes schema by ID
   * @param {string} id - identifier of a schema
   */
  deleteSchema(id: string): Promise<void> {
    return api.grok_Sticky_DeleteSchema(id);
  }
}