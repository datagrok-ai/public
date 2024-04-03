import {toDart, toJs} from "./wrappers";
import {Schema, EntityType, EntityProperty} from "./entities";
import {IDartApi} from "./api/grok_api.g";
import { Column, DataFrame } from "./dataframe";

const api: IDartApi = <any>window;

export class StickyMeta {

  async getSchemas(): Promise<Schema[]> {
    return toJs(await api.grok_Sticky_GetSchemas());
  }

  setAllValues(schema: Schema, keys: Column, values: DataFrame): Promise<void> {
    return api.grok_Sticky_SetAllValues(toDart(schema), toDart(keys), toDart(values));
  }

  getAllValues(schema: Schema, keys: Column): Promise<DataFrame> {
    return toJs(api.grok_Sticky_GetAllValues(toDart(schema), toDart(keys)));
  }

  createSchema(name: string, types: Map<String, String>[], properties: Map<String, String>[]): Promise<Schema> {
    var schema = Schema.create(name);
    let entityTypes: Array<EntityType> = [];
    for (var typeIdx in types) {
      var typeDesc = types[typeIdx];
      entityTypes.push(EntityType.create(typeDesc['name'], typeDesc['matchBy']));
    }
    let entityProperties: Array<EntityProperty> = [];
    for (var propIdx in properties) {
      var propDesc = properties[propIdx];
      entityProperties.push(EntityProperty.create(propDesc['name'], propDesc['type']));
    }
    schema.entityTypes = entityTypes;
    schema.properties = entityProperties;
    return this.saveSchema(schema).then(() => schema);
  }

  saveSchema(schema: Schema): Promise<void> {
    return api.grok_Sticky_SaveSchema(toDart(schema));
  }

  deleteSchema(id: string): Promise<void> {
    return api.grok_Sticky_DeleteSchema(id);
  }
}