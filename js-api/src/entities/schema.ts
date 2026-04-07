/**
 * Schema, EntityType, and HistoryEntry classes.
 * @module entities/schema
 */

import {toDart, toJs} from "../wrappers";
import {IDartApi} from "../api/grok_api.g";
import {EntityProperty} from "./property";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export class HistoryEntry {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  };

  get object(): object { return toJs(api.grok_HistoryEntry_Get_Object(this.dart)); }
  get time(): object { return toJs(api.grok_HistoryEntry_Get_Time(this.dart)); }
}


export class EntityType {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  };

  static create(name: string, matching: string): EntityType {
    return toJs(api.grok_EntityType_Create(toDart(name), toDart(matching)));
  }

  get name(): string { return toJs(api.grok_EntityType_Get_Name(this.dart)); }
  set name(s: string) { api.grok_EntityType_Set_Name(this.dart, toDart(s)); }

  get matching(): string { return toJs(api.grok_EntityType_Get_Matching(this.dart)); }
  set matching(s: string) { api.grok_EntityType_Set_Matching(this.dart, toDart(s)); }
}


/** Represents dynamic property schema, associated with the entity type. */
export class Schema {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  static create(name: string): Schema {
    return toJs(api.grok_Schema_Create(toDart(name)));
  }

  get name(): string { return api.grok_Schema_Get_Name(this.dart); }

  /** Schema properties */
  get properties(): EntityProperty[] { return toJs(api.grok_Schema_Get_Properties(this.dart)); }
  set properties(p: EntityProperty[]) { api.grok_Schema_Set_Properties(this.dart, p); }

  /** Entity types associated with this schema. */
  get entityTypes(): EntityType[] { return toJs(api.grok_Schema_Get_EntityTypes(this.dart)); }
  set entityTypes(et: EntityType[]) { api.grok_Schema_Set_EntityTypes(this.dart, et); }
}
