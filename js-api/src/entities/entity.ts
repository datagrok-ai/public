/**
 * Base Entity class for all database-stored objects.
 * @module entities/entity
 */

import {toJs} from "../wrappers";
import {IDartApi} from "../api/grok_api.g";
import dayjs from "dayjs";

declare var grok: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

/** @class
 * Base class for system objects stored in the database in a structured manner.
 * Contains base properties: id, name and path
 * */
export class Entity {

  public dart: any;

  /** @constructs Entity*/
  constructor(dart: any) {
    this.dart = dart;
  }

  /** Entity ID (GUID)
   *  @type {string} */
  get id(): string { return api.grok_Entity_Get_Id(this.dart); }
  set id(x: string) { api.grok_Entity_Set_Id(this.dart, x); }

  /** Generates new {@link id} for this entity. */
  newId(): void { api.grok_Entity_New_Id(this.dart); }

  /** Entity friendly name
   *  @type {string} */
  get friendlyName(): string { return api.grok_Entity_Get_FriendlyName(this.dart); }
  set friendlyName(x: string) { api.grok_Entity_Set_FriendlyName(this.dart, x); }

  /** Entity short name */
  get name(): string { return api.grok_Entity_Get_Name(this.dart); }
  set name(x: string) { api.grok_Entity_Set_Name(this.dart, x); }

  /** Entity full-qualified name */
  get nqName(): string { return api.grok_Entity_Get_nqName(this.dart); }

  /** Entity path */
  get path(): string { return api.grok_Entity_Path(this.dart); }

  /** Time when entity was created **/
  get createdOn(): dayjs.Dayjs { return dayjs(api.grok_Entity_Get_CreatedOn(this.dart)); }

  /** Time when entity was updated **/
  get updatedOn(): dayjs.Dayjs | null {
    const d = api.grok_Entity_Get_UpdatedOn(this.dart);
    return d ? dayjs(d) : null;
  }

  /** Who created entity **/
  get author(): any {
    // Import User dynamically to avoid circular dependency
    const {User} = require('./user');
    return new User(api.grok_Entity_Get_Author(this.dart));
  }

  /** Entity type name **/
  get entityType(): string { return api.grok_Entity_Get_EntityType(this.dart); }

  /** Gets entity properties */
  getProperties(): Promise<{[index: string]: any}> {
    return api.grok_EntitiesDataSource_GetProperties(grok.dapi.entities.dart, this.dart);
  }

  /** Sets entity properties */
  setProperties(props: {[index: string]: any}): Promise<any> {
    return api.grok_EntitiesDataSource_SetProperties(grok.dapi.entities.dart, this.dart, props);
  }

  /** Returns a string representing the object */
  toString(): string { return api.grok_Object_ToString(this.dart); }

  hasTag(tag: string): boolean { return api.grok_Entity_Has_Tag(this.dart, tag); }

  /** Adds a specified tag */
  tag(tag: string): boolean { return api.grok_Entity_Tag(this.dart, tag); }

  /** Removes a specified tag */
  unTag(tag: string): boolean { return api.grok_Entity_UnTag(this.dart, tag); }
}
