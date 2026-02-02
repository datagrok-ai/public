/**
 * Project class.
 * @module entities/project
 */

import {toJs} from "../wrappers";
import {MapProxy} from "../proxies";
import {DataFrame} from "../dataframe";
import {IDartApi} from "../api/grok_api.g";
import {Entity} from "./entity";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export interface ProjectOpenOptions {
  closeAll: boolean;
  openViews: 'all' | 'saved' | 'none';
}

/** Represents a project */
export class Project extends Entity {
  constructor(dart: any) {
    super(dart);

    this.options = new MapProxy(api.grok_Project_Get_Options(this.dart), 'options');
  }

  public options: any;

  static create(): Project {return toJs(api.grok_Project_From_Id(null)); };

  get pictureUrl(): string {
    return api.grok_PictureMixin_Get_PictureUrl(this.dart);
  }

  get path(): string {
    return api.grok_Project_Get_Path(this.dart);
  }

  get isOnServer(): string {
    return api.grok_Project_Get_IsOnServer(this.dart);
  }

  get isLocal(): string {
    return api.grok_Project_Get_IsLocal(this.dart);
  }

  /** Project description */
  get description(): string {
    return api.grok_Project_Description(this.dart);
  }

  /** Project changes flag */
  get isDirty(): boolean {
    return api.grok_Project_IsDirty(this.dart);
  }

  /** Project is empty flag */
  get isEmpty(): boolean {
    return api.grok_Project_IsEmpty(this.dart);
  }

  get isDashboard(): boolean {
    return api.grok_Project_IsDashboard(this.dart);
  }

  get isPackage(): boolean {
    return api.grok_Project_IsPackage(this.dart);
  }

  toMarkup(): string {
    return api.grok_Project_ToMarkup(this.dart);
  }

  /** Opens the project in workspace */
  open(options?: ProjectOpenOptions): Promise<Project> {
    let openViews = options?.openViews ?? 'all';
    return api.grok_Project_Open(this.dart, options?.closeAll ?? false, openViews == 'all' || openViews == 'saved', openViews == 'all');
  }

  get links(): Entity[] {
    return toJs(api.grok_Project_GetRelations(this.dart, true));
  }

  get children(): Entity[] {
    return toJs(api.grok_Project_GetRelations(this.dart, false));
  }

  addLink(entity: Entity | DataFrame): void {
    if (entity instanceof DataFrame)
      entity = entity.getTableInfo();
    api.grok_Project_AddRelation(this.dart, entity.dart, true);
  }

  addChild(entity: Entity |DataFrame): void {
    if (entity instanceof DataFrame)
      entity = entity.getTableInfo();
    api.grok_Project_AddRelation(this.dart, entity.dart, false);
  }

  removeLink(entity: Entity): void {
    api.grok_Project_RemoveRelation(this.dart, entity.dart);
  }

  removeChild(entity: Entity): void {
    api.grok_Project_RemoveRelation(this.dart, entity.dart);
  }

}
