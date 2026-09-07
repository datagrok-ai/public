/**
 * Project class.
 * @module entities/project
 */

import {toJs} from "../wrappers";
import {MapProxy} from "../proxies";
import {DataFrame} from "../dataframe";
import {IDartApi} from "../api/grok_api.g";
import {Entity} from "./entity";
import {TableView} from "../views/view";

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
    this.meta = new MapProxy(api.grok_Project_Get_Meta(this.dart), 'meta') as any;
  }

  /** Project options as key-value pairs. */
  public options: any;

  /** Project metadata (`metaParams`) — user-facing key:value bag, e.g. `demoPath`. */
  public meta: {[key: string]: string};

  /** Creates a new, unsaved project. */
  static create(): Project {return toJs(api.grok_Project_From_Id(null)); };

  /** Opens the "Save project" dialog for the given tables only (no workspace scan), with
   * `views[i]` or `layouts[i]` as the layout. Pass `project` to re-publish into it. */
  static showSaveDialog(options: {tables: DataFrame[], views?: (TableView | null)[],
      layouts?: (string | null)[], name?: string, description?: string,
      project?: Project | string}): Promise<Project | null> {
    return (api as any).grok_Project_OpenSaveDialog(
      options.tables.map((t) => t.dart),
      (options.views ?? []).map((v) => v == null ? null : v.dart),
      options.layouts ?? [],
      options.name ?? '', options.description ?? '',
      typeof options.project === 'string' ? options.project : options.project?.id ?? '');
  }

  /** URL of the project picture. */
  get pictureUrl(): string {
    return api.grok_PictureMixin_Get_PictureUrl(this.dart);
  }

  get path(): string {
    return api.grok_Project_Get_Path(this.dart);
  }

  /** Whether the project has been saved to the server. */
  get isOnServer(): string {
    return api.grok_Project_Get_IsOnServer(this.dart);
  }

  /** Whether the project exists only in this session. */
  get isLocal(): string {
    return api.grok_Project_Get_IsLocal(this.dart);
  }

  /** Project description */
  get description(): string {
    return api.grok_Project_Description(this.dart);
  }

  set description(s: string) {
    api.grok_Project_Set_Description(this.dart, s);
  }

  /** Project changes flag */
  get isDirty(): boolean {
    return api.grok_Project_IsDirty(this.dart);
  }

  /** Project is empty flag */
  get isEmpty(): boolean {
    return api.grok_Project_IsEmpty(this.dart);
  }

  /** Whether the project is a dashboard. */
  get isDashboard(): boolean {
    return api.grok_Project_IsDashboard(this.dart);
  }

  /** Whether the project belongs to a package. */
  get isPackage(): boolean {
    return api.grok_Project_IsPackage(this.dart);
  }

  /** True for Spaces — hierarchical containers with their own storage. See {@link grok.dapi.spaces}. */
  get isSpace(): boolean {
    return api.grok_Project_Get_IsSpace(this.dart);
  }

  /** The markup that references this project. */
  toMarkup(): string {
    return api.grok_Project_ToMarkup(this.dart);
  }

  /** Opens the project in workspace */
  open(options?: ProjectOpenOptions): Promise<Project> {
    let openViews = options?.openViews ?? 'all';
    return api.grok_Project_Open(this.dart, options?.closeAll ?? false, openViews == 'all' || openViews == 'saved', openViews == 'all');
  }

  /** Closes the project: removes its tables and views from the workspace. */
  close(): void {
    api.grok_Project_Close(this.dart);
  }

  /** Entities linked to the project (shared, not owned). */
  get links(): Entity[] {
    return toJs(api.grok_Project_GetRelations(this.dart, true));
  }

  /** Entities owned by the project. */
  get children(): Entity[] {
    return toJs(api.grok_Project_GetRelations(this.dart, false));
  }

  /** Links an entity (or the table info of a DataFrame) to the project. */
  addLink(entity: Entity | DataFrame): void {
    if (entity instanceof DataFrame)
      entity = entity.getTableInfo();
    api.grok_Project_AddRelation(this.dart, entity.dart, true);
  }

  /** Adds an entity (or the table info of a DataFrame) as a child of the project. */
  addChild(entity: Entity |DataFrame): void {
    if (entity instanceof DataFrame)
      entity = entity.getTableInfo();
    api.grok_Project_AddRelation(this.dart, entity.dart, false);
  }

  /** Removes a link. */
  removeLink(entity: Entity): void {
    api.grok_Project_RemoveRelation(this.dart, entity.dart);
  }

  /** Removes a child. */
  removeChild(entity: Entity): void {
    api.grok_Project_RemoveRelation(this.dart, entity.dart);
  }

  /** Adds a table to this project and opens it as a view */
  addTableView(table: DataFrame): TableView {
    return toJs(api.grok_Project_AddTableView(this.dart, table.dart));
  }
}
