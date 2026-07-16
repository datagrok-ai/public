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

  public options: any;

  /** Project metadata (`metaParams`) — user-facing key:value bag, e.g. `demoPath`. */
  public meta: {[key: string]: string};

  static create(): Project {return toJs(api.grok_Project_From_Id(null)); };

  /** Opens the platform "Save project" dialog for an externally assembled set of
   * tables — the standard dashboard-publishing UI (per-table data-sync toggles,
   * creation-script dependencies, share link, upload) without scanning the
   * workspace: only the passed tables are offered. `views[i]`, when given, is
   * the (possibly detached) TableView showing `tables[i]`; its layout is saved
   * with the project and reopens with the table. When `views[i]` is absent,
   * `layouts[i]` (a `ViewLayout.viewState` string) ships the layout without a
   * live view. Data sync defaults on for tables whose `df` carries a creation
   * script (`.script` tag). Pass `project` (a previously saved project or its
   * id) to re-publish into the SAME project — its table/view children are
   * replaced by the fresh ones instead of creating a new project per save.
   * @returns the saved project, or null when the dialog is cancelled. */
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

  /** True for Spaces — hierarchical containers with their own storage. See {@link grok.dapi.spaces}. */
  get isSpace(): boolean {
    return api.grok_Project_Get_IsSpace(this.dart);
  }

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

  /** Adds a table to this project and opens it as a view */
  addTableView(table: DataFrame): TableView {
    return toJs(api.grok_Project_AddTableView(this.dart, table.dart));
  }
}
