/**
 * ViewLayout and ViewInfo classes.
 * @module entities/view-layout
 */

import {toJs} from "../wrappers";
import {IDartApi} from "../api/grok_api.g";
import {View} from "../views/view";
import {Entity} from "./entity";
import {ColumnInfo, TableInfo} from "./table-info";

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


export class ViewLayout extends Entity {


  constructor(dart: any) {
    super(dart);
  }

  static fromJson(json: string): ViewLayout {
    return toJs(api.grok_ViewLayout_FromJson(json));
  }

  static fromViewState(state: string): ViewLayout {
    return toJs(api.grok_ViewLayout_FromViewState(state));
  }

  get viewState(): string {
    return api.grok_ViewLayout_Get_ViewState(this.dart);
  }

  set viewState(state: string) {
    api.grok_ViewLayout_Set_ViewState(this.dart, state);
  }

  getUserDataValue(key: string): string {
    return api.grok_ViewLayout_Get_UserDataValue(this.dart, key);
  }

  setUserDataValue(key: string, value: string) {
    return api.grok_ViewLayout_Set_UserDataValue(this.dart, key, value);
  }

  toJson(): string {
    return api.grok_ViewLayout_ToJson(this.dart);
  }

  get columns(): ColumnInfo[] {
    return toJs(api.grok_ViewLayout_Get_Columns(this.dart));
  }

  /** Names of the columns the layout's viewers are bound to, read from its state the way the platform reads a
   * viewer's settings: every `…ColumnName` and `…ColumnNames` setting, the Filters viewer's filters and the inner
   * viewers of a trellis plot (a viewer whose package is not loaded, from its stored settings). Each name comes once,
   * whatever its case. The grid binds none, since it shows whatever columns the table has. Known for any
   * layout, one made with {@link fromViewState} too, unlike {@link columns}, which only a layout saved from a table
   * view carries. Reading it creates no viewer.
   * @example
   * const layout = DG.ViewLayout.fromViewState(grok.shell.tv.saveLayout().viewState);
   * grok.shell.info(layout.boundColumnNames.join(', '));
   * @see {@link https://public.datagrok.ai/js/samples/ui/views/layout-bound-columns} */
  get boundColumnNames(): string[] {
    return toJs(api.grok_ViewLayout_Get_BoundColumnNames(this.dart));
  }
}

export class ViewInfo extends Entity {


  constructor(dart: any) {
    super(dart);
  }

  static fromJson(json: string): ViewInfo {
    return new ViewInfo(api.grok_ViewInfo_FromJson(json));
  }

  static fromViewState(state: string): ViewInfo {
    return new ViewInfo(api.grok_ViewInfo_FromViewState(state));
  }

  get table() : TableInfo {
    return toJs(api.grok_ViewInfo_Get_Table(this.dart));
  }

  /** Only defined within the context of the OnViewLayoutXXX events */
  get view(): View {
    return toJs(api.grok_ViewInfo_Get_View(this.dart));
  }

  get viewState(): string {
    return api.grok_ViewInfo_Get_ViewState(this.dart);
  }

  set viewState(state: string) {
    api.grok_ViewInfo_Set_ViewState(this.dart, state);
  }

  getUserDataValue(key: string): string {
    return api.grok_ViewInfo_Get_UserDataValue(this.dart, key);
  }

  setUserDataValue(key: string, value: string) {
    return api.grok_ViewInfo_Set_UserDataValue(this.dart, key, value);
  }

  toJson(): string {
    return api.grok_ViewInfo_ToJson(this.dart);
  }
}
