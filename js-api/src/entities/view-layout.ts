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

  /** @constructs ViewLayout */
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

}

export class ViewInfo extends Entity {

  /** @constructs ViewInfo */
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
