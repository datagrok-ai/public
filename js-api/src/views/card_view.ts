import {View} from "./view";
import {ObjectHandler} from "../../ui";
import {toJs} from "../wrappers";
import {IDartApi} from "../api/grok_api.g";
const api: IDartApi = <any>window;

export enum RENDER_MODE {
  BRIEF = "Brief",
  CARD = "Card",
  GRID = "Grid",
}
/** Base view for working with a collection of objects that reside on the server.
 *  Typically, results are filtered by applying AND operation between two
 *  filters: {@link permanentFilter} (which is set programmatically and is not visible)
 *  and {@link searchValue} entered by the user.
 *
 *  More details on the smart search syntax: {@link https://datagrok.ai/help/datagrok/smart-search}
 *
 */
export class CardView extends View {
  declare dart: any;

  constructor(dart: any) {
    super(dart);
  }

  /** Creates a new CardView object with the specified options. */
  static create(options?: any): CardView {
    return new CardView(api.grok_CardView_Create(options));
  }

  /**
   * User-specified {@link https://datagrok.ai/help/datagrok/smart-search | filter expression}.
   * @type {string} */
  get searchValue(): string { return api.grok_CardView_Get_SearchValue(this.dart); }
  set searchValue(s: string) { api.grok_CardView_Set_SearchValue(this.dart, s); }

  /** Object handler (instructions how to render, drag-and-drop, etc) */
  get meta(): ObjectHandler { return api.grok_CardView_Get_Meta(this.dart); }
  set meta(s: ObjectHandler) { api.grok_CardView_Set_Meta(this.dart, s); }

  /** Semantic type of the items. */
  get objectType(): string { return api.grok_CardView_Get_Type(this.dart); }
  set objectType(s: string) { api.grok_CardView_Set_Type(this.dart, s); }

  get searchFields(): string[] { return toJs(api.grok_CardView_Get_SearchFields(this.dart)); }
  set searchFields(s: string[]) { api.grok_CardView_Set_SearchFields(this.dart, s); }

  /** Programmatically defined invisible
   * {@link https://datagrok.ai/help/datagrok/smart-search | filter expression}.
   *  @type {string} */
  get permanentFilter(): string { return api.grok_CardView_Get_PermanentFilter(this.dart); }
  set permanentFilter(s: string) { api.grok_CardView_Set_PermanentFilter(this.dart, s); }

  /** Category filter properties list */
  get categoryFilters(): { [property: string]: string } { return toJs(api.grok_CardView_Get_CategoryFilters(this.dart)); }
  set categoryFilters(ff: { [property: string]: string }) { api.grok_CardView_Set_CategoryFilters(this.dart, ff); }

  /** Text filter properties list */
  get filters(): { [property: string]: string } { return toJs(api.grok_CardView_Get_Filters(this.dart)); }
  set filters(ff: { [property: string]: string }) { api.grok_CardView_Set_Filters(this.dart, ff); }

  /** Grouping properties list */
  get hierarchy(): string[] { return toJs(api.grok_CardView_Get_Hierarchy(this.dart)); }
  set hierarchy(s: string[]) { api.grok_CardView_Set_Hierarchy(this.dart, s); }

  /** All possible grouping properties list */
  get hierarchyProperties(): { [property: string]: string } { return toJs(api.grok_CardView_Get_HierarchyProperties(this.dart)); }
  set hierarchyProperties(s: { [property: string]: string }) { api.grok_CardView_Set_HierarchyProperties(this.dart, s); }

  /** Grouping mode on */
  get showTree(): boolean { return api.grok_CardView_Get_ShowTree(this.dart); }
  set showTree(s: boolean) { api.grok_CardView_Set_ShowTree(this.dart, s); }

  /** Render mode */
  get renderMode(): RENDER_MODE { return api.grok_CardView_Get_RenderMode(this.dart); }
  set renderMode(s: RENDER_MODE) { api.grok_CardView_Set_RenderMode(this.dart, s); }

  refresh() {api.grok_CardView_Refresh(this.dart);}
  repaint() {api.grok_CardView_Repaint(this.dart);}
}

export class CustomCardView extends CardView {
  constructor(options: any) {
    super(api.grok_CardView_Create(options));
  }
}

/** Projects view */
export class ProjectsView extends CardView {
  /** @constructs ProjectsView */
  constructor(dart: any) {
    super(dart);
  }

  static create(params: object): ProjectsView {
    return new ProjectsView(api.grok_ProjectsView(params));
  }
}


/** Scripts view */
export class ScriptsView extends CardView {
  /** @constructs ProjectsView */
  constructor(dart: any) {
    super(dart);
  }

  static create(params: object): ScriptsView {
    return new ScriptsView(api.grok_ScriptsView(params));
  }
}
