import {View} from "./view";
import {ObjectHandler} from "../../ui";
import {toJs} from "../wrappers";
let api = <any>window;

/** Base view for working with a collection of objects that reside on the server.
 *  Typically, results are filtered by applying AND operation between two
 *  filters: {@link permanentFilter} (which is set programmatically and is not visible)
 *  and {@link searchValue} entered by the user.
 *
 *  More details on the smart search syntax: {@link https://datagrok.ai/help/overview/smart-search}
 *
 */
export class CardView extends View {
  dart: any;

  constructor(dart: any) {
    super(dart);
  }

  /** Creates a new CardView object with the specifed options. */
  static create(options?: any): CardView {
    return new CardView(api.grok_CardView_Create(options));
  }

  /**
   * User-specified {@link https://datagrok.ai/help/overview/smart-search | filter expression}.
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
   * {@link https://datagrok.ai/help/overview/smart-search | filter expression}.
   *  @type {string} */
  get permanentFilter(): string { return api.grok_CardView_Get_PermanentFilter(this.dart); }
  set permanentFilter(s: string) { api.grok_CardView_Set_PermanentFilter(this.dart, s); }

  get categoryFilters(): object { return toJs(api.grok_CardView_Get_CategoryFilters(this.dart)); }
  set categoryFilters(ff: object) { api.grok_CardView_Set_CategoryFilters(this.dart, ff); }

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
