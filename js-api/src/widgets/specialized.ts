/**
 * Specialized widgets: FilesWidget, FunctionsWidget, Favorites, VisualDbQueryEditor, Breadcrumbs, DropDown.
 * @module widgets/specialized
 */

import {toDart, toJs} from "../wrappers";
import * as rxjs from "rxjs";
import {fromEvent, Observable} from "rxjs";
import {filter, map} from 'rxjs/operators';
import {__obs, StreamSubscription} from "../events";
import {Func, TableInfo, TableQuery} from "../entities";
import {IDartApi} from "../api/grok_api.g";
import {Grid} from "../grid";
import {DartWidget, TypedEventArgs, Widget} from "./base";
import {TagEditor} from "./tree";
import {fileShares} from "./types";
import {Menu} from "./menu";

import '../../css/breadcrumbs.css';
import '../../css/drop-down.css';

declare let ui: any;
const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;


/** File browser widget */
export class FilesWidget extends DartWidget {
  /** Creates a [FilesWidget] and opens a directory, if [path] is specified.
   * [path] accepts a full-qualified name (see [Entity.nqName]). */
  static create(params: {path?: string, dataSourceFilter?: fileShares[]} = {}): FilesWidget {
    return toJs(api.grok_FilesWidget(params));
  }
}

export class Breadcrumbs {
  path: string[];
  root: HTMLDivElement;

  constructor(path: string[]) {
    this.root = ui.div();
    this.path = path;

    this.root = ui.divH(path.map((element) => ui.div(
      ui.link(element, () => {}, '', 'ui-breadcrumbs-text-element'), 'ui-breadcrumbs-element')), 'ui-breadcrumbs');

    const rootElements = this.root.getElementsByClassName('ui-breadcrumbs-element');
    for (let i = 0; i < rootElements.length - 1; i++)
      rootElements[i].after(ui.iconFA('chevron-right'));
  }


  get onPathClick(): Observable<string[]> {
    const pathElements = this.root.getElementsByClassName('ui-breadcrumbs-text-element');

    return fromEvent<MouseEvent>(pathElements, 'click')
      .pipe(
        map((event) => {
          const currentElement = (event.target as HTMLElement).innerText;
          return this.path.slice(0, this.path.indexOf(currentElement) + 1);
        })
      );
  }
}


/** Menu items definition: functions for items, objects for nested groups.
 * @example
 * {
 *   'Add': () => {},
 *   'Edit': () => {},
 *   'More': { 'Option 1': () => {}, 'Option 2': () => {} }
 * }
 */
export type DropDownMenuItems = {[key: string]: (() => void) | DropDownMenuItems};

/** Builder function for full Menu control (separators, headers, etc.)
 * @example
 * (menu) => {
 *   menu.item('Add', () => {});
 *   menu.separator();
 *   menu.item('Delete', () => {});
 * }
 */
export type DropDownMenuBuilder = (menu: Menu) => void;

/** Options for dropdown with list items */
export interface DropDownListOptions {
  onItemClick?: (item: string, index: number) => void;
}

/** Options for dropdown */
export interface DropDownOptions {
  onExpand?: () => void;
  onCollapse?: () => void;
  /** For custom dropdowns: close when an interactive element is clicked (default: true) */
  closeOnClick?: boolean;
}

export class DropDown extends Widget {
  private _isExpanded: boolean;
  private _currentMenu: Menu | null;
  private _menuCloseSub: rxjs.Subscription | null;
  private _options?: DropDownOptions;
  private _contentContainer?: HTMLElement;

  constructor(label: string | Element, options?: DropDownOptions) {
    const labelElement = typeof label === 'string' ? ui.span([label]) : label;
    super(ui.div([labelElement], 'ui-drop-down-root'));
    this._isExpanded = false;
    this._currentMenu = null;
    this._menuCloseSub = null;
    this._options = options;
  }

  /** Removes the widget and cleans up resources */
  detach(): void {
    this.removeMenuSubscription();
    if (this._contentContainer?.parentNode)
      this._contentContainer.parentNode.removeChild(this._contentContainer);
    super.detach();
  }

  removeMenuSubscription(): void {
    this._menuCloseSub?.unsubscribe();
    this._menuCloseSub = null;
  }

  get isExpanded(): boolean {
    return this._isExpanded;
  }

  expand(): void {
    if (this._isExpanded)
      return;
    this._isExpanded = true;
    this.root.classList.add('ui-drop-down-root-expanded');
    if (this._options?.onExpand)
      this._options.onExpand();
  }

  collapse(hideMenu: boolean = true): void {
    if (!this._isExpanded)
      return;
    this._isExpanded = false;
    this.removeMenuSubscription();
    this.root.classList.remove('ui-drop-down-root-expanded');
    if (hideMenu && this._currentMenu) {
      const menu = this._currentMenu;
      this._currentMenu = null;
      menu.hide();
    }
    else
      this._currentMenu = null;
    if (this._options?.onCollapse)
      this._options.onCollapse();
  }

  private _initMenuDropDown(items: DropDownMenuItems | DropDownMenuBuilder): void {
    const isBuilder = typeof items === 'function';
    const populateMenu = (menu: Menu, menuItems: DropDownMenuItems) => {
      for (const key of Object.keys(menuItems)) {
        const value = menuItems[key];
        if (typeof value === 'function') {
          menu.item(key, () => {
            value();
            this.collapse();
          });
        }
        else if (typeof value === 'object') {
          populateMenu(menu.group(key), value);
          menu.endGroup();
        }
      }
    };

    this.root.addEventListener('click', (e) => {
      if (e.button !== 0)
        return;

      if (this._isExpanded)
        this.collapse();
      else {
        this.expand();
        this._currentMenu = Menu.popup();
        if (isBuilder)
          items(this._currentMenu);
        else
          populateMenu(this._currentMenu, items);
        this._menuCloseSub = this._currentMenu.onClose.subscribe(() => this.collapse(false));

        const rect = this.root.getBoundingClientRect();
        this._currentMenu.show({x: rect.left, y: rect.bottom});
      }
    });
  }

  private _initListDropDown(items: string[], listOptions?: DropDownListOptions): void {
    const menuItems: DropDownMenuItems = {};
    items.forEach((item, index) => {
      menuItems[item] = () => {
        if (listOptions?.onItemClick)
          listOptions.onItemClick(item, index);
      };
    });
    this._initMenuDropDown(menuItems);
  }

  private _initCustomDropDown(createElement: () => HTMLElement): void {
    const contentContainer = ui.div([], 'ui-drop-down-content');
    this._contentContainer = contentContainer;
    const closeOnClick = this._options?.closeOnClick !== false;

    const hideDropDown = () => {
      this.collapse();
      ui.empty(contentContainer);
      contentContainer.classList.remove('ui-drop-down-content-visible');
    };
    const showDropDown = () => {
      this.expand();
      const element = createElement();
      contentContainer.append(element);
      const rect = this.root.getBoundingClientRect();
      contentContainer.style.left = `${rect.left}px`;
      contentContainer.style.top = `${rect.bottom + 2}px`;
      contentContainer.classList.add('ui-drop-down-content-visible');
    };
    document.body.appendChild(contentContainer);

    let justOpened = false;
    this.root.addEventListener('click', (e) => {
      if (e.button !== 0)
        return;
      if (this._isExpanded)
        hideDropDown();
      else {
        showDropDown();
        justOpened = true;
        setTimeout(() => justOpened = false, 0);
      }
    });

    if (closeOnClick) {
      contentContainer.addEventListener('mouseup', (e) => {
        const target = e.target as HTMLElement;
        if (target.matches('h1, h2, h3, h4, h5, h6') || (target.matches('span, div') && !target.closest('button, a, .ui-btn, .d4-link-action')))
          return;
        // The delay to let the button handlers execute first
        setTimeout(() => hideDropDown(), 10);
      }, true);
    }

    this.sub(fromEvent<MouseEvent>(document, 'click').subscribe((event) => {
      if (!this._isExpanded || justOpened)
        return;
      if (this.root.contains(event.target as Node) || contentContainer.contains(event.target as Node))
        return;
      hideDropDown();
    }));
  }

  /** Creates a dropdown with menu items. This is the recommended way to create dropdowns.
   * @param label - Text or element to display as the dropdown trigger
   * @param items - Menu items object, or a builder function for full Menu access
   * @param options - Optional dropdown options
   * @example
   * // Simple object syntax
   * DropDown.menu('Actions', {
   *   'Add': () => grok.shell.info('add'),
   *   'Edit': () => grok.shell.info('edit'),
   * });
   *
   * // Builder function for separators, headers, etc.
   * DropDown.menu('Actions', (menu) => {
   *   menu.item('Add', () => grok.shell.info('add'));
   *   menu.separator();
   *   menu.item('Delete', () => grok.shell.info('delete'));
   * });
   */
  static menu(label: string | Element, items: DropDownMenuItems | DropDownMenuBuilder, options?: DropDownOptions): DropDown {
    const dropdown = new DropDown(label, options);
    dropdown._initMenuDropDown(items);
    return dropdown;
  }

  /** Creates a dropdown with a list of selectable items.
   * @param label - Text or element to display as the dropdown trigger
   * @param items - Array of string items to display
   * @param listOptions - Options including onItemClick callback
   * @param options - Optional dropdown options
   */
  static list(label: string | Element, items: string[], listOptions?: DropDownListOptions, options?: DropDownOptions): DropDown {
    const dropdown = new DropDown(label, options);
    dropdown._initListDropDown(items, listOptions);
    return dropdown;
  }

  /** Creates a dropdown with custom content.
   * @param label - Text or element to display as the dropdown trigger
   * @param createElement - Function that creates the dropdown content element
   * @param options - Optional dropdown options
   */
  static custom(label: string | Element, createElement: () => HTMLElement, options?: DropDownOptions): DropDown {
    const dropdown = new DropDown(label, options);
    dropdown._initCustomDropDown(createElement);
    return dropdown;
  }
}


export class FunctionsWidget extends DartWidget {
  constructor(dart: any) {
    super(dart);
  }

    /** Observes platform events with the specified eventId. */
    onEvent(eventId: string | null = null): rxjs.Observable<any> {
      if (eventId !== null)
        return __obs(eventId, this.dart);

      let dartStream = api.grok_Viewer_Get_EventBus_Events(this.dart);
      return rxjs.fromEventPattern(
        function (handler) {
          return api.grok_Stream_Listen(dartStream, function (x: any) {
            handler(new TypedEventArgs(x));
          });
        },
        function (handler, dart) {
          new StreamSubscription(dart).cancel();
        }
      );
    }

  get onActionClicked(): rxjs.Observable<Func> { return this.onEvent('d4-action-click'); }
  get onActionPlusIconClicked(): rxjs.Observable<Func> { return this.onEvent('d4-action-plus-icon-click'); }
}

export class Favorites extends DartWidget {
  constructor(dart: any) {
    super(dart);
  }

  static async add(entity: any, group?: any): Promise<void> {
    await api.grok_Favorites_Add(entity.dart, group?.dart);
  }

  static async remove(entity: any, group?: any): Promise<void> {
    await api.grok_Favorites_Remove(entity.dart, group?.dart);
  }
}

export class VisualDbQueryEditor extends DartWidget {
  constructor(dart: any) {
    super(dart);
  }

  static fromDbTable(table: TableInfo): VisualDbQueryEditor {
    return api.grok_VisualDbQueryEditor_FromDbTable(toDart(table));
  }

  static fromQuery(query: TableQuery): VisualDbQueryEditor {
    return api.grok_VisualDbQueryEditor_FromQuery(toDart(query));
  }

  get grid(): Grid {
    return api.grok_VisualDbQueryEditor_Get_Grid(this.dart);
  }

  get query(): TableQuery {
    return api.grok_VisualDbQueryEditor_Get_Query(this.dart);
  }

  get inputSchemas(): {[key: string]: TableInfo[]} {
    return api.grok_VisualDbQueryEditor_Get_InputSchemas(this.dart);
  }

  get pivotTag(): TagEditor {
    return api.grok_VisualDbQueryEditor_Get_PivotTag(this.dart);
  }

  get havingTag(): TagEditor {
    return api.grok_VisualDbQueryEditor_Get_HavingTag(this.dart);
  }

  get orderTag(): TagEditor {
    return api.grok_VisualDbQueryEditor_Get_OrderTag(this.dart);
  }

  get whereTag(): TagEditor {
    return api.grok_VisualDbQueryEditor_Get_WhereTag(this.dart);
  }

  get groupByTag(): TagEditor {
    return api.grok_VisualDbQueryEditor_Get_GroupByTag(this.dart);
  }

  get aggregateTag(): TagEditor {
    return api.grok_VisualDbQueryEditor_Get_AggrTag(this.dart);
  }

  get mainTag(): TagEditor {
    return api.grok_VisualDbQueryEditor_Get_MainTag(this.dart);
  }

  getTableInfoName(table: TableInfo): string {
    return api.grok_VisualDbQueryEditor_Get_TableInfoName(this.dart, toDart(table));
  }

  getTableInfoByName(name: string): TableInfo {
    return api.grok_VisualDbQueryEditor_Get_TableInfoByName(this.dart, name);
  }

  refreshQuery(): void {
    api.grok_VisualDbQueryEditor_RefreshQuery(this.dart);
  }

  isInit(): Promise<void> {
    return api.grok_VisualDbQueryEditor_IsInit(this.dart);
  }

  get onChanged(): Observable<any> { return api.grok_VisualDbQueryEditor_OnChanged(this.dart); }

  set showAddToWorkspaceBtn(show: boolean) {
    api.grok_VisualDbQueryEditor_Set_ShowAddToWorkspaceBtn(this.dart, show);
  }

  setSingleColumnMode(column: string): Promise<void> {
    return api.grok_VisualDbQueryEditor_Set_SingleColumnMode(this.dart, column);
  }
}
