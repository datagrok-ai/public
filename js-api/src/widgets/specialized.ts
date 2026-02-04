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


export class DropDown extends Widget {
  private _dropDownElement: HTMLDivElement;
  private _isMouseOverElement: boolean;
  isExpanded: boolean;
  createElement: () => HTMLElement;

  constructor(label: string | Element, createElement: () => HTMLElement) {
    const dropDownElement = ui.div(ui.div([], 'ui-drop-down-content'), 'ui-combo-drop-down-fixed');
    super(ui.div([dropDownElement, label], 'ui-drop-down-root'));

    this._isMouseOverElement = false;
    this.isExpanded = false;
    this.createElement = createElement;
    this._dropDownElement = dropDownElement;
    this._dropDownElement.style.visibility = 'hidden';
    this._initEventListeners();
  }

  private _initEventListeners() {
    this.root.addEventListener('mousedown', (e) => {
      // check if the button is LMB
      if (e.button !== 0 || this._isMouseOverElement)
        return;
      this._setExpandedState(this.isExpanded);
    });

    this._dropDownElement.addEventListener('mouseover', () => {
      this._isMouseOverElement = true;
    }, false);

    this._dropDownElement.addEventListener('mouseleave', () => {
      this._isMouseOverElement = false;
    }, false);

    this.sub(fromEvent<MouseEvent>(document, 'click').subscribe((event) => {
      if (this.root.contains(event.target as Node))
        return;
      if (!this.isExpanded)
        return;
      this._setExpandedState(this.isExpanded);
    }));
  }

  private _setExpandedState(isExpanded: boolean) {
    this.isExpanded = !isExpanded;
    if (isExpanded) {
      this.root.classList.remove('ui-drop-down-root-expanded');
      ui.empty(this._dropDownElement);
      this._dropDownElement.style.visibility = 'hidden';
      return;
    }

    this.root.classList.add('ui-drop-down-root-expanded');
    const element = this.createElement();
    element.classList.add('ui-drop-down-content');    // this is not right - we should not change it
    this._dropDownElement.append(element);
    this._dropDownElement.style.visibility = 'visible';
  }


  get onExpand(): Observable<MouseEvent> {
    return fromEvent<MouseEvent>(this.root, 'click').pipe(filter(() => !this._isMouseOverElement));
  }

  get onElementClick(): Observable<MouseEvent> {
    return fromEvent<MouseEvent>(this._dropDownElement, 'click');
  }

  // TODO: add list constructor to DropDown
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
