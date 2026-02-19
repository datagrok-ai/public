import * as rxjs from 'rxjs';
import {Observable} from 'rxjs';
import * as rxjsOperators from 'rxjs/operators';
import {filter} from 'rxjs/operators';
import {toJs} from './wrappers';
import {FileInfo, Package} from './entities';
import {Accordion, Dialog, InputBase, TreeViewNode, Widget} from "./widgets";
import {View} from './views/view';
import {ViewInfo} from './entities';
import {Viewer} from "./viewer";
import {Column, DataFrame} from "./dataframe";
import {GridCell} from "./grid";
import {IDartApi} from "./api/grok_api.g";
import {LogMessage} from './logger';
import {EVENT_TYPE} from './const';

const api: IDartApi = (typeof window !== 'undefined' ? window : global.window) as any;

export function debounce<T>(observable: rxjs.Observable<T>, milliseconds: number = 100): rxjs.Observable<T> {
  return observable.pipe(rxjsOperators.debounceTime(milliseconds));
}

/**
 * Creates an RxJS Observable from a Dart event stream.
 *
 * **Note on double-underscore prefix**: This function uses `__` to indicate it is a
 * low-level internal utility. While exported for advanced use cases (e.g., custom event
 * handling in packages), prefer using the typed event getters on the {@link Events} class
 * (e.g., `grok.events.onTableAdded`) for type safety and discoverability.
 *
 * @typeParam T - Type of event data emitted
 * @param eventId - Event identifier (use {@link EVENT_TYPE} constants)
 * @param object - Optional object to scope events to (null for global events)
 * @returns RxJS Observable that emits event data
 *
 * @example
 * // Prefer typed getters:
 * grok.events.onTableAdded.subscribe(e => console.log(e.args.dataFrame));
 *
 * // Direct usage (advanced):
 * __obs(EVENT_TYPE.TABLE_ADDED).subscribe(e => console.log(e));
 *
 * @internal
 */
export function __obs<T = any>(eventId: string, object: any = null): Observable<T> {
  if (object == null) {
    return rxjs.fromEventPattern(
      function (handler) {
        return api.grok_OnEvent(eventId, function (x: any) {
          const jso = toJs(x);
          handler(jso);
        });
      },
      function (handler, dart) {
        new StreamSubscription(dart).cancel();
      }
    );
  }
  else {
    return rxjs.fromEventPattern(
      function (handler) {
        return api.grok_OnObjectEvent(object, eventId, function (x: any) {
          handler(toJs(x));
        });
      },
      function (handler, dart) {
        new StreamSubscription(dart).cancel();
      }
    );
  }
}

/**
 * Converts Dart stream to rxjs.Observable.
 * */
export function observeStream(dartStream: any): Observable<any> {
  return rxjs.fromEventPattern(
    function (handler) {
      return api.grok_Stream_Listen(dartStream, function (x: any) {
        handler(toJs(x));
      });
    },
    function (handler, dart) {
      new StreamSubscription(dart).cancel();
    }
  );
}


/** Global platform events. */
export class Events {
  private customEventBus: EventBus;

  constructor() {
    this.customEventBus = new EventBus();
  }

  /** Observes platform events with the specified eventId.
   * Sample: {@link https://public.datagrok.ai/js/samples/ui/ui-events} */
  onEvent(eventId: string): rxjs.Observable<any> {
    return __obs(eventId);
  }

  /** Observes custom events with the specified eventId.
   * Sample: {@link https://public.datagrok.ai/js/samples/events/custom-events} */
  onCustomEvent(eventId: string): rxjs.Observable<any> {
    return this.customEventBus.onEvent(eventId);
  }

  /** Observes events with the specified eventId.
   * To see which events are getting fired, use the Inspector tool. 
   * Open it (Alt+I), go to the "Client Log" tab, and perform the action that you want to 
   * intercept. In the panel, you will see one or more of the events,
   *  click on them to inspect event parameters. To simplify the development process, 
   * we also generate JavaScript code for handling this particular event, 
   * copy-paste it from the context panel into your code if needed.
   * Sample: {@link https://public.datagrok.ai/js/samples/events/custom-events}
   * @param {string} eventId - such as 'd4-current-view-changed'
   * @param args - event arguments*/
  fireCustomEvent(eventId: string, args: any): void { this.customEventBus.fire(eventId, args); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/events/viewer-events} */
  get onContextMenu(): rxjs.Observable<any> { return __obs(EVENT_TYPE.CONTEXT_MENU); }

  get onAIGenerationAbortRequest(): rxjs.Observable<any> { return __obs(EVENT_TYPE.AI_GENERATION_ABORT); }

  get onAIPanelToggleRequest(): rxjs.Observable<Widget> { return __obs(EVENT_TYPE.AI_PANEL_TOGGLE); }

  get onContextMenuClosed(): rxjs.Observable<any> { return __obs(EVENT_TYPE.CONTEXT_MENU_CLOSED); }

  get onCurrentViewChanged(): rxjs.Observable<any> { return __obs(EVENT_TYPE.CURRENT_VIEW_CHANGED); }

  get onCurrentViewChanging(): rxjs.Observable<EventData<ViewArgs>> { return __obs(EVENT_TYPE.CURRENT_VIEW_CHANGING); }

  get onCurrentObjectChanged(): rxjs.Observable<EventData<EventArgs>> { return __obs(EVENT_TYPE.CURRENT_OBJECT_CHANGED); }

  get onCurrentCellChanged(): rxjs.Observable<any> { return __obs(EVENT_TYPE.CURRENT_CELL_CHANGED); }

  get onInputCreated(): rxjs.Observable<InputBase> { return __obs(EVENT_TYPE.INPUT_CREATED); }

  get onDialogShown(): rxjs.Observable<Dialog> { return __obs(EVENT_TYPE.DIALOG_SHOWN); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/events/global-events} */
  get onTableAdded(): rxjs.Observable<EventData<DataFrameArgs>> { return __obs(EVENT_TYPE.TABLE_ADDED); }

  get onTableRemoved(): rxjs.Observable<EventData<DataFrameArgs>> { return __obs(EVENT_TYPE.TABLE_REMOVED); }

  get onQueryStarted(): rxjs.Observable<any> { return __obs(EVENT_TYPE.QUERY_STARTED); }

  get onQueryFinished(): rxjs.Observable<any> { return __obs(EVENT_TYPE.QUERY_FINISHED); }

  get onViewChanged(): rxjs.Observable<any> { return __obs(EVENT_TYPE.VIEW_CHANGED); }

  get onViewChanging(): rxjs.Observable<any> { return __obs(EVENT_TYPE.VIEW_CHANGING); }

  get onViewAdded(): rxjs.Observable<View> { return __obs(EVENT_TYPE.VIEW_ADDED); }

  get onViewAdding(): rxjs.Observable<View> { return __obs(EVENT_TYPE.VIEW_ADDING); }

  get onViewRemoved(): rxjs.Observable<View> { return __obs(EVENT_TYPE.VIEW_REMOVED); }

  get onViewRemoving(): rxjs.Observable<EventData<ViewArgs>> { return __obs(EVENT_TYPE.VIEW_REMOVING); }

  get onViewRenamed(): rxjs.Observable<View> { return __obs(EVENT_TYPE.VIEW_RENAMED); }

  get onResetFilterRequest(): rxjs.Observable<any> { return __obs(EVENT_TYPE.RESET_FILTER_REQUEST); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/events/layout-events} */
  get onViewLayoutGenerated(): rxjs.Observable<ViewInfo> { return __obs(EVENT_TYPE.VIEW_LAYOUT_GENERATED); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/events/layout-events} */
  get onViewLayoutApplying(): rxjs.Observable<ViewInfo> { return __obs(EVENT_TYPE.VIEW_LAYOUT_APPLYING); }

  /** Sample: {@link https://public.datagrok.ai/js/samples/events/layout-events} */
  get onViewLayoutApplied(): rxjs.Observable<ViewInfo> { return __obs(EVENT_TYPE.VIEW_LAYOUT_APPLIED); }

  /** File in the file share has been edited and saved by the user. */
  get onFileEdited(): rxjs.Observable<FileInfo> { return __obs(EVENT_TYPE.FILE_EDITED); }

  get onCurrentProjectChanged(): rxjs.Observable<any> { return __obs(EVENT_TYPE.CURRENT_PROJECT_CHANGED); }

  get onProjectUploaded(): rxjs.Observable<any> { return __obs(EVENT_TYPE.PROJECT_SAVED); }

  get onProjectSaved(): rxjs.Observable<any> { return __obs(EVENT_TYPE.PROJECT_SAVED); }

  get onProjectSaving(): rxjs.Observable<any> { return __obs(EVENT_TYPE.PROJECT_SAVING); }

  get onProjectOpened(): rxjs.Observable<any> { return __obs(EVENT_TYPE.PROJECT_OPENED); }

  get onProjectClosing(): rxjs.Observable<any> { return __obs(EVENT_TYPE.PROJECT_CLOSING); }

  get onProjectClosed(): rxjs.Observable<any> { return __obs(EVENT_TYPE.PROJECT_CLOSED); }

  get onProjectModified(): rxjs.Observable<any> { return __obs(EVENT_TYPE.PROJECT_MODIFIED); }

  get onTooltipRequest(): rxjs.Observable<any> { return __obs(EVENT_TYPE.TOOLTIP_REQUEST); }

  get onTooltipShown(): rxjs.Observable<any> { return __obs(EVENT_TYPE.TOOLTIP_SHOWN); }

  get onTooltipClosed(): rxjs.Observable<any> { return __obs(EVENT_TYPE.TOOLTIP_CLOSED); }

  get onViewerAdded(): rxjs.Observable<EventData<ViewerArgs>> { return __obs(EVENT_TYPE.VIEWER_ADDED); }

  get onViewerClosed(): rxjs.Observable<EventData<ViewerArgs>> { return __obs(EVENT_TYPE.VIEWER_CLOSED); }

  get onFormCreating(): rxjs.Observable<EventData<ColumnsArgs>> { return __obs(EVENT_TYPE.FORM_CREATING); }

  /** You can use it to dynamically add panes for the context panel */
  get onAccordionConstructed(): rxjs.Observable<Accordion> { return __obs(EVENT_TYPE.ACCORDION_CONSTRUCTED); }

  /** Occurs when a package is successfully loaded. */
  get onPackageLoaded(): rxjs.Observable<Package> { return __obs(EVENT_TYPE.PACKAGE_LOADED); }

  /** You can use it to override the default implementation of file import. */
  get onFileImportRequest(): rxjs.Observable<EventData<FileImportArgs>> { return __obs(EVENT_TYPE.FILE_IMPORT_REQUEST); }

  get onGridCellLinkClicked(): rxjs.Observable<EventData<GridCellArgs>> {return __obs(EVENT_TYPE.GRID_CELL_LINK_CLICKED); }

  get onBrowseNodeCreated(): rxjs.Observable<TreeViewNode> {
    return __obs<TreeViewNode>(EVENT_TYPE.TREE_VIEW_NODE_ADDED).pipe(filter(n => n.rootNode.tag == 'Browse' ));
  }

  get onLog(): Observable<LogMessage> { return api.grok_Logger_OnLog(); }

  get onServerMessage(): Observable<IServerMessageEventArgs> { return __obs(EVENT_TYPE.SERVER_MESSAGE); }
}


/** Subscription to an event stream. Call [cancel] to stop listening. */
export class StreamSubscription {
  private dart: any;
  constructor(dart: any) {
    this.dart = dart;
  }

  unsubscribe(): void { this.cancel(); }

  cancel(): void { api.grok_Subscription_Cancel(this.dart); }
}

/** @see Event arguments. {@link args} contains event details.
 *  Sample: {@link https://public.datagrok.ai/js/samples/events/global-events}*/
export class EventData<TArgs = any> {
  public dart: any;

  constructor(dart: any) {
    this.dart = dart;
  }

  /** @type {UIEvent} */
  get causedBy(): UIEvent {
    return api.grok_EventData_Get_CausedBy(this.dart);
  }

  /** Sender - typically, but not always, an HTML Element */
  get sender(): any {
    return toJs(api.grok_EventData_Get_Sender(this.dart));
  }

  /** Whether the default event handling is prevented. See also {@link preventDefault} */
  get isDefaultPrevented(): boolean {
    return api.grok_EventData_Get_IsDefaultPrevented(this.dart);
  }

  /** Prevents default handling. See also {@link isDefaultPrevented}.
   * Sample: {@link https://public.datagrok.ai/js/samples/events/prevented-event} */
  preventDefault(): void {
    api.grok_EventData_PreventDefault(this.dart);
  }

  /** Event details. */
  get args(): TArgs {
    let x = api.grok_EventData_Get_Args(this.dart);
    let result: any = {};
    for (const property in x)
      if (x.hasOwnProperty(property))
        result[property] = toJs(x[property]);
    return result;
  }
}

/** Central event hub. */
export class EventBus {
  private _streams: Map<any, any>;

  constructor() {
    this._streams = new Map();
  }

  onEvent(type: string): rxjs.Observable<any> {
    let subject = this._getSubject(type);

    return new rxjs.Observable(function subscribe(observer) {
      subject.subscribe({
        next: (v: any) => observer.next(v),
        error: (err: any) => observer.error(err),
        complete: () => observer.complete()
      });
    });
  }

  _getSubject(type: string): any {
    if (!this._streams.has(type)) {
      let s = new rxjs.Subject();
      this._streams.set(type, s);
      return s;
    }

    return this._streams.get(type);
  }

  fire(type: string, data: any): void {
    let subject = this._getSubject(type);
    subject.next(data);
  }
}

/**
 * Wraps a Dart stream subscription as a JavaScript StreamSubscription.
 * @param dart - Dart subscription handle
 * @returns StreamSubscription wrapper
 * @internal
 */
export function _sub(dart: any): StreamSubscription {
  return new StreamSubscription(dart);
}

export interface MapChangeArgs<K, V> {
  source: any;
  change: string;
  key: K;
  value: V;
}

export interface ViewerArgs {
  viewer: Viewer;
}

export type FileImportArgs = {
  file: File,
  tables: DataFrame[]
  preventDefault: () => void
}

export interface ViewArgs {
  view: View;
}

export class ColumnsArgs extends EventData {
  get columns(): Column[] {
    return toJs(api.grok_ColumnsArgs_Get_Columns(this.dart));
  }

  set columns(list: Column[]) {
    api.grok_ColumnsArgs_Set_Columns(this.dart, list);
  }
}

export interface GridCellArgs {
  gridCell: GridCell;
  link: string;
}

export interface DataFrameArgs {
  dataFrame: DataFrame;
}

export interface InputArgs {
  input: InputBase;
}

export interface EventArgs {
  value: any;
  trigger: HTMLElement;
  event: UIEvent;
}

export type JsonPrimitive = string | number | boolean | null;
export type JsonObject = { [key: string]: JsonValue };
export type JsonArray = JsonValue[];
export type JsonValue = JsonPrimitive | JsonObject | JsonArray;

export interface IServerMessageEventArgs {
  eventType: string;
  message: any;
}
