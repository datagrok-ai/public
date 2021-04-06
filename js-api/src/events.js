import * as rxjs from 'rxjs';
import * as rxjsOperators from 'rxjs/operators';
import { toJs } from './wrappers';
let api = window;
export function debounce(observable, milliseconds = 100) {
  return observable.pipe(rxjsOperators.debounceTime(milliseconds));
}
export function __obs(eventId, object = null) {
  console.log(`${eventId} initialized.`);
  if (object == null) {
    let observable = rxjs.fromEventPattern(function (handler) {
      return api.grok_OnEvent(eventId, function (x) {
        handler(toJs(x));
      });
    }, function (handler, d) {
      new StreamSubscription(d).cancel();
    });
    return observable;
  }
  else {
    let o2 = rxjs.fromEventPattern(function (handler) {
      return api.grok_OnObjectEvent(object, eventId, function (x) {
        handler(toJs(x));
      });
    }, function (handler, d) {
      new StreamSubscription(d).cancel();
    });
    return o2;
  }
}
/**
 * Converts Dart stream to rxjs.Observable.
 * @param {Object} dartStream
 * @returns {rxjs.Observable}
 * */
export function observeStream(dartStream) {
  let observable = rxjs.fromEventPattern(function (handler) {
    return api.grok_Stream_Listen(dartStream, function (x) {
      handler(toJs(x));
    });
  }, function (handler, d) {
    new StreamSubscription(d).cancel();
  });
  return observable;
}
/** Global platform events. */
export class Events {
  constructor() {
    this.customEventBus = new EventBus();
  }
  /** Observes platform events with the specified eventId.
     * @returns {rxjs.Observable} */
  onEvent(eventId) {
    return __obs(eventId);
  }
  /** Observes custom events with the specified eventId.
     * {@link https://public.datagrok.ai/js/samples/events/custom-events}
     * @returns {rxjs.Observable} */
  onCustomEvent(eventId) {
    return this.customEventBus.onEvent(eventId);
  }
  /** Observes events with the specified eventId.
     * {@link https://public.datagrok.ai/js/samples/events/custom-events}
     * @param {string} eventId
     * @param args - event arguments*/
  fireCustomEvent(eventId, args) {
    return this.customEventBus.fire(eventId, args);
  }
  /** @returns {rxjs.Observable} */ get onContextMenu() {
    return __obs('d4-context-menu');
  }
  /** @returns {rxjs.Observable} */ get onContextMenuClosed() {
    return __obs('d4-menu-closed');
  }
  /** @returns {rxjs.Observable} */ get onCurrentViewChanged() {
    return __obs('d4-current-view-changed');
  }
  /** @returns {rxjs.Observable} */ get onCurrentCellChanged() {
    return __obs('d4-current-cell-changed');
  }
  /** @returns {rxjs.Observable} */ get onTableAdded() {
    return __obs('d4-table-added');
  }
  /** @returns {rxjs.Observable} */ get onTableRemoved() {
    return __obs('d4-table-removed');
  }
  /** @returns {rxjs.Observable} */ get onQueryStarted() {
    return __obs('d4-query-started');
  }
  /** @returns {rxjs.Observable} */ get onQueryFinished() {
    return __obs('d4-query-finished');
  }
  /** @returns {rxjs.Observable} */ get onViewChanged() {
    return __obs('grok-view-changed');
  }
  /** @returns {rxjs.Observable} */ get onViewAdded() {
    return __obs('grok-view-added');
  }
  /** @returns {rxjs.Observable} */ get onViewRemoved() {
    return __obs('grok-view-removed');
  }
  /** @returns {rxjs.Observable} */ get onViewRenamed() {
    return __obs('grok-view-renamed');
  }
  /** @returns {rxjs.Observable} */ get onViewLayoutGenerated() {
    return __obs('d4-view-layout-generated');
  }
  /** @returns {rxjs.Observable} */ get onViewLayoutApplying() {
    return __obs('d4-view-layout-applying');
  }
  /** @returns {rxjs.Observable} */ get onViewLayoutApplied() {
    return __obs('d4-view-layout-applied');
  }
  /** @returns {rxjs.Observable} */ get onCurrentProjectChanged() {
    return __obs('grok-current-project-changed');
  }
  /** @returns {rxjs.Observable} */ get onProjectUploaded() {
    return __obs('grok-project-uploaded');
  }
  /** @returns {rxjs.Observable} */ get onProjectSaved() {
    return __obs('grok-project-saved');
  }
  /** @returns {rxjs.Observable} */ get onProjectOpened() {
    return __obs('grok-project-opened');
  }
  /** @returns {rxjs.Observable} */ get onProjectClosed() {
    return __obs('grok-project-closed');
  }
  /** @returns {rxjs.Observable} */ get onProjectModified() {
    return __obs('grok-project-modified');
  }
  /** @returns {rxjs.Observable} */ get onTooltipRequest() {
    return __obs('d4-tooltip-request');
  }
  /** @returns {rxjs.Observable} */ get onTooltipShown() {
    return __obs('d4-tooltip-shown');
  }
  /** @returns {rxjs.Observable} */ get onTooltipClosed() {
    return __obs('d4-tooltip-closed');
  }
  /** @returns {rxjs.Observable} */ get onViewerAdded() {
    return __obs('d4-viewer-added');
  }
  /** @returns {rxjs.Observable} */ get onViewerClosed() {
    return __obs('d4-viewer-closed');
  }
}
/*

export class Stream {
  private d: any;
  constructor(d: any) {
    this.d = d;
  }

  listen(onData: any) {
    return new StreamSubscription(api.grok_Stream_Listen(this.d, onData));
  }

  toObservable() {
    let observable = rxjs.fromEventPattern(
      function (handler) {
        return api.grok_OnEvent(eventId, function (x) {
          handler(w ? toJs(x) : x);
        });
      },
      function (handler, streamSubscription) {
        streamSubscription.cancel();
      }
    );
    return observable;
  }
}

*/
/** Subscription to an event stream. Call [cancel] to stop listening. */
export class StreamSubscription {
  constructor(d) {
    this.d = d;
  }
  unsubscribe() {
    this.cancel();
  }
  cancel() {
    api.grok_Subscription_Cancel(this.d);
  }
}
/** Event arguments. {@see args} contains event details.
 *  Sample: {@link https://public.datagrok.ai/js/samples/events/global-events}*/
export class EventData {
  constructor(d) {
    this.d = d;
  }
  /** @type {UIEvent} */
  get causedBy() {
    return api.grok_EventData_Get_CausedBy(this.d);
  }
  /** Whether the default event handling is prevented. See also {@link preventDefault}
     * @returns {boolean} */
  get isDefaultPrevented() {
    return api.grok_EventData_Get_IsDefaultPrevented(this.d);
  }
  /** Prevents default handling. See also {@link isDefaultPrevented} */
  preventDefault() {
    api.grok_EventData_PreventDefault(this.d);
  }
  /** Event details. */
  get args() {
    let x = api.grok_EventData_Get_Args(this.d);
    let result = {};
    for (const property in x)
      if (x.hasOwnProperty(property))
        result[property] = toJs(x[property]);
    return result;
  }
}
/** Central event hub. */
export class EventBus {
  constructor() {
    this._streams = new Map();
  }
  onEvent(type) {
    let subject = this._getSubject(type);
    return new rxjs.Observable(function subscribe(observer) {
      subject.subscribe({
        next: (v) => observer.next(v),
        error: (err) => observer.error(err),
        complete: () => observer.complete()
      });
    });
  }
  _getSubject(type) {
    if (!this._streams.has(type)) {
      let s = new rxjs.Subject();
      this._streams.set(type, s);
      return s;
    }
    return this._streams.get(type);
  }
  fire(type, data) {
    let subject = this._getSubject(type);
    subject.next(data);
  }
}
export function _sub(d) {
  return new StreamSubscription(d);
}
