import * as rxjs from "rxjs";
import * as rxjsOperators from "rxjs/operators";
import {toJs} from "./wrappers";

export function debounce(observable, milliseconds = 100) {
    return observable.pipe(rxjsOperators.debounceTime(milliseconds));
}

export function __obs(eventId, object = null) {
    console.log(`${eventId} initialized.`);

    if (object == null) {
        let observable = rxjs.fromEventPattern(
            function(handler) { return grok_OnEvent(eventId, function (x) { handler(toJs(x)); }); },
            function(handler, d) { new StreamSubscription(d).cancel(); }
        );
        return observable;
    }
    else {
        let o2 = rxjs.fromEventPattern(
            function(handler) { return grok_OnObjectEvent(object, eventId, function (x) { handler(toJs(x)); }); },
            function(handler, d) { new StreamSubscription(d).cancel(); }
        );
        return o2;
    }
}


export function observeStream(dartStream) {
    let observable = rxjs.fromEventPattern(
        function(handler) { return grok_Stream_Listen(dartStream, function (x) { handler(toJs(x)); }); },
        function(handler, d) { new StreamSubscription(d).cancel(); }
    );
    return observable;
}


/** Global platform events. */
export class Events {

    constructor() {
        this.customEventBus = new EventBus();
    }

    /** Observes platform events with the specified eventId.
     * @returns {Observable} */
    onEvent(eventId) { return __obs(eventId); }

    /** Observes custom events with the specified eventId.
     * @returns {Observable} */
    onCustomEvent(eventId) {
        return this.customEventBus.onEvent(eventId);
    }

    /** Observes events with the specified eventId.
     * @param {string} eventId
     * @param args - event arguments*/
    fireCustomEvent(eventId, args) {
        return this.customEventBus.fire(eventId, args);
    }

    /** @returns {Observable} */ get onCurrentViewChanged () { return __obs('d4-current-view-changed'); }

    /** @returns {Observable} */ get onCurrentCellChanged () { return __obs('d4-current-cell-changed'); }
    /** @returns {Observable} */ get onTableAdded () { return __obs('d4-table-added'); }
    /** @returns {Observable} */ get onTableRemoved () { return __obs('d4-table-removed'); }
    /** @returns {Observable} */ get onQueryStarted () { return __obs('d4-query-started'); }
    /** @returns {Observable} */ get onQueryFinished () { return __obs('d4-query-finished'); }

    /** @returns {Observable} */ get onViewChanged () { return __obs('grok-view-changed'); }
    /** @returns {Observable} */ get onViewAdded () { return __obs('grok-view-added'); }
    /** @returns {Observable} */ get onViewRemoved () { return __obs('grok-view-removed'); }
    /** @returns {Observable} */ get onViewRenamed () { return __obs('grok-view-renamed'); }

    /** @returns {Observable} */ get onCurrentProjectChanged () { return __obs('grok-current-project-changed'); }
    /** @returns {Observable} */ get onProjectUploaded () { return __obs('grok-project-uploaded'); }
    /** @returns {Observable} */ get onProjectSaved () { return __obs('grok-project-saved'); }
    /** @returns {Observable} */ get onProjectOpened () { return __obs('grok-project-opened'); }
    /** @returns {Observable} */ get onProjectClosed () { return __obs('grok-project-closed'); }
    /** @returns {Observable} */ get onProjectModified () { return __obs('grok-project-modified'); }
}


export class Stream {
    constructor(d) { this.d = d; }

    listen(onData) {
        return new StreamSubscription(grok_Stream_Listen(this.d, onData));
    }

    toObservable() {
        let observable = rxjs.fromEventPattern(
            function(handler) { return grok_OnEvent(eventId, function (x) { handler(w ? toJs(x) : x); }); },
            function(handler, streamSubscription) { streamSubscription.cancel(); }
        );
        return observable;
    }
}


/** Subscription to an event stream. Call [cancel] to stop listening. */
export class StreamSubscription {
    constructor(d) { this.d = d; }

    cancel() { grok_Subscription_Cancel(this.d); }
}

export class EventData {
    constructor(d) {
        this.d = d;
    }

    get causedBy() { return grok_EventData_Get_CausedBy(this.d); }
    get isDefaultPrevented() { return grok_EventData_Get_IsDefaultPrevented(this.d); }
    preventDefault() { grok_EventData_PreventDefault(this.d); }

    get args() {
        let x = grok_EventData_Get_Args(this.d);
        let result = {};
        for (const property in x)
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

export function _sub(d) { return new StreamSubscription(d); }
