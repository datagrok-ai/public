import * as rxjs from "rxjs";
import * as rxjsOperators from "rxjs/operators";
import {toDart, toJs} from "./wrappers";

export function debounce(observable, milliseconds = 100) {
    return observable.pipe(rxjsOperators.debounceTime(milliseconds));
}

export function __obs(eventId, object = null) {
    console.log(`${eventId} initialized.`);

    if (object == null) {
        let observable = rxjs.fromEventPattern(
            function(handler) { return grok_OnEvent(eventId, function (x, w) { handler(w ? toJs(x) : x); }); },
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


export class Events {
    onEvent(eventId) { return __obs(eventId); }
    fireEvent(eventId, arg) {
        return grok_FireEvent(eventId, arg);
    }

    get onCurrentViewChanged () { return __obs('d4-current-view-changed'); }
    get onCurrentCellChanged () { return __obs('d4-current-cell-changed'); }
    get onTableAdded () { return __obs('d4-table-added'); }
    get onTableRemoved () { return __obs('d4-table-removed'); }
    get onQueryStarted () { return __obs('d4-query-started'); }
    get onQueryFinished () { return __obs('d4-query-finished'); }

    get onViewChanged () { return __obs('grok-view-changed'); }
    get onViewAdded () { return __obs('grok-view-added'); }
    get onViewRemoved () { return __obs('grok-view-removed'); }
    get onViewRenamed () { return __obs('grok-view-renamed'); }

    get onCurrentProjectChanged () { return __obs('grok-current-project-changed'); }
    get onProjectUploaded () { return __obs('grok-project-uploaded'); }
    get onProjectSaved () { return __obs('grok-project-saved'); }
    get onProjectOpened () { return __obs('grok-project-opened'); }
    get onProjectClosed () { return __obs('grok-project-closed'); }
    get onProjectModified () { return __obs('grok-project-modified'); }
}


export class Stream {
    constructor(d) { this.d = d; }

    listen(onData) {
        return new StreamSubscription(grok_Stream_Listen(this.d, onData));
    }

    toObservable() {
        let observable = rxjs.fromEventPattern(
            function(handler) { return grok_OnEvent(eventId, function (x, w) { handler(w ? toJs(x) : x); }); },
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

}

export function _sub(d) { return new StreamSubscription(d); }
