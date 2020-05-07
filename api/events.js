function __obs(eventId) {
    console.log(`$eventId initialized.`)
    let observable = rxjs.fromEventPattern(
        function(handler) { return grok_OnEvent(eventId, function (x) { handler(_toJs(x)); }); },
        function(handler, streamSubscription) { streamSubscription.cancel(); }
    );
    return observable;
}

class GrokEvents {
    on(eventId) { return __obs(eventId); }

    get onViewChanged() { return __obs('grok-view-changed'); }
    get onCurrentViewChanged() { return __obs('d4-current-view-changed'); }
    get onTableAdded() { return __obs('d4-table-added'); }
    get onTableRemoved() { return __obs('d4-table-removed'); }

}


/** Subscription to an event stream. Call [cancel] to stop listening. */
class StreamSubscription {
    constructor(d) { this.d = d; }

    cancel() { grok_Subscription_Cancel(this.d); }
}

class EventData {
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
            result[property] = _toJs(x[property]);
        return result;
    }

}

/** Central event hub. */
class EventBus {

}

function _sub(d) { return new StreamSubscription(d); }

