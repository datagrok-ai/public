﻿import {Observable} from "rxjs";

type DartStream = any;

export function debounce<T>(observable: Observable<T>, milliseconds?: number): Observable<T>;

export type ObsReturn = ReturnType<typeof __obs>;

//TODO: any!!
export function __obs(eventId: string, object?: any | null): Observable<any>;


export function observeStream(dartStream: DartStream): Observable<any>;


/** Global platform events. */
export class Events {


    /** @returns {Observable} */ get onContextMenu(): Observable<any>;
    
    /** @returns {Observable} */ get onCurrentViewChanged(): Observable<any>;

    /** @returns {Observable} */ get onCurrentCellChanged(): Observable<any>;

    /** @returns {Observable} */ get onTableAdded(): Observable<any>;

    /** @returns {Observable} */ get onTableRemoved(): Observable<any>;

    /** @returns {Observable} */ get onQueryStarted(): Observable<any>;

    /** @returns {Observable} */ get onQueryFinished(): Observable<any>;

    /** @returns {Observable} */ get onViewChanged(): Observable<any>;

    /** @returns {Observable} */ get onViewAdded(): Observable<any>;

    /** @returns {Observable} */ get onViewRemoved(): Observable<any>;

    /** @returns {Observable} */ get onViewRenamed(): Observable<any>;

    /** @returns {Observable} */ get onCurrentProjectChanged(): Observable<any>;

    /** @returns {Observable} */ get onProjectUploaded(): Observable<any>;

    /** @returns {Observable} */ get onProjectSaved(): Observable<any>;

    /** @returns {Observable} */ get onProjectOpened(): Observable<any>;

    /** @returns {Observable} */ get onProjectClosed(): Observable<any>;

    /** @returns {Observable} */ get onProjectModified(): Observable<any>;

    /** Observes platform events with the specified eventId.
     * @returns {Observable} */
    onEvent(eventId: string): Observable<any>

    /** Observes custom events with the specified eventId.
     * {@link https://public.datagrok.ai/js/samples/events/custom-events}
     * @returns {Observable} */
    onCustomEvent(eventId: string): Observable<any>

    /** Observes events with the specified eventId.
     * {@link https://public.datagrok.ai/js/samples/events/custom-events}
     * @param {string} eventId
     * @param args - event arguments*/
    fireCustomEvent(eventId: string, args: any): void;
}


export class Stream {

    listen(onData: Function): StreamSubscription

    toObservable(): Observable<any>
}


/** Subscription to an event stream. Call [cancel] to stop listening. */
export class StreamSubscription {
    constructor(d: DartStream);

    cancel(): void;
}

/** Event arguments. {@see args} contains event details.
 *  Sample: {@link https://public.datagrok.ai/js/samples/events/global-events}*/
export class EventData {

    /** @type {UIEvent} */
    get causedBy(): UIEvent

    /** Whether the default event handling is prevented. See also {@link preventDefault}
     * @returns {boolean} */
    get isDefaultPrevented(): boolean

    /** Event details. */
    get args(): { [key: string]: object };

    /** Prevents default handling. See also {@link isDefaultPrevented} */
    preventDefault(): void
}

/** Central event hub. */
export class EventBus {

    onEvent(type: string): Observable<any>

    /*

    _getSubject(type) {
        if (!this._streams.has(type)) {
            let s = new rxjs.Subject();
            this._streams.set(type, s);
            return s;
        }

        return this._streams.get(type);
    }*/


    fire(type: string, data: any): void;

}

export function _sub(d: DartStream): StreamSubscription;