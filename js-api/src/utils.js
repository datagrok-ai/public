import {Balloon} from "./widgets";
import * as rxjs from 'rxjs';

export function _jsThen(promise, f) {
    promise.then(f);
}

export function _toJson(x) {
    return x === null ? null : JSON.stringify(x);
}

export function *range(length) {
    for (let i = 0; i < length; i++)
        yield i;
}

/** Returns an 'identity' array where the element in idx-th position is equals to idx. */
export function identity(length) {
    let res = new Array(length);
    for (let i = 0; i < length; i++)
        res[i] = i;
    return res;
}

/*window.onerror = function (message, url, lineNumber, columnNumber, errorObject) {
    return grok_Error(message, url, lineNumber, columnNumber, errorObject);
};*/

export function time(s, f) {
    let start = new Date();
    let result = f();
    let stop = new Date();
    console.log(`${s}: ${stop - start}ms`);
    Balloon.info(`${s}: ${stop - start}ms`);
    return result;
}

/** @returns {rxjs.Observable} */
export function _onSizeChanged(element) {
    return rxjs.Observable.create(function(observer) {
        const resizeObserver = new ResizeObserver(observerEntries => {
            // trigger a new item on the stream when resizes happen
            for (const entry of observerEntries) {
                observer.next(entry);
            }
        });

        // start listening for resize events
        resizeObserver.observe(element);

        // cancel resize observer on cancelation
        return () => resizeObserver.disconnect();
    });
}